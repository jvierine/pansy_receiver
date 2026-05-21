#!/usr/bin/env python3
"""Restart pansy_uhd_rx when raw-voltage output or phasecal goes bad."""
import datetime as dt
import math
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np


def load_env(path):
    values = {}
    if not path.exists():
        return values
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key] = value.strip().strip('"')
    return values


def config():
    cfg_path = Path(os.environ.get("PANSY_RECEIVER_CONFIG", Path.home() / ".config/pansy-receiver/pansy-receiver.env"))
    file_cfg = load_env(cfg_path)

    def get(name, default):
        return os.environ.get(name, file_cfg.get(name, default))

    return {
        "repo": Path(get("PANSY_RECEIVER_REPO", Path.home() / "src/git/pansy_receiver")),
        "service": get("PANSY_UHD_RX_SERVICE", "pansy-uhd-rx.service"),
        "raw_root": Path(get("PANSY_UHD_RX_OUTDIR", "/media/archive")),
        "channels": [ch.strip() for ch in get("PANSY_WATCH_CHANNELS", "ch000,ch001,ch002,ch003,ch004,ch005,ch006,ch007").split(",") if ch.strip()],
        "stale_seconds": float(get("PANSY_STALE_SECONDS", "60")),
        "check_seconds": float(get("PANSY_CHECK_SECONDS", "10")),
        "cooldown_seconds": float(get("PANSY_RESTART_COOLDOWN_SECONDS", "120")),
        "stop_timeout_seconds": float(get("PANSY_STOP_TIMEOUT_SECONDS", "45")),
        "phase_root": Path(get("PANSY_PHASE_ROOT", "/media/analysis/metadata/phase")),
        "phase_max_jump_deg": float(get("PANSY_PHASE_MAX_JUMP_DEG", "30")),
        "phase_lookback_days": float(get("PANSY_PHASE_LOOKBACK_DAYS", "30")),
        "phase_min_baseline": int(get("PANSY_PHASE_MIN_BASELINE", "10")),
        "phase_check_timeout_seconds": float(get("PANSY_PHASE_CHECK_TIMEOUT_SECONDS", "1200")),
        "phase_max_restarts": int(get("PANSY_PHASE_MAX_RESTARTS", "3")),
    }


def log(msg):
    now = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    print(f"{now} {msg}", flush=True)


def run(cmd, check=False):
    log("+ " + " ".join(cmd))
    return subprocess.run(cmd, check=check, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


def pids_for(name):
    out = subprocess.run(["pgrep", "-x", name], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if out.returncode != 0:
        return []
    return [int(pid) for pid in out.stdout.split()]


def wait_no_pansy_uhd(timeout):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if not pids_for("pansy_uhd_rx"):
            return True
        time.sleep(1)
    return not pids_for("pansy_uhd_rx")


def clean_restart(service, stop_timeout, reason):
    log(f"Restart requested: {reason}")
    run(["systemctl", "--user", "stop", service])
    if not wait_no_pansy_uhd(stop_timeout):
        log("pansy_uhd_rx still alive after systemd stop; sending SIGINT")
        for pid in pids_for("pansy_uhd_rx"):
            try:
                os.kill(pid, signal.SIGINT)
            except ProcessLookupError:
                pass
        if not wait_no_pansy_uhd(15):
            log("pansy_uhd_rx still alive after SIGINT; sending SIGTERM")
            for pid in pids_for("pansy_uhd_rx"):
                try:
                    os.kill(pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            wait_no_pansy_uhd(10)
    result = run(["systemctl", "--user", "start", service])
    if result.returncode != 0:
        log(result.stdout.strip())
        raise RuntimeError(f"failed to start {service}")


def scan_h5_mtime(root, max_recent_dirs=6):
    if not root.exists():
        return None, None

    try:
        entries = list(root.iterdir())
    except OSError:
        return None, None

    # Digital RF sample directories are timestamped, e.g.
    # 2026-05-21T10-00-00. Ignore channel metadata directories here; their
    # HDF5 files are not evidence that raw voltage samples are still arriving.
    dirs = [p for p in entries if p.is_dir() and p.name[:1].isdigit()]
    dirs.sort(key=lambda p: (p.name, p.stat().st_mtime), reverse=True)

    latest_mtime = None
    latest_path = None
    for path in (p for p in entries if p.is_file()):
        if not path.name.endswith((".h5", ".hdf5")):
            continue
        try:
            mtime = path.stat().st_mtime
        except OSError:
            continue
        if latest_mtime is None or mtime > latest_mtime:
            latest_mtime = mtime
            latest_path = path

    for base in dirs[:max_recent_dirs]:
        try:
            walker = os.walk(base)
            for dirpath, _, filenames in walker:
                for filename in filenames:
                    if not filename.endswith((".h5", ".hdf5")):
                        continue
                    path = Path(dirpath) / filename
                    try:
                        mtime = path.stat().st_mtime
                    except OSError:
                        continue
                    if latest_mtime is None or mtime > latest_mtime:
                        latest_mtime = mtime
                        latest_path = path
        except OSError:
            continue

    return latest_mtime, latest_path


def stale_channels(raw_root, channels, stale_seconds):
    now = time.time()
    stale = []
    ages = {}
    for channel in channels:
        mtime, path = scan_h5_mtime(raw_root / channel)
        if mtime is None:
            stale.append((channel, math.inf, None))
            ages[channel] = math.inf
            continue
        age = now - mtime
        ages[channel] = age
        if age > stale_seconds:
            stale.append((channel, age, path))
    return stale, ages


def sample_to_datetime(sample):
    return dt.datetime.fromtimestamp(int(sample) / 1e6, tz=dt.timezone.utc)


def read_phase_file(path):
    records = []
    try:
        with h5py.File(path, "r") as h5:
            for key in h5.keys():
                if not key.isdigit():
                    continue
                group = h5[key]
                if "xphase" not in group:
                    continue
                xphase = np.asarray(group["xphase"])
                if xphase.size == 0:
                    continue
                phase = np.full(8, np.nan, dtype=np.float64)
                n = min(8, xphase.size)
                phase[:n] = np.angle(xphase[:n])
                records.append((sample_to_datetime(key), phase))
    except OSError:
        return []
    return records


def phase_files(phase_root, lookback_days):
    if not phase_root.exists():
        return []
    cutoff_day = (dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=lookback_days + 2)).date().isoformat()
    files = []
    for day_dir in phase_root.glob("20*"):
        if not day_dir.is_dir() or day_dir.name[:10] < cutoff_day:
            continue
        files.extend(day_dir.glob("txphase@*.h5"))
    return sorted(files)


def circular_median(phases):
    return np.angle(np.nanmedian(np.exp(1j * phases), axis=0))


def angular_diff_deg(a, b):
    return np.rad2deg(np.angle(np.exp(1j * (a - b))))


def latest_phase_status(phase_root, lookback_days, min_baseline):
    records = []
    for path in phase_files(phase_root, lookback_days):
        records.extend(read_phase_file(path))

    if not records:
        return None

    records.sort(key=lambda item: item[0])
    latest_time, latest_phase = records[-1]
    cutoff = latest_time - dt.timedelta(days=lookback_days)
    baseline = np.asarray([phase for t, phase in records[:-1] if t >= cutoff])
    if baseline.shape[0] < min_baseline:
        return {
            "latest_time": latest_time,
            "latest_phase": latest_phase,
            "ok": None,
            "reason": f"only {baseline.shape[0]} baseline phase samples",
        }

    median = circular_median(baseline)
    diff = angular_diff_deg(latest_phase, median)
    max_abs = float(np.nanmax(np.abs(diff)))
    channel = int(np.nanargmax(np.abs(diff))) + 1
    return {
        "latest_time": latest_time,
        "latest_phase": latest_phase,
        "max_abs_deg": max_abs,
        "channel": channel,
        "ok": max_abs <= 30.0,
        "reason": f"latest phase differs from {lookback_days:g}-day median by {max_abs:.1f} deg on TX {channel}",
    }


def main():
    cfg = config()
    threshold = cfg["phase_max_jump_deg"]
    log(
        "watching "
        f"{cfg['raw_root']} channels={','.join(cfg['channels'])} "
        f"stale>{cfg['stale_seconds']:.0f}s phase_jump>{threshold:.1f}deg"
    )

    last_restart = 0.0
    pending_phase_restart_time = None
    phase_restart_count = 0

    while True:
        try:
            stale, ages = stale_channels(cfg["raw_root"], cfg["channels"], cfg["stale_seconds"])
            if stale and time.time() - last_restart >= cfg["cooldown_seconds"]:
                reason = ", ".join(
                    f"{ch} age={'missing' if math.isinf(age) else f'{age:.1f}s'}"
                    for ch, age, _ in stale
                )
                clean_restart(cfg["service"], cfg["stop_timeout_seconds"], f"stale raw-voltage data: {reason}")
                last_restart = time.time()
                pending_phase_restart_time = dt.datetime.now(dt.timezone.utc)
                phase_restart_count = 0
            elif stale:
                log("stale data detected during restart cooldown: " + ", ".join(ch for ch, _, _ in stale))

            if pending_phase_restart_time is not None:
                status = latest_phase_status(
                    cfg["phase_root"],
                    cfg["phase_lookback_days"],
                    cfg["phase_min_baseline"],
                )
                now = dt.datetime.now(dt.timezone.utc)
                if status is None:
                    if (now - pending_phase_restart_time).total_seconds() > cfg["phase_check_timeout_seconds"]:
                        log("phase check timed out: no phase metadata found")
                        pending_phase_restart_time = None
                elif status["latest_time"] <= pending_phase_restart_time:
                    if (now - pending_phase_restart_time).total_seconds() > cfg["phase_check_timeout_seconds"]:
                        log("phase check timed out: no new phase metadata after restart")
                        pending_phase_restart_time = None
                elif status["ok"] is None:
                    log("phase check skipped: " + status["reason"])
                    pending_phase_restart_time = None
                else:
                    max_abs = float(status["max_abs_deg"])
                    ok = max_abs <= threshold
                    if ok:
                        log("phase check ok: " + status["reason"])
                        pending_phase_restart_time = None
                        phase_restart_count = 0
                    elif phase_restart_count < cfg["phase_max_restarts"]:
                        phase_restart_count += 1
                        clean_restart(
                            cfg["service"],
                            cfg["stop_timeout_seconds"],
                            f"bad phasecal after restart: {status['reason']}",
                        )
                        last_restart = time.time()
                        pending_phase_restart_time = dt.datetime.now(dt.timezone.utc)
                    else:
                        log("phase check failed but restart limit reached: " + status["reason"])
                        pending_phase_restart_time = None

        except Exception as exc:
            log(f"watchdog error: {exc}")

        time.sleep(cfg["check_seconds"])


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(0)
