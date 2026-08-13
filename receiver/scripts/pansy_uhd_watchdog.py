#!/usr/bin/env python3
"""Restart pansy_uhd_rx when raw-voltage output or phasecal goes bad."""
import datetime as dt
import json
import math
import os
import signal
import subprocess
import sys
import time
import warnings
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
        "restart_settle_seconds": float(get("PANSY_RESTART_SETTLE_SECONDS", "15")),
        "restart_state_path": Path(
            get(
                "PANSY_RESTART_STATE_PATH",
                Path.home() / ".local/state/pansy-receiver/receiver_restart.json",
            )
        ),
        "phase_root": Path(get("PANSY_PHASE_ROOT", "/media/analysis/metadata/phase")),
        "tx_metadata_root": Path(get("PANSY_TX_METADATA_ROOT", "/media/analysis/metadata/tx")),
        "phase_max_jump_deg": float(get("PANSY_PHASE_MAX_JUMP_DEG", "30")),
        "phase_lookback_days": float(get("PANSY_PHASE_LOOKBACK_DAYS", "30")),
        "phase_min_baseline": int(get("PANSY_PHASE_MIN_BASELINE", "10")),
        "phase_check_timeout_seconds": float(get("PANSY_PHASE_CHECK_TIMEOUT_SECONDS", "1200")),
        "phase_max_restarts": int(get("PANSY_PHASE_MAX_RESTARTS", "3")),
        "phase_check_seconds": float(get("PANSY_PHASE_CHECK_SECONDS", "60")),
        "phase_post_restart_wait_seconds": float(
            get("PANSY_PHASE_POST_RESTART_WAIT_SECONDS", "90")
        ),
        "phase_raw_search_samples": int(get("PANSY_PHASE_RAW_SEARCH_SAMPLES", "500000")),
        "phase_raw_pulse_samples": int(get("PANSY_PHASE_RAW_PULSE_SAMPLES", "100")),
        "phase_tx_lookback_seconds": float(get("PANSY_PHASE_TX_LOOKBACK_SECONDS", "600")),
        "phase_min_valid_channels": int(get("PANSY_PHASE_MIN_VALID_CHANNELS", "7")),
        "phase_min_mean_amplitude": float(get("PANSY_PHASE_MIN_MEAN_AMPLITUDE", "1e7")),
        "phase_reference_deg": get("PANSY_PHASE_REFERENCE_DEG", ""),
    }


def log(msg):
    now = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    print(f"{now} {msg}", flush=True)


def run(cmd, check=False):
    log("+ " + " ".join(cmd))
    return subprocess.run(cmd, check=check, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


def write_restart_state(path, timestamp, reason):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(
            {
                "receiver_started_utc": timestamp.isoformat(),
                "reason": reason,
            },
            indent=2,
        )
        + "\n"
    )
    temporary.replace(path)


def read_restart_state(path):
    try:
        value = json.loads(path.read_text())["receiver_started_utc"]
        timestamp = dt.datetime.fromisoformat(value)
        if timestamp.tzinfo is None:
            timestamp = timestamp.replace(tzinfo=dt.timezone.utc)
        return timestamp.astimezone(dt.timezone.utc)
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None


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


def clean_restart(service, stop_timeout, settle_seconds, state_path, reason):
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
    if pids_for("pansy_uhd_rx"):
        raise RuntimeError("refusing to start a second receiver while pansy_uhd_rx is still alive")
    if settle_seconds > 0:
        log(f"all receiver processes stopped; waiting {settle_seconds:g}s for UHD devices to settle")
        time.sleep(settle_seconds)
    result = run(["systemctl", "--user", "start", service])
    if result.returncode != 0:
        log(result.stdout.strip())
        raise RuntimeError(f"failed to start {service}")
    started_at = dt.datetime.now(dt.timezone.utc)
    write_restart_state(state_path, started_at, reason)
    log(f"all receivers started; restart recorded at {started_at.isoformat()}")
    return started_at


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


def historical_phase_reference(phase_root, lookback_days, min_baseline):
    records = []
    for path in phase_files(phase_root, lookback_days):
        records.extend(read_phase_file(path))
    if not records:
        return None, 0
    records.sort(key=lambda item: item[0])
    cutoff = records[-1][0] - dt.timedelta(days=lookback_days)
    baseline = np.asarray([phase for timestamp, phase in records if timestamp >= cutoff])
    if baseline.shape[0] < min_baseline:
        return None, int(baseline.shape[0])
    return circular_median(baseline), int(baseline.shape[0])


def classify_phase_vector(phase, amplitude, reference, threshold_deg, min_valid_channels):
    phase = np.asarray(phase, dtype=np.float64)
    amplitude = np.asarray(amplitude, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    n = min(phase.size, amplitude.size, reference.size)
    valid = (
        np.isfinite(phase[:n])
        & np.isfinite(amplitude[:n])
        & (amplitude[:n] > 0.0)
        & np.isfinite(reference[:n])
    )
    valid_count = int(np.count_nonzero(valid))
    if valid_count < min_valid_channels:
        return {
            "ok": None,
            "valid_channel_count": valid_count,
            "reason": f"only {valid_count} valid raw TX phase channels",
        }
    diff = angular_diff_deg(phase[:n], reference[:n])
    max_abs = float(np.nanmax(np.where(valid, np.abs(diff), np.nan)))
    channel = int(np.nanargmax(np.where(valid, np.abs(diff), np.nan)))
    return {
        "ok": max_abs <= threshold_deg,
        "valid_channel_count": valid_count,
        "max_abs_deg": max_abs,
        "channel": channel,
        "phase_deg": np.rad2deg(phase[:n]),
        "diff_deg": diff,
        "reason": f"raw TX phase differs from reference by {max_abs:.1f} deg on {channel}",
    }


def metadata_mode_id(value):
    values = np.asarray(value).reshape(-1)
    if values.size == 0:
        return None
    try:
        return int(values[0])
    except (TypeError, ValueError):
        return None


def latest_mesosphere_phase_status(
    raw_root,
    tx_metadata_root,
    channels,
    reference,
    threshold_deg,
    tx_lookback_seconds,
    pulse_samples,
    min_valid_channels,
    min_mean_amplitude,
    after_time=None,
):
    import digital_rf as drf

    reader = drf.DigitalRFReader(str(raw_root))
    bounds = [reader.get_bounds(channel) for channel in channels]
    start_bound = max(bound[0] for bound in bounds)
    raw_end = min(bound[1] for bound in bounds) - pulse_samples - 1
    if raw_end <= start_bound:
        return {"ok": None, "reason": "insufficient common raw-voltage samples"}

    try:
        metadata = drf.DigitalMetadataReader(str(tx_metadata_root))
        tx_start, tx_end = metadata.get_bounds()
    except Exception as exc:
        return {"ok": None, "reason": f"TX metadata unavailable ({exc})"}

    search_start = max(tx_start, tx_end - int(tx_lookback_seconds * 1e6))
    if after_time is not None:
        search_start = max(search_start, int(after_time.timestamp() * 1e6) + 1)
    if tx_end < search_start:
        return {"ok": None, "reason": "no post-restart mesosphere TX metadata yet"}

    try:
        records = metadata.read(search_start, tx_end, "id")
    except Exception as exc:
        return {"ok": None, "reason": f"could not read TX metadata ({exc})"}
    candidates = sorted(
        int(sample)
        for sample, mode_id in records.items()
        if metadata_mode_id(mode_id) == 1
        and start_bound <= int(sample) <= raw_end
    )
    if not candidates:
        qualifier = " post-restart" if after_time is not None else ""
        return {"ok": None, "reason": f"no{qualifier} mesosphere-mode TX pulse available"}
    pulse_start = candidates[-1]

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="The read_vector_c81d method is deprecated.*")
        reference_voltage = reader.read_vector_c81d(pulse_start, pulse_samples, channels[0])
        xphase = np.full(len(channels), np.nan + 1j * np.nan, dtype=np.complex128)
        for index, channel in enumerate(channels):
            voltage = reader.read_vector_c81d(pulse_start, pulse_samples, channel)
            xphase[index] = np.mean(reference_voltage * np.conj(voltage))

    amplitude = np.abs(xphase)
    mean_amplitude = float(np.nanmean(amplitude))
    if not np.isfinite(mean_amplitude) or mean_amplitude < min_mean_amplitude:
        return {
            "ok": None,
            "sample_idx": pulse_start,
            "latest_time": sample_to_datetime(pulse_start),
            "mean_amplitude": mean_amplitude,
            "reason": f"raw TX pulse amplitude too small ({mean_amplitude:.3g})",
        }

    status = classify_phase_vector(
        np.angle(xphase),
        amplitude,
        reference,
        threshold_deg,
        min_valid_channels,
    )
    status.update(
        {
            "sample_idx": pulse_start,
            "latest_time": sample_to_datetime(pulse_start),
            "mean_amplitude": mean_amplitude,
        }
    )
    return status


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
    last_receiver_restart = read_restart_state(cfg["restart_state_path"])
    phase_restart_count = 0
    last_phase_check = 0.0
    last_phase_sample = None
    configured_reference = np.fromstring(cfg["phase_reference_deg"], sep=",")
    if configured_reference.size == len(cfg["channels"]):
        phase_reference = np.deg2rad(configured_reference)
        reference_count = 0
        reference_description = "configured known-good raw TX phase"
    else:
        phase_reference, reference_count = historical_phase_reference(
            cfg["phase_root"],
            cfg["phase_lookback_days"],
            cfg["phase_min_baseline"],
        )
        reference_description = f"{reference_count} historical samples"
    if phase_reference is None:
        log(f"raw phase watchdog disabled: only {reference_count} historical phase samples")
    else:
        log(
            f"raw phase reference from {reference_description}: "
            + " ".join(f"{value:.1f}" for value in np.rad2deg(phase_reference))
            + " deg"
        )

    while True:
        try:
            stale, ages = stale_channels(cfg["raw_root"], cfg["channels"], cfg["stale_seconds"])
            if stale and time.time() - last_restart >= cfg["cooldown_seconds"]:
                reason = ", ".join(
                    f"{ch} age={'missing' if math.isinf(age) else f'{age:.1f}s'}"
                    for ch, age, _ in stale
                )
                last_receiver_restart = clean_restart(
                    cfg["service"],
                    cfg["stop_timeout_seconds"],
                    cfg["restart_settle_seconds"],
                    cfg["restart_state_path"],
                    f"stale raw-voltage data: {reason}",
                )
                last_restart = time.time()
            elif stale:
                log("stale data detected during restart cooldown: " + ", ".join(ch for ch, _, _ in stale))

            if phase_reference is not None and time.time() - last_phase_check >= cfg["phase_check_seconds"]:
                last_phase_check = time.time()
                now = dt.datetime.now(dt.timezone.utc)
                if last_receiver_restart is not None:
                    wait_until = last_receiver_restart + dt.timedelta(
                        seconds=cfg["phase_post_restart_wait_seconds"]
                    )
                    if now < wait_until:
                        remaining = (wait_until - now).total_seconds()
                        log(f"raw phase check waiting {remaining:.0f}s for post-restart TX pulses")
                        time.sleep(cfg["check_seconds"])
                        continue
                status = latest_mesosphere_phase_status(
                    cfg["raw_root"],
                    cfg["tx_metadata_root"],
                    cfg["channels"],
                    phase_reference,
                    threshold,
                    cfg["phase_tx_lookback_seconds"],
                    cfg["phase_raw_pulse_samples"],
                    cfg["phase_min_valid_channels"],
                    cfg["phase_min_mean_amplitude"],
                    after_time=last_receiver_restart,
                )
                sample_idx = status.get("sample_idx")
                if (
                    last_receiver_restart is not None
                    and status.get("latest_time") is not None
                    and status["latest_time"] <= last_receiver_restart
                ):
                    log("raw phase check skipped: newest TX pulse predates receiver restart")
                elif sample_idx is not None and sample_idx == last_phase_sample:
                    log("raw phase check skipped: no new raw TX pulse")
                elif status["ok"] is None:
                    log("raw phase check skipped: " + status["reason"])
                elif status["ok"]:
                    log("raw phase check ok: " + status["reason"])
                    phase_restart_count = 0
                elif time.time() - last_restart < cfg["cooldown_seconds"]:
                    log("bad raw TX phase detected during restart cooldown: " + status["reason"])
                elif phase_restart_count < cfg["phase_max_restarts"]:
                    phase_restart_count += 1
                    last_receiver_restart = clean_restart(
                        cfg["service"],
                        cfg["stop_timeout_seconds"],
                        cfg["restart_settle_seconds"],
                        cfg["restart_state_path"],
                        f"bad raw TX phase ({phase_restart_count}/{cfg['phase_max_restarts']}): {status['reason']}",
                    )
                    last_restart = time.time()
                else:
                    log("raw phase check failed but restart limit reached: " + status["reason"])
                last_phase_sample = sample_idx

        except Exception as exc:
            log(f"watchdog error: {exc}")

        time.sleep(cfg["check_seconds"])


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(0)
