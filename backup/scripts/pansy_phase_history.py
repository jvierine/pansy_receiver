#!/usr/bin/env python3
import datetime as dt
import json
import os
import subprocess
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

warnings.filterwarnings("ignore", category=DeprecationWarning)

import h5py
import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt


WARNING_THRESHOLD_DEG = 30.0


def env(name, default):
    return os.environ.get(name, default)


def load_env_file(path):
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


def day_from_dir(path):
    return path.parent.name[:10]


def sample_to_datetime64(sample):
    return np.datetime64(int(sample), "us")


def read_phase_file(path):
    times = []
    phases = []
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
                times.append(sample_to_datetime64(key))
                phases.append(np.angle(xphase[:8]))
    except OSError:
        return np.array([], dtype="datetime64[us]"), np.empty((0, 8), dtype=np.float32)

    if not times:
        return np.array([], dtype="datetime64[us]"), np.empty((0, 8), dtype=np.float32)

    out = np.full((len(phases), 8), np.nan, dtype=np.float32)
    for i, phase in enumerate(phases):
        n = min(8, phase.size)
        out[i, :n] = phase[:n]

    return np.asarray(times, dtype="datetime64[us]"), out


def read_phase_file_record(path):
    return read_phase_file(Path(path))


def files_by_directory_day(phase_root):
    days = {}
    for path in phase_root.glob("20*/txphase@*.h5"):
        days.setdefault(day_from_dir(path), []).append(str(path))
    return {day: sorted(paths) for day, paths in days.items()}


def load_cache(path):
    if not path.exists():
        return {"days": {}}
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return {"days": {}}


def save_cache(path, cache):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(cache, indent=2, sort_keys=True))
    tmp.replace(path)


def latest_cached_day(cache):
    days = cache.get("days", {})
    return max(days) if days else None


def rescan_start_day(cache, rescan_days):
    latest = latest_cached_day(cache)
    if latest is None:
        return None
    latest_date = dt.date.fromisoformat(latest)
    return (latest_date - dt.timedelta(days=max(rescan_days - 1, 0))).isoformat()


def count_day(day, paths, workers, day_cache_dir):
    times = []
    phases = []

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(read_phase_file_record, path) for path in paths]
        for future in as_completed(futures):
            t, p = future.result()
            if t.size:
                times.append(t)
                phases.append(p)

    if times:
        t_all = np.concatenate(times)
        p_all = np.concatenate(phases, axis=0)
        order = np.argsort(t_all)
        t_all = t_all[order]
        p_all = p_all[order]
    else:
        t_all = np.array([], dtype="datetime64[us]")
        p_all = np.empty((0, 8), dtype=np.float32)

    day_cache_dir.mkdir(parents=True, exist_ok=True)
    out = day_cache_dir / f"{day}.npz"
    np.savez_compressed(out, time=t_all.astype("datetime64[us]").astype("int64"), phase=p_all)
    return out, int(t_all.size)


def load_all_cached_days(day_cache_dir, days):
    times = []
    phases = []
    for day in sorted(days):
        path = day_cache_dir / f"{day}.npz"
        if not path.exists():
            continue
        try:
            data = np.load(path)
            t = data["time"].astype("datetime64[us]")
            p = data["phase"]
        except Exception:
            continue
        if t.size:
            times.append(t)
            phases.append(p)

    if not times:
        return np.array([], dtype="datetime64[us]"), np.empty((0, 8), dtype=np.float32)

    return np.concatenate(times), np.concatenate(phases, axis=0)


def circular_median(phases):
    return np.angle(np.nanmedian(np.exp(1j * phases), axis=0))


def angular_diff_deg(a, b):
    return np.rad2deg(np.angle(np.exp(1j * (a - b))))


def phase_warning(times, phases):
    if times.size == 0:
        return None

    order = np.argsort(times)
    times = times[order]
    phases = phases[order]
    latest_phase = phases[-1]
    latest_day = times[-1].astype("datetime64[D]")
    start = latest_day - np.timedelta64(30, "D")
    mask = (times.astype("datetime64[D]") >= start) & (times.astype("datetime64[D]") < latest_day)

    if np.count_nonzero(mask) < 2:
        return None

    median = circular_median(phases[mask])
    diff = angular_diff_deg(latest_phase, median)
    max_abs = float(np.nanmax(np.abs(diff)))
    if max_abs <= WARNING_THRESHOLD_DEG:
        return None

    channel = int(np.nanargmax(np.abs(diff))) + 1
    return f"WARNING: latest phasecal differs from previous 30-day median by {max_abs:.1f} deg on TX {channel}"


def plot_phase_history(times, phases, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12, 5.6), constrained_layout=True)
    warning = phase_warning(times, phases)

    if times.size:
        x = times.astype("datetime64[ms]").astype(object)
        for i in range(min(8, phases.shape[1])):
            ax.scatter(
                x,
                np.rad2deg(phases[:, i]),
                s=2,
                alpha=0.55,
                linewidths=0,
                rasterized=True,
                label=f"TX {i + 1}",
            )
        ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=4, maxticks=9))
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
        ax.set_ylim(-190, 190)
    else:
        ax.text(0.5, 0.5, "No phase calibration data found", ha="center", va="center", transform=ax.transAxes)

    title = f"PANSY transmitter phase calibration: {times.size:,} samples"
    if warning:
        title = f"{title}\n{warning}"
    ax.set_title(title)
    ax.set_xlabel("Date (UTC)")
    ax.set_ylabel("Transmitter phase (deg)")
    ax.grid(True, color="#cbd5e1", linewidth=0.7, alpha=0.8)
    ax.legend(loc="upper right", ncol=4, fontsize=8, frameon=False)

    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def deploy_plot(output_path, web_host, web_root):
    if not output_path.exists() or not web_host or not web_root:
        return
    subprocess.run(
        ["rsync", "-az", str(output_path), f"{web_host}:{web_root}/"],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def main():
    config_path = Path(env("PANSY_BACKUP_CONFIG", str(Path.home() / ".config/pansy-backup/pansy-backup.env")))
    config = load_env_file(config_path)
    local_root = Path(config.get("LOCAL_METADATA_ROOT", "/mnt/data/juha/pansy/metadata"))
    state_dir = Path(config.get("STATE_DIR", "/mnt/data/juha/pansy/backup_state"))
    web_dir = Path(config.get("WEB_BUILD_DIR", "/mnt/data/juha/pansy/web"))
    workers = int(config.get("PHASE_HISTORY_WORKERS", env("PHASE_HISTORY_WORKERS", "8")))
    rescan_days = int(config.get("PHASE_HISTORY_RESCAN_DAYS", env("PHASE_HISTORY_RESCAN_DAYS", "3")))
    web_host = config.get("WEB_HOST", "")
    web_root = config.get("WEB_ROOT", "")

    phase_root = local_root / "phase"
    cache_dir = state_dir / "phase_history"
    day_cache_dir = cache_dir / "days"
    cache_path = cache_dir / "phase_history_cache.json"
    output_path = web_dir / "phase_history.png"

    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = load_cache(cache_path)
    days = cache.setdefault("days", {})
    day_files = files_by_directory_day(phase_root)
    rescan_from = rescan_start_day(cache, rescan_days)

    for day in sorted(day_files):
        if rescan_from is not None and day < rescan_from and day in days:
            continue

        npz_path, sample_count = count_day(day, day_files[day], workers, day_cache_dir)
        days[day] = {
            "sample_count": sample_count,
            "file_count": len(day_files[day]),
            "cache_file": str(npz_path),
            "scanned_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        }
        cache["generated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        save_cache(cache_path, cache)

        times, phases = load_all_cached_days(day_cache_dir, days)
        plot_phase_history(times, phases, output_path)
        deploy_plot(output_path, web_host, web_root)

    cache["generated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    save_cache(cache_path, cache)
    times, phases = load_all_cached_days(day_cache_dir, days)
    plot_phase_history(times, phases, output_path)
    deploy_plot(output_path, web_host, web_root)


if __name__ == "__main__":
    main()
