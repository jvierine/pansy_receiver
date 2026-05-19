#!/usr/bin/env python3
import datetime as dt
import json
import os
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import h5py

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt


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


def day_from_sample(sample):
    timestamp = int(sample) / 1e6
    return dt.datetime.fromtimestamp(timestamp, tz=dt.timezone.utc).date().isoformat()


def day_from_dir(path):
    return path.parent.name[:10]


def count_file(path):
    counts = Counter()
    with h5py.File(path, "r") as h5:
        for key in h5.keys():
            if key.isdigit():
                counts[day_from_sample(key)] += 1
    return dict(counts)


def count_file_record(path):
    try:
        return count_file(Path(path))
    except OSError:
        return {}


def load_cache(path):
    if not path.exists():
        return {"days": {}}
    try:
        cache = json.loads(path.read_text())
    except json.JSONDecodeError:
        return {"days": {}}

    # Migrate from the first cache format, which only stored aggregate counts.
    if "days" not in cache:
        cache["days"] = {
            day: {
                "count": count,
                "cached_from_legacy": True,
            }
            for day, count in cache.get("daily_counts", {}).items()
        }
    return cache


def save_cache(path, cache):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(cache, indent=2, sort_keys=True))
    tmp.replace(path)


def files_by_directory_day(detections_root):
    days = defaultdict(list)
    for path in detections_root.glob("20*/det@*.h5"):
        days[day_from_dir(path)].append(str(path))
    return {day: sorted(paths) for day, paths in days.items()}


def aggregate_daily_counts(cache):
    return {
        day: int(record.get("count", 0))
        for day, record in sorted(cache.get("days", {}).items())
    }


def latest_cached_day(cache):
    daily_counts = cache.get("daily_counts", {})
    if daily_counts:
        return max(daily_counts)

    days = cache.get("days", {})
    if days:
        return max(days)

    return None


def rescan_start_day(cache, rescan_days):
    latest = latest_cached_day(cache)
    if latest is None:
        return None
    latest_date = dt.date.fromisoformat(latest)
    return (latest_date - dt.timedelta(days=max(rescan_days - 1, 0))).isoformat()


def count_day(day, paths, workers):
    counts = Counter()
    if not paths:
        return 0

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(count_file_record, path) for path in paths]
        for future in as_completed(futures):
            counts.update(future.result())

    return int(counts.get(day, 0))


def plot_counts(daily_counts, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    dates = [dt.date.fromisoformat(day) for day in sorted(daily_counts)]
    counts = [daily_counts[day.isoformat()] for day in dates]
    total = sum(counts)

    fig, ax = plt.subplots(figsize=(11, 4.8), constrained_layout=True)
    if dates:
        ax.bar(dates, counts, width=0.85, color="#0f766e", edgecolor="#134e4a", linewidth=0.4)
        ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=4, maxticks=9))
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
        ax.set_xlim(dates[0] - dt.timedelta(days=1), dates[-1] + dt.timedelta(days=1))
    else:
        ax.text(0.5, 0.5, "No detections found", ha="center", va="center", transform=ax.transAxes)

    ax.set_title(f"PANSY meteor head echo detections: {total:,} total")
    ax.set_xlabel("Date (UTC)")
    ax.set_ylabel("Daily detections")
    ax.grid(axis="y", color="#cbd5e1", linewidth=0.7, alpha=0.8)
    ax.set_axisbelow(True)

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
    workers = int(config.get("DETECTION_HISTORY_WORKERS", env("DETECTION_HISTORY_WORKERS", "8")))
    rescan_days = int(config.get("DETECTION_HISTORY_RESCAN_DAYS", env("DETECTION_HISTORY_RESCAN_DAYS", "3")))
    web_host = config.get("WEB_HOST", "")
    web_root = config.get("WEB_ROOT", "")

    detections_root = local_root / "detections"
    cache_path = state_dir / "detection_history_cache.json"
    output_path = web_dir / "detections_daily.png"

    state_dir.mkdir(parents=True, exist_ok=True)
    cache = load_cache(cache_path)
    day_files = files_by_directory_day(detections_root)
    days = cache.setdefault("days", {})
    rescan_from = rescan_start_day(cache, rescan_days)

    for day in sorted(day_files):
        # Historical days are frozen once cached. The last few plotted days,
        # and all later days, are rescanned because rsync may still be adding
        # files there.
        if rescan_from is not None and day < rescan_from and day in days:
            continue

        count = count_day(day, day_files[day], workers)
        days[day] = {
            "count": count,
            "file_count": len(day_files[day]),
            "scanned_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        }
        cache["generated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        cache["daily_counts"] = aggregate_daily_counts(cache)
        save_cache(cache_path, cache)
        plot_counts(cache["daily_counts"], output_path)
        deploy_plot(output_path, web_host, web_root)

    # Still refresh the plot if every day was served from cache.
    cache["generated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    cache["daily_counts"] = aggregate_daily_counts(cache)
    save_cache(cache_path, cache)
    plot_counts(cache["daily_counts"], output_path)
    deploy_plot(output_path, web_host, web_root)


if __name__ == "__main__":
    main()
