#!/usr/bin/env python3
"""Build and publish current PANSY radiant monitor products."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import subprocess
import time
from pathlib import Path

import h5py
import numpy as np

from collect_daily_radiants_from_orbit_metadata import collect, daily_counts
from plot_fitted_radiant_distribution import arrays_to_rows, plot_radiants_with_options
from plot_fitted_radiant_healpix import healpix_histogram, plot_healpix_with_showers, read_radiants, write_h5
from shower_radiant_overlay import active_showers


DEFAULT_ORBIT_METADATA_DIR = Path("/mnt/data/juha/pansy/metadata/orbit")
DEFAULT_OUTPUT_DIR = Path("/mnt/data/juha/pansy/events/current_distribution_plots")
DEFAULT_RSYNC_DEST = "j@4.235.86.214:/var/www/html/pansy/"


def utc_day_from_epoch(epoch_unix: float) -> str:
    return dt.datetime.fromtimestamp(float(epoch_unix), tz=dt.timezone.utc).strftime("%Y-%m-%d")


def write_radiant_h5(path: Path, rows: np.ndarray, files_read: int, files_skipped: int) -> np.ndarray:
    path.parent.mkdir(parents=True, exist_ok=True)
    counts = daily_counts(rows)
    with h5py.File(path, "w") as h:
        h.attrs["script"] = Path(__file__).name
        h.attrs["source"] = "compact hourly orbit metadata events tables"
        h.attrs["radiant_definition"] = "DASST zenith-attraction corrected orbit radiant"
        h.attrs["longitude_convention"] = "lambda_radiant - lambda_sun wrapped to [0, 360) deg; signed copy in [-180, 180)"
        h.attrs["files_read"] = int(files_read)
        h.attrs["files_skipped"] = int(files_skipped)
        h.create_dataset("radiants", data=rows)
        h.create_dataset("daily_counts", data=counts)
    return counts


def write_status_json(
    path: Path,
    rows: np.ndarray,
    counts: np.ndarray,
    files_read: int,
    files_skipped: int,
    scatter_frames: int,
    histogram_frames: int,
) -> dict:
    latest_day = counts["utc_day"][-1].decode() if len(counts) else None
    latest_day_count = int(counts["count"][-1]) if len(counts) else 0
    first_day = counts["utc_day"][0].decode() if len(counts) else None
    latest_epoch = float(np.nanmax(rows["epoch_unix"])) if len(rows) else None
    status = {
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "first_analyzed_day": first_day,
        "latest_analyzed_day": latest_day,
        "latest_radiant_utc": utc_day_from_epoch(latest_epoch) if latest_epoch is not None else None,
        "latest_day_radiants": latest_day_count,
        "total_radiants": int(len(rows)),
        "daily_bins": int(len(counts)),
        "orbit_files_read": int(files_read),
        "orbit_files_skipped": int(files_skipped),
        "scatter_frames": int(scatter_frames),
        "histogram_frames": int(histogram_frames),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
    return status


def run_animation(input_h5: Path, output_dir: Path, bins: int, window_days: int) -> tuple[int, int]:
    from animate_radiant_distribution import make_animation, read_radiants as read_animation_radiants

    arr = read_animation_radiants(input_h5)
    arr = arr[np.argsort(arr["epoch_unix"])]
    frame_dir = output_dir / "radiant_animation_frames"
    scatter_frames = make_animation(
        arr,
        output_dir / "current_radiant_scatter_daily.gif",
        frame_dir / "scatter",
        "scatter",
        None,
        5.0,
        False,
        window_days,
        bins,
    )
    histogram_frames = make_animation(
        arr,
        output_dir / "current_radiant_histogram_daily.gif",
        frame_dir / "histogram",
        "histogram",
        None,
        5.0,
        False,
        window_days,
        bins,
    )
    return scatter_frames, histogram_frames


def build_products(args) -> dict:
    output_dir = args.output_dir
    rows, files_read, files_skipped = collect(args.orbit_metadata_dir)
    radiant_h5 = output_dir / "current_fitted_radiant_distribution.h5"
    counts = write_radiant_h5(radiant_h5, rows, files_read, files_skipped)
    write_radiant_h5(output_dir / "current_daily_fitted_radiants.h5", rows, files_read, files_skipped)

    plot_radiants_with_options(
        arrays_to_rows(rows),
        output_dir / "current_fitted_radiant_distribution.png",
        show_shower_overlays=True,
        shower_catalog=None,
        shower_solar_longitude=None,
        shower_date=None,
        shower_peak_tolerance_deg=5.0,
    )

    lon, lat, speed = read_radiants(radiant_h5)
    healpix_count, mean_speed = healpix_histogram(lon, lat, speed, args.nside, args.min_count_for_mean_speed)
    write_h5(
        output_dir / "current_fitted_radiant_distribution_healpix.h5",
        radiant_h5,
        args.nside,
        healpix_count,
        mean_speed,
        len(lon),
    )
    showers, query = active_showers(rows, peak_tolerance_deg=5.0)
    plot_healpix_with_showers(
        healpix_count,
        output_dir / "current_fitted_radiant_distribution_healpix.png",
        len(lon),
        showers,
        query,
        rows=rows,
    )

    scatter_frames, histogram_frames = run_animation(radiant_h5, output_dir, args.bins, args.window_days)
    status = write_status_json(
        output_dir / "radiant_monitor.json",
        rows,
        counts,
        files_read,
        files_skipped,
        scatter_frames,
        histogram_frames,
    )
    return status


def rsync_products(output_dir: Path, dest: str, bwlimit: str) -> None:
    files = [
        "current_fitted_radiant_distribution.png",
        "current_fitted_radiant_distribution_healpix.png",
        "current_radiant_scatter_daily.gif",
        "current_radiant_histogram_daily.gif",
        "radiant_monitor.json",
    ]
    cmd = ["rsync", "-avz", "--bwlimit", str(bwlimit), *[str(output_dir / name) for name in files], dest]
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-metadata-dir", type=Path, default=DEFAULT_ORBIT_METADATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--rsync-dest", default=DEFAULT_RSYNC_DEST)
    parser.add_argument("--rsync-bwlimit", default="1000")
    parser.add_argument("--nside", type=int, default=64)
    parser.add_argument("--bins", type=int, default=144)
    parser.add_argument("--window-days", type=int, default=3)
    parser.add_argument("--min-count-for-mean-speed", type=int, default=3)
    parser.add_argument("--loop", action="store_true")
    parser.add_argument("--interval-s", type=float, default=1800.0)
    parser.add_argument("--no-rsync", action="store_true")
    args = parser.parse_args()

    while True:
        status = build_products(args)
        if not args.no_rsync:
            rsync_products(args.output_dir, args.rsync_dest, args.rsync_bwlimit)
        print(
            "radiant_monitor_publish",
            status.get("latest_analyzed_day"),
            f"radiants={status.get('total_radiants')}",
            flush=True,
        )
        if not args.loop:
            break
        time.sleep(max(60.0, float(args.interval_s)))


if __name__ == "__main__":
    main()
