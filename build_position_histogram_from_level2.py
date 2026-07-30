#!/usr/bin/env python3
"""Aggregate released Level 2 positions on the interferometry direction grid."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

from plot_meteor_position_histograms import (
    TX_BEAM_COUNT,
    interferometry_grid_offsets_deg,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level2-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--start-date", default="2025-03-15")
    parser.add_argument("--grid-n", type=int, default=501)
    parser.add_argument("--chunk-rows", type=int, default=1_000_000)
    return parser.parse_args()


def start_sample(date_text: str) -> int:
    value = datetime.strptime(date_text, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    return int(value.timestamp() * 1_000_000)


def add_chunk(measurements, slc, minimum_sample, all_counts, beam_counts):
    sample = np.asarray(measurements["time_sample_idx"][slc], dtype=np.int64)
    keep = sample >= minimum_sample
    keep &= np.asarray(measurements["selection_keep"][slc], dtype=bool)
    if not np.any(keep):
        return 0

    east = np.asarray(measurements["east_km"][slc], dtype=np.float64)[keep]
    north = np.asarray(measurements["north_km"][slc], dtype=np.float64)[keep]
    up = np.asarray(measurements["up_km"][slc], dtype=np.float64)[keep]
    beam = np.asarray(measurements["beam_id"][slc], dtype=np.int64)[keep]
    radius = np.sqrt(east**2 + north**2 + up**2)
    valid = (
        np.isfinite(radius)
        & (radius > 0.0)
        & np.isfinite(east)
        & np.isfinite(north)
        & np.isfinite(up)
        & (beam >= 0)
        & (beam < TX_BEAM_COUNT)
    )
    if not np.any(valid):
        return 0

    # The interferometry convention points along wave propagation, from the
    # meteor toward the array, hence the minus signs relative to ENU position.
    u = -east[valid] / radius[valid]
    v = -north[valid] / radius[valid]
    beam = beam[valid]
    grid_n = all_counts.shape[0]
    col = np.rint((u + 1.0) * 0.5 * (grid_n - 1)).astype(np.int64)
    row = np.rint((v + 1.0) * 0.5 * (grid_n - 1)).astype(np.int64)
    valid = (
        (u * u + v * v <= 1.0)
        & (row >= 0)
        & (row < grid_n)
        & (col >= 0)
        & (col < grid_n)
    )
    row = row[valid]
    col = col[valid]
    beam = beam[valid]
    np.add.at(all_counts, (row, col), 1)
    np.add.at(beam_counts, (beam, row, col), 1)
    return len(row)


def main():
    args = parse_args()
    paths = sorted(args.level2_dir.glob("pansy_level2_*.h5"))
    if not paths:
        raise FileNotFoundError(f"no Level 2 monthly files in {args.level2_dir}")

    minimum_sample = start_sample(args.start_date)
    east_grid, north_grid = interferometry_grid_offsets_deg(args.grid_n)[2:4]
    monthly_counts = []
    monthly_beam_counts = []
    month_labels = []
    total = 0

    for path in paths:
        month = path.stem.removeprefix("pansy_level2_")
        all_counts = np.zeros((args.grid_n, args.grid_n), dtype=np.int32)
        beam_counts = np.zeros(
            (TX_BEAM_COUNT, args.grid_n, args.grid_n), dtype=np.int32
        )
        with h5py.File(path, "r") as handle:
            measurements = handle["measurements"]
            n_rows = len(measurements["time_sample_idx"])
            month_total = 0
            for first in range(0, n_rows, args.chunk_rows):
                last = min(first + args.chunk_rows, n_rows)
                month_total += add_chunk(
                    measurements,
                    slice(first, last),
                    minimum_sample,
                    all_counts,
                    beam_counts,
                )
        if month_total == 0:
            continue
        label = max(f"{month}-01", args.start_date) if month == args.start_date[:7] else f"{month}-01"
        month_labels.append(label)
        monthly_counts.append(all_counts)
        monthly_beam_counts.append(beam_counts)
        total += month_total
        print(f"{month}: {month_total:,} retained pulse positions", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["schema"] = "pansy_level2_position_histogram_stack_v1"
        handle.attrs["source_level2_dir"] = str(args.level2_dir)
        handle.attrs["start_date"] = args.start_date
        handle.attrs["grid_n"] = args.grid_n
        handle.attrs["total_retained_measurements"] = total
        handle.create_dataset("utc_day", data=np.asarray(month_labels, dtype="S10"))
        handle.create_dataset(
            "east_grid_deg", data=east_grid.astype(np.float32), compression="gzip"
        )
        handle.create_dataset(
            "north_grid_deg", data=north_grid.astype(np.float32), compression="gzip"
        )
        handle.create_dataset(
            "all_counts",
            data=np.asarray(monthly_counts, dtype=np.int32),
            compression="gzip",
            shuffle=True,
        )
        handle.create_dataset(
            "beam_counts",
            data=np.asarray(monthly_beam_counts, dtype=np.int32),
            compression="gzip",
            shuffle=True,
        )
    temporary.replace(args.output)
    print(f"{total:,} retained pulse positions")
    print(args.output)


if __name__ == "__main__":
    main()
