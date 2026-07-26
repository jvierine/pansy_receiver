#!/usr/bin/env python3
"""Build an HDF5 index of selected orbit solutions beyond the TX-beam angle limit."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import numpy as np


ROW_DTYPE = np.dtype(
    [
        ("sample_idx", "<i8"),
        ("tx_beam_weighted_mean_deg", "<f4"),
        ("initial_detection_height_km", "<f4"),
        ("selected_hypothesis", "S8"),
    ]
)


def iter_orbit_files(orbit_dir: Path):
    yield from sorted(orbit_dir.glob("20*/*.h5"))


def selected_rows(path: Path):
    try:
        with h5py.File(path, "r") as handle:
            if "events" not in handle or "aliases" not in handle:
                return []
            events = handle["events"][:]
            aliases = handle["aliases"][:]
    except OSError:
        return []

    selected = {
        int(event["sample_idx"]): (
            event["selected_hypothesis"],
            float(event["initial_detection_height_km"]),
        )
        for event in events
    }
    rows = []
    for alias in aliases:
        sample_idx = int(alias["sample_idx"])
        event = selected.get(sample_idx)
        if event is None or alias["hypothesis_label"] != event[0]:
            continue
        rows.append(
            (
                sample_idx,
                float(alias["tx_beam_weighted_mean_deg"]),
                event[1],
                event[0],
            )
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--orbit-dir",
        type=Path,
        default=Path("/mnt/data/juha/pansy/metadata/orbit"),
    )
    parser.add_argument("--output-h5", type=Path, required=True)
    parser.add_argument("--max-angle-deg", type=float, default=15.0)
    args = parser.parse_args()

    selected_count = 0
    affected = []
    files = list(iter_orbit_files(args.orbit_dir))
    for file_number, path in enumerate(files, start=1):
        rows = selected_rows(path)
        selected_count += len(rows)
        affected.extend(
            row
            for row in rows
            if np.isfinite(row[1]) and row[1] > float(args.max_angle_deg)
        )
        if file_number % 1000 == 0:
            print(
                f"progress files={file_number}/{len(files)} "
                f"selected={selected_count} affected={len(affected)}",
                flush=True,
            )

    table = np.asarray(affected, dtype=ROW_DTYPE)
    if len(table):
        _, unique_index = np.unique(table["sample_idx"], return_index=True)
        table = table[np.sort(unique_index)]
        table.sort(order="sample_idx")

    args.output_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output_h5, "w") as handle:
        handle.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        handle.attrs["source_program"] = Path(__file__).name
        handle.attrs["orbit_dir"] = str(args.orbit_dir)
        handle.attrs["max_tx_beam_weighted_mean_angle_deg"] = float(
            args.max_angle_deg
        )
        handle.attrs["orbit_files_read"] = int(len(files))
        handle.attrs["selected_solution_count"] = int(selected_count)
        handle.attrs["affected_event_count"] = int(len(table))
        handle.create_dataset(
            "sample_idx",
            data=table["sample_idx"],
            compression="gzip",
            shuffle=True,
        )
        handle.create_dataset("table", data=table, compression="gzip", shuffle=True)

    print(
        f"tx_beam_angle_reprocess_index selected={selected_count} "
        f"affected={len(table)} max_angle_deg={args.max_angle_deg:g}",
        flush=True,
    )
    print(args.output_h5)


if __name__ == "__main__":
    main()
