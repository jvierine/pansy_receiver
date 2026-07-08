#!/usr/bin/env python3
"""Save guarded cut-metadata noise measurements to an HDF5 sidecar."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

import paper_plot_noise as ppn


def save_hdf5(data: dict[str, np.ndarray | str | int], output: Path, source: Path, chunk_seconds: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        handle.attrs["day"] = str(data["day"])
        handle.attrs["source"] = str(source)
        handle.attrs["guard_samples"] = int(data["guard_samples"])
        handle.attrs["chunk_seconds"] = int(chunk_seconds)
        handle.attrs["power_estimator"] = "mean"
        handle.attrs["description"] = (
            "Per-pulse guarded cut-padding raw noise power. "
            "noise_power is mean(I^2+Q^2) over receiver channels and retained pre/post padding samples."
        )
        handle.create_dataset("times_s", data=np.asarray(data["times_s"], dtype=np.float64), compression="gzip", shuffle=True)
        handle.create_dataset("beam_id", data=np.asarray(data["beam_id"], dtype=np.int16), compression="gzip", shuffle=True)
        handle.create_dataset("noise_power", data=np.asarray(data["noise_power"], dtype=np.float64), compression="gzip", shuffle=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--day", default="2025-05-10")
    parser.add_argument("--metadata", type=Path, default=Path("/mnt/data/juha/pansy/metadata/cut"))
    parser.add_argument("--output", type=Path, default=Path("figs/cut_noise_measurements_2025-05-10_guard25_mean.h5"))
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--chunk-seconds", type=int, default=60)
    parser.add_argument("--guard-samples", type=int, default=25)
    args = parser.parse_args()

    data = ppn.read_cut_day(
        args.metadata,
        args.day,
        workers=args.workers,
        chunk_seconds=args.chunk_seconds,
        guard_samples=args.guard_samples,
    )
    save_hdf5(data, args.output, args.metadata, args.chunk_seconds)
    print(f"output {args.output}")
    print(f"samples {len(data['noise_power'])}")
    print(f"guard_samples {int(data['guard_samples'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
