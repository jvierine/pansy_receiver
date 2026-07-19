#!/usr/bin/env python3
"""Wait for finite adaptive mass profiles and render a reproducible example batch."""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import os
import subprocess
import time
from pathlib import Path

import h5py
import numpy as np


def finite_adaptive_sample(path: Path) -> int | None:
    try:
        with h5py.File(path, "r") as handle:
            if handle.attrs.get("profile_strategy", "") != "adaptive_log_radius":
                return None
            lower_bounded = bool(handle["result/ci95_lower_bounded"][()])
            upper_bounded = bool(handle["result/ci95_upper_bounded"][()])
            lower = float(handle["result/ci95_lower_radius_um"][()])
            upper = float(handle["result/ci95_upper_radius_um"][()])
            if not (lower_bounded and upper_bounded and np.isfinite(lower) and np.isfinite(upper)):
                return None
            return int(handle.attrs["sample_idx"])
    except (KeyError, OSError, ValueError):
        return None


def completed_samples(profile_dir: Path) -> list[int]:
    samples = []
    for path in sorted(profile_dir.glob("mass_profile_phase_aware_*.h5")):
        sample = finite_adaptive_sample(path)
        if sample is not None:
            samples.append(sample)
    return samples


def write_manifest(path: Path, samples: list[int], profile_dir: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle["sample_idx"] = np.asarray(samples, dtype=np.int64)
        handle.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        handle.attrs["selection"] = "first completed adaptive profiles with finite lower and upper 95% profile bounds"
        handle.attrs["profile_dir"] = str(profile_dir)
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--count", type=int, default=50)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--poll-seconds", type=float, default=60.0)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = args.output_dir / "finite_adaptive_examples.h5"

    while True:
        samples = completed_samples(args.profile_dir)
        print(
            dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            f"finite_profiles={len(samples)}/{args.count}",
            flush=True,
        )
        if len(samples) >= args.count:
            samples = samples[: args.count]
            break
        time.sleep(args.poll_seconds)
    write_manifest(manifest, samples, args.profile_dir)

    def render(worker_index: int) -> None:
        command = [
            os.environ.get("PANSY_PROFILE_PYTHON", "python3"),
            str(args.repo / "render_pansy_true_r0_batch.py"),
            "--sample-index-h5",
            str(manifest),
            "--limit",
            str(args.count),
            "--base",
            str(args.base),
            "--output-dir",
            str(args.output_dir),
            "--profile-input-dir",
            str(args.profile_dir),
            "--worker-index",
            str(worker_index),
            "--worker-count",
            str(args.workers),
        ]
        log_path = args.output_dir / f"render_worker_{worker_index:03d}.log"
        with log_path.open("w") as log:
            subprocess.run(command, cwd=args.repo, stdout=log, stderr=subprocess.STDOUT, check=True)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(render, worker) for worker in range(args.workers)]
        for future in concurrent.futures.as_completed(futures):
            future.result()
    (args.output_dir / "COMPLETE").write_text(
        dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds") + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
