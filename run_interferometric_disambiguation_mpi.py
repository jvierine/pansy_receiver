#!/usr/bin/env python3
"""MPI batch runner for PANSY interferometric disambiguation event plots."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import h5py
from mpi4py import MPI


def discover_cut_sample_indices(cut_dir: Path) -> list[int]:
    sample_indices: list[int] = []
    for path in sorted(cut_dir.glob("*/*.h5")):
        try:
            with h5py.File(path, "r") as handle:
                for key in handle.keys():
                    try:
                        sample_indices.append(int(key))
                    except ValueError:
                        continue
        except OSError:
            continue
    return sorted(set(sample_indices))


def run_one(sample_idx: int, args, rank: int) -> tuple[bool, float, Path]:
    output_dir = Path(args.output_dir)
    summary = output_dir / f"pansy_interferometer_disambiguation_summary_{sample_idx}.png"
    diagnostics_h5 = output_dir / f"pansy_disambiguation_diagnostics_{sample_idx}.h5"
    log_path = output_dir / "logs" / f"{sample_idx}.rank{rank:03d}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if args.skip_existing and summary.exists() and diagnostics_h5.exists():
        return True, 0.0, log_path

    cmd = [
        sys.executable,
        str(Path(__file__).with_name("plot_interferometric_disambiguation.py")),
        "--sample-idx",
        str(sample_idx),
        "--cut-dir",
        str(args.cut_dir),
        "--output-dir",
        str(output_dir),
        "--grid-n",
        str(args.grid_n),
        "--coherence-threshold",
        str(args.coherence_threshold),
        "--max-peaks-per-pulse",
        str(args.max_peaks_per_pulse),
        "--snr-threshold",
        str(args.snr_threshold),
        "--overview-only",
        "--orbit-samples",
        str(args.orbit_samples),
        "--orbit-metadata-dir",
        str(args.orbit_metadata_dir),
    ]
    if args.run_dasst:
        cmd.append("--run-dasst")
    if args.recompute_cut_observables:
        cmd.append("--recompute-cut-observables")
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    t0 = time.time()
    with log_path.open("w") as log:
        log.write(" ".join(cmd) + "\n\n")
        log.flush()
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
    return proc.returncode == 0, time.time() - t0, log_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Batch PANSY disambiguation overview plots with MPI.")
    parser.add_argument("--cut-dir", type=Path, default=Path("data/metadata/cut"))
    parser.add_argument("--output-dir", type=Path, default=Path("test_plots"))
    parser.add_argument("--grid-n", type=int, default=501)
    parser.add_argument("--coherence-threshold", type=float, default=0.80)
    parser.add_argument("--max-peaks-per-pulse", type=int, default=32)
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--orbit-samples", type=int, default=20)
    parser.add_argument("--orbit-metadata-dir", type=Path, default=Path("data/metadata/orbit"))
    parser.add_argument("--run-dasst", action="store_true", help="Run DASST for the winning hypothesis when the local DASST module is available.")
    parser.add_argument("--recompute-cut-observables", action="store_true", help="Recompute full per-pulse range-Doppler observables instead of using cached cut detections.")
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        sample_indices = discover_cut_sample_indices(args.cut_dir)
        if args.limit is not None:
            sample_indices = sample_indices[: args.limit]
        print(f"discovered_events {len(sample_indices)}")
        print(f"mpi_ranks {size}")
    else:
        sample_indices = None

    sample_indices = comm.bcast(sample_indices, root=0)
    my_indices = sample_indices[rank::size]
    ok_count = 0
    fail_count = 0
    for n_done, sample_idx in enumerate(my_indices, start=1):
        ok, elapsed, log_path = run_one(sample_idx, args, rank)
        if ok:
            ok_count += 1
            status = "ok"
        else:
            fail_count += 1
            status = "failed"
        print(
            f"rank {rank:03d} {status} sample_idx {sample_idx} "
            f"elapsed_s {elapsed:.1f} progress {n_done}/{len(my_indices)} log {log_path}",
            flush=True,
        )

    totals = comm.gather((ok_count, fail_count), root=0)
    if rank == 0:
        total_ok = sum(x[0] for x in totals)
        total_fail = sum(x[1] for x in totals)
        print(f"batch_complete ok {total_ok} failed {total_fail}")
        if total_fail:
            sys.exit(1)


if __name__ == "__main__":
    main()
