#!/usr/bin/env python3
"""MPI batch runner for PANSY interferometric disambiguation event plots."""

from __future__ import annotations

import argparse
import datetime as dt
import os
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np
from mpi4py import MPI


def discover_cut_sample_indices(cut_dir: Path, day: str | None = None) -> list[int]:
    sample_indices: list[int] = []
    pattern = f"{day}T*/*.h5" if day else "*/*.h5"
    for path in sorted(cut_dir.glob(pattern)):
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


def discover_cut_days(cut_dir: Path) -> list[str]:
    days = set()
    for path in cut_dir.iterdir():
        if not path.is_dir():
            continue
        name = path.name
        if len(name) >= 10 and name[4] == "-" and name[7] == "-" and "T" in name:
            days.add(name[:10])
    return sorted(days)


def load_sample_index_h5(path: Path) -> list[int]:
    with h5py.File(path, "r") as handle:
        if "sample_idx" not in handle:
            raise KeyError(f"{path} does not contain sample_idx")
        return [int(value) for value in np.asarray(handle["sample_idx"], dtype=np.int64)]


def run_one(sample_idx: int, args, rank: int) -> tuple[bool, float, Path]:
    output_dir = Path(args.output_dir)
    if args.daily_output_dirs:
        day = dt.datetime.fromtimestamp(sample_idx / 1e6, tz=dt.timezone.utc).strftime("%Y-%m-%d")
        output_dir = output_dir / day
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
    if args.tx_phase_quality_h5 is not None:
        cmd.extend(["--tx-phase-quality-h5", str(args.tx_phase_quality_h5)])
    if args.require_good_tx_phase:
        cmd.append("--require-good-tx-phase")
    if args.tx_phase_max_age_s is not None:
        cmd.extend(["--tx-phase-max-age-s", str(args.tx_phase_max_age_s)])
    if args.compute_missing_tx_phase:
        cmd.append("--compute-missing-tx-phase")
    if args.tx_phase_raw_voltage_dir is not None:
        cmd.extend(["--tx-phase-raw-voltage-dir", str(args.tx_phase_raw_voltage_dir)])
    if args.tx_phase_tx_metadata_dir is not None:
        cmd.extend(["--tx-phase-tx-metadata-dir", str(args.tx_phase_tx_metadata_dir)])
    cmd.extend(["--tx-phase-fallback-search-radius-s", str(args.tx_phase_fallback_search_radius_s)])
    cmd.extend(["--tx-phase-fallback-samples", str(args.tx_phase_fallback_samples)])
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
    parser.add_argument("--daily-output-dirs", action="store_true", help="Write event products under output-dir/YYYY-MM-DD based on sample time.")
    parser.add_argument("--grid-n", type=int, default=501)
    parser.add_argument("--coherence-threshold", type=float, default=0.80)
    parser.add_argument("--max-peaks-per-pulse", type=int, default=32)
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--sample-idx", type=int, nargs="*", help="Explicit sample indices to process instead of discovering all cuts.")
    parser.add_argument("--sample-index-h5", type=Path, default=None, help="HDF5 event index with sample_idx dataset.")
    parser.add_argument("--daily-input-batches", action="store_true", help="Discover and process cut events one UTC day at a time.")
    parser.add_argument("--day", type=str, nargs="*", help="Restrict daily input batches to these YYYY-MM-DD UTC days.")
    parser.add_argument("--start-day", type=str, default=None, help="First YYYY-MM-DD UTC day for daily input batches.")
    parser.add_argument("--end-day", type=str, default=None, help="Last YYYY-MM-DD UTC day for daily input batches.")
    parser.add_argument("--orbit-samples", type=int, default=10)
    parser.add_argument("--orbit-metadata-dir", type=Path, default=Path("data/metadata/orbit"))
    parser.add_argument("--run-dasst", action="store_true", help="Run DASST for the winning hypothesis when the local DASST module is available.")
    parser.add_argument("--recompute-cut-observables", action="store_true", help="Recompute full per-pulse range-Doppler observables instead of using cached cut detections.")
    parser.add_argument("--tx-phase-quality-h5", type=Path, default=None, help="HDF5 TX cross-phase quality table produced by tx_phase_quality.py.")
    parser.add_argument("--require-good-tx-phase", action="store_true", help="Reject events whose nearest TX cross-phase sample is missing or misaligned.")
    parser.add_argument("--tx-phase-max-age-s", type=float, default=None, help="Override maximum nearest TX phase sample age in seconds.")
    parser.add_argument("--compute-missing-tx-phase", action="store_true", help="Compute TX cross-phase from raw voltage when no nearby txphase metadata exists.")
    parser.add_argument("--tx-phase-raw-voltage-dir", type=Path, default=None, help="Digital RF raw voltage directory for missing TX phase fallback.")
    parser.add_argument("--tx-phase-tx-metadata-dir", type=Path, default=None, help="TX metadata directory for missing TX phase fallback.")
    parser.add_argument("--tx-phase-fallback-search-radius-s", type=float, default=2.0)
    parser.add_argument("--tx-phase-fallback-samples", type=int, default=120)
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    def run_sample_indices(sample_indices: list[int], label: str) -> tuple[int, int]:
        sample_indices = comm.bcast(sample_indices if rank == 0 else None, root=0)
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
                f"rank {rank:03d} {status} batch {label} sample_idx {sample_idx} "
                f"elapsed_s {elapsed:.1f} progress {n_done}/{len(my_indices)} log {log_path}",
                flush=True,
            )
        totals = comm.gather((ok_count, fail_count), root=0)
        if rank == 0:
            total_ok = sum(x[0] for x in totals)
            total_fail = sum(x[1] for x in totals)
            print(f"batch_complete {label} ok {total_ok} failed {total_fail}", flush=True)
            return total_ok, total_fail
        return 0, 0

    if rank == 0:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        args.orbit_metadata_dir.mkdir(parents=True, exist_ok=True)
        print(f"mpi_ranks {size}", flush=True)

    total_ok = 0
    total_fail = 0
    if args.daily_input_batches and not args.sample_idx and args.sample_index_h5 is None:
        if rank == 0:
            days = discover_cut_days(args.cut_dir)
            if args.day:
                wanted = set(args.day)
                days = [day for day in days if day in wanted]
            if args.start_day is not None:
                days = [day for day in days if day >= args.start_day]
            if args.end_day is not None:
                days = [day for day in days if day <= args.end_day]
            print(f"discovered_days {len(days)}", flush=True)
        else:
            days = None
        days = comm.bcast(days, root=0)
        remaining_limit = args.limit
        for day in days:
            if rank == 0:
                if remaining_limit is not None and remaining_limit <= 0:
                    sample_indices = []
                else:
                    sample_indices = discover_cut_sample_indices(args.cut_dir, day=day)
                if remaining_limit is not None and sample_indices:
                    sample_indices = sample_indices[:remaining_limit]
                    remaining_limit -= len(sample_indices)
                print(f"day_discovered {day} events {len(sample_indices)}", flush=True)
            else:
                sample_indices = None
            day_ok, day_fail = run_sample_indices(sample_indices, day)
            total_ok += day_ok
            total_fail += day_fail
    else:
        if rank == 0:
            if args.sample_idx:
                sample_indices = sorted(set(int(x) for x in args.sample_idx))
            elif args.sample_index_h5 is not None:
                sample_indices = load_sample_index_h5(args.sample_index_h5)
            else:
                sample_indices = discover_cut_sample_indices(args.cut_dir)
            if args.limit is not None and not args.sample_idx:
                sample_indices = sample_indices[: args.limit]
            print(f"discovered_events {len(sample_indices)}", flush=True)
        else:
            sample_indices = None
        total_ok, total_fail = run_sample_indices(sample_indices, "all")

    if rank == 0:
        print(f"all_batches_complete ok {total_ok} failed {total_fail}", flush=True)
        if total_fail:
            sys.exit(1)


if __name__ == "__main__":
    main()
