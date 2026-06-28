#!/usr/bin/env python3
"""MPI backfill of compact orbit metadata from existing candidate-state HDF5 files."""

from __future__ import annotations

import argparse
import datetime as dt
import subprocess
import sys
import time
from pathlib import Path

from mpi4py import MPI

import orbit_metadata_table


STATE_PREFIX = "pansy_candidate_orbit_states_"
STATE_SUFFIX = ".h5"


def sample_idx_from_state_path(path: Path) -> int:
    name = path.name
    if not (name.startswith(STATE_PREFIX) and name.endswith(STATE_SUFFIX)):
        raise ValueError(f"not a candidate-state path: {path}")
    return int(name[len(STATE_PREFIX) : -len(STATE_SUFFIX)])


def day_from_sample_idx(sample_idx: int) -> str:
    return dt.datetime.fromtimestamp(sample_idx / 1e6, tz=dt.timezone.utc).strftime("%Y-%m-%d")


def discover_state_paths(events_dir: Path, start_day: str | None, end_day: str | None) -> list[Path]:
    paths: list[Path] = []
    for day_dir in sorted(events_dir.glob("2025-*")):
        if not day_dir.is_dir():
            continue
        day = day_dir.name[:10]
        if start_day is not None and day < start_day:
            continue
        if end_day is not None and day > end_day:
            continue
        paths.extend(sorted(day_dir.glob(f"{STATE_PREFIX}*{STATE_SUFFIX}")))
    return paths


def output_h5_for_state(state_path: Path) -> Path:
    sample_idx = sample_idx_from_state_path(state_path)
    return state_path.with_name(f"pansy_candidate_orbits_dasst_{sample_idx}.h5")


def run_one(state_path: Path, args, rank: int) -> tuple[str, float, int, str]:
    sample_idx = sample_idx_from_state_path(state_path)
    t0 = time.time()
    if args.skip_existing and orbit_metadata_table.has_sample(args.orbit_metadata_dir, sample_idx):
        return "skip", 0.0, sample_idx, ""
    cmd = [
        sys.executable,
        str(Path(__file__).with_name("dasst_orbits_from_candidate_states.py")),
        str(state_path),
        "--output-h5",
        str(output_h5_for_state(state_path)),
        "--metadata-dir",
        str(args.orbit_metadata_dir),
    ]
    log_path = state_path.parent / "logs" / f"{sample_idx}.orbit_backfill.rank{rank:03d}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        log.write(" ".join(cmd) + "\n\n")
        log.flush()
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=Path(__file__).resolve().parent)
    status = "ok" if proc.returncode == 0 and orbit_metadata_table.has_sample(args.orbit_metadata_dir, sample_idx) else "failed"
    return status, time.time() - t0, sample_idx, str(log_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events-dir", type=Path, default=Path("test_plots"))
    parser.add_argument("--orbit-metadata-dir", type=Path, default=Path("data/metadata/orbit"))
    parser.add_argument("--start-day", type=str, default=None)
    parser.add_argument("--end-day", type=str, default=None)
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        state_paths = discover_state_paths(args.events_dir, args.start_day, args.end_day)
        if args.limit is not None:
            state_paths = state_paths[: args.limit]
        print(f"orbit_backfill_discovered {len(state_paths)} state_files ranks {size}", flush=True)
    else:
        state_paths = None
    state_paths = comm.bcast(state_paths, root=0)

    counts = {"ok": 0, "skip": 0, "failed": 0}
    my_paths = state_paths[rank::size]
    for n_done, state_path in enumerate(my_paths, start=1):
        status, elapsed, sample_idx, log_path = run_one(state_path, args, rank)
        counts[status] += 1
        print(
            f"rank {rank:03d} {status} sample_idx {sample_idx} day {day_from_sample_idx(sample_idx)} "
            f"elapsed_s {elapsed:.1f} progress {n_done}/{len(my_paths)} log {log_path}",
            flush=True,
        )

    gathered = comm.gather(counts, root=0)
    if rank == 0:
        total = {"ok": 0, "skip": 0, "failed": 0}
        for item in gathered:
            for key in total:
                total[key] += item[key]
        print(f"orbit_backfill_complete ok {total['ok']} skip {total['skip']} failed {total['failed']}", flush=True)
        if total["failed"]:
            sys.exit(1)


if __name__ == "__main__":
    main()
