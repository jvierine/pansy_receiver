#!/usr/bin/env python3
"""Backfill robust-fit path masks into compact orbit metadata tables."""

from __future__ import annotations

import argparse
import datetime as dt
import multiprocessing as mp
import os
import time
from pathlib import Path

import h5py
import numpy as np

import orbit_metadata_table as omt


def candidate_state_path(events_dir: Path, sample_idx: int) -> Path:
    epoch = dt.datetime.fromtimestamp(sample_idx / 1_000_000.0, tz=dt.timezone.utc)
    return events_dir / epoch.strftime("%Y-%m-%d") / f"pansy_candidate_orbit_states_{sample_idx}.h5"


def winning_group(handle: h5py.File, label: bytes | str) -> h5py.Group | None:
    name = label.decode("ascii", "ignore") if isinstance(label, bytes) else str(label)
    if name in handle and isinstance(handle[name], h5py.Group):
        return handle[name]
    rank_zero = [group for group in handle.values() if isinstance(group, h5py.Group) and int(group.attrs.get("combined_rank", -1)) == 0]
    return rank_zero[0] if len(rank_zero) == 1 else None


def load_selection_keep(path: Path, label: bytes | str, expected_length: int) -> np.ndarray | None:
    try:
        with h5py.File(path, "r") as handle:
            group = winning_group(handle, label)
            if group is None or "path_selection_keep" not in group:
                return None
            keep = np.asarray(group["path_selection_keep"][()], dtype=bool).ravel()
    except OSError:
        return None
    return keep if len(keep) == expected_length else None


def rewrite_file(path: Path, events_dir: Path) -> tuple[int, int, int, int]:
    with h5py.File(path, "r") as src:
        attrs = dict(src.attrs)
        events = omt._coerce_structured(src["events"][()] if "events" in src else np.zeros(0, omt.EVENT_DTYPE), omt.EVENT_DTYPE)
        aliases = omt._coerce_structured(src["aliases"][()] if "aliases" in src else np.zeros(0, omt.ALIAS_DTYPE), omt.ALIAS_DTYPE)
        paths = omt._coerce_structured(src["paths"][()] if "paths" in src else np.zeros(0, omt.PATH_DTYPE), omt.PATH_DTYPE)

    paths["selection_keep"] = False
    labels = {int(row["sample_idx"]): row["selected_hypothesis"] for row in events}
    found = 0
    missing = 0
    kept = 0
    for sample_idx in np.unique(paths["sample_idx"]):
        row_idx = np.flatnonzero(paths["sample_idx"] == sample_idx)
        mask = load_selection_keep(
            candidate_state_path(events_dir, int(sample_idx)),
            labels.get(int(sample_idx), b""),
            len(row_idx),
        )
        if mask is None:
            missing += 1
            continue
        paths["selection_keep"][row_idx] = mask
        found += 1
        kept += int(np.count_nonzero(mask))

    tmp = path.with_suffix(path.suffix + ".selection_keep_tmp")
    if tmp.exists():
        tmp.unlink()
    with h5py.File(tmp, "w") as dst:
        for name, value in attrs.items():
            dst.attrs[name] = value
        dst.attrs["schema"] = "pansy_orbit_table_v3"
        dst.attrs["path_selection_keep_source"] = "winning candidate-state robust-fit mask"
        dst.attrs["path_selection_keep_events_found"] = found
        dst.attrs["path_selection_keep_events_missing"] = missing
        dst.create_dataset("events", data=events)
        dst.create_dataset("aliases", data=aliases)
        dst.create_dataset("paths", data=paths)
    with h5py.File(tmp, "r"):
        pass
    os.replace(tmp, path)
    return len(events), found, missing, kept


def task(args: tuple[int, int, str, str]) -> tuple[int, int, str, int, int, int, int]:
    index, total, path_str, events_dir_str = args
    n_events, found, missing, kept = rewrite_file(Path(path_str), Path(events_dir_str))
    return index, total, path_str, n_events, found, missing, kept


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-metadata-dir", type=Path, required=True)
    parser.add_argument("--events-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--skip-recent-minutes", type=float, default=60.0)
    args = parser.parse_args()

    files = sorted(args.orbit_metadata_dir.glob("**/orbit@*.h5"))
    if args.skip_recent_minutes > 0:
        cutoff = time.time() - args.skip_recent_minutes * 60.0
        files = [path for path in files if path.stat().st_mtime < cutoff]
    if args.limit is not None:
        files = files[: args.limit]
    tasks = [(i, len(files), str(path), str(args.events_dir)) for i, path in enumerate(files, 1)]

    found_total = 0
    missing_total = 0
    kept_total = 0
    if args.workers <= 1:
        iterator = map(task, tasks)
        pool = None
    else:
        pool = mp.get_context("spawn").Pool(args.workers)
        iterator = pool.imap_unordered(task, tasks, chunksize=1)
    try:
        for index, total, path, n_events, found, missing, kept in iterator:
            found_total += found
            missing_total += missing
            kept_total += kept
            print(
                f"backfill_path_keep {index}/{total} events={n_events} found={found} "
                f"missing={missing} kept_points={kept} {path}",
                flush=True,
            )
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    print(
        f"backfill_path_keep_done files={len(files)} found={found_total} "
        f"missing={missing_total} kept_points={kept_total}",
        flush=True,
    )


if __name__ == "__main__":
    main()
