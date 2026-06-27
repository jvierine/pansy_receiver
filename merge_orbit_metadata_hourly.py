#!/usr/bin/env python3
"""Merge existing PANSY orbit minute files into hourly table files."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from pathlib import Path

import h5py
import numpy as np

import orbit_metadata_table as omt
from tableize_orbit_metadata import payload_from_group


def read_any_orbit_file(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    events = []
    aliases = []
    paths = []
    with h5py.File(path, "r") as h:
        if "events" in h:
            return (
                h["events"][()].astype(omt.EVENT_DTYPE, copy=False),
                h["aliases"][()].astype(omt.ALIAS_DTYPE, copy=False) if "aliases" in h else np.zeros(0, omt.ALIAS_DTYPE),
                h["paths"][()].astype(omt.PATH_DTYPE, copy=False) if "paths" in h else np.zeros(0, omt.PATH_DTYPE),
            )
        for name, obj in h.items():
            if not isinstance(obj, h5py.Group):
                continue
            payload = payload_from_group(obj)
            sample_idx = int(payload.get("sample_idx", int(name)))
            event, alias_rows, path_rows = omt.payload_to_rows(sample_idx, payload)
            events.append(event)
            aliases.append(alias_rows)
            paths.append(path_rows)
    return (
        np.concatenate(events) if events else np.zeros(0, omt.EVENT_DTYPE),
        np.concatenate(aliases) if aliases else np.zeros(0, omt.ALIAS_DTYPE),
        np.concatenate(paths) if paths else np.zeros(0, omt.PATH_DTYPE),
    )


def merge_hour(task):
    i, n, hour_dir_str = task
    hour_dir = Path(hour_dir_str)
    files = sorted(hour_dir.glob("orbit@*.h5"))
    events = []
    aliases = []
    paths = []
    before = sum(path.stat().st_size for path in files)
    for path in files:
        e, a, p = read_any_orbit_file(path)
        events.append(e)
        aliases.append(a)
        paths.append(p)
    event_table = np.concatenate(events) if events else np.zeros(0, omt.EVENT_DTYPE)
    alias_table = np.concatenate(aliases) if aliases else np.zeros(0, omt.ALIAS_DTYPE)
    path_table = np.concatenate(paths) if paths else np.zeros(0, omt.PATH_DTYPE)
    if len(event_table):
        _, unique_idx = np.unique(event_table["sample_idx"], return_index=True)
        event_table = event_table[np.sort(unique_idx)]
        event_table = event_table[np.argsort(event_table["sample_idx"])]
        keep_samples = np.isin(alias_table["sample_idx"], event_table["sample_idx"]) if len(alias_table) else []
        alias_table = alias_table[keep_samples] if len(alias_table) else alias_table
        keep_samples = np.isin(path_table["sample_idx"], event_table["sample_idx"]) if len(path_table) else []
        path_table = path_table[keep_samples] if len(path_table) else path_table
    hour_s = int(hour_dir.name.replace("-", ":").replace("T", " ").split(":")[0] and 0)
    if len(event_table):
        hour_s = int((int(event_table["sample_idx"][0]) // 1_000_000 // 3600) * 3600)
    else:
        first = files[0].stem.split("@")[-1]
        hour_s = (int(first) // 3600) * 3600
    out = hour_dir / f"orbit@{hour_s}.h5"
    tmp = out.with_suffix(out.suffix + ".hour_tmp")
    if tmp.exists():
        tmp.unlink()
    with h5py.File(tmp, "w") as h:
        h.attrs["schema"] = "pansy_orbit_table_v1"
        h.attrs["file_cadence_seconds"] = 3600
        h.create_dataset("events", data=event_table)
        h.create_dataset("aliases", data=alias_table)
        h.create_dataset("paths", data=path_table)
    with h5py.File(tmp, "r"):
        pass
    os.replace(tmp, out)
    for path in files:
        if path != out:
            path.unlink()
    after = out.stat().st_size
    return i, n, before, after, len(files), len(event_table), str(out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orbit_metadata_dir", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()
    hour_dirs = sorted(path for path in args.orbit_metadata_dir.iterdir() if path.is_dir() and list(path.glob("orbit@*.h5")))
    if args.limit is not None:
        hour_dirs = hour_dirs[: args.limit]
    tasks = [(i, len(hour_dirs), str(path)) for i, path in enumerate(hour_dirs, 1)]
    total_before = 0
    total_after = 0
    total_events = 0
    if args.workers <= 1:
        iterator = map(merge_hour, tasks)
    else:
        pool = mp.get_context("spawn").Pool(args.workers)
        iterator = pool.imap_unordered(merge_hour, tasks, chunksize=1)
    try:
        for i, n, before, after, n_files, n_events, out in iterator:
            total_before += before
            total_after += after
            total_events += n_events
            print(
                f"merge_orbit_hour {i}/{n} files={n_files} events={n_events} "
                f"before={before} after={after} {out}",
                flush=True,
            )
    finally:
        if args.workers > 1:
            pool.close()
            pool.join()
    print(f"merge_orbit_hour_done hours={len(hour_dirs)} events={total_events} before={total_before} after={total_after}", flush=True)


if __name__ == "__main__":
    main()
