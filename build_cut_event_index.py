#!/usr/bin/env python3
"""Build a compact HDF5 event index from PANSY cut metadata."""

from __future__ import annotations

import argparse
import datetime as dt
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import h5py
import numpy as np


def _read_keys(path: str) -> list[int]:
    try:
        with h5py.File(path, "r") as handle:
            keys = [int(key) for key in handle.keys() if key.isdigit()]
    except OSError:
        keys = []
    return keys


def _read_key_chunk(paths: list[str]) -> np.ndarray:
    keys: list[int] = []
    for path in paths:
        keys.extend(_read_keys(path))
    return np.asarray(keys, dtype=np.int64)


def build_index(cut_dir: Path, output_h5: Path, workers: int = 16) -> int:
    files = sorted(str(path) for path in cut_dir.glob("20*/*.h5"))
    workers = max(1, int(workers))
    chunk_size = max(64, int(np.ceil(len(files) / (workers * 16)))) if files else 64
    file_chunks = [files[i : i + chunk_size] for i in range(0, len(files), chunk_size)]
    chunks: list[np.ndarray] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_read_key_chunk, chunk) for chunk in file_chunks]
        for n_done, future in enumerate(as_completed(futures), start=1):
            keys = future.result()
            if keys.size:
                chunks.append(keys)
            if n_done % 25 == 0 or n_done == len(file_chunks):
                indexed_files = min(n_done * chunk_size, len(files))
                print(
                    f"indexed_files {indexed_files}/{len(files)} "
                    f"chunks {n_done}/{len(file_chunks)} events_so_far {sum(len(c) for c in chunks)}",
                    flush=True,
                )
    if chunks:
        sample_idx = np.unique(np.concatenate(chunks).astype(np.int64))
    else:
        sample_idx = np.empty(0, dtype=np.int64)
    output_h5.parent.mkdir(parents=True, exist_ok=True)
    tmp = output_h5.with_suffix(output_h5.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()
    with h5py.File(tmp, "w") as h5:
        h5.attrs["schema_version"] = "pansy_cut_event_index_v1"
        h5.attrs["source_program"] = "build_cut_event_index.py"
        h5.attrs["cut_dir"] = str(cut_dir)
        h5.attrs["source_file_count"] = int(len(files))
        h5.attrs["source_chunk_count"] = int(len(file_chunks))
        h5.attrs["source_chunk_size"] = int(chunk_size)
        h5.attrs["event_count"] = int(sample_idx.size)
        h5.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        h5.create_dataset("sample_idx", data=sample_idx, compression="gzip", shuffle=True)
    tmp.replace(output_h5)
    return int(sample_idx.size)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cut-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    count = build_index(args.cut_dir, args.output, workers=args.workers)
    print(f"wrote {args.output}")
    print(f"event_count {count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
