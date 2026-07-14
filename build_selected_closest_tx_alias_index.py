#!/usr/bin/env python3
"""Build an HDF5 index of events whose selected alias is closest to the TX beam.

This uses per-event disambiguation diagnostics, not compact orbit metadata,
because the compact orbit alias table only retains a subset of aliases.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
from pathlib import Path

import h5py
import numpy as np


RESULT_DTYPE = np.dtype(
    [
        ("sample_idx", "<i8"),
        ("selected_hypothesis", "S8"),
        ("closest_hypothesis", "S8"),
        ("selected_tx_beam_center_distance_dc", "<f4"),
        ("closest_tx_beam_center_distance_dc", "<f4"),
        ("match", "?"),
    ]
)


def diagnostic_paths(events_dir: Path):
    yield from sorted(events_dir.glob("**/pansy_disambiguation_diagnostics_*.h5"))


def read_one(path: Path):
    try:
        with h5py.File(path, "r") as h:
            sample_idx = int(h.attrs.get("sample_idx", -1))
            selected = h.attrs.get("selected_hypothesis", b"")
            if isinstance(selected, str):
                selected = selected.encode("utf-8", "ignore")
            selected = np.asarray(selected, dtype="S8").item()
            if sample_idx < 0 or selected == b"" or "hypotheses" not in h:
                return None
            labels = []
            tx_distance = []
            for label, group in h["hypotheses"].items():
                value = float(group.attrs.get("tx_beam_snr_weighted_mean_dc", np.nan))
                if np.isfinite(value):
                    labels.append(label.encode("utf-8", "ignore"))
                    tx_distance.append(value)
            if not labels:
                return None
            labels = np.asarray(labels, dtype="S8")
            tx_distance = np.asarray(tx_distance, dtype=np.float64)
            closest_i = int(np.argmin(tx_distance))
            selected_i = np.flatnonzero(labels == selected)
            selected_tx = float(tx_distance[selected_i[0]]) if len(selected_i) else np.nan
            closest = labels[closest_i]
            return (
                sample_idx,
                selected,
                closest,
                selected_tx,
                float(tx_distance[closest_i]),
                bool(selected == closest),
            )
    except OSError:
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events-dir", type=Path, default=Path("/mnt/data/juha/pansy/events"))
    parser.add_argument("--output-h5", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=10000)
    args = parser.parse_args()

    rows = []
    paths = list(diagnostic_paths(args.events_dir))
    files_read = len(paths)
    if args.workers <= 1:
        iterator = map(read_one, paths)
        for n_done, row in enumerate(iterator, start=1):
            if row is not None:
                rows.append(row)
            if args.progress_every > 0 and n_done % int(args.progress_every) == 0:
                print(f"progress {n_done}/{files_read} rows={len(rows)}", flush=True)
    else:
        with mp.Pool(processes=int(args.workers)) as pool:
            for n_done, row in enumerate(pool.imap_unordered(read_one, paths, chunksize=128), start=1):
                if row is not None:
                    rows.append(row)
                if args.progress_every > 0 and n_done % int(args.progress_every) == 0:
                    print(f"progress {n_done}/{files_read} rows={len(rows)}", flush=True)

    table = np.zeros(len(rows), dtype=RESULT_DTYPE)
    for i, row in enumerate(rows):
        table[i] = row
    if len(table):
        table = table[np.argsort(table["sample_idx"])]
    matched = table[table["match"]]
    mismatched = table[~table["match"]]

    args.output_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output_h5, "w") as h:
        h.attrs["script"] = Path(__file__).name
        h.attrs["events_dir"] = str(args.events_dir)
        h.attrs["diagnostic_files_read"] = int(files_read)
        h.attrs["events_with_selected_and_tx"] = int(len(table))
        h.attrs["matched_event_count"] = int(len(matched))
        h.attrs["mismatched_event_count"] = int(len(mismatched))
        h.create_dataset("table", data=table, compression="gzip", shuffle=True)
        h.create_dataset("sample_idx", data=matched["sample_idx"], compression="gzip", shuffle=True)
        h.create_dataset("matched_sample_idx", data=matched["sample_idx"], compression="gzip", shuffle=True)
        h.create_dataset("mismatched_sample_idx", data=mismatched["sample_idx"], compression="gzip", shuffle=True)

    print(
        f"selected_closest_tx_alias_index files_read={files_read} "
        f"events={len(table)} matched={len(matched)} mismatched={len(mismatched)}"
    )
    print(args.output_h5)


if __name__ == "__main__":
    main()
