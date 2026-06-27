#!/usr/bin/env python3
"""Convert existing group-style orbit metadata files to compact table files."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from pathlib import Path

import h5py
import numpy as np

import orbit_metadata_table as omt


def dataset(group: h5py.Group, name: str, default):
    if name in group:
        return group[name][()]
    return default


def payload_from_packed_group(group: h5py.Group) -> dict:
    event = group["event"][()]
    payload = {name: event[name] for name in event.dtype.names}
    for name in [
        "initial_state_gcrs_m_mps",
        "fit_parameters",
        "fit_parameter_covariance",
        "ceplecha_parameters",
        "ceplecha_parameter_std",
        "ceplecha_parameter_covariance",
        "kepler",
        "kepler_std",
        "kepler_covariance",
    ]:
        payload[name] = dataset(group, name, np.asarray([]))
    if "aliases" in group:
        aliases = group["aliases"][()]
        payload.update(
            {
                "alias_hypothesis_labels": aliases["hypothesis_label"],
                "alias_candidate_numbers": aliases["candidate_number"],
                "alias_combined_ranks": aliases["combined_rank"],
                "alias_combined_scores": aliases["combined_score"],
                "alias_selection_model_types": aliases["selection_model_type"],
                "alias_plausibility_models": aliases["plausibility_model"],
                "alias_plausibility_redchi": aliases["plausibility_redchi"],
                "alias_tx_beam_snr_weighted_mean_dc": aliases["tx_beam_snr_weighted_mean_dc"],
                "alias_tx_beam_snr_weighted_rms_dc": aliases["tx_beam_snr_weighted_rms_dc"],
                "alias_tx_beam_weighted_mean_deg": aliases["tx_beam_weighted_mean_deg"],
                "alias_tx_beam_weighted_rms_deg": aliases["tx_beam_weighted_rms_deg"],
                "alias_tx_lobe_snr_weighted_mean_dc": aliases["tx_lobe_snr_weighted_mean_dc"],
                "alias_tx_lobe_snr_weighted_rms_dc": aliases["tx_lobe_snr_weighted_rms_dc"],
                "alias_tx_lobe_p90_dc": aliases["tx_lobe_p90_dc"],
                "alias_kepler": aliases["kepler"],
                "alias_kepler_std": aliases["kepler_std"],
                "alias_frac_e_gt_1": aliases["frac_e_gt_1"],
                "alias_interstellar_nominal": aliases["interstellar_nominal"],
            }
        )
    if "path" in group:
        path = group["path"][()]
        payload.update(
            {
                "path_t_rel_s": path["t_rel_s"],
                "path_position_enu_km": path["position_enu_km"],
                "path_snr": path["snr"],
                "path_beam_id": path["beam_id"],
            }
        )
    return payload


def payload_from_group(group: h5py.Group) -> dict:
    if "event" in group:
        return payload_from_packed_group(group)
    return {name: obj[()] for name, obj in group.items() if isinstance(obj, h5py.Dataset)}


def tableize_file(path: Path) -> tuple[int, int, int]:
    before = path.stat().st_size
    events = []
    aliases = []
    paths = []
    with h5py.File(path, "r") as src:
        if "events" in src:
            return before, before, int(len(src["events"]))
        for name, obj in src.items():
            if not isinstance(obj, h5py.Group):
                continue
            payload = payload_from_group(obj)
            sample_idx = int(payload.get("sample_idx", int(name)))
            event, alias, path_rows = omt.payload_to_rows(sample_idx, payload)
            events.append(event)
            aliases.append(alias)
            paths.append(path_rows)
    event_table = np.concatenate(events) if events else np.zeros(0, dtype=omt.EVENT_DTYPE)
    alias_table = np.concatenate(aliases) if aliases else np.zeros(0, dtype=omt.ALIAS_DTYPE)
    path_table = np.concatenate(paths) if paths else np.zeros(0, dtype=omt.PATH_DTYPE)
    order = np.argsort(event_table["sample_idx"]) if len(event_table) else []
    event_table = event_table[order]
    tmp = path.with_suffix(path.suffix + ".table_tmp")
    if tmp.exists():
        tmp.unlink()
    with h5py.File(tmp, "w") as dst:
        dst.attrs["schema"] = "pansy_orbit_table_v1"
        dst.create_dataset("events", data=event_table)
        dst.create_dataset("aliases", data=alias_table)
        dst.create_dataset("paths", data=path_table)
    with h5py.File(tmp, "r"):
        pass
    os.replace(tmp, path)
    return before, path.stat().st_size, len(event_table)


def task(args):
    i, n, path = args
    before, after, n_events = tableize_file(Path(path))
    return i, n, before, after, n_events, path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orbit_metadata_dir", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()
    files = sorted(args.orbit_metadata_dir.glob("**/orbit@*.h5"))
    if args.limit is not None:
        files = files[: args.limit]
    tasks = [(i, len(files), str(path)) for i, path in enumerate(files, 1)]
    total_before = sum(Path(path).stat().st_size for _i, _n, path in tasks)
    total_after = 0
    total_events = 0
    if args.workers <= 1:
        iterator = map(task, tasks)
    else:
        pool = mp.get_context("spawn").Pool(args.workers)
        iterator = pool.imap_unordered(task, tasks, chunksize=8)
    try:
        for i, n, before, after, n_events, path in iterator:
            total_after += after
            total_events += n_events
            print(f"tableize_orbit {i}/{n} before={before} after={after} events={n_events} {path}", flush=True)
    finally:
        if args.workers > 1:
            pool.close()
            pool.join()
    print(f"tableize_orbit_done files={len(files)} events={total_events} before={total_before} after={total_after}", flush=True)


if __name__ == "__main__":
    main()
