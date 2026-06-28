#!/usr/bin/env python3
"""Rewrite PANSY orbit metadata files with only catalogue-scale fields."""

from __future__ import annotations

import argparse
import fcntl
import multiprocessing as mp
import os
import shutil
import subprocess
import time
from pathlib import Path

import h5py
import numpy as np

import orbit_metadata_table as omt


RETAIN_DATASETS = {
    "sample_idx",
    "selected_hypothesis",
    "candidate_number",
    "combined_rank",
    "combined_score",
    "selection_model_type",
    "orbit_solution_type",
    "initial_state_gcrs_m_mps",
    "fit_parameters",
    "fit_parameter_covariance",
    "ceplecha_initial_radius_m",
    "ceplecha_initial_radius_std_m",
    "ceplecha_initial_mass_kg",
    "ceplecha_initial_mass_std_kg",
    "ceplecha_log10_radius_std",
    "ceplecha_covariance_available",
    "ceplecha_reduced_chi2",
    "ceplecha_bic",
    "ceplecha_n",
    "ceplecha_dof",
    "ceplecha_parameters",
    "ceplecha_parameter_std",
    "ceplecha_parameter_covariance",
    "log10_beta_kg_m2",
    "sigma_log10_beta",
    "initial_detection_height_km",
    "kepler",
    "kepler_std",
    "kepler_covariance",
    "n_uncertainty_samples",
    "frac_e_gt_1",
    "radiant_ra_deg",
    "radiant_dec_deg",
    "radiant_speed_km_s",
    "v_g_km_s",
    "radiant_ecliptic_lon_deg",
    "radiant_ecliptic_lat_deg",
    "radiant_sun_ecliptic_lon_deg",
    "radiant_sun_ecliptic_lat_deg",
    "path_t_rel_s",
    "path_position_enu_km",
    "path_doppler_mps",
    "path_snr",
    "path_beam_id",
    "alias_hypothesis_labels",
    "alias_candidate_numbers",
    "alias_combined_ranks",
    "alias_combined_scores",
    "alias_selection_model_types",
    "alias_plausibility_models",
    "alias_plausibility_redchi",
    "alias_tx_beam_snr_weighted_mean_dc",
    "alias_tx_beam_snr_weighted_rms_dc",
    "alias_tx_beam_weighted_mean_deg",
    "alias_tx_beam_weighted_rms_deg",
    "alias_tx_lobe_snr_weighted_mean_dc",
    "alias_tx_lobe_snr_weighted_rms_dc",
    "alias_tx_lobe_p90_dc",
    "alias_kepler",
    "alias_kepler_std",
    "alias_frac_e_gt_1",
    "alias_interstellar_nominal",
    "all_aliases_interstellar_nominal",
    "n_aliases_orbit_tested",
}


def copy_attrs(src, dst) -> None:
    for key, value in src.attrs.items():
        dst.attrs[key] = value


def copy_retained_group(src: h5py.Group, dst: h5py.Group) -> tuple[int, int]:
    kept = 0
    dropped = 0
    copy_attrs(src, dst)
    for name, obj in src.items():
        if isinstance(obj, h5py.Dataset):
            if name in RETAIN_DATASETS:
                src.copy(name, dst)
                kept += 1
            else:
                dropped += 1
        elif isinstance(obj, h5py.Group):
            child = dst.create_group(name)
            child_kept, child_dropped = copy_retained_group(obj, child)
            kept += child_kept
            dropped += child_dropped
    return kept, dropped


def rewrite_table_file(src: h5py.File, dst: h5py.File) -> tuple[int, int]:
    events = omt._coerce_structured(src["events"][()] if "events" in src else np.zeros(0, omt.EVENT_DTYPE), omt.EVENT_DTYPE)
    aliases = omt._coerce_structured(src["aliases"][()] if "aliases" in src else np.zeros(0, omt.ALIAS_DTYPE), omt.ALIAS_DTYPE)
    paths = omt._coerce_structured(src["paths"][()] if "paths" in src else np.zeros(0, omt.PATH_DTYPE), omt.PATH_DTYPE)
    dst.attrs["schema"] = "pansy_orbit_table_v2"
    dst.create_dataset("events", data=events)
    dst.create_dataset("aliases", data=aliases)
    dst.create_dataset("paths", data=paths)
    return 3, 0


def compact_file(path: Path, repack: bool, lock_path: Path | None = None) -> tuple[int, int, int, int]:
    tmp = path.with_suffix(path.suffix + ".compact_tmp")
    repacked = path.with_suffix(path.suffix + ".repack_tmp")
    for extra in (tmp, repacked):
        if extra.exists():
            extra.unlink()
    before = path.stat().st_size
    kept = 0
    dropped = 0
    lock_fh = None
    try:
        if lock_path is not None:
            lock_path.parent.mkdir(parents=True, exist_ok=True)
            lock_fh = lock_path.open("w")
            fcntl.flock(lock_fh, fcntl.LOCK_EX)
        with h5py.File(path, "r") as src, h5py.File(tmp, "w") as dst:
            copy_attrs(src, dst)
            if "events" in src:
                kept, dropped = rewrite_table_file(src, dst)
            else:
                for name, obj in src.items():
                    if isinstance(obj, h5py.Group):
                        group = dst.create_group(name)
                        group_kept, group_dropped = copy_retained_group(obj, group)
                        kept += group_kept
                        dropped += group_dropped
                    elif isinstance(obj, h5py.Dataset):
                        if name in RETAIN_DATASETS:
                            src.copy(name, dst)
                            kept += 1
                        else:
                            dropped += 1
        final_tmp = tmp
        if repack:
            h5repack = shutil.which("h5repack")
            if h5repack is not None:
                subprocess.run([h5repack, str(tmp), str(repacked)], check=True)
                tmp.unlink()
                final_tmp = repacked
        with h5py.File(final_tmp, "r"):
            pass
        os.replace(final_tmp, path)
    finally:
        if lock_fh is not None:
            fcntl.flock(lock_fh, fcntl.LOCK_UN)
            lock_fh.close()
    after = path.stat().st_size
    return before, after, kept, dropped


def compact_task(task: tuple[int, int, str, bool, str | None]) -> tuple[int, int, int, int, int, int, str]:
    index, total, path_str, repack, lock_path_str = task
    lock_path = Path(lock_path_str) if lock_path_str else None
    before, after, kept, dropped = compact_file(Path(path_str), repack=repack, lock_path=lock_path)
    return index, total, before, after, kept, dropped, path_str


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orbit_metadata_dir", type=Path)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--no-repack", action="store_true")
    parser.add_argument("--use-writer-lock", action="store_true", help="Serialize rewrites with orbit_metadata_table writers.")
    parser.add_argument("--skip-recent-minutes", type=float, default=0.0, help="Skip files modified in the last N minutes.")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    files = sorted(args.orbit_metadata_dir.glob("**/orbit@*.h5"))
    if args.skip_recent_minutes > 0.0:
        cutoff = time.time() - args.skip_recent_minutes * 60.0
        files = [path for path in files if path.stat().st_mtime < cutoff]
    if args.limit is not None:
        files = files[: args.limit]
    total_before = 0
    total_after = 0
    total_kept = 0
    total_dropped = 0
    for path in files:
        before = path.stat().st_size
        total_before += before
        if args.dry_run:
            print(f"dry_run {path} before={before}")
    if args.dry_run:
        print(f"compact_orbit_done files={len(files)} before={total_before} after=0 kept=0 dropped=0", flush=True)
        return

    lock_path = args.orbit_metadata_dir / ".orbit_metadata.lock" if args.use_writer_lock else None
    tasks = [(i, len(files), str(path), not args.no_repack, str(lock_path) if lock_path is not None else None) for i, path in enumerate(files, start=1)]
    if args.workers <= 1:
        iterator = map(compact_task, tasks)
    else:
        pool = mp.get_context("spawn").Pool(processes=args.workers)
        iterator = pool.imap_unordered(compact_task, tasks, chunksize=4)
    try:
        for i, n_files, after_before, after, kept, dropped, path_str in iterator:
            total_after += after
            total_kept += kept
            total_dropped += dropped
            print(
                f"compact_orbit {i}/{n_files} before={after_before} after={after} "
                f"kept={kept} dropped={dropped} {path_str}",
                flush=True,
            )
    finally:
        if args.workers > 1:
            pool.close()
            pool.join()
    print(
        f"compact_orbit_done files={len(files)} before={total_before} after={total_after} "
        f"kept={total_kept} dropped={total_dropped}",
        flush=True,
    )


if __name__ == "__main__":
    main()
