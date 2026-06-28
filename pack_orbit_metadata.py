#!/usr/bin/env python3
"""Pack compact PANSY orbit metadata into a small number of HDF5 datasets."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np


EVENT_DTYPE = np.dtype(
    [
        ("sample_idx", "<i8"),
        ("candidate_number", "<i4"),
        ("combined_rank", "<i4"),
        ("combined_score", "<f4"),
        ("log10_beta_kg_m2", "<f4"),
        ("sigma_log10_beta", "<f4"),
        ("initial_detection_height_km", "<f4"),
        ("ceplecha_initial_radius_m", "<f4"),
        ("ceplecha_initial_radius_std_m", "<f4"),
        ("ceplecha_initial_mass_kg", "<f4"),
        ("ceplecha_initial_mass_std_kg", "<f4"),
        ("ceplecha_log10_radius_std", "<f4"),
        ("ceplecha_covariance_available", "?"),
        ("ceplecha_reduced_chi2", "<f4"),
        ("ceplecha_bic", "<f4"),
        ("ceplecha_n", "<i4"),
        ("ceplecha_dof", "<i4"),
        ("n_uncertainty_samples", "<i4"),
        ("frac_e_gt_1", "<f4"),
        ("radiant_ra_deg", "<f4"),
        ("radiant_dec_deg", "<f4"),
        ("radiant_speed_km_s", "<f4"),
        ("v_g_km_s", "<f4"),
        ("radiant_ecliptic_lon_deg", "<f4"),
        ("radiant_ecliptic_lat_deg", "<f4"),
        ("radiant_sun_ecliptic_lon_deg", "<f4"),
        ("radiant_sun_ecliptic_lat_deg", "<f4"),
        ("all_aliases_interstellar_nominal", "?"),
        ("n_aliases_orbit_tested", "<i4"),
        ("selected_hypothesis", "S8"),
        ("selection_model_type", "S32"),
        ("orbit_solution_type", "S32"),
    ]
)

ALIAS_DTYPE = np.dtype(
    [
        ("hypothesis_label", "S8"),
        ("candidate_number", "<i4"),
        ("combined_rank", "<i4"),
        ("combined_score", "<f4"),
        ("selection_model_type", "S32"),
        ("plausibility_model", "S32"),
        ("plausibility_redchi", "<f4"),
        ("tx_beam_snr_weighted_mean_dc", "<f4"),
        ("tx_beam_snr_weighted_rms_dc", "<f4"),
        ("tx_beam_weighted_mean_deg", "<f4"),
        ("tx_beam_weighted_rms_deg", "<f4"),
        ("tx_lobe_snr_weighted_mean_dc", "<f4"),
        ("tx_lobe_snr_weighted_rms_dc", "<f4"),
        ("tx_lobe_p90_dc", "<f4"),
        ("kepler", "<f4", (7,)),
        ("kepler_std", "<f4", (7,)),
        ("frac_e_gt_1", "<f4"),
        ("interstellar_nominal", "?"),
    ]
)

PATH_DTYPE = np.dtype(
    [
        ("t_rel_s", "<f4"),
        ("position_enu_km", "<f4", (3,)),
        ("doppler_mps", "<f4"),
        ("snr", "<f4"),
        ("beam_id", "<i1"),
    ]
)


def dataset(group: h5py.Group, name: str, default):
    if name in group:
        return group[name][()]
    return default


def scalar(group: h5py.Group, name: str, default=np.nan):
    value = dataset(group, name, default)
    arr = np.asarray(value)
    return arr.item() if arr.shape == () else value


def fixed_string(group: h5py.Group, name: str, length: int):
    value = scalar(group, name, b"")
    if isinstance(value, str):
        value = value.encode("utf-8", "ignore")
    return np.asarray(value, dtype=f"S{length}").item()


def coerce_structured(arr: np.ndarray, dtype: np.dtype) -> np.ndarray:
    arr = np.asarray(arr)
    out = np.zeros(arr.shape, dtype=dtype)
    if arr.dtype.names is None:
        return out
    for name in dtype.names:
        if name in arr.dtype.names:
            out[name] = arr[name]
        elif dtype[name].kind == "f":
            out[name] = np.nan
    return out


def make_event(group: h5py.Group) -> np.ndarray:
    row = np.zeros((), dtype=EVENT_DTYPE)
    for name in EVENT_DTYPE.names:
        if EVENT_DTYPE[name].kind == "S":
            row[name] = fixed_string(group, name, EVENT_DTYPE[name].itemsize)
        elif EVENT_DTYPE[name].kind == "b":
            row[name] = bool(scalar(group, name, False))
        elif EVENT_DTYPE[name].kind in "iu":
            row[name] = int(scalar(group, name, 0))
        else:
            row[name] = float(scalar(group, name, np.nan))
    return row


def make_aliases(group: h5py.Group) -> np.ndarray:
    labels = np.asarray(dataset(group, "alias_hypothesis_labels", np.asarray([], dtype="S8")), dtype="S8")
    labels = np.atleast_1d(labels)
    if labels.shape == (1,) and labels[0] == b"":
        labels = labels[:0]
    n = len(labels)
    arr = np.zeros(n, dtype=ALIAS_DTYPE)
    arr["hypothesis_label"] = labels
    fields = [
        ("candidate_number", "alias_candidate_numbers", 0),
        ("combined_rank", "alias_combined_ranks", 0),
        ("combined_score", "alias_combined_scores", np.nan),
        ("selection_model_type", "alias_selection_model_types", b""),
        ("plausibility_model", "alias_plausibility_models", b""),
        ("plausibility_redchi", "alias_plausibility_redchi", np.nan),
        ("tx_beam_snr_weighted_mean_dc", "alias_tx_beam_snr_weighted_mean_dc", np.nan),
        ("tx_beam_snr_weighted_rms_dc", "alias_tx_beam_snr_weighted_rms_dc", np.nan),
        ("tx_beam_weighted_mean_deg", "alias_tx_beam_weighted_mean_deg", np.nan),
        ("tx_beam_weighted_rms_deg", "alias_tx_beam_weighted_rms_deg", np.nan),
        ("tx_lobe_snr_weighted_mean_dc", "alias_tx_lobe_snr_weighted_mean_dc", np.nan),
        ("tx_lobe_snr_weighted_rms_dc", "alias_tx_lobe_snr_weighted_rms_dc", np.nan),
        ("tx_lobe_p90_dc", "alias_tx_lobe_p90_dc", np.nan),
        ("frac_e_gt_1", "alias_frac_e_gt_1", np.nan),
        ("interstellar_nominal", "alias_interstellar_nominal", False),
    ]
    for dst, src, default in fields:
        vals = np.atleast_1d(np.asarray(dataset(group, src, np.full(n, default))))
        if len(vals) == n:
            arr[dst] = vals
    for dst, src in [("kepler", "alias_kepler"), ("kepler_std", "alias_kepler_std")]:
        vals = np.asarray(dataset(group, src, np.full((n, 7), np.nan)), dtype=np.float64)
        if vals.shape == (n, 7):
            arr[dst] = vals
    return arr


def make_path(group: h5py.Group) -> np.ndarray:
    t = np.asarray(dataset(group, "path_t_rel_s", np.asarray([], dtype=np.float64)), dtype=np.float64)
    pos = np.asarray(dataset(group, "path_position_enu_km", np.empty((0, 3))), dtype=np.float64)
    snr = np.asarray(dataset(group, "path_snr", np.asarray([], dtype=np.float64)), dtype=np.float64)
    dop = np.asarray(dataset(group, "path_doppler_mps", np.asarray([], dtype=np.float64)), dtype=np.float64)
    beam = np.asarray(dataset(group, "path_beam_id", np.asarray([], dtype=np.int64)), dtype=np.int64)
    if len(dop) == 0 and len(t):
        dop = np.full(len(t), np.nan, dtype=np.float64)
    n = min(len(t), len(pos), len(dop), len(snr), len(beam))
    arr = np.zeros(n, dtype=PATH_DTYPE)
    if n:
        arr["t_rel_s"] = t[:n]
        arr["position_enu_km"] = pos[:n]
        arr["doppler_mps"] = dop[:n]
        arr["snr"] = snr[:n]
        arr["beam_id"] = beam[:n]
    return arr


def copy_packed_group(src: h5py.Group, dst: h5py.Group) -> None:
    for key, value in src.attrs.items():
        dst.attrs[key] = value
    if "event" in src and "aliases" in src and "path" in src:
        dst.create_dataset("event", data=coerce_structured(src["event"][()], EVENT_DTYPE))
        dst.create_dataset("aliases", data=coerce_structured(src["aliases"][()], ALIAS_DTYPE))
        dst.create_dataset("path", data=coerce_structured(src["path"][()], PATH_DTYPE))
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
            if name in src:
                dst.create_dataset(name, data=np.asarray(src[name][()], dtype=np.float32))
        return
    dst.create_dataset("event", data=make_event(src))
    for name, default in [
        ("initial_state_gcrs_m_mps", np.full(6, np.nan)),
        ("fit_parameters", np.asarray([], dtype=np.float64)),
        ("fit_parameter_covariance", np.empty((0, 0))),
        ("ceplecha_parameters", np.asarray([], dtype=np.float64)),
        ("ceplecha_parameter_std", np.asarray([], dtype=np.float64)),
        ("ceplecha_parameter_covariance", np.empty((0, 0))),
        ("kepler", np.full(7, np.nan)),
        ("kepler_std", np.full(7, np.nan)),
        ("kepler_covariance", np.full((7, 7), np.nan)),
    ]:
        dst.create_dataset(name, data=np.asarray(dataset(src, name, default), dtype=np.float32))
    dst.create_dataset("aliases", data=make_aliases(src))
    dst.create_dataset("path", data=make_path(src))


def pack_file(path: Path, repack: bool) -> tuple[int, int]:
    tmp = path.with_suffix(path.suffix + ".pack_tmp")
    repacked = path.with_suffix(path.suffix + ".pack_repack_tmp")
    for extra in (tmp, repacked):
        if extra.exists():
            extra.unlink()
    before = path.stat().st_size
    with h5py.File(path, "r") as src, h5py.File(tmp, "w") as dst:
        for key, value in src.attrs.items():
            dst.attrs[key] = value
        for name, obj in src.items():
            if isinstance(obj, h5py.Group):
                copy_packed_group(obj, dst.create_group(name))
    final_tmp = tmp
    if repack and shutil.which("h5repack"):
        subprocess.run(["h5repack", str(tmp), str(repacked)], check=True)
        tmp.unlink()
        final_tmp = repacked
    with h5py.File(final_tmp, "r"):
        pass
    os.replace(final_tmp, path)
    return before, path.stat().st_size


def task(args: tuple[int, int, str, bool]):
    i, n, path, repack = args
    before, after = pack_file(Path(path), repack)
    return i, n, before, after, path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orbit_metadata_dir", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--no-repack", action="store_true")
    args = parser.parse_args()

    files = sorted(args.orbit_metadata_dir.glob("**/orbit@*.h5"))
    if args.limit is not None:
        files = files[: args.limit]
    tasks = [(i, len(files), str(path), not args.no_repack) for i, path in enumerate(files, 1)]
    total_before = sum(Path(path).stat().st_size for _i, _n, path, _r in tasks)
    total_after = 0
    if args.workers <= 1:
        iterator = map(task, tasks)
    else:
        pool = mp.get_context("spawn").Pool(args.workers)
        iterator = pool.imap_unordered(task, tasks, chunksize=4)
    try:
        for i, n, before, after, path in iterator:
            total_after += after
            print(f"pack_orbit {i}/{n} before={before} after={after} {path}", flush=True)
    finally:
        if args.workers > 1:
            pool.close()
            pool.join()
    print(f"pack_orbit_done files={len(files)} before={total_before} after={total_after}", flush=True)


if __name__ == "__main__":
    main()
