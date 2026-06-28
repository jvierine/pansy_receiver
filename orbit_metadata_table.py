#!/usr/bin/env python3
"""Compact table storage for PANSY orbit catalogue metadata."""

from __future__ import annotations

import datetime as dt
import fcntl
import os
from pathlib import Path

import h5py
import numpy as np


EVENT_DTYPE = np.dtype(
    [
        ("sample_idx", "<i8"),
        ("candidate_number", "<i8"),
        ("combined_rank", "<i8"),
        ("combined_score", "<f8"),
        ("log10_beta_kg_m2", "<f8"),
        ("sigma_log10_beta", "<f8"),
        ("initial_detection_height_km", "<f8"),
        ("ceplecha_initial_radius_m", "<f8"),
        ("ceplecha_initial_radius_std_m", "<f8"),
        ("ceplecha_initial_mass_kg", "<f8"),
        ("ceplecha_initial_mass_std_kg", "<f8"),
        ("ceplecha_log10_radius_std", "<f8"),
        ("ceplecha_covariance_available", "?"),
        ("ceplecha_reduced_chi2", "<f8"),
        ("ceplecha_bic", "<f8"),
        ("ceplecha_n", "<i8"),
        ("ceplecha_dof", "<i8"),
        ("n_uncertainty_samples", "<i8"),
        ("frac_e_gt_1", "<f8"),
        ("radiant_ra_deg", "<f8"),
        ("radiant_dec_deg", "<f8"),
        ("radiant_speed_km_s", "<f8"),
        ("v_g_km_s", "<f8"),
        ("radiant_ecliptic_lon_deg", "<f8"),
        ("radiant_ecliptic_lat_deg", "<f8"),
        ("radiant_sun_ecliptic_lon_deg", "<f8"),
        ("radiant_sun_ecliptic_lat_deg", "<f8"),
        ("all_aliases_interstellar_nominal", "?"),
        ("n_aliases_orbit_tested", "<i8"),
        ("initial_state_gcrs_m_mps", "<f8", (6,)),
        ("fit_parameters", "<f8", (7,)),
        ("fit_parameter_covariance", "<f8", (7, 7)),
        ("ceplecha_parameters", "<f8", (7,)),
        ("ceplecha_parameter_std", "<f8", (7,)),
        ("ceplecha_parameter_covariance", "<f8", (7, 7)),
        ("kepler", "<f8", (7,)),
        ("kepler_std", "<f8", (7,)),
        ("kepler_covariance", "<f8", (7, 7)),
        ("selected_hypothesis", "S8"),
        ("selection_model_type", "S32"),
        ("orbit_solution_type", "S32"),
    ]
)

ALIAS_DTYPE = np.dtype(
    [
        ("sample_idx", "<i8"),
        ("hypothesis_label", "S8"),
        ("candidate_number", "<i8"),
        ("combined_rank", "<i8"),
        ("combined_score", "<f8"),
        ("selection_model_type", "S32"),
        ("plausibility_model", "S32"),
        ("plausibility_redchi", "<f8"),
        ("tx_beam_snr_weighted_mean_dc", "<f8"),
        ("tx_beam_snr_weighted_rms_dc", "<f8"),
        ("tx_beam_weighted_mean_deg", "<f8"),
        ("tx_beam_weighted_rms_deg", "<f8"),
        ("tx_lobe_snr_weighted_mean_dc", "<f8"),
        ("tx_lobe_snr_weighted_rms_dc", "<f8"),
        ("tx_lobe_p90_dc", "<f8"),
        ("kepler", "<f8", (7,)),
        ("kepler_std", "<f8", (7,)),
        ("frac_e_gt_1", "<f8"),
        ("interstellar_nominal", "?"),
    ]
)

PATH_DTYPE = np.dtype(
    [
        ("sample_idx", "<i8"),
        ("t_rel_s", "<f8"),
        ("position_enu_km", "<f8", (3,)),
        ("snr", "<f8"),
        ("beam_id", "<i8"),
    ]
)


def _bytes(value, dtype: str):
    if isinstance(value, str):
        value = value.encode("utf-8", "ignore")
    arr = np.asarray(value, dtype=dtype)
    return arr.item() if arr.shape == () else arr


def _fixed(values, n: int, fill=np.nan) -> np.ndarray:
    out = np.full(n, fill, dtype=np.float64)
    arr = np.asarray(values, dtype=np.float64).ravel()
    out[: min(n, len(arr))] = arr[:n]
    return out


def _fixed_matrix(values, n: int, fill=np.nan) -> np.ndarray:
    out = np.full((n, n), fill, dtype=np.float64)
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim == 2:
        nr = min(n, arr.shape[0])
        nc = min(n, arr.shape[1])
        out[:nr, :nc] = arr[:nr, :nc]
    return out


def table_path(root: Path, sample_idx: int) -> Path:
    epoch_s = int(sample_idx // 1_000_000)
    hour_s = (epoch_s // 3600) * 3600
    hour = dt.datetime.fromtimestamp(hour_s, tz=dt.timezone.utc)
    subdir = root / hour.strftime("%Y-%m-%dT%H-00-00")
    return subdir / f"orbit@{hour_s}.h5"


def payload_to_rows(sample_idx: int, payload: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    event = np.zeros(1, dtype=EVENT_DTYPE)
    row = event[0]
    for name in EVENT_DTYPE.names:
        if name in {
            "initial_state_gcrs_m_mps",
            "fit_parameters",
            "ceplecha_parameters",
            "ceplecha_parameter_std",
            "kepler",
            "kepler_std",
        }:
            row[name] = _fixed(payload.get(name, []), EVENT_DTYPE[name].shape[0])
        elif name in {"fit_parameter_covariance", "ceplecha_parameter_covariance", "kepler_covariance"}:
            row[name] = _fixed_matrix(payload.get(name, np.empty((0, 0))), EVENT_DTYPE[name].shape[0])
        elif EVENT_DTYPE[name].kind == "S":
            row[name] = _bytes(payload.get(name, b""), EVENT_DTYPE[name].str)
        elif EVENT_DTYPE[name].kind == "b":
            row[name] = bool(payload.get(name, False))
        elif EVENT_DTYPE[name].kind in "iu":
            row[name] = int(payload.get(name, 0))
        else:
            row[name] = float(payload.get(name, np.nan))
    row["sample_idx"] = int(sample_idx)

    labels = np.asarray(payload.get("alias_hypothesis_labels", np.asarray([], dtype="S8")), dtype="S8")
    labels = np.atleast_1d(labels)
    if labels.shape == (1,) and labels[0] == b"":
        labels = labels[:0]
    aliases = np.zeros(len(labels), dtype=ALIAS_DTYPE)
    aliases["sample_idx"] = int(sample_idx)
    aliases["hypothesis_label"] = labels
    alias_map = [
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
    for dst, src, default in alias_map:
        vals = np.atleast_1d(np.asarray(payload.get(src, np.full(len(labels), default))))
        if len(vals) == len(labels):
            aliases[dst] = vals
    for dst, src in [("kepler", "alias_kepler"), ("kepler_std", "alias_kepler_std")]:
        vals = np.asarray(payload.get(src, np.full((len(labels), 7), np.nan)), dtype=np.float64)
        if vals.shape == (len(labels), 7):
            aliases[dst] = vals

    t = np.asarray(payload.get("path_t_rel_s", np.asarray([], dtype=np.float64)), dtype=np.float64)
    pos = np.asarray(payload.get("path_position_enu_km", np.empty((0, 3))), dtype=np.float64)
    snr = np.asarray(payload.get("path_snr", np.asarray([], dtype=np.float64)), dtype=np.float64)
    beam = np.asarray(payload.get("path_beam_id", np.asarray([], dtype=np.int64)), dtype=np.int64)
    n_path = min(len(t), len(pos), len(snr), len(beam))
    paths = np.zeros(n_path, dtype=PATH_DTYPE)
    if n_path:
        paths["sample_idx"] = int(sample_idx)
        paths["t_rel_s"] = t[:n_path]
        paths["position_enu_km"] = pos[:n_path]
        paths["snr"] = snr[:n_path]
        paths["beam_id"] = beam[:n_path]
    return event, aliases, paths


def _read_table(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not path.exists():
        return np.zeros(0, EVENT_DTYPE), np.zeros(0, ALIAS_DTYPE), np.zeros(0, PATH_DTYPE)
    with h5py.File(path, "r") as h:
        if "events" not in h:
            return np.zeros(0, EVENT_DTYPE), np.zeros(0, ALIAS_DTYPE), np.zeros(0, PATH_DTYPE)
        events = h["events"][()]
        aliases = h["aliases"][()] if "aliases" in h else np.zeros(0, ALIAS_DTYPE)
        paths = h["paths"][()] if "paths" in h else np.zeros(0, PATH_DTYPE)
    return events.astype(EVENT_DTYPE, copy=False), aliases.astype(ALIAS_DTYPE, copy=False), paths.astype(PATH_DTYPE, copy=False)


def write_payload(root: Path, sample_idx: int, payload: dict) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    lock_path = root / ".orbit_metadata.lock"
    out = table_path(root, sample_idx)
    out.parent.mkdir(parents=True, exist_ok=True)
    event, aliases, paths = payload_to_rows(sample_idx, payload)
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        old_events, old_aliases, old_paths = _read_table(out)
        keep_events = old_events["sample_idx"] != sample_idx if len(old_events) else np.asarray([], dtype=bool)
        keep_aliases = old_aliases["sample_idx"] != sample_idx if len(old_aliases) else np.asarray([], dtype=bool)
        keep_paths = old_paths["sample_idx"] != sample_idx if len(old_paths) else np.asarray([], dtype=bool)
        new_events = np.concatenate([old_events[keep_events], event])
        new_aliases = np.concatenate([old_aliases[keep_aliases], aliases])
        new_paths = np.concatenate([old_paths[keep_paths], paths])
        order = np.argsort(new_events["sample_idx"])
        new_events = new_events[order]
        tmp = out.with_suffix(out.suffix + ".tmp")
        if tmp.exists():
            tmp.unlink()
        with h5py.File(tmp, "w") as h:
            h.attrs["schema"] = "pansy_orbit_table_v1"
            h.create_dataset("events", data=new_events)
            h.create_dataset("aliases", data=new_aliases)
            h.create_dataset("paths", data=new_paths)
        with h5py.File(tmp, "r"):
            pass
        os.replace(tmp, out)
        fcntl.flock(lock, fcntl.LOCK_UN)
    return out


def has_sample(root: Path, sample_idx: int) -> bool:
    out = table_path(root, sample_idx)
    if not out.exists():
        return False
    try:
        events, _aliases, _paths = _read_table(out)
    except OSError:
        return False
    return bool(len(events) and np.any(events["sample_idx"] == int(sample_idx)))


def delete_sample(root: Path, sample_idx: int) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    lock_path = root / ".orbit_metadata.lock"
    out = table_path(root, sample_idx)
    if not out.exists():
        return out
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        old_events, old_aliases, old_paths = _read_table(out)
        keep_events = old_events["sample_idx"] != sample_idx if len(old_events) else np.asarray([], dtype=bool)
        keep_aliases = old_aliases["sample_idx"] != sample_idx if len(old_aliases) else np.asarray([], dtype=bool)
        keep_paths = old_paths["sample_idx"] != sample_idx if len(old_paths) else np.asarray([], dtype=bool)
        tmp = out.with_suffix(out.suffix + ".tmp")
        if tmp.exists():
            tmp.unlink()
        with h5py.File(tmp, "w") as h:
            h.attrs["schema"] = "pansy_orbit_table_v1"
            h.create_dataset("events", data=old_events[keep_events])
            h.create_dataset("aliases", data=old_aliases[keep_aliases])
            h.create_dataset("paths", data=old_paths[keep_paths])
        with h5py.File(tmp, "r"):
            pass
        os.replace(tmp, out)
        fcntl.flock(lock, fcntl.LOCK_UN)
    return out
