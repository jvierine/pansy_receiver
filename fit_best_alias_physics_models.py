#!/usr/bin/env python3
"""Fit best-alias meteor physics models and write browseable Digital RF metadata."""

from __future__ import annotations

import argparse
import concurrent.futures as futures
import fcntl
import sys
from pathlib import Path

import digital_rf as drf
import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

import pansy_ballistic as pbal

try:
    from meteor_trajectory_models import integrate_ceplecha
except Exception:
    _models_src = Path(__file__).resolve().parent.parent / "meteor_trajectory_models" / "src"
    if _models_src.exists():
        sys.path.insert(0, str(_models_src))
    from meteor_trajectory_models import integrate_ceplecha


WRITER_ARGS = {
    "subdirectory_cadence_seconds": 3600,
    "file_cadence_seconds": 60,
    "samples_per_second_numerator": 1000000,
    "samples_per_second_denominator": 1,
    "file_name": "meteor_physics",
}

FIXED_AM_MIN_M2_KG = 1e-8
FIXED_AM_MAX_M2_KG = 1e4
DEFAULT_SIGMA_POS_KM = 0.5
DEFAULT_SIGMA_DOP_KM_S = 1.0
SHRINKING_RADIUS_MIN_M = 1e-7
SHRINKING_RADIUS_MAX_M = 1e-2
SHRINKING_METEOROID_DENSITY_KG_M3 = 3000.0
SHRINKING_ABLATION_SIGMA_KG_J = 1e-8
SHRINKING_SAMPLE_DT_S = 5e-4
SHRINKING_RADIUS_START_GRID_M = np.array([3e-6, 7e-6, 15e-6, 30e-6, 80e-6, 300e-6, 1e-3], dtype=np.float64)


def metadata_writer(path: Path) -> drf.DigitalMetadataWriter:
    path.mkdir(parents=True, exist_ok=True)
    return drf.DigitalMetadataWriter(
        str(path),
        WRITER_ARGS["subdirectory_cadence_seconds"],
        WRITER_ARGS["file_cadence_seconds"],
        WRITER_ARGS["samples_per_second_numerator"],
        WRITER_ARGS["samples_per_second_denominator"],
        WRITER_ARGS["file_name"],
    )


def write_metadata(path: Path, sample_idx: int, payload: dict, reanalyze: bool = False) -> None:
    if reanalyze:
        pbal.delete_metadata_key(path, sample_idx, WRITER_ARGS["file_name"])
    writer = metadata_writer(path)
    lock_path = path / ".meteor_physics_metadata.lock"
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            writer.write(int(sample_idx), payload)
        except ValueError as exc:
            if "name already exists" not in str(exc):
                raise
            writer = metadata_writer(path)
            writer.write(int(sample_idx), payload)
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def load_best_alias(diagnostics_h5: Path) -> dict:
    with h5py.File(diagnostics_h5, "r") as h:
        sample_idx = int(h.attrs["sample_idx"])
        sample_epoch_unix = float(h.attrs["sample_epoch_unix"])
        label = h.attrs.get("selected_hypothesis", "")
        if not label:
            ranked = []
            for name, grp in h["hypotheses"].items():
                ranked.append((int(grp.attrs.get("combined_rank", 999999)), name))
            label = sorted(ranked)[0][1]
        grp = h["hypotheses"][label]
        if "position_enu_km" in grp:
            points = grp["position_enu_km"][()]
        elif "fit_t_s" in grp and "direction_cosines_uvw" in grp:
            points = grp["direction_cosines_uvw"][()] * grp["range_km"][()][:, None]
        else:
            raise RuntimeError(f"best hypothesis {label} does not contain fitted positions")
        t = grp["t_rel_s"][()]
        doppler_km_s = grp["doppler_mps"][()] / 1e3
        snr = grp["snr"][()] if "snr" in grp else np.full(len(t), np.nan)
        return {
            "sample_idx": sample_idx,
            "sample_epoch_unix": sample_epoch_unix,
            "label": str(label),
            "hypothesis_id": int(grp.attrs.get("hypothesis_id", -1)),
            "candidate_number": int(grp.attrs.get("candidate_number", -1)),
            "combined_rank": int(grp.attrs.get("combined_rank", -1)),
            "combined_score": float(grp.attrs.get("combined_score", np.nan)),
            "t_s": np.asarray(t, dtype=np.float64),
            "points_km": np.asarray(points, dtype=np.float64),
            "doppler_km_s": np.asarray(doppler_km_s, dtype=np.float64),
            "snr": np.asarray(snr, dtype=np.float64),
            "diagnostics_h5": str(diagnostics_h5),
        }


def linear_initial_guess(t_s: np.ndarray, points_km: np.ndarray) -> np.ndarray:
    tau = np.asarray(t_s, dtype=np.float64) - float(np.min(t_s))
    design = np.column_stack([np.ones_like(tau), tau])
    coeff = np.linalg.lstsq(design, np.asarray(points_km) * 1e3, rcond=None)[0]
    return np.concatenate([coeff[0], coeff[1]])


def constant_velocity_model(params: np.ndarray, t_s: np.ndarray):
    tau = np.asarray(t_s, dtype=np.float64) - float(np.min(t_s))
    pos_m = params[:3][None, :] + tau[:, None] * params[3:6][None, :]
    vel_mps = np.repeat(params[3:6][None, :], len(tau), axis=0)
    return pos_m / 1e3, vel_mps / 1e3


def fixed_am_acceleration(state: np.ndarray, rho_of_alt_m, cd_a_over_m_m2_kg: float) -> np.ndarray:
    pos = state[:3]
    vel = state[3:]
    alt_m = max(float(pos[2]) + pbal.PANSY_ALT_M, 0.0)
    rho = float(rho_of_alt_m(alt_m))
    speed = float(np.linalg.norm(vel))
    if speed <= 0.0:
        drag = np.zeros(3, dtype=np.float64)
    else:
        drag = -float(cd_a_over_m_m2_kg) * rho * speed * vel
    gravity = np.array([0.0, 0.0, -9.81 * (pbal.EARTH_RADIUS_M / (pbal.EARTH_RADIUS_M + alt_m)) ** 2.0])
    return gravity + drag


def fixed_am_rk4_step(state: np.ndarray, dt_s: float, rho_of_alt_m, cd_a_over_m_m2_kg: float) -> np.ndarray:
    def deriv(y):
        return np.concatenate([y[3:], fixed_am_acceleration(y, rho_of_alt_m, cd_a_over_m_m2_kg)])

    k1 = deriv(state)
    k2 = deriv(state + 0.5 * dt_s * k1)
    k3 = deriv(state + 0.5 * dt_s * k2)
    k4 = deriv(state + dt_s * k3)
    return state + (dt_s / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def fixed_am_model(params: np.ndarray, t_s: np.ndarray, rho_of_alt_m, dt_max_s: float = 0.002):
    tau = np.asarray(t_s, dtype=np.float64) - float(np.min(t_s))
    order = np.argsort(tau)
    cd_a_over_m = float(np.clip(10.0 ** params[6], FIXED_AM_MIN_M2_KG, FIXED_AM_MAX_M2_KG))
    state = np.concatenate([params[:3], params[3:6]]).astype(np.float64)
    positions = np.empty((len(tau), 3), dtype=np.float64)
    velocities = np.empty((len(tau), 3), dtype=np.float64)
    t_prev = 0.0
    for idx in order:
        target = float(tau[idx])
        while t_prev + 1e-12 < target:
            dt_step = min(dt_max_s, target - t_prev)
            state = fixed_am_rk4_step(state, dt_step, rho_of_alt_m, cd_a_over_m)
            t_prev += dt_step
        positions[idx] = state[:3]
        velocities[idx] = state[3:6]
    return positions / 1e3, velocities / 1e3, cd_a_over_m


def propagate_shrinking_radius_model(params: np.ndarray, t_s: np.ndarray, rho_of_alt_m):
    tau = np.asarray(t_s, dtype=np.float64) - float(np.min(t_s))
    order = np.argsort(tau)
    tau_sorted = tau[order]
    t1 = max(float(np.max(tau_sorted)) if len(tau_sorted) else 0.0, SHRINKING_SAMPLE_DT_S)
    radius0_m = float(np.clip(10.0 ** params[6], SHRINKING_RADIUS_MIN_M, SHRINKING_RADIUS_MAX_M))

    result = integrate_ceplecha(
        params[:3],
        params[3:6],
        radius0_m,
        rho_of_alt_m,
        meteoroid_density_kg_m3=SHRINKING_METEOROID_DENSITY_KG_M3,
        ablation_sigma_kg_j=SHRINKING_ABLATION_SIGMA_KG_J,
        t_span_s=(0.0, t1),
        sample_dt_s=min(SHRINKING_SAMPLE_DT_S, max(t1 / 5.0, 1e-6)),
        height_function=lambda r: max(float(r[2]) + pbal.PANSY_ALT_M, 0.0),
    )
    if result.time_s.size < 2:
        raise RuntimeError(f"shrinking-radius integration returned too few samples: {result.message}")
    pos_sorted = np.column_stack([np.interp(tau_sorted, result.time_s, result.position_m[:, dim]) for dim in range(3)])
    vel_sorted = np.column_stack([np.interp(tau_sorted, result.time_s, result.velocity_mps[:, dim]) for dim in range(3)])
    radius_sorted = np.interp(tau_sorted, result.time_s, result.radius_m)
    mass_sorted = np.interp(tau_sorted, result.time_s, result.mass_kg)
    pos = np.empty_like(pos_sorted)
    vel = np.empty_like(vel_sorted)
    radius = np.empty_like(radius_sorted)
    mass = np.empty_like(mass_sorted)
    pos[order] = pos_sorted
    vel[order] = vel_sorted
    radius[order] = radius_sorted
    mass[order] = mass_sorted
    return pos / 1e3, vel / 1e3, radius, mass, bool(result.success), str(result.message)


def predicted_doppler(model_km: np.ndarray, velocity_km_s: np.ndarray) -> np.ndarray:
    rng = np.linalg.norm(model_km, axis=1)
    return np.sum(model_km * velocity_km_s, axis=1) / np.maximum(rng, 1e-9)


def fit_model(name, t_s, points_km, doppler_km_s, rho_of_alt_m, sigma_pos_km, sigma_dop_km_s, p0_list):
    t_s = np.asarray(t_s, dtype=np.float64)
    points_km = np.asarray(points_km, dtype=np.float64)
    doppler_km_s = np.asarray(doppler_km_s, dtype=np.float64)
    finite = np.isfinite(t_s) & np.all(np.isfinite(points_km), axis=1) & np.isfinite(doppler_km_s)

    if name == "fixed_velocity":
        n_params = 6
        bounds = (
            np.array([-np.inf, -np.inf, 20e3, -90e3, -90e3, -90e3]),
            np.array([np.inf, np.inf, 220e3, 90e3, 90e3, 90e3]),
        )
        x_scale = np.array([1e5, 1e5, 1e5, 7e4, 7e4, 7e4])

        def model(x):
            pos, vel = constant_velocity_model(x, t_s)
            return pos, vel, {}

    elif name == "fixed_am":
        n_params = 7
        bounds = (
            np.array([-np.inf, -np.inf, 20e3, -90e3, -90e3, -90e3, np.log10(FIXED_AM_MIN_M2_KG)]),
            np.array([np.inf, np.inf, 220e3, 90e3, 90e3, 90e3, np.log10(FIXED_AM_MAX_M2_KG)]),
        )
        x_scale = np.array([1e5, 1e5, 1e5, 7e4, 7e4, 7e4, 1.0])

        def model(x):
            pos, vel, cd_a_over_m = fixed_am_model(x, t_s, rho_of_alt_m)
            return pos, vel, {"cd_a_over_m_m2_kg": cd_a_over_m}

    elif name == "shrinking_radius":
        n_params = 7
        bounds = (
            np.array([-np.inf, -np.inf, 20e3, -90e3, -90e3, -90e3, np.log10(SHRINKING_RADIUS_MIN_M)]),
            np.array([np.inf, np.inf, 220e3, 90e3, 90e3, 90e3, np.log10(SHRINKING_RADIUS_MAX_M)]),
        )
        x_scale = np.array([1e5, 1e5, 1e5, 7e4, 7e4, 7e4, 1.0])

        def model(x):
            pos, vel, radius, mass, success, message = propagate_shrinking_radius_model(x, t_s, rho_of_alt_m)
            return pos, vel, {
                "radius_m": radius,
                "mass_kg": mass,
                "integrator_success": success,
                "integrator_message": message,
            }

    else:
        raise ValueError(name)

    def residual(x, keep):
        pos, vel, _extra = model(x)
        pred_dop = predicted_doppler(pos, vel)
        pos_res = ((points_km[keep] - pos[keep]) / sigma_pos_km).ravel()
        dop_res = (doppler_km_s[keep] - pred_dop[keep]) / sigma_dop_km_s
        return np.concatenate([pos_res, dop_res])

    best = None
    failures = 0
    for p0 in p0_list:
        try:
            result0 = opt.least_squares(
                lambda x: residual(x, finite),
                np.asarray(p0, dtype=np.float64),
                bounds=bounds,
                x_scale=x_scale,
                loss="soft_l1",
                f_scale=1.0,
                max_nfev=300,
            )
        except Exception:
            failures += 1
            continue
        if best is None or result0.cost < best.cost:
            best = result0
    if best is None:
        raise RuntimeError(f"{name}: all starts failed")

    pos, vel, extra = model(best.x)
    pred_dop = predicted_doppler(pos, vel)
    norm_per_point = np.sqrt(
        np.sum(((points_km - pos) / sigma_pos_km) ** 2, axis=1)
        + ((doppler_km_s - pred_dop) / sigma_dop_km_s) ** 2
    )
    keep = finite & (norm_per_point < 3.5)
    if np.count_nonzero(keep) < max(n_params + 2, 10):
        keep = finite.copy()

    result = opt.least_squares(
        lambda x: residual(x, keep),
        best.x,
        bounds=bounds,
        x_scale=x_scale,
        loss="linear",
        max_nfev=300,
    )
    pos, vel, extra = model(result.x)
    pred_dop = predicted_doppler(pos, vel)
    pos_res = points_km - pos
    dop_res = doppler_km_s - pred_dop
    raw = np.concatenate([(pos_res[keep] / sigma_pos_km).ravel(), ((dop_res[keep] / sigma_dop_km_s).ravel())])
    chi2 = float(np.sum(raw**2.0))
    n_res = int(len(raw))
    dof = max(1, n_res - n_params)
    cov, std, cov_ok, cov_dof, cov_res_var = pbal.covariance_from_lsq(result, n_res)
    return {
        "name": name,
        "params": result.x,
        "parameter_std": std,
        "parameter_covariance": cov,
        "covariance_available": bool(cov_ok),
        "covariance_degrees_of_freedom": int(cov_dof),
        "covariance_residual_variance": float(cov_res_var),
        "model_enu_km": pos,
        "velocity_km_s": vel,
        "pred_doppler_km_s": pred_dop,
        "pos_res_km": pos_res,
        "dop_res_km_s": dop_res,
        "keep_mask": keep.astype(np.int8),
        "n_points": int(np.count_nonzero(keep)),
        "n_params": int(n_params),
        "chi2": chi2,
        "dof": int(dof),
        "reduced_chi2": chi2 / dof,
        "aic": chi2 + 2.0 * n_params,
        "bic": chi2 + n_params * np.log(max(1, n_res)),
        "pos_rms_km": float(np.sqrt(np.mean(np.sum(pos_res[keep] ** 2.0, axis=1)))),
        "dop_rms_km_s": float(np.sqrt(np.mean(dop_res[keep] ** 2.0))),
        "multistart_failures": int(failures),
        **extra,
    }


def load_dasst_best(dasst_h5: Path | None, label: str) -> dict:
    out = {
        "dasst_available": False,
        "dasst_kepler": np.full(7, np.nan),
        "dasst_kepler_std": np.full(7, np.nan),
        "dasst_kepler_covariance": np.full((7, 7), np.nan),
        "dasst_frac_e_gt_1": np.nan,
    }
    if dasst_h5 is None or not dasst_h5.exists():
        return out
    with h5py.File(dasst_h5, "r") as h:
        key = label if label in h else None
        if key is None:
            ranked = [(int(g.attrs.get("combined_rank", 999999)), name) for name, g in h.items()]
            if not ranked:
                return out
            key = sorted(ranked)[0][1]
        grp = h[key]
        out["dasst_available"] = True
        out["dasst_hypothesis"] = str(key)
        out["dasst_kepler"] = grp["kepler"][()] if "kepler" in grp else out["dasst_kepler"]
        out["dasst_kepler_std"] = grp["kepler_std"][()] if "kepler_std" in grp else out["dasst_kepler_std"]
        out["dasst_kepler_covariance"] = grp["kepler_covariance"][()] if "kepler_covariance" in grp else out["dasst_kepler_covariance"]
        out["dasst_frac_e_gt_1"] = float(grp.attrs.get("frac_e_gt_1", np.nan))
    return out


def add_fit_payload(payload: dict, prefix: str, fit: dict):
    params = np.asarray(fit["params"], dtype=np.float64)
    payload[f"{prefix}_params"] = params
    payload[f"{prefix}_parameter_std"] = np.asarray(fit["parameter_std"], dtype=np.float64)
    payload[f"{prefix}_parameter_covariance"] = np.asarray(fit["parameter_covariance"], dtype=np.float64)
    payload[f"{prefix}_covariance_available"] = bool(fit["covariance_available"])
    payload[f"{prefix}_n_points"] = int(fit["n_points"])
    payload[f"{prefix}_n_params"] = int(fit["n_params"])
    payload[f"{prefix}_chi2"] = float(fit["chi2"])
    payload[f"{prefix}_dof"] = int(fit["dof"])
    payload[f"{prefix}_reduced_chi2"] = float(fit["reduced_chi2"])
    payload[f"{prefix}_aic"] = float(fit["aic"])
    payload[f"{prefix}_bic"] = float(fit["bic"])
    payload[f"{prefix}_pos_rms_km"] = float(fit["pos_rms_km"])
    payload[f"{prefix}_dop_rms_km_s"] = float(fit["dop_rms_km_s"])
    payload[f"{prefix}_initial_height_km"] = float(params[2] / 1e3)
    payload[f"{prefix}_initial_velocity_mps"] = np.asarray(params[3:6], dtype=np.float64)
    payload[f"{prefix}_initial_speed_km_s"] = float(np.linalg.norm(params[3:6]) / 1e3)
    payload[f"{prefix}_keep_mask"] = np.asarray(fit["keep_mask"], dtype=np.int8)
    payload[f"{prefix}_model_enu_km"] = np.asarray(fit["model_enu_km"], dtype=np.float64)
    payload[f"{prefix}_velocity_km_s"] = np.asarray(fit["velocity_km_s"], dtype=np.float64)
    payload[f"{prefix}_pred_doppler_km_s"] = np.asarray(fit["pred_doppler_km_s"], dtype=np.float64)
    payload[f"{prefix}_pos_res_km"] = np.asarray(fit["pos_res_km"], dtype=np.float64)
    payload[f"{prefix}_dop_res_km_s"] = np.asarray(fit["dop_res_km_s"], dtype=np.float64)

    if prefix == "fixed_am":
        cd_a_over_m = float(fit["cd_a_over_m_m2_kg"])
        payload["fixed_am_cd_a_over_m_m2_kg"] = cd_a_over_m
        payload["fixed_am_log10_cd_a_over_m_std"] = float(fit["parameter_std"][6])
    if prefix == "shrinking_radius":
        radius = np.asarray(fit["radius_m"], dtype=np.float64)
        mass = np.asarray(fit["mass_kg"], dtype=np.float64)
        payload["shrinking_radius_initial_radius_m"] = float(radius[0])
        payload["shrinking_radius_final_radius_m"] = float(radius[-1])
        payload["shrinking_radius_initial_mass_kg"] = float(mass[0])
        payload["shrinking_radius_final_mass_kg"] = float(mass[-1])
        payload["shrinking_radius_log10_radius_std"] = float(fit["parameter_std"][6])
        payload["shrinking_radius_radius_m"] = radius
        payload["shrinking_radius_mass_kg"] = mass
        payload["shrinking_radius_integrator_success"] = bool(fit["integrator_success"])
        payload["shrinking_radius_integrator_message"] = str(fit["integrator_message"])


def best_model_name(fits: dict[str, dict]) -> str:
    return min(fits, key=lambda name: float(fits[name]["bic"]))


def model_label(name: str) -> str:
    return {
        "fixed_velocity": "fixed velocity",
        "fixed_am": r"fixed $C_dA/m$",
        "shrinking_radius": "shrinking radius",
    }[name]


def finite_sigma(value, scale=1.0):
    value = float(value)
    return value * scale if np.isfinite(value) else np.nan


def plot_event_fit(output: Path, obs: dict, fits: dict[str, dict], payload: dict) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    best_name = str(payload["best_physics_model"])
    colors = {
        "fixed_velocity": "tab:blue",
        "fixed_am": "tab:green",
        "shrinking_radius": "tab:orange",
    }
    names = ["fixed_velocity", "fixed_am", "shrinking_radius"]
    pts = np.asarray(obs["points_km"], dtype=np.float64)
    t = np.asarray(obs["t_s"], dtype=np.float64)
    order = np.argsort(t)
    snr = np.asarray(obs["snr"], dtype=np.float64)
    snr_color = np.log10(np.maximum(snr, 1e-12))
    best_fit = fits[best_name]
    best_model = np.asarray(best_fit["model_enu_km"], dtype=np.float64)
    best_doppler = np.asarray(best_fit["pred_doppler_km_s"], dtype=np.float64)
    doppler = np.asarray(obs["doppler_km_s"], dtype=np.float64)
    residual = doppler - best_doppler
    tau = t - float(np.nanmin(t))
    along_direction = np.asarray(best_fit["params"][3:6], dtype=np.float64)
    along_norm = np.linalg.norm(along_direction)
    if along_norm > 0.0:
        along_unit = along_direction / along_norm
    else:
        along_unit = np.array([1.0, 0.0, 0.0])
    along_obs = pts @ along_unit / 1e3
    along_model = best_model @ along_unit / 1e3

    fig, axes = plt.subplots(3, 3, figsize=(16.0, 14.0), constrained_layout=True)

    ax = axes[0, 0]
    sc = ax.scatter(pts[:, 0], pts[:, 1], s=16, c=snr_color, cmap="viridis", edgecolors="none", label="measurements")
    ax.plot(best_model[order, 0], best_model[order, 1], color="black", lw=2.1, label=model_label(best_name))
    ax.set_xlabel("East (km)")
    ax.set_ylabel("North (km)")
    ax.set_title("1a. Best-alias horizontal trajectory")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, loc="best")

    ax = axes[0, 1]
    ax.scatter(along_obs[order], pts[order, 2], s=16, c=snr_color[order], cmap="viridis", edgecolors="none")
    ax.plot(along_model[order], best_model[order, 2], color="black", lw=2.1)
    ax.set_xlabel("Along-track coordinate (km)")
    ax.set_ylabel("Up (km)")
    ax.set_title("1b. Vertical trajectory")
    ax.grid(True, alpha=0.25)

    ax = axes[0, 2]
    ax.scatter(tau[order], residual[order] * 1e3, s=16, c=snr_color[order], cmap="viridis", edgecolors="none")
    ax.axhline(0.0, color="black", lw=1.0, alpha=0.8)
    ax.set_xlabel("Time since first echo (s)")
    ax.set_ylabel("Doppler residual (m/s)")
    ax.set_title(f"1c. Best-model residuals; RMS={best_fit['doppler_rms_km_s'] * 1e3:.1f} m/s")
    ax.grid(True, alpha=0.25)

    ax = axes[1, 0]
    labels = ["East", "North", "Up"]
    for dim, label in enumerate(labels):
        ax.plot(tau[order], pts[order, dim] - best_model[order, dim], ".", ms=4.0, label=label)
    ax.axhline(0.0, color="black", lw=1.0, alpha=0.65)
    ax.set_xlabel("Time since first echo (s)")
    ax.set_ylabel("Position residual (km)")
    ax.set_title(f"2a. Position residuals; RMS={best_fit['position_rms_km']:.2f} km")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, ncols=3, loc="best")

    ax = axes[1, 1]
    ax.axis("off")
    p = np.asarray(best_fit["params"], dtype=np.float64)
    ps = np.asarray(best_fit["parameter_std"], dtype=np.float64)
    fixed_am = fits["fixed_am"]
    shrink = fits["shrinking_radius"]
    lines = [
        f"Sample                 {obs['sample_idx']}",
        f"Best alias             {obs['label']}  candidate {obs['candidate_number']}",
        f"Best physics model     {model_label(best_name)}",
        "",
        "Best-model state at first echo",
        f"  h0                   {p[2] / 1e3:8.2f} +/- {finite_sigma(ps[2], 1e-3):.2f} km",
        f"  |v0|                 {np.linalg.norm(p[3:6]) / 1e3:8.2f} km/s",
        f"  vE,vN,vU             {p[3]/1e3:8.2f} {p[4]/1e3:8.2f} {p[5]/1e3:8.2f} km/s",
        "",
        "Secondary physical parameters",
        f"  fixed CdA/m          {fixed_am['cd_a_over_m_m2_kg']:8.3g} m2/kg",
        f"  sigma log10(CdA/m)   {fixed_am['parameter_std'][6]:8.3g}",
        f"  shrink r0            {shrink['radius_m'][0] * 1e6:8.2f} um",
        f"  shrink m0            {shrink['mass_kg'][0]:8.3g} kg",
        f"  sigma log10(r0)      {shrink['parameter_std'][6]:8.3g}",
    ]
    ax.text(0.02, 0.98, "\n".join(lines), ha="left", va="top", transform=ax.transAxes, family="monospace", fontsize=9.2)
    ax.set_title("2b. Fit parameters and 1-sigma uncertainties")

    ax = axes[1, 2]
    bic = np.asarray([fits[name]["bic"] for name in names], dtype=np.float64)
    delta = bic - np.nanmin(bic)
    y = np.arange(len(names))
    bars = ax.barh(y, delta, color=[colors[name] for name in names], alpha=0.85)
    for bar, name, value in zip(bars, names, delta):
        if name == best_name:
            bar.set_edgecolor("black")
            bar.set_linewidth(1.8)
        ax.text(max(value, 0.02), bar.get_y() + bar.get_height() / 2.0, f" {value:.1f}", va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([model_label(name) for name in names])
    ax.set_xlabel(r"$\Delta$BIC relative to best")
    ax.set_title("2c. Model support")
    ax.grid(True, axis="x", alpha=0.25)
    ax.invert_yaxis()

    ax = axes[2, 0]
    ax.scatter(tau[order], doppler[order], s=17, c="black", alpha=0.55, label="measured")
    for name in names:
        fit = fits[name]
        ax.plot(
            tau[order],
            fit["pred_doppler_km_s"][order],
            color="black" if name == best_name else colors[name],
            lw=2.4 if name == best_name else 1.0,
            alpha=0.98 if name == best_name else 0.35,
            label=f"{model_label(name)}  chi2nu={fit['reduced_chi2']:.2f}",
        )
    ax.set_xlabel("Time since first echo (s)")
    ax.set_ylabel("Doppler/range rate (km/s)")
    ax.set_title(f"3a. Doppler fit: BIC-best {model_label(best_name)}")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7.5, loc="best")

    ax = axes[2, 1]
    rows = []
    for name in names:
        fit = fits[name]
        rows.append(
            [
                model_label(name),
                f"{fit['chi2']:.1f}",
                f"{fit['reduced_chi2']:.2f}",
                f"{fit['position_rms_km']:.2f}",
                f"{fit['doppler_rms_km_s'] * 1e3:.1f}",
                f"{fit['bic']:.1f}",
            ]
        )
    ax.axis("off")
    table = ax.table(
        cellText=rows,
        colLabels=["model", "chi2", "chi2nu", "pos km", "dop m/s", "BIC"],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.2)
    table.scale(1.0, 1.45)
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor("0.9")
        elif rows[row - 1][0] == model_label(best_name):
            cell.set_edgecolor("black")
            cell.set_linewidth(1.4)
    ax.set_title("3b. Fit quality")

    ax = axes[2, 2]
    ax.axis("off")
    orbit_lines = ["DASST zenithal-attraction corrected"]
    if bool(payload.get("dasst_available", False)):
        kep = np.asarray(payload["dasst_kepler"], dtype=np.float64)
        kstd = np.asarray(payload["dasst_kepler_std"], dtype=np.float64)
        orbit_lines.extend(
            [
                f"a = {kep[0]:.3g} +/- {kstd[0]:.2g} AU",
                f"e = {kep[1]:.3f} +/- {kstd[1]:.2g}",
                f"i = {kep[2]:.2f} +/- {kstd[2]:.2g} deg",
                f"q = {kep[6]:.3g} +/- {kstd[6]:.2g} AU",
                f"P(e>1) = {float(payload.get('dasst_frac_e_gt_1', np.nan)):.3f}",
            ]
        )
    else:
        orbit_lines.append("No matching DASST result found")
    ax.text(0.03, 0.96, "\n".join(orbit_lines), ha="left", va="top", transform=ax.transAxes, family="monospace", fontsize=10)
    ax.set_title("3c. Orbit metadata")

    fig.suptitle(f"PANSY best-alias physics fits: {obs['sample_idx']}")
    fig.savefig(output, dpi=200)
    plt.close(fig)


def fit_event(
    diagnostics_h5: Path,
    dasst_h5: Path | None,
    metadata_dir: Path,
    plot_dir: Path | None,
    sigma_pos_km: float,
    sigma_dop_km_s: float,
    reanalyze: bool = False,
) -> dict:
    obs = load_best_alias(diagnostics_h5)
    rho_of_alt_m, msis_meta = pbal.density_interpolator(obs["sample_epoch_unix"])
    p0_linear = linear_initial_guess(obs["t_s"], obs["points_km"])

    fixed_velocity = fit_model(
        "fixed_velocity",
        obs["t_s"],
        obs["points_km"],
        obs["doppler_km_s"],
        rho_of_alt_m,
        sigma_pos_km,
        sigma_dop_km_s,
        [p0_linear],
    )
    fixed_am_starts = [np.concatenate([fixed_velocity["params"][:6], [log_b]]) for log_b in (-4.0, -2.0, 0.0, 1.0, 2.0)]
    fixed_am = fit_model(
        "fixed_am",
        obs["t_s"],
        obs["points_km"],
        obs["doppler_km_s"],
        rho_of_alt_m,
        sigma_pos_km,
        sigma_dop_km_s,
        fixed_am_starts,
    )
    fixed_am_radius_guess = 3.0 / (4.0 * SHRINKING_METEOROID_DENSITY_KG_M3 * max(float(fixed_am["cd_a_over_m_m2_kg"]), 1e-30))
    radius_starts = sorted(set(float(np.clip(r, SHRINKING_RADIUS_MIN_M, SHRINKING_RADIUS_MAX_M)) for r in [fixed_am_radius_guess, *SHRINKING_RADIUS_START_GRID_M]))
    shrinking_starts = [np.concatenate([fixed_am["params"][:6], [np.log10(radius)]]) for radius in radius_starts]
    shrinking_radius = fit_model(
        "shrinking_radius",
        obs["t_s"],
        obs["points_km"],
        obs["doppler_km_s"],
        rho_of_alt_m,
        sigma_pos_km,
        sigma_dop_km_s,
        shrinking_starts,
    )
    fits = {
        "fixed_velocity": fixed_velocity,
        "fixed_am": fixed_am,
        "shrinking_radius": shrinking_radius,
    }
    best_name = best_model_name(fits)

    payload = {
        "schema_version": "pansy_best_alias_meteor_physics_v1",
        "source_program": "fit_best_alias_physics_models.py",
        "sample_idx": int(obs["sample_idx"]),
        "sample_epoch_unix": float(obs["sample_epoch_unix"]),
        "diagnostics_h5": obs["diagnostics_h5"],
        "best_hypothesis": obs["label"],
        "hypothesis_id": int(obs["hypothesis_id"]),
        "candidate_number": int(obs["candidate_number"]),
        "combined_rank": int(obs["combined_rank"]),
        "combined_score": float(obs["combined_score"]),
        "t_s": obs["t_s"],
        "observed_position_enu_km": obs["points_km"],
        "observed_doppler_km_s": obs["doppler_km_s"],
        "observed_snr": obs["snr"],
        "best_physics_model": best_name,
        "best_physics_model_criterion": "minimum_bic",
        "sigma_pos_km": float(sigma_pos_km),
        "sigma_dop_km_s": float(sigma_dop_km_s),
        "meteoroid_density_kg_m3": float(SHRINKING_METEOROID_DENSITY_KG_M3),
        "shrinking_ablation_sigma_kg_j": float(SHRINKING_ABLATION_SIGMA_KG_J),
        "msis_alt_grid_km": msis_meta["msis_alt_grid_km"],
        "msis_density_kg_m3": msis_meta["msis_density_kg_m3"],
    }
    add_fit_payload(payload, "fixed_velocity", fixed_velocity)
    add_fit_payload(payload, "fixed_am", fixed_am)
    add_fit_payload(payload, "shrinking_radius", shrinking_radius)
    dasst = load_dasst_best(dasst_h5, obs["label"])
    payload.update(dasst)
    payload["dasst_kepler_names"] = np.asarray(["a_au", "e", "i_deg", "raan_deg", "argp_deg", "nu_deg", "q_au"], dtype="S16")

    write_metadata(metadata_dir, obs["sample_idx"], payload, reanalyze=reanalyze)
    if plot_dir is not None:
        plot_event_fit(plot_dir / f"meteor_physics_fit_{obs['sample_idx']}.png", obs, fits, payload)
    return {
        "sample_idx": int(obs["sample_idx"]),
        "best_model": best_name,
        "fixed_velocity_bic": float(fixed_velocity["bic"]),
        "fixed_am_bic": float(fixed_am["bic"]),
        "shrinking_radius_bic": float(shrinking_radius["bic"]),
        "metadata_dir": str(metadata_dir),
    }


def dasst_path_for(diagnostics_h5: Path, dasst_dir: Path | None) -> Path | None:
    if dasst_dir is None:
        return None
    sample = diagnostics_h5.stem.split("_")[-1]
    path = dasst_dir / f"pansy_candidate_orbits_dasst_{sample}.h5"
    return path if path.exists() else None


def run_one_from_args(args) -> dict:
    return fit_event(
        args.diagnostics_h5,
        args.dasst_h5,
        args.metadata_dir,
        args.plot_dir,
        args.sigma_pos_km,
        args.sigma_dop_km_s,
        reanalyze=args.reanalyze,
    )


def run_batch_item(item):
    diagnostics_h5, dasst_dir, metadata_dir, plot_dir, sigma_pos_km, sigma_dop_km_s, reanalyze = item
    return fit_event(
        diagnostics_h5,
        dasst_path_for(diagnostics_h5, dasst_dir),
        metadata_dir,
        plot_dir,
        sigma_pos_km,
        sigma_dop_km_s,
        reanalyze=reanalyze,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Fit best-alias meteor physics models and write Digital RF metadata.")
    parser.add_argument("diagnostics_h5", type=Path, nargs="?")
    parser.add_argument("--dasst-h5", type=Path, default=None)
    parser.add_argument("--dasst-dir", type=Path, default=Path("test_plots"))
    parser.add_argument("--metadata-dir", type=Path, default=Path("data/metadata/meteor_physics"))
    parser.add_argument("--plot-dir", type=Path, default=None)
    parser.add_argument("--sigma-pos-km", type=float, default=DEFAULT_SIGMA_POS_KM)
    parser.add_argument("--sigma-dop-km-s", type=float, default=DEFAULT_SIGMA_DOP_KM_S)
    parser.add_argument("--reanalyze", action="store_true")
    parser.add_argument("--batch-glob", type=str, default=None)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--max-events", type=int, default=None)
    args = parser.parse_args()

    if args.batch_glob is not None:
        paths = sorted(Path().glob(args.batch_glob))
        if args.max_events is not None:
            paths = paths[: args.max_events]
        if not paths:
            raise RuntimeError(f"no diagnostics matched {args.batch_glob}")
        args.metadata_dir.mkdir(parents=True, exist_ok=True)
        if args.plot_dir is not None:
            args.plot_dir.mkdir(parents=True, exist_ok=True)
        work = [
            (path, args.dasst_dir, args.metadata_dir, args.plot_dir, args.sigma_pos_km, args.sigma_dop_km_s, args.reanalyze)
            for path in paths
        ]
        ok = 0
        failed = 0
        if args.jobs <= 1:
            iterator = map(run_batch_item, work)
            for result in iterator:
                ok += 1
                print(
                    f"batch_ok {result['sample_idx']} best={result['best_model']} "
                    f"bic_v={result['fixed_velocity_bic']:.1f} bic_am={result['fixed_am_bic']:.1f} "
                    f"bic_shrink={result['shrinking_radius_bic']:.1f}",
                    flush=True,
                )
        else:
            with futures.ProcessPoolExecutor(max_workers=args.jobs) as ex:
                futs = {ex.submit(run_batch_item, item): item[0] for item in work}
                for fut in futures.as_completed(futs):
                    path = futs[fut]
                    try:
                        result = fut.result()
                    except Exception as exc:
                        failed += 1
                        print(f"batch_failed {path} {type(exc).__name__}: {exc}", flush=True)
                        continue
                    ok += 1
                    print(
                        f"batch_ok {result['sample_idx']} best={result['best_model']} "
                        f"bic_v={result['fixed_velocity_bic']:.1f} bic_am={result['fixed_am_bic']:.1f} "
                        f"bic_shrink={result['shrinking_radius_bic']:.1f}",
                        flush=True,
                    )
        print(f"batch_complete ok={ok} failed={failed} total={len(paths)}")
        return 0 if failed == 0 else 1

    if args.diagnostics_h5 is None:
        raise RuntimeError("provide diagnostics_h5 or --batch-glob")

    result = run_one_from_args(args)

    for prefix in ["fixed_velocity", "fixed_am", "shrinking_radius"]:
        print(
            prefix,
            f"bic={result[prefix + '_bic']:.3f}",
        )
    print(f"best_model {result['best_model']}")
    print(f"metadata_write {args.metadata_dir} {result['sample_idx']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
