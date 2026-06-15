#!/usr/bin/env python3
"""Robust ballistic fits for PANSY simple_meteor_fit events."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import digital_rf as drf
import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as so
from astropy import units as u
from astropy.coordinates import (
    CartesianDifferential,
    CartesianRepresentation,
    EarthLocation,
    GCRS,
    ITRS,
)
from astropy.time import Time

import pansy_config as pc


try:
    from pymsis import msis, utils as msis_utils
except Exception:  # pragma: no cover - handled at runtime with clear message
    msis = None
    msis_utils = None


EARTH_RADIUS_M = 6371.0e3
DRAG_FACTOR = 3.0 / 8.0
ALT_GRID_KM = np.linspace(10.0, 200.0, 300)
PANSY_ALT_M = 0.0

BALLISTIC_WRITER_ARGS = {
    "subdirectory_cadence_seconds": 3600,
    "file_cadence_seconds": 60,
    "samples_per_second_numerator": 1000000,
    "samples_per_second_denominator": 1,
    "file_name": "ballistic",
}


def metadata_writer(path: str | Path) -> drf.DigitalMetadataWriter:
    Path(path).mkdir(parents=True, exist_ok=True)
    return drf.DigitalMetadataWriter(
        str(path),
        BALLISTIC_WRITER_ARGS["subdirectory_cadence_seconds"],
        BALLISTIC_WRITER_ARGS["file_cadence_seconds"],
        BALLISTIC_WRITER_ARGS["samples_per_second_numerator"],
        BALLISTIC_WRITER_ARGS["samples_per_second_denominator"],
        BALLISTIC_WRITER_ARGS["file_name"],
    )


def metadata_file_for_sample(channel_dir: str | Path, sample_idx: int, file_name: str) -> Path:
    sample_sec = int(sample_idx // 1_000_000)
    hour_sec = (sample_sec // 3600) * 3600
    file_sec = (sample_sec // 60) * 60
    subdir = dt.datetime.fromtimestamp(hour_sec, tz=dt.timezone.utc).strftime("%Y-%m-%dT%H-00-00")
    return Path(channel_dir) / subdir / f"{file_name}@{file_sec}.h5"


def delete_metadata_key(channel_dir: str | Path, sample_idx: int, file_name: str) -> bool:
    """Delete one Digital Metadata sample group if it exists.

    Digital RF metadata files used by this project are HDF5 files with a top-level
    group named by the integer sample index. Reanalysis needs this because the
    writer refuses to overwrite an existing sample key.
    """
    path = metadata_file_for_sample(channel_dir, sample_idx, file_name)
    if not path.exists():
        return False
    key = str(int(sample_idx))
    deleted = False
    with h5py.File(path, "a") as handle:
        if key in handle:
            del handle[key]
            deleted = True
        empty = len(handle.keys()) == 0
    if empty:
        path.unlink()
    return deleted


def enu_basis_ecef(lat_deg: float = pc.lat, lon_deg: float = pc.lon) -> np.ndarray:
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    return np.column_stack([east, north, up])


def enu_state_to_gcrs(state_enu: np.ndarray, unix_time_s: float) -> np.ndarray:
    """Convert an ENU position/velocity state at PANSY to GCRS meters and m/s."""
    state_enu = np.asarray(state_enu, dtype=np.float64)
    location = EarthLocation(lat=pc.lat * u.deg, lon=pc.lon * u.deg, height=PANSY_ALT_M * u.m)
    obstime = Time(float(unix_time_s), format="unix", scale="utc")
    origin = location.get_itrs(obstime=obstime).cartesian.xyz.to_value(u.m)
    rot = enu_basis_ecef()
    pos_itrs = origin + rot @ state_enu[:3]
    vel_itrs = rot @ state_enu[3:]
    rep = CartesianRepresentation(
        pos_itrs[0] * u.m,
        pos_itrs[1] * u.m,
        pos_itrs[2] * u.m,
        differentials=CartesianDifferential(
            vel_itrs[0] * u.m / u.s,
            vel_itrs[1] * u.m / u.s,
            vel_itrs[2] * u.m / u.s,
        ),
    )
    gcrs = ITRS(rep, obstime=obstime).transform_to(GCRS(obstime=obstime))
    pos = gcrs.cartesian.without_differentials().xyz.to_value(u.m)
    vel = gcrs.cartesian.differentials["s"].d_xyz.to_value(u.m / u.s)
    return np.concatenate([pos, vel])


def unix_us(value: str) -> int:
    when = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    return int(when.timestamp() * 1_000_000)


def density_interpolator(unix_time_s: float, lat_deg: float = pc.lat, lon_deg: float = pc.lon):
    if msis is None:
        raise RuntimeError("pymsis is required for ballistic fitting. Install requirements.txt.")
    date0 = np.datetime64(int(unix_time_s * 1e9), "ns")
    try:
        data = msis.run([date0], [float(lon_deg)], [float(lat_deg)], ALT_GRID_KM, geomagnetic_activity=1)
    except Exception:
        msis_utils.download_f107_ap()
        data = msis.run([date0], [float(lon_deg)], [float(lat_deg)], ALT_GRID_KM, geomagnetic_activity=1)
    rho = np.asarray(data[0, 0, 0, :, 0], dtype=np.float64)

    def rho_of_alt_m(alt_m):
        alt_km = np.clip(np.asarray(alt_m, dtype=np.float64) / 1e3, ALT_GRID_KM[0], ALT_GRID_KM[-1])
        return np.interp(alt_km, ALT_GRID_KM, rho)

    return rho_of_alt_m, {
        "msis_date_utc": str(date0),
        "msis_lat_deg": float(lat_deg),
        "msis_lon_deg": float(lon_deg),
        "msis_alt_grid_km": ALT_GRID_KM,
        "msis_density_kg_m3": rho,
    }


def initial_guess(t_rel_s: np.ndarray, measured_m: np.ndarray, log10_beta: float = 0.0) -> np.ndarray:
    design = np.column_stack([np.ones_like(t_rel_s), t_rel_s])
    coeff = np.linalg.lstsq(design, measured_m, rcond=None)[0]
    return np.concatenate([coeff[0], coeff[1], [float(log10_beta)]])


def acceleration(state: np.ndarray, rho_of_alt_m, beta_kg_m2: float) -> np.ndarray:
    pos = state[:3]
    vel = state[3:]
    alt_m = max(float(pos[2]) + PANSY_ALT_M, 0.0)
    rho = float(rho_of_alt_m(alt_m))
    speed = float(np.linalg.norm(vel))
    if speed <= 0.0:
        return np.zeros(3, dtype=np.float64)
    drag = -DRAG_FACTOR * rho / max(beta_kg_m2, 1e-12) * speed * vel
    gravity = np.array([0.0, 0.0, -9.81 * (EARTH_RADIUS_M / (EARTH_RADIUS_M + alt_m)) ** 2.0])
    return gravity + drag


def rk4_step(state: np.ndarray, dt_s: float, rho_of_alt_m, beta_kg_m2: float) -> np.ndarray:
    def deriv(y):
        return np.concatenate([y[3:], acceleration(y, rho_of_alt_m, beta_kg_m2)])

    k1 = deriv(state)
    k2 = deriv(state + 0.5 * dt_s * k1)
    k3 = deriv(state + 0.5 * dt_s * k2)
    k4 = deriv(state + dt_s * k3)
    return state + (dt_s / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def propagate(params: np.ndarray, t_rel_s: np.ndarray, rho_of_alt_m, dt_max_s: float = 0.002):
    beta = float(np.clip(10.0 ** params[6], 1e-4, 1e6))
    state = np.concatenate([params[:3], params[3:6]]).astype(np.float64)
    positions = []
    velocities = []
    t_prev = 0.0
    for t in np.asarray(t_rel_s, dtype=np.float64):
        while t_prev + 1e-12 < t:
            dt_step = min(dt_max_s, float(t - t_prev))
            state = rk4_step(state, dt_step, rho_of_alt_m, beta)
            t_prev += dt_step
        positions.append(state[:3].copy())
        velocities.append(state[3:].copy())
    return np.asarray(positions), np.asarray(velocities), beta


def reverse_to_above_atmosphere(
    params: np.ndarray,
    rho_of_alt_m,
    max_alt_m: float = 140e3,
    max_t_s: float = 20.0,
    dt_s: float = 0.002,
):
    beta = float(np.clip(10.0 ** params[6], 1e-4, 1e6))
    state = np.concatenate([params[:3], params[3:6]]).astype(np.float64)
    times = [0.0]
    states = [state.copy()]
    t = 0.0
    while t > -max_t_s and state[2] < max_alt_m:
        state = rk4_step(state, -abs(dt_s), rho_of_alt_m, beta)
        t -= abs(dt_s)
        times.append(t)
        states.append(state.copy())
    states = np.asarray(states)
    times = np.asarray(times)
    return times, states[:, :3], states[:, 3:], beta


def covariance_from_lsq(result, n_residuals: int):
    n_params = int(result.x.size)
    dof = int(n_residuals - n_params)
    cov = np.full((n_params, n_params), np.nan, dtype=np.float64)
    if dof <= 0:
        return cov, np.full(n_params, np.nan), False, dof, np.nan
    try:
        residual_variance = float(2.0 * result.cost / dof)
        cov = np.linalg.pinv(result.jac.T @ result.jac) * residual_variance
        std = np.sqrt(np.maximum(np.diag(cov), 0.0))
        ok = bool(np.all(np.isfinite(std)))
    except np.linalg.LinAlgError:
        residual_variance = np.nan
        std = np.full(n_params, np.nan)
        ok = False
    return cov, std, ok, dof, residual_variance


def fit_ballistic_event(record: dict, sample_idx: int, robust_f_scale: float = 2.0):
    txidx = np.asarray(record["txidx"], dtype=np.float64)
    t_rel_s = txidx / 1e6 - txidx[0] / 1e6
    measured_m = np.column_stack([record["ew"], record["ns"], record["up"]]).astype(np.float64) * 1e3
    snr = np.asarray(record["snr"], dtype=np.float64)
    sigma_m = 80.0 + 600.0 / np.sqrt(np.maximum(snr, 1.0))
    rho_of_alt_m, msis_meta = density_interpolator(sample_idx / 1e6)

    p0 = initial_guess(t_rel_s, measured_m)

    def residual(params, keep=None):
        idx = np.ones(len(t_rel_s), dtype=bool) if keep is None else keep
        pred, _vel, _beta = propagate(params, t_rel_s[idx], rho_of_alt_m)
        return ((pred - measured_m[idx]) / sigma_m[idx, None]).ravel()

    fit1 = so.least_squares(
        residual,
        p0,
        bounds=(
            np.array([-np.inf, -np.inf, 40e3, -90e3, -90e3, -90e3, -4.0]),
            np.array([np.inf, np.inf, 180e3, 90e3, 90e3, 90e3, 6.0]),
        ),
        x_scale=np.array([1e5, 1e5, 1e5, 7e4, 7e4, 7e4, 1.0]),
        loss="soft_l1",
        f_scale=robust_f_scale,
        max_nfev=250,
    )

    pred1, _v1, _b1 = propagate(fit1.x, t_rel_s, rho_of_alt_m)
    norm_res1 = np.linalg.norm((pred1 - measured_m) / sigma_m[:, None], axis=1) / np.sqrt(3.0)
    keep = norm_res1 < 3.0
    if np.count_nonzero(keep) < max(7, len(keep) // 2):
        keep[:] = True

    fit2 = so.least_squares(
        lambda p: residual(p, keep),
        fit1.x,
        bounds=(
            np.array([-np.inf, -np.inf, 40e3, -90e3, -90e3, -90e3, -4.0]),
            np.array([np.inf, np.inf, 180e3, 90e3, 90e3, 90e3, 6.0]),
        ),
        x_scale=np.array([1e5, 1e5, 1e5, 7e4, 7e4, 7e4, 1.0]),
        loss="linear",
        max_nfev=250,
    )

    pred, vel, beta = propagate(fit2.x, t_rel_s, rho_of_alt_m)
    residuals_m = pred - measured_m
    normalized = residuals_m / sigma_m[:, None]
    cov, std, cov_ok, dof, residual_variance = covariance_from_lsq(fit2, len(residual(fit2.x, keep)))
    rev_t, rev_pos, rev_vel, _beta = reverse_to_above_atmosphere(fit2.x, rho_of_alt_m)
    initial_state_enu = fit2.x[:6]
    above_state_enu = np.concatenate([rev_pos[-1], rev_vel[-1]])
    initial_state_gcrs = enu_state_to_gcrs(initial_state_enu, sample_idx / 1e6)
    above_epoch_unix = sample_idx / 1e6 + float(rev_t[-1])
    above_state_gcrs = enu_state_to_gcrs(above_state_enu, above_epoch_unix)

    return {
        "sample_idx": int(sample_idx),
        "t0_unix": float(sample_idx / 1e6),
        "txidx": txidx,
        "t_rel_s": t_rel_s,
        "measured_enu_m": measured_m,
        "model_enu_m": pred,
        "model_velocity_enu_mps": vel,
        "residual_enu_m": residuals_m,
        "normalized_residual_enu": normalized,
        "sigma_m": sigma_m,
        "snr": snr,
        "mfs": np.asarray(record["mfs"], dtype=np.float64),
        "inlier_mask": keep.astype(np.int8),
        "params": fit2.x,
        "param_names": np.asarray(["east_m", "north_m", "up_m", "ve_mps", "vn_mps", "vu_mps", "log10_beta_kg_m2"]),
        "parameter_covariance": cov,
        "parameter_std": std,
        "covariance_available": bool(cov_ok),
        "covariance_degrees_of_freedom": int(dof),
        "covariance_residual_variance": float(residual_variance),
        "ballistic_coefficient_kg_m2": float(beta),
        "reference_frame": np.asarray("GCRS", dtype="S8"),
        "local_fit_frame": np.asarray("PANSY ENU", dtype="S16"),
        "initial_state_enu": initial_state_enu,
        "initial_state_gcrs_m_mps": initial_state_gcrs,
        "initial_speed_mps": float(np.linalg.norm(fit2.x[3:6])),
        "above_atmosphere_dt_s": float(rev_t[-1]),
        "above_atmosphere_epoch_unix": float(above_epoch_unix),
        "above_atmosphere_state_enu": above_state_enu,
        "above_atmosphere_state_gcrs_m_mps": above_state_gcrs,
        "above_atmosphere_speed_mps": float(np.linalg.norm(rev_vel[-1])),
        "above_atmosphere_alt_m": float(rev_pos[-1, 2]),
        "speed_increase_to_above_atmosphere_mps": float(np.linalg.norm(rev_vel[-1]) - np.linalg.norm(fit2.x[3:6])),
        "reverse_time_s": rev_t,
        "reverse_position_enu_m": rev_pos,
        "reverse_velocity_enu_mps": rev_vel,
        "msis_alt_grid_km": msis_meta["msis_alt_grid_km"],
        "msis_density_kg_m3": msis_meta["msis_density_kg_m3"],
        "msis_lat_deg": msis_meta["msis_lat_deg"],
        "msis_lon_deg": msis_meta["msis_lon_deg"],
        "rms_residual_m": float(np.sqrt(np.mean(residuals_m[keep] ** 2.0))),
        "optimizer_success": bool(fit2.success),
        "optimizer_cost": float(fit2.cost),
        "optimizer_nfev": int(fit2.nfev),
    }


def write_fit(writer: drf.DigitalMetadataWriter, fit: dict) -> None:
    payload = {k: v for k, v in fit.items() if k != "param_names"}
    payload["param_names"] = np.asarray(fit["param_names"], dtype="S32")
    writer.write(int(fit["sample_idx"]), payload)


def plot_fit(fit: dict, output_dir: str | Path) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    t = fit["t_rel_s"]
    labels = ["East", "North", "Up"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    inlier = fit["inlier_mask"].astype(bool)
    sigma95 = 1.96 * fit["sigma_m"]
    fit_lo = None
    fit_hi = None
    cov = np.asarray(fit["parameter_covariance"], dtype=np.float64)
    if cov.shape == (7, 7) and np.all(np.isfinite(cov)):
        try:
            rng = np.random.default_rng(0)
            rho_of_alt_m, _meta = density_interpolator(float(fit["t0_unix"]))
            samples = rng.multivariate_normal(np.asarray(fit["params"], dtype=np.float64), cov, size=120)
            models = []
            for sample in samples:
                sample[6] = np.clip(sample[6], -4.0, 6.0)
                model, _vel, _beta = propagate(sample, t, rho_of_alt_m)
                models.append(model)
            models = np.asarray(models)
            fit_lo = np.nanpercentile(models, 2.5, axis=0)
            fit_hi = np.nanpercentile(models, 97.5, axis=0)
        except Exception:
            fit_lo = None
            fit_hi = None
    for comp, ax in enumerate(axes.flat[:3]):
        ax.errorbar(t[inlier], fit["measured_enu_m"][inlier, comp] / 1e3, yerr=sigma95[inlier] / 1e3, fmt=".", color="tab:blue", label="inlier")
        if np.any(~inlier):
            ax.errorbar(t[~inlier], fit["measured_enu_m"][~inlier, comp] / 1e3, yerr=sigma95[~inlier] / 1e3, fmt="x", color="tab:red", label="outlier")
        if fit_lo is not None and fit_hi is not None:
            ax.fill_between(t, fit_lo[:, comp] / 1e3, fit_hi[:, comp] / 1e3, color="tab:green", alpha=0.20, label="95% fit")
        ax.plot(t, fit["model_enu_m"][:, comp] / 1e3, color="black", lw=1.5, label="ballistic fit")
        ax.set_xlabel("Time since first detection (s)")
        ax.set_ylabel(f"{labels[comp]} (km)")
        ax.grid(True, alpha=0.3)
    axes.flat[0].legend(loc="best", fontsize=8)
    axes.flat[3].plot(fit["reverse_position_enu_m"][:, 2] / 1e3, np.linalg.norm(fit["reverse_velocity_enu_mps"], axis=1) / 1e3, color="tab:green")
    axes.flat[3].scatter([fit["initial_state_enu"][2] / 1e3], [fit["initial_speed_mps"] / 1e3], color="tab:blue", label="first detection")
    axes.flat[3].scatter([fit["above_atmosphere_alt_m"] / 1e3], [fit["above_atmosphere_speed_mps"] / 1e3], color="tab:orange", label="above atmosphere")
    axes.flat[3].set_xlabel("Altitude (km)")
    axes.flat[3].set_ylabel("Speed (km/s)")
    axes.flat[3].grid(True, alpha=0.3)
    axes.flat[3].legend(fontsize=8)
    fig.suptitle(
        f"PANSY ballistic fit {dt.datetime.fromtimestamp(fit['t0_unix'], tz=dt.timezone.utc).isoformat()}\n"
        f"beta={fit['ballistic_coefficient_kg_m2']:.3g} kg/m2, "
        f"delta-v_above={fit['speed_increase_to_above_atmosphere_mps'] / 1e3:.2f} km/s, "
        f"rms={fit['rms_residual_m']:.0f} m"
    )
    path = output_dir / f"ballistic_fit_{fit['sample_idx']}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def select_nice_event(simple_dir: str | Path, start_us: int, end_us: int) -> tuple[int, dict]:
    dm = drf.DigitalMetadataReader(str(simple_dir))
    records = dm.read(start_us, end_us)
    if not records:
        raise RuntimeError("No simple_meteor_fit records found in requested interval")
    candidates = []
    for key, rec in records.items():
        n_points = len(rec["txidx"])
        speed = float(np.linalg.norm(rec["v0"]))
        max_snr = float(np.nanmax(rec["snr"]))
        std = np.asarray(rec.get("std", [np.nan, np.nan, np.nan]), dtype=float)
        if n_points >= 12 and 10.0 <= speed <= 75.0:
            score = max_snr + 2.0 * n_points - 2.0 * np.nanmean(std) * 1e3
            candidates.append((score, key, rec))
    if not candidates:
        key = max(records, key=lambda k: len(records[k]["txidx"]))
        return int(key), records[key]
    _score, key, rec = max(candidates, key=lambda item: item[0])
    return int(key), rec


def fit_one(
    simple_dir: str | Path,
    output_dir: str | Path,
    plot_dir: str | Path,
    start: str,
    end: str,
    reanalyze: bool = False,
):
    key, rec = select_nice_event(simple_dir, unix_us(start), unix_us(end))
    fit = fit_ballistic_event(rec, key)
    if reanalyze:
        delete_metadata_key(output_dir, key, BALLISTIC_WRITER_ARGS["file_name"])
    writer = metadata_writer(output_dir)
    write_fit(writer, fit)
    plot_path = plot_fit(fit, plot_dir)
    return fit, plot_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Fit robust PANSY ballistic trajectories from simple_meteor_fit metadata.")
    parser.add_argument("--simple-dir", type=Path, default=Path("data/metadata/simple_meteor_fit"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/metadata/ballistic_fit"))
    parser.add_argument("--plot-dir", type=Path, default=Path("plots/ballistic_fit"))
    parser.add_argument("--start", default="2025-05-06T00:00:00+00:00")
    parser.add_argument("--end", default="2025-05-07T00:00:00+00:00")
    parser.add_argument("--reanalyze", action="store_true", help="Delete an existing ballistic_fit key before writing.")
    args = parser.parse_args()
    fit, plot_path = fit_one(args.simple_dir, args.output_dir, args.plot_dir, args.start, args.end, reanalyze=args.reanalyze)
    print(f"wrote ballistic_fit metadata for {fit['sample_idx']}")
    print(f"plot: {plot_path}")
    print(f"speed increase to above atmosphere: {fit['speed_increase_to_above_atmosphere_mps'] / 1e3:.3f} km/s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
