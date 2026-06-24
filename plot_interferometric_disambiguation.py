#!/usr/bin/env python3
"""Make memo figures for PANSY interferometric alias disambiguation."""

from __future__ import annotations

import argparse
import itertools
import subprocess
import sys
import time
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py
import numpy as np
import scipy.constants as sc
import scipy.ndimage as ndi
import scipy.optimize as opt
import scipy.stats as st

import pansy_interferometry as pint
import pansy_config as pc
import pansy_modes as pmm
import pansy_ballistic as pbal
import pansy_orbit as porb
import tx_phase_quality as txpq
from interferometer_alias_diagnostics import RangeDopplerSearch, amp_scale, load_cut, recompute_cut_observables


def plot_antenna_positions(output: Path):
    """Plot PANSY antenna positions and the receiver modules used here."""
    output.parent.mkdir(parents=True, exist_ok=True)
    ready = []
    not_ready = []
    rftx = None
    for serial, ant in pc.antenna.items():
        pos = [ant["x"], ant["y"], ant["z"]]
        if serial == "RFTX":
            rftx = pos
            continue
        try:
            is_ready = int(ant["ready"]) == 1
        except ValueError:
            is_ready = str(ant["ready"]).strip().lower() in {"true", "yes", "ready"}
        if is_ready:
            ready.append(pos)
        else:
            not_ready.append(pos)

    ready = np.asarray(ready, dtype=np.float64)
    not_ready = np.asarray(not_ready, dtype=np.float64) if not_ready else np.empty((0, 3))
    rx_centers = []
    rx_labels = []
    for i, name in enumerate(pc.connections):
        if name == "RFTX" or name not in pc.module_center:
            continue
        rx_centers.append(pc.module_center[name])
        rx_labels.append(str(i))
    rx_centers = np.asarray(rx_centers, dtype=np.float64)

    fig, ax = plt.subplots(figsize=(7.0, 6.8), constrained_layout=True)
    if len(not_ready):
        ax.scatter(not_ready[:, 0], not_ready[:, 1], s=12, color="0.75", label="not ready antennas")
    ax.scatter(ready[:, 0], ready[:, 1], s=16, color="tab:blue", label="ready antennas")
    if rftx is not None:
        ax.scatter([rftx[0]], [rftx[1]], s=120, marker="*", color="black", label="RF reference")
    if len(rx_centers):
        ax.scatter(rx_centers[:, 0], rx_centers[:, 1], s=80, color="tab:red", label="receiver module centers")
        for label, pos in zip(rx_labels, rx_centers):
            ax.text(
                pos[0],
                pos[1],
                label,
                ha="center",
                va="center",
                fontsize=9,
                bbox={"facecolor": "white", "edgecolor": "tab:red", "alpha": 0.8, "boxstyle": "circle"},
            )
    ax.set_xlabel("Array x coordinate (m)")
    ax.set_ylabel("Array y coordinate (m)")
    ax.set_title("PANSY antenna positions and meteor receiver modules")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_tx_array_beam_patterns(output: Path, grid_n=501, model="module_incoherent"):
    """Plot the five transmit array-factor beam patterns over the visible sky."""
    output.parent.mkdir(parents=True, exist_ok=True)
    u, v, w, valid = horizon_grid(grid_n)
    beam_vecs = tx_beam_unit_vectors()
    gain_maps = precompute_tx_array_gain_maps(u, v, w, valid, beam_vecs=beam_vecs, model=model)
    beam_names = ["Zenith", "North", "East", "South", "West"]

    fig, axes = plt.subplots(1, len(beam_vecs), figsize=(16.0, 3.8), constrained_layout=True)
    for beam_i, ax in enumerate(axes):
        gain = gain_maps[beam_i]
        im = ax.pcolormesh(u, v, gain, shading="auto", cmap="viridis", vmin=-20.0, vmax=0.0)
        ax.contour(u, v, gain, levels=[-12.0, -6.0, -3.0], colors="white", linewidths=0.7)
        ax.plot(beam_vecs[beam_i, 0], beam_vecs[beam_i, 1], "rx", ms=7, mew=1.5)
        ax.add_patch(plt.Circle((0.0, 0.0), 1.0, color="black", fill=False, linewidth=0.9))
        ax.set_title(f"{beam_i}: {beam_names[beam_i] if beam_i < len(beam_names) else 'beam'}")
        ax.set_xlabel("u/east direction cosine")
        if beam_i == 0:
            ax.set_ylabel("v/north direction cosine")
        ax.set_aspect("equal")
        ax.set_xlim(-1.02, 1.02)
        ax.set_ylim(-1.02, 1.02)
        ax.grid(True, alpha=0.18)
    fig.colorbar(im, ax=axes, label="Relative TX array-factor gain (dB)")
    fig.suptitle(f"PANSY transmit beam patterns ({model.replace('_', ' ')})")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def horizon_grid(n: int):
    """Return a square u/v grid covering all directions above the horizon."""
    axis = np.linspace(-1.0, 1.0, n)
    u, v = np.meshgrid(axis, axis)
    rho2 = u**2 + v**2
    valid = rho2 <= 1.0
    w = np.full_like(u, np.nan)
    w[valid] = -np.sqrt(1.0 - rho2[valid])
    return u, v, w, valid


def tx_beam_unit_vectors():
    """Transmit beam pointing directions in the same u/v/w convention."""
    mode = pmm.get_m_mode()
    vecs = []
    for az_deg, za_deg in mode["beam_pos_az_za"]:
        el_deg = 90.0 - za_deg
        p_h = np.cos(np.deg2rad(el_deg))
        w = -np.sin(np.deg2rad(el_deg))
        v = p_h * np.cos(-np.deg2rad(az_deg))
        u = -p_h * np.sin(-np.deg2rad(az_deg))
        vecs.append([u, v, w])
    return np.asarray(vecs, dtype=np.float64)


TX_BEAM_SHORT_NAMES = np.asarray(["Z", "N", "E", "S", "W"])


def tx_beam_proximity_label(track, prefix="TX"):
    """Compact plot label for the active commanded-beam proximity metric."""
    mean_dc = track.get("tx_beam_snr_weighted_mean_dc", np.nan)
    tx_term = track.get("combined_tx_term", np.nan)
    if not np.isfinite(mean_dc):
        return ""
    if np.isfinite(mean_dc) and np.isfinite(tx_term):
        return f"{prefix} d={mean_dc:.3f} T={tx_term:.1f}"
    return f"{prefix} d={mean_dc:.3f}"


def tx_beam_center_projection_points(track, candidates):
    """Return all commanded beam-center EW/NS points, marking active beams."""
    if "idx" not in track:
        return []
    beam_vecs = tx_beam_unit_vectors()
    rows = [candidates[i] for i in track["idx"]]
    if not rows:
        return []
    beam_id = np.asarray([r["beam_id"] for r in rows], dtype=np.int64)
    ranges = np.asarray([r["range_km"] for r in rows], dtype=np.float64)
    snr = np.asarray([r["snr"] for r in rows], dtype=np.float64)
    distances = np.asarray(track.get("tx_beam_center_distance_dc", np.full(len(rows), np.nan)), dtype=np.float64)
    all_weight = np.maximum(snr, 1e-6)
    fallback_range_km = float(np.sum(all_weight * ranges) / np.sum(all_weight))
    out = []
    for bid in range(len(beam_vecs)):
        if bid < 0 or bid >= len(beam_vecs):
            continue
        mask = beam_id == bid
        active = bool(np.any(mask))
        if active:
            weight = np.maximum(snr[mask], 1e-6)
            mean_range_km = float(np.sum(weight * ranges[mask]) / np.sum(weight))
            d = distances[mask]
            good = np.isfinite(d)
            mean_d = float(np.sum(weight[good] * d[good]) / np.sum(weight[good])) if np.any(good) else np.nan
            total_weight = float(np.sum(weight))
        else:
            mean_range_km = fallback_range_km
            mean_d = np.nan
            total_weight = 0.0
        ew_ns = mean_range_km * beam_vecs[bid, :2]
        out.append(
            {
                "beam_id": int(bid),
                "label": str(TX_BEAM_SHORT_NAMES[bid]),
                "east_km": float(ew_ns[0]),
                "north_km": float(ew_ns[1]),
                "radius_10deg_km": float(mean_range_km * np.sin(np.deg2rad(10.0))),
                "mean_dc": mean_d,
                "weight": total_weight,
                "active": active,
            }
        )
    return out


def draw_tx_beam_projection_marker(ax, beam, color="0.45", fontsize=6.5):
    """Draw commanded TX beam center and its 10 degree angular footprint in EW/NS."""
    beam_color = f"C{int(beam['beam_id']) % 10}"
    if beam.get("active", False):
        ax.add_patch(
            plt.Circle(
                (beam["east_km"], beam["north_km"]),
                beam["radius_10deg_km"],
                facecolor="0.82",
                edgecolor="none",
                alpha=0.2,
                zorder=0,
            )
        )
    ax.text(
        beam["east_km"],
        beam["north_km"],
        beam["label"],
        fontsize=fontsize + 1.5,
        color=beam_color,
        ha="center",
        va="center",
        fontweight="bold",
        alpha=0.95,
        zorder=0.2,
    )


def scatter_points_by_tx_beam(ax, points, beam_id, mask=None, marker=".", size=3.0, alpha=0.65):
    """Scatter ENU points with categorical colors for transmit beam id."""
    points = np.asarray(points, dtype=np.float64)
    beam_id = np.asarray(beam_id, dtype=np.int64)
    if mask is None:
        mask = np.ones(len(points), dtype=bool)
    else:
        mask = np.asarray(mask, dtype=bool)
    for bid in range(len(TX_BEAM_SHORT_NAMES)):
        use = mask & (beam_id == bid)
        if np.any(use):
            ax.plot(
                points[use, 0],
                points[use, 1],
                marker,
                color=f"C{bid}",
                ms=size,
                alpha=alpha,
                zorder=2.5,
            )


def winning_track_for_rti(tracks):
    scored = [t for t in tracks if "combined_rank" in t and ("selection_pred_doppler_km_s" in t or "ballistic_pred_doppler_km_s" in t)]
    if scored:
        return sorted(scored, key=lambda t: t["combined_rank"])[0]
    ballistic = [t for t in tracks if "ballistic_pred_doppler_km_s" in t]
    if ballistic:
        return sorted(ballistic, key=lambda t: t.get("ballistic_rank", 999))[0]
    return None


def circular_phase_residual(observed_phase_rad, model_phase_rad):
    return np.angle(np.exp(1j * observed_phase_rad) * np.exp(-1j * model_phase_rad))


def ballistic_covariance_band(track, t, event_epoch_unix, observable):
    """Linearized 95% band for a scalar observable of the ballistic state."""
    model_type = track.get("selection_model_type", "msis_drag")
    params = np.asarray(track.get("selection_params", track.get("ballistic_params", [])), dtype=np.float64)
    cov = np.asarray(track.get("selection_parameter_covariance", track.get("ballistic_parameter_covariance", [])), dtype=np.float64)
    if model_type == "fixed_velocity":
        expected_shape = (6,)
    else:
        expected_shape = (7,)
    if params.shape != expected_shape or cov.shape != (len(params), len(params)) or not np.all(np.isfinite(cov)):
        return None
    rho_of_alt_m = None
    if model_type != "fixed_velocity":
        try:
            rho_of_alt_m, _meta = pbal.density_interpolator(event_epoch_unix)
        except Exception:
            return None
    steps = np.array([10.0, 10.0, 10.0, 10.0, 10.0, 10.0] + ([] if model_type == "fixed_velocity" else [1e-4]), dtype=np.float64)
    steps = np.maximum(steps, np.abs(params) * 1e-5)
    jac = np.empty((len(t), len(params)), dtype=np.float64)

    def evaluate(p):
        if model_type == "fixed_velocity":
            tau = np.asarray(t, dtype=np.float64) - float(np.min(t))
            pos = (p[:3][None, :] + tau[:, None] * p[3:6][None, :]) / 1e3
            vel = np.repeat((p[3:6] / 1e3)[None, :], len(tau), axis=0)
            return pos, vel
        pos, vel, _beta, _ = propagate_drag_model(p, t, rho_of_alt_m)
        return pos, vel

    for i, step in enumerate(steps):
        p_hi = params.copy()
        p_lo = params.copy()
        p_hi[i] += step
        p_lo[i] -= step
        try:
            hi_pos, hi_vel = evaluate(p_hi)
            lo_pos, lo_vel = evaluate(p_lo)
        except Exception:
            return None
        if observable == "range":
            hi_val = np.linalg.norm(hi_pos, axis=1)
            lo_val = np.linalg.norm(lo_pos, axis=1)
        elif observable == "doppler":
            hi_rng = np.linalg.norm(hi_pos, axis=1)
            lo_rng = np.linalg.norm(lo_pos, axis=1)
            hi_val = np.sum(hi_pos * hi_vel, axis=1) / np.maximum(hi_rng, 1e-9)
            lo_val = np.sum(lo_pos * lo_vel, axis=1) / np.maximum(lo_rng, 1e-9)
        else:
            raise ValueError(f"unknown observable {observable}")
        jac[:, i] = (hi_val - lo_val) / (2.0 * step)
    sigma2 = np.einsum("ij,jk,ik->i", jac, cov, jac)
    sigma = np.sqrt(np.maximum(sigma2, 0.0))
    sigma[~np.isfinite(sigma)] = np.nan
    return 1.96 * sigma


def range_uncertainty_from_ballistic_covariance(track, t, event_epoch_unix):
    return ballistic_covariance_band(track, t, event_epoch_unix, "range")


def doppler_uncertainty_from_ballistic_covariance(track, t, event_epoch_unix):
    return ballistic_covariance_band(track, t, event_epoch_unix, "doppler")


KEPLER_LABELS = ("a", "e", "i", "raan", "argp", "nu", "q")
KEPLER_UNITS = ("AU", "", "deg", "deg", "deg", "deg", "AU")


def read_dasst_orbit_group(orbit_h5_path: Path | None, hypothesis: str):
    if orbit_h5_path is None or not Path(orbit_h5_path).exists():
        return None
    try:
        with h5py.File(orbit_h5_path, "r") as h:
            if hypothesis not in h:
                return None
            grp = h[hypothesis]
            return {
                "source": "DASST orbit",
                "kepler": grp["kepler"][()],
                "kepler_std": grp["kepler_std"][()] if "kepler_std" in grp else np.full(7, np.nan),
                "kepler_samples": grp["kepler_samples"][()] if "kepler_samples" in grp else np.empty((0, 7)),
                "frac_e_gt_1": float(grp.attrs.get("frac_e_gt_1", np.nan)),
            }
    except Exception:
        return None


def candidate_orbit_result(track):
    orbit = track.get("candidate_orbit") if track is not None else None
    if not orbit:
        return None
    kep = np.asarray(orbit.get("kepler", np.full(7, np.nan)), dtype=np.float64)
    std = np.asarray(orbit.get("kepler_std", np.full(7, np.nan)), dtype=np.float64)
    if kep.shape != (7,):
        return None
    return {
        "source": "preliminary orbit",
        "kepler": kep,
        "kepler_std": std if std.shape == (7,) else np.full(7, np.nan),
        "kepler_samples": np.empty((0, 7)),
        "frac_e_gt_1": np.nan,
    }


def orbit_result_for_track(orbit_h5_path: Path | None, track):
    hypothesis = hypothesis_label(track) if track is not None else "H03"
    return read_dasst_orbit_group(orbit_h5_path, hypothesis) or candidate_orbit_result(track)


def format_kepler_lines(result):
    if result is None:
        return []
    kep = np.asarray(result["kepler"], dtype=np.float64)
    std = np.asarray(result.get("kepler_std", np.full(7, np.nan)), dtype=np.float64)
    lines = [str(result.get("source", "orbit"))]
    for i, (label, unit) in enumerate(zip(KEPLER_LABELS, KEPLER_UNITS)):
        value = kep[i] if i < len(kep) else np.nan
        sigma = std[i] if i < len(std) else np.nan
        unit_text = f" {unit}" if unit else ""
        if np.isfinite(sigma):
            lines.append(f"  {label:<5} {value:12.5g} +/- {sigma:.2g}{unit_text}")
        else:
            lines.append(f"  {label:<5} {value:12.5g}{unit_text}")
    frac_e_gt_1 = float(result.get("frac_e_gt_1", np.nan))
    if np.isfinite(frac_e_gt_1):
        lines.append(f"  P(e>1) {frac_e_gt_1:12.3f}")
    return lines


def summary_orbit_text(orbit_h5_path: Path | None, track):
    result = orbit_result_for_track(orbit_h5_path, track)
    return "\n".join(format_kepler_lines(result))


def plot_best_range_fit_panel(ax, winner, competitors, event_epoch_unix, annotate_panel):
    if winner is None or ("selection_model" not in winner and "ballistic_model" not in winner):
        annotate_panel(ax, "No trajectory winner available")
        ax.set_title("1b. Range fit")
        ax.set_xlabel("Time since cut start (s)")
        ax.set_ylabel("Range (km)")
        return
    for track in competitors:
        if track is winner or "ballistic_t" not in track:
            continue
        t_alt = np.asarray(track["ballistic_t"], dtype=np.float64)
        model_alt = np.asarray(track.get("selection_model", track.get("ballistic_model")), dtype=np.float64)
        if len(t_alt) == 0 or len(model_alt) != len(t_alt):
            continue
        order_alt = np.argsort(t_alt)
        ax.plot(
            t_alt[order_alt],
            np.linalg.norm(model_alt, axis=1)[order_alt],
            color="0.70",
            lw=1.0,
            ls="--" if track.get("combined_reject", False) else "-",
            alpha=0.72,
            zorder=1,
        )
        if len(order_alt):
            ax.text(
                t_alt[order_alt][-1],
                np.linalg.norm(model_alt, axis=1)[order_alt][-1],
                hypothesis_label(track),
                fontsize=6.5,
                color="0.55",
                ha="left",
                va="center",
                clip_on=True,
            )
    t = np.asarray(winner["ballistic_t"], dtype=np.float64)
    order = np.argsort(t)
    points = np.asarray(winner["fit_points"], dtype=np.float64)
    model = np.asarray(winner.get("selection_model", winner.get("ballistic_model")), dtype=np.float64)
    keep = np.asarray(winner.get("selection_keep", winner.get("ballistic_keep", np.ones(len(t), dtype=bool))), dtype=bool)
    measured_range = np.linalg.norm(points, axis=1)
    model_range = np.linalg.norm(model, axis=1)
    band_95 = range_uncertainty_from_ballistic_covariance(winner, t, event_epoch_unix)

    ax.plot(t[order], measured_range[order], ".", ms=8.0, color="0.25", alpha=0.55, label="measurements", zorder=100)
    if np.any(~keep):
        out = order[~keep[order]]
        ax.plot(t[out], measured_range[out], "x", ms=9.6, color="tab:red", alpha=0.85, label="clipped", zorder=101)
    ax.plot(t[order], model_range[order], "-", lw=2.0, color="black", label=f"{hypothesis_label(winner)} {winner.get('selection_model_label', 'best fit')}", zorder=5)
    if band_95 is not None and np.any(np.isfinite(band_95)):
        lo = model_range - band_95
        hi = model_range + band_95
        ax.fill_between(t[order], lo[order], hi[order], color="black", alpha=0.15, linewidth=0.0, label="95% fit region", zorder=2)
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Range (km)")
    ax.set_title("1b. Range measurements and best fit")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7, loc="best")


def plot_best_fit_summary_panel(ax, winner, orbit_h5_path, annotate_panel):
    ax.axis("off")
    if winner is None or ("selection_params" not in winner and "ballistic_params" not in winner):
        annotate_panel(ax, "No trajectory fit summary available")
        ax.set_title("2b. Fit summary")
        return
    params = np.asarray(winner.get("selection_params", winner.get("ballistic_params", np.full(7, np.nan))), dtype=np.float64)
    std = np.asarray(winner.get("selection_parameter_std", winner.get("ballistic_parameter_std", np.full(len(params), np.nan))), dtype=np.float64)
    points = np.asarray(winner.get("fit_points", np.empty((0, 3))), dtype=np.float64)
    model = np.asarray(winner.get("selection_model", winner.get("ballistic_model", np.empty((0, 3)))), dtype=np.float64)
    keep = np.asarray(winner.get("selection_keep", winner.get("ballistic_keep", np.ones(len(points), dtype=bool))), dtype=bool)
    range_res_m = np.full(len(points), np.nan)
    if len(points) and len(model) == len(points):
        range_res_m = (np.linalg.norm(points, axis=1) - np.linalg.norm(model, axis=1)) * 1e3
    kept_range = range_res_m[keep] if len(range_res_m) == len(keep) else range_res_m
    kept_range = kept_range[np.isfinite(kept_range)]
    range_rms_m = float(np.sqrt(np.mean(kept_range**2))) if len(kept_range) else np.nan
    beta = float(winner.get("ballistic_coefficient_kg_m2", np.nan)) if winner.get("selection_model_type") != "fixed_velocity" else np.nan
    beta_sigma = np.log(10.0) * beta * std[6] if std.shape == (7,) and np.isfinite(beta) else np.nan
    speed_km_s = np.linalg.norm(params[3:6]) / 1e3
    speed_sigma_km_s = np.linalg.norm(std[3:6]) / 1e3 if std.shape == (7,) else np.nan
    lines = [
        f"Alias                 {hypothesis_label(winner)}  candidate {candidate_number(winner)}",
        f"Alias ranking fit     {winner.get('selection_model_label', winner.get('ballistic_model_type', 'MSIS drag'))}",
        f"Trajectory model      {winner.get('physics_best_model_label', winner.get('physics_best_model', 'not fit'))}",
        f"Selection rule        TX center if chi2nu <= {winner.get('combined_good_fit_redchi_threshold', np.nan):.3g}",
        f"Combined rank/score   {winner.get('combined_rank', -1)} / {winner.get('combined_score', np.nan):.3g}",
        f"N kept / clipped      {winner.get('selection_n', winner.get('ballistic_n', 0))} / {winner.get('selection_n_outliers', winner.get('ballistic_n_outliers', 0))}",
        "",
        f"chi2 / dof            {winner.get('selection_chi2', winner.get('ballistic_chi2', np.nan)):.1f} / {winner.get('selection_dof', winner.get('ballistic_dof', 0))}",
        f"reduced chi2          {winner.get('selection_reduced_chi2', winner.get('ballistic_reduced_chi2', np.nan)):.3g}",
        f"position RMS          {winner.get('selection_pos_rms_km', winner.get('ballistic_pos_rms_km', np.nan)):.3g} km",
        f"range RMS             {range_rms_m:.3g} m",
        f"Doppler RMS           {winner.get('selection_dop_rms_km_s', winner.get('ballistic_dop_rms_km_s', np.nan)) * 1e3:.3g} m/s",
        "",
        f"h0                    {params[2] / 1e3:.3f} +/- {std[2] / 1e3:.3f} km",
        f"|v0|                  {speed_km_s:.3f} +/- {speed_sigma_km_s:.3f} km/s",
        f"vE,vN,vU              {params[3]/1e3:.2f} {params[4]/1e3:.2f} {params[5]/1e3:.2f} km/s",
        f"beta                  {beta:.3g} +/- {beta_sigma:.2g} kg/m2",
    ]
    if "physics_bic_fixed_velocity" in winner:
        lines.extend(
            [
                "",
                "Trajectory model BIC",
                f"  fixed velocity       {winner.get('physics_bic_fixed_velocity', np.nan):.3g}",
                f"  fixed CdA/m          {winner.get('physics_bic_fixed_am', np.nan):.3g}",
                f"  shrinking radius     {winner.get('physics_bic_shrinking_radius', np.nan):.3g}",
            ]
        )
    elif "physics_model_error" in winner:
        lines.extend(["", f"Trajectory model fit failed: {winner['physics_model_error']}"])
    orbit_line = summary_orbit_text(orbit_h5_path, winner)
    if orbit_line:
        lines.extend(["", orbit_line.replace("\n", "\n")])
    ax.text(0.02, 0.98, "\n".join(lines), ha="left", va="top", transform=ax.transAxes, family="monospace", fontsize=8.5)
    ax.set_title("2b. Best-fit summary")


def doppler_series_for_rti(obs, tracks, fit_t0_s=None):
    """Best available per-pulse Doppler for Doppler-sliced RTI."""
    doppler_mps = np.asarray(obs["doppler_mps"], dtype=np.float64).copy()
    track = winning_track_for_rti(tracks)
    if track is None:
        return doppler_mps, "measured peak Doppler"
    t_obs = np.asarray(obs["tx_idx"], dtype=np.float64) / 1e6
    if fit_t0_s is None:
        fit_t0_s = float(t_obs[0])
    t_rel = t_obs - float(fit_t0_s)
    t_fit = np.asarray(track.get("ballistic_t", []), dtype=np.float64)
    pred = np.asarray(track.get("ballistic_pred_doppler_km_s", []), dtype=np.float64) * 1e3
    pred = np.asarray(track.get("selection_pred_doppler_km_s", pred / 1e3), dtype=np.float64) * 1e3
    good = np.isfinite(t_fit) & np.isfinite(pred)
    if np.count_nonzero(good) >= 2:
        order = np.argsort(t_fit[good])
        doppler_mps = np.interp(t_rel, t_fit[good][order], pred[good][order], left=pred[good][order][0], right=pred[good][order][-1])
        return doppler_mps, f"{hypothesis_label(track)} {track.get('selection_model_label', 'trajectory')} Doppler"
    if np.count_nonzero(good) == 1:
        doppler_mps[:] = pred[good][0]
        return doppler_mps, f"{hypothesis_label(track)} {track.get('selection_model_label', 'trajectory')} Doppler"
    return doppler_mps, "measured peak Doppler"


def pulse_pair_phase_from_cut(obs, cut, interp=1):
    """Build same-beam pulse-pair phase from matched-filter bins.

    The already-determined range/Doppler peak of the earlier pulse selects the
    bin.  The later same-beam pulse is evaluated at that exact bin, using its
    own transmit pulse, so the phase observable is repeatable and tied to the
    initial within-pulse Doppler pick.
    """
    z_rx = np.asarray(cut["zrx_echoes_re"], dtype=np.complex64) + 1j * np.asarray(cut["zrx_echoes_im"], dtype=np.complex64)
    z_tx = np.asarray(cut["ztx_pulses_re"], dtype=np.complex64) + 1j * np.asarray(cut["ztx_pulses_im"], dtype=np.complex64)
    scale = amp_scale()
    for i in range(z_rx.shape[1]):
        z_rx[:, i, :] *= scale[i]

    rds = RangeDopplerSearch(
        txlen=z_tx.shape[1],
        echolen=z_rx.shape[2],
        n_channels=z_rx.shape[1],
        interp=interp,
    )
    zf_by_pulse = []
    for pulse in range(z_rx.shape[0]):
        txi = np.repeat(z_tx[pulse, :], interp)
        zf_channels = []
        for chi in range(z_rx.shape[1]):
            zr = np.repeat(z_rx[pulse, chi, :], interp)
            z = zr[rds.idx_mat] * txi[None, :]
            zd = rds.decim(z)
            zf = np.fft.fftshift(np.fft.fft(zd, rds.fftlen, axis=1), axes=1)
            zf_channels.append(zf)
        zf_by_pulse.append(zf_channels)

    tx_idx = np.asarray(obs["tx_idx"], dtype=np.float64)
    beam_id = np.asarray(obs["beam_id"], dtype=np.int64)
    doppler_mps = np.asarray(obs["doppler_mps"], dtype=np.float64)
    range_km = np.asarray(obs["range_km"], dtype=np.float64)
    snr = np.asarray(obs["snr"], dtype=np.float64)
    delays = np.asarray(cut["delays"], dtype=np.float64)
    drg_km = sc.c / 1e6 / 2.0 / 1e3 / interp
    peak_rg = np.rint(range_km / drg_km - delays).astype(np.int64)
    peak_di = np.asarray([int(np.argmin(np.abs(rds.dopv - dop))) for dop in doppler_mps], dtype=np.int64)
    tx_time_s = tx_idx / 1e6

    zpp = np.full(len(tx_idx), np.nan + 1j * np.nan, dtype=np.complex64)
    phase_dt_s = np.full(len(tx_idx), np.nan, dtype=np.float64)
    coarse_phase_rad = np.full(len(tx_idx), np.nan, dtype=np.float64)
    pair_snr = np.full(len(tx_idx), np.nan, dtype=np.float64)
    prev_t_rel = np.full(len(tx_idx), np.nan, dtype=np.float64)
    prev_pulse = np.full(len(tx_idx), -1, dtype=np.int64)

    for beam in np.unique(beam_id):
        idx = np.where(beam_id == beam)[0]
        idx = idx[np.argsort(tx_time_s[idx])]
        for prev_i, cur_i in zip(idx[:-1], idx[1:]):
            dt_s = tx_time_s[cur_i] - tx_time_s[prev_i]
            if not np.isfinite(dt_s) or dt_s <= 0.0:
                continue
            rg = peak_rg[prev_i]
            di = peak_di[prev_i]
            if rg < 0 or rg >= rds.n_rg or di < 0 or di >= rds.fftlen:
                continue
            cross = 0.0j
            for chi in range(z_rx.shape[1]):
                cross += zf_by_pulse[prev_i][chi][rg, di] * np.conj(zf_by_pulse[cur_i][chi][rg, di])
            if not np.isfinite(cross.real) or not np.isfinite(cross.imag) or np.abs(cross) == 0.0:
                continue
            zpp[cur_i] = cross / np.abs(cross)
            phase_dt_s[cur_i] = dt_s
            coarse_hz = 2.0 * doppler_mps[prev_i] / pc.wavelength
            coarse_phase_rad[cur_i] = np.angle(np.exp(-1j * 2.0 * np.pi * coarse_hz * dt_s))
            pair_snr[cur_i] = min(snr[prev_i], snr[cur_i])
            prev_t_rel[cur_i] = tx_time_s[prev_i] - tx_time_s[0]
            prev_pulse[cur_i] = prev_i
    return zpp, phase_dt_s, coarse_phase_rad, pair_snr, prev_t_rel, prev_pulse


def add_pulse_pair_phase_observables(obs, cut):
    """Add wrapped same-beam pulse-to-pulse phase samples to cut observables."""
    zpp, phase_dt_s, coarse_phase_rad, pair_snr, prev_t_rel, prev_pulse = pulse_pair_phase_from_cut(obs, cut, interp=1)
    obs = dict(obs)
    obs["zpp"] = zpp
    obs["zpp_phase_rad"] = np.angle(zpp)
    obs["zpp_phase_dt_s"] = phase_dt_s
    obs["zpp_doppler_mps"] = -obs["zpp_phase_rad"] * pc.wavelength / (4.0 * np.pi * phase_dt_s)
    obs["zpp_coarse_phase_rad"] = coarse_phase_rad
    obs["zpp_coarse_doppler_mps"] = -coarse_phase_rad * pc.wavelength / (4.0 * np.pi * phase_dt_s)
    obs["zpp_pair_snr"] = pair_snr
    obs["zpp_prev_t_rel"] = prev_t_rel
    obs["zpp_prev_pulse"] = prev_pulse
    return obs


def cached_cut_observables(cut: dict, interp: int = 1):
    """Use cached cut detections and compute only interferometric phase at those bins."""
    required = {"c_tx_idx", "c_beam_idx", "c_range_km", "c_doppler", "c_snr"}
    if not required.issubset(cut):
        raise KeyError("cut does not contain cached detection observables")

    z_rx = np.asarray(cut["zrx_echoes_re"], dtype=np.complex64) + 1j * np.asarray(cut["zrx_echoes_im"], dtype=np.complex64)
    z_tx = np.asarray(cut["ztx_pulses_re"], dtype=np.complex64) + 1j * np.asarray(cut["ztx_pulses_im"], dtype=np.complex64)
    scale = amp_scale()
    for i in range(z_rx.shape[1]):
        z_rx[:, i, :] *= scale[i]

    tx_idx_all = np.asarray(cut["tx_idx"], dtype=np.int64)
    delays_all = np.asarray(cut["delays"], dtype=np.int64)
    tx_idx = np.asarray(cut["c_tx_idx"], dtype=np.float64)
    cached_beam = np.asarray(cut["c_beam_idx"], dtype=np.int64)
    range_km = np.asarray(cut["c_range_km"], dtype=np.float64)
    doppler_mps = np.asarray(cut["c_doppler"], dtype=np.float64)
    snr = np.asarray(cut["c_snr"], dtype=np.float64)

    pulse_lookup = {int(sample): i for i, sample in enumerate(tx_idx_all)}
    pulse_i = np.asarray([pulse_lookup.get(int(sample), -1) for sample in tx_idx], dtype=np.int64)
    ok = pulse_i >= 0
    tx_idx = tx_idx[ok]
    beam_id = np.asarray(cut["beam_id"], dtype=np.int64)[pulse_i[ok]]
    if np.any((beam_id < 0) | (beam_id >= 5)):
        beam_id = cached_beam[ok]
    range_km = range_km[ok]
    doppler_mps = doppler_mps[ok]
    snr = snr[ok]
    pulse_i = pulse_i[ok]

    rds = RangeDopplerSearch(
        txlen=z_tx.shape[1],
        echolen=z_rx.shape[2],
        interp=interp,
        n_channels=z_rx.shape[1],
    )
    drg_km = sc.c / 1e6 / 2.0 / 1e3 / interp
    xct = np.full((len(tx_idx), len(rds.ch_pairs)), np.nan + 1j * np.nan, dtype=np.complex64)
    for out_i, pi in enumerate(pulse_i):
        rg = int(np.rint(range_km[out_i] / drg_km - delays_all[pi]))
        if rg < 0 or rg >= rds.n_rg:
            continue
        di = int(np.argmin(np.abs(rds.dopv - doppler_mps[out_i])))
        txi = np.repeat(z_tx[pi, :], interp)
        zf = []
        row_idx = rg + np.arange(rds.txlen, dtype=np.int64)
        for chi in range(z_rx.shape[1]):
            zr = np.repeat(z_rx[pi, chi, :], interp)
            z = zr[row_idx] * txi
            new_width = int(rds.txlen / rds.fdec)
            zd = np.zeros(new_width, dtype=np.complex64)
            base = np.arange(new_width, dtype=np.int64) * rds.fdec
            for dec_i in range(rds.fdec):
                zd += z[base + dec_i]
            spec = np.fft.fftshift(np.fft.fft(zd, rds.fftlen))
            zf.append(spec[di])
        for pair_i, (a, b) in enumerate(rds.ch_pairs):
            xct[out_i, pair_i] = zf[a] * np.conj(zf[b])

    rti_snr = np.zeros((len(tx_idx), 1600), dtype=np.float32)
    dti_mps = np.zeros((len(tx_idx), 1600), dtype=np.float32)
    full_range_grid_km = np.arange(1600, dtype=np.float64) * (sc.c / 1e6 / 2.0 / 1e3)
    full_bin = np.rint(range_km / (sc.c / 1e6 / 2.0 / 1e3)).astype(np.int64)
    for out_i, bin_i, val_snr, val_dop in zip(range(len(tx_idx)), full_bin, snr, doppler_mps):
        if 0 <= bin_i < rti_snr.shape[1]:
            rti_snr[out_i, bin_i] = val_snr
            dti_mps[out_i, bin_i] = val_dop

    return {
        "tx_idx": tx_idx,
        "beam_id": beam_id,
        "range_km": range_km,
        "doppler_mps": doppler_mps,
        "snr": snr,
        "xc": xct,
        "rti_snr": rti_snr,
        "dti_mps": dti_mps,
        "range_grid_km": full_range_grid_km,
        "cached_peak_rti": True,
        "radar_frequency_hz": float(pc.freq),
    }


def matched_filter_rti_at_doppler(cut, doppler_mps, interp=1):
    """Compute RTI by slicing the range-Doppler matched filter at a specified Doppler per pulse."""
    z_rx = np.asarray(cut["zrx_echoes_re"], dtype=np.complex64) + 1j * np.asarray(cut["zrx_echoes_im"], dtype=np.complex64)
    z_tx = np.asarray(cut["ztx_pulses_re"], dtype=np.complex64) + 1j * np.asarray(cut["ztx_pulses_im"], dtype=np.complex64)
    scale = amp_scale()
    for i in range(z_rx.shape[1]):
        z_rx[:, i, :] *= scale[i]
    rds = RangeDopplerSearch(txlen=z_tx.shape[1], echolen=z_rx.shape[2], interp=interp, n_channels=z_rx.shape[1])
    delays = np.asarray(cut["delays"], dtype=np.int64)
    if len(doppler_mps) != z_rx.shape[0]:
        raise ValueError(f"doppler series has {len(doppler_mps)} samples, cut has {z_rx.shape[0]} pulses")
    rti_snr = np.zeros((z_rx.shape[0], 1600), dtype=np.float32)
    ridx = np.arange(rds.n_rg, dtype=int)
    doppler_mps = np.asarray(doppler_mps, dtype=np.float64)
    for ti in range(z_rx.shape[0]):
        mf, _pprof, _peak_dopv, noise_floor, _xc, _rgmax = rds.mf(z_tx[ti, :], z_rx[ti, :, :])
        dop_idx = int(np.argmin(np.abs(rds.dopv - doppler_mps[ti])))
        prof = mf[:, dop_idx]
        valid = ridx + delays[ti]
        good = (valid >= 0) & (valid < rti_snr.shape[1])
        rti_snr[ti, valid[good]] = (prof[good] - noise_floor) / max(float(noise_floor), 1e-12)
    return rti_snr


def plot_matched_filter_rti_panel(ax, obs_rti, cut, tracks, fit_t0_s=None, obs_overlay=None):
    """Panel 1c: range-time-intensity at the best Doppler model."""
    if bool(obs_rti.get("cached_peak_rti", False)):
        raw_tx_idx = np.asarray(cut["tx_idx"], dtype=np.float64)
        t0 = float(raw_tx_idx[0] / 1e6) if fit_t0_s is None else float(fit_t0_s)
        t_rel = raw_tx_idx / 1e6 - t0
        cached_tx_s = np.asarray(obs_rti["tx_idx"], dtype=np.float64) / 1e6
        cached_dop = np.asarray(obs_rti["doppler_mps"], dtype=np.float64)
        if len(cached_tx_s) and len(cached_dop):
            order = np.argsort(cached_tx_s)
            raw_dop = np.interp(raw_tx_idx / 1e6, cached_tx_s[order], cached_dop[order], left=cached_dop[order][0], right=cached_dop[order][-1])
        else:
            raw_dop = np.zeros_like(raw_tx_idx, dtype=np.float64)
        raw_obs = {
            "tx_idx": raw_tx_idx,
            "doppler_mps": raw_dop,
            "range_grid_km": np.asarray(obs_rti["range_grid_km"], dtype=np.float64),
        }
        doppler_mps, label = doppler_series_for_rti(raw_obs, tracks, fit_t0_s=t0)
        try:
            rti = matched_filter_rti_at_doppler(cut, doppler_mps, interp=1)
            title = f"1c. Range-Doppler RTI at {label}"
        except Exception as exc:
            t0 = float(obs_rti["tx_idx"][0] / 1e6) if fit_t0_s is None else float(fit_t0_s)
            t_rel = np.asarray(obs_rti["tx_idx"], dtype=np.float64) / 1e6 - t0
            rti = np.asarray(obs_rti.get("rti_snr", np.empty((0, 0))), dtype=np.float64)
            title = f"1c. Range-Doppler RTI fallback ({type(exc).__name__})"
    else:
        t0 = float(obs_rti["tx_idx"][0] / 1e6) if fit_t0_s is None else float(fit_t0_s)
        t_rel = np.asarray(obs_rti["tx_idx"], dtype=np.float64) / 1e6 - t0
        doppler_mps, label = doppler_series_for_rti(obs_rti, tracks, fit_t0_s=t0)
        try:
            rti = matched_filter_rti_at_doppler(cut, doppler_mps, interp=1)
            title = f"1c. Range-Doppler RTI at {label}"
        except Exception as exc:
            rti = np.asarray(obs_rti.get("rti_snr", np.empty((0, 0))), dtype=np.float64)
            title = f"1c. Range-Doppler RTI fallback ({type(exc).__name__})"
    rti_db = 10.0 * np.log10(np.maximum(rti, 1e-6))
    finite = rti_db[np.isfinite(rti_db)]
    if finite.size:
        vmax = max(3.0, float(np.nanpercentile(finite, 99.8)))
        vmin = max(-18.0, vmax - 32.0)
    else:
        vmin = -18.0
        vmax = 20.0
    im = ax.pcolormesh(
        t_rel,
        obs_rti["range_grid_km"],
        rti_db.T,
        shading="auto",
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
    )
    overlay = obs_rti if obs_overlay is None else obs_overlay
    if len(overlay["range_km"]):
        overlay_t = np.asarray(overlay["tx_idx"], dtype=np.float64) / 1e6 - t0
        ax.scatter(
            overlay_t,
            np.asarray(overlay["range_km"], dtype=np.float64),
            s=5.0,
            c="black",
            marker="o",
            linewidths=0.0,
            alpha=0.85,
            zorder=6,
        )
    if len(overlay["range_km"]):
        pad = 8.0
        ax.set_ylim(max(0.0, float(np.nanmin(overlay["range_km"]) - pad)), float(np.nanmax(overlay["range_km"]) + pad))
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Range (km)")
    ax.set_title(title)
    ax.grid(False)
    return im


def winning_track_for_view(tracks):
    scored = [t for t in tracks if "combined_rank" in t and "fit_points" in t]
    if scored:
        return sorted(scored, key=lambda t: t["combined_rank"])[0]
    ballistic = [t for t in tracks if "ballistic_rank" in t and "fit_points" in t]
    if ballistic:
        return sorted(ballistic, key=lambda t: t["ballistic_rank"])[0]
    return None


def final_trajectory_tracks_for_view(tracks):
    """Tracks with a final selected trajectory fit suitable for event summary panels."""
    selected = [
        t
        for t in tracks
        if t.get("reason") == "kept"
        and "combined_rank" in t
        and "fit_points" in t
        and "selection_model" in t
        and "selection_pred_doppler_km_s" in t
        and not bool(t.get("combined_reject", False))
    ]
    if not selected:
        selected = [
            t
            for t in tracks
            if t.get("reason") == "kept"
            and "combined_rank" in t
            and "fit_points" in t
            and "selection_model" in t
            and "selection_pred_doppler_km_s" in t
        ]
    return sorted(selected, key=lambda t: t.get("combined_rank", 999))


def straight_line_sumsq_km2(track):
    """Sum squared distances to the line joining first/last kept fit points."""
    points = np.asarray(track.get("fit_points", np.empty((0, 3))), dtype=np.float64)
    if len(points) < 2:
        return np.nan
    keep = np.asarray(
        track.get(
            "selection_keep",
            track.get("fixed_velocity_keep", track.get("ballistic_keep", np.ones(len(points), dtype=bool))),
        ),
        dtype=bool,
    )
    if len(keep) != len(points):
        keep = np.ones(len(points), dtype=bool)
    times = np.asarray(track.get("ballistic_t", track.get("fit_t", np.arange(len(points)))), dtype=np.float64)
    if len(times) != len(points):
        times = np.arange(len(points), dtype=np.float64)
    idx = np.flatnonzero(keep & np.all(np.isfinite(points), axis=1) & np.isfinite(times))
    if len(idx) < 2:
        return np.nan
    idx = idx[np.argsort(times[idx])]
    p0 = points[idx[0]]
    p1 = points[idx[-1]]
    direction = p1 - p0
    norm = float(np.linalg.norm(direction))
    if not np.isfinite(norm) or norm <= 0.0:
        return np.nan
    direction = direction / norm
    residual = points[idx] - (p0 + np.outer((points[idx] - p0) @ direction, direction))
    return float(np.sum(np.linalg.norm(residual, axis=1) ** 2))


def set_winning_path_view(ax, tracks, candidates, pad_km=8.0):
    """Constrain EW/NS view to the winning path and its active TX beam centers."""
    winner = winning_track_for_view(tracks)
    if winner is None:
        return
    pts = [np.asarray(winner["fit_points"][:, :2], dtype=np.float64)]
    if "selection_model" in winner:
        pts.append(np.asarray(winner["selection_model"][:, :2], dtype=np.float64))
    elif "ballistic_model" in winner:
        pts.append(np.asarray(winner["ballistic_model"][:, :2], dtype=np.float64))
    active_beams = [
        [beam["east_km"], beam["north_km"]]
        for beam in tx_beam_center_projection_points(winner, candidates)
        if beam.get("active", False)
    ]
    if active_beams:
        pts.append(np.asarray(active_beams, dtype=np.float64))
    xy = np.vstack([p for p in pts if p.size])
    if not len(xy):
        return
    xmin, ymin = np.nanmin(xy, axis=0)
    xmax, ymax = np.nanmax(xy, axis=0)
    span = max(float(xmax - xmin), float(ymax - ymin), 2.0)
    pad = max(pad_km, 0.18 * span)
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    half = 0.5 * span + pad
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)


def tx_beam_decision_note(tracks):
    """Explain whether TX beam proximity affected the final choice."""
    scored = [
        t
        for t in tracks
        if "combined_rank" in t
        and "selection_reduced_chi2" in t
        and np.isfinite(t.get("selection_reduced_chi2", np.nan))
    ]
    if len(scored) < 2:
        return ""
    by_combined = sorted(scored, key=lambda t: t["combined_rank"])
    by_ballistic = sorted(scored, key=lambda t: t["selection_reduced_chi2"])
    winner = by_combined[0]
    best_ballistic = by_ballistic[0]
    second_ballistic = by_ballistic[1]
    close_ballistic = (
        second_ballistic["selection_reduced_chi2"]
        <= best_ballistic["selection_reduced_chi2"] + max(1.0, 0.35 * best_ballistic["selection_reduced_chi2"])
    )
    if winner is not best_ballistic:
        return (
            "TX proximity changed winner:\n"
            f"combined {hypothesis_label(winner)} vs trajectory {hypothesis_label(best_ballistic)}\n"
            f"{hypothesis_label(winner)} TX d={winner.get('tx_beam_snr_weighted_mean_dc', np.nan):.3f}, "
            f"{hypothesis_label(best_ballistic)} d={best_ballistic.get('tx_beam_snr_weighted_mean_dc', np.nan):.3f}"
        )
    if close_ballistic:
        return (
            "Trajectory fits are close;\n"
            "TX proximity helps choose:\n"
            f"{hypothesis_label(winner)} d={winner.get('tx_beam_snr_weighted_mean_dc', np.nan):.3f}, "
            f"{hypothesis_label(second_ballistic)} d={second_ballistic.get('tx_beam_snr_weighted_mean_dc', np.nan):.3f}"
        )
    return (
        "Winner primarily trajectory;\n"
        f"TX d={winner.get('tx_beam_snr_weighted_mean_dc', np.nan):.3f}"
    )


def tx_array_positions():
    """Return ready transmit antenna positions for array-factor diagnostics."""
    pos = []
    for serial, ant in pc.antenna.items():
        if serial == "RFTX":
            continue
        try:
            ready = int(ant["ready"]) == 1
        except ValueError:
            ready = str(ant["ready"]).strip().lower() in {"true", "yes", "ready"}
        if ready:
            pos.append([ant["x"], ant["y"], ant["z"]])
    if not pos:
        raise RuntimeError("no ready transmit antennas found")
    return np.asarray(pos, dtype=np.float64)


def tx_array_module_positions():
    """Ready transmit antenna positions grouped by physical module."""
    groups = []
    for name in sorted(pc.modules):
        pos = []
        for ant in pc.antenna.values():
            if ant["name"] != name or ant["serial"] == "RFTX":
                continue
            try:
                ready = int(ant["ready"]) == 1
            except ValueError:
                ready = str(ant["ready"]).strip().lower() in {"true", "yes", "ready"}
            if ready:
                pos.append([ant["x"], ant["y"], ant["z"]])
        if pos:
            groups.append(np.asarray(pos, dtype=np.float64))
    if not groups:
        raise RuntimeError("no ready transmit antenna modules found")
    return groups


def tx_module_center_positions():
    """One coherent array element at each PANSY antenna module center."""
    pos = []
    for name, center in sorted(pc.module_center.items()):
        if name == "RFTX":
            continue
        pos.append(np.asarray(center, dtype=np.float64))
    if not pos:
        raise RuntimeError("no transmit module centers found")
    pos = np.asarray(pos, dtype=np.float64)
    return pos - np.mean(pos, axis=0, keepdims=True)


def tx_array_gain_db(uvw, beam_id, tx_pos=None, beam_vecs=None, model="module_center_coherent"):
    """Relative transmit array-factor power gain in dB for each trial direction.

    ``model='module_center_coherent'`` uses one coherent array element at each
    PANSY module center.  This is the beam pattern used for the transmit-lobe
    consistency diagnostic.
    """
    uvw = np.asarray(uvw, dtype=np.float64)
    beam_id = np.asarray(beam_id, dtype=np.int64)
    if beam_vecs is None:
        beam_vecs = tx_beam_unit_vectors()
    k0 = 2.0 * np.pi / pc.wavelength
    gain = np.full(len(uvw), np.nan, dtype=np.float64)
    if model == "module_center_coherent":
        centers = tx_module_center_positions() if tx_pos is None else tx_pos
        for i, direction in enumerate(uvw):
            if not np.all(np.isfinite(direction)):
                continue
            steer = beam_vecs[beam_id[i]]
            phase = k0 * (centers @ (direction - steer))
            af = np.abs(np.mean(np.exp(1j * phase)))
            gain[i] = 20.0 * np.log10(max(af, 1e-8))
        return gain

    if model == "module_incoherent":
        modules = tx_array_module_positions() if tx_pos is None else tx_pos
        for i, direction in enumerate(uvw):
            if not np.all(np.isfinite(direction)):
                continue
            steer = beam_vecs[beam_id[i]]
            powers = []
            for pos in modules:
                phase = k0 * (pos @ (direction - steer))
                powers.append(np.abs(np.mean(np.exp(1j * phase))) ** 2)
            gain[i] = 10.0 * np.log10(max(float(np.mean(powers)), 1e-16))
        return gain

    if tx_pos is None:
        tx_pos = tx_array_positions()
    for i, direction in enumerate(uvw):
        if not np.all(np.isfinite(direction)):
            continue
        steer = beam_vecs[beam_id[i]]
        phase = k0 * (tx_pos @ (direction - steer))
        af = np.abs(np.mean(np.exp(1j * phase)))
        gain[i] = 20.0 * np.log10(max(af, 1e-8))
    return gain


def precompute_tx_array_gain_maps(u, v, w, valid, tx_pos=None, beam_vecs=None, model="module_center_coherent"):
    """Precompute TX array-factor gain maps for all transmit beams on a sky grid."""
    if tx_pos is None and model == "module_center_coherent":
        tx_pos = tx_module_center_positions()
    if tx_pos is None and model == "element_coherent":
        tx_pos = tx_array_positions()
    if beam_vecs is None:
        beam_vecs = tx_beam_unit_vectors()
    uvw = np.column_stack([u[valid], v[valid], w[valid]])
    gain_maps = np.full((len(beam_vecs),) + u.shape, np.nan, dtype=np.float32)
    k0 = 2.0 * np.pi / pc.wavelength
    if model == "module_center_coherent":
        for beam_i, steer in enumerate(beam_vecs):
            phase = k0 * (tx_pos @ (uvw - steer).T)
            af = np.abs(np.mean(np.exp(1j * phase), axis=0))
            flat = np.full(u.shape, np.nan, dtype=np.float32)
            flat[valid] = (20.0 * np.log10(np.maximum(af, 1e-8))).astype(np.float32)
            gain_maps[beam_i] = flat
        return gain_maps

    if model == "module_incoherent":
        modules = tx_array_module_positions() if tx_pos is None else tx_pos
        for beam_i, steer in enumerate(beam_vecs):
            power = np.zeros(len(uvw), dtype=np.float64)
            for pos in modules:
                phase = k0 * (pos @ (uvw - steer).T)
                power += np.abs(np.mean(np.exp(1j * phase), axis=0)) ** 2
            power /= len(modules)
            flat = np.full(u.shape, np.nan, dtype=np.float32)
            flat[valid] = (10.0 * np.log10(np.maximum(power, 1e-16))).astype(np.float32)
            gain_maps[beam_i] = flat
        return gain_maps

    for beam_i, steer in enumerate(beam_vecs):
        phase = k0 * (tx_pos @ (uvw - steer).T)
        af = np.abs(np.mean(np.exp(1j * phase), axis=0))
        flat = np.full(u.shape, np.nan, dtype=np.float32)
        flat[valid] = (20.0 * np.log10(np.maximum(af, 1e-8))).astype(np.float32)
        gain_maps[beam_i] = flat
    return gain_maps


def tx_grating_lobe_centroids_from_gain_maps(
    gain_maps,
    u,
    v,
    valid,
    gaussian_width_dc=0.1,
    min_separation_dc=0.25,
    threshold_rel=0.12,
):
    """Find smoothed TX grating-lobe centers on the full horizon in u/v space."""
    du = float(np.nanmedian(np.diff(u[0, :])))
    dv = float(np.nanmedian(np.diff(v[:, 0])))
    pixel = max(abs(du), abs(dv))
    sigma_pix = max(1.0, gaussian_width_dc / pixel)
    sep_pix = max(3, int(round(min_separation_dc / pixel)))
    if sep_pix % 2 == 0:
        sep_pix += 1

    centroids = []
    smoothed_maps = []
    for beam_i in range(gain_maps.shape[0]):
        power = 10.0 ** (np.asarray(gain_maps[beam_i], dtype=np.float64) / 10.0)
        power = np.where(valid & np.isfinite(power), power, 0.0)
        smooth = ndi.gaussian_filter(power, sigma=sigma_pix, mode="constant", cval=0.0)
        smooth = np.where(valid, smooth, np.nan)
        maxval = np.nanmax(smooth)
        smoothed_maps.append(smooth / maxval if np.isfinite(maxval) and maxval > 0.0 else smooth)
        if not np.isfinite(maxval) or maxval <= 0.0:
            centroids.append(np.empty((0, 3), dtype=np.float64))
            continue
        local = ndi.maximum_filter(np.nan_to_num(smooth, nan=-np.inf), size=sep_pix, mode="constant", cval=-np.inf)
        mask = valid & np.isfinite(smooth) & (smooth == local) & (smooth >= threshold_rel * maxval)
        ii, jj = np.where(mask)
        if len(ii) == 0:
            centroids.append(np.empty((0, 3), dtype=np.float64))
            continue
        rel = smooth[ii, jj] / maxval
        order = np.argsort(rel)[::-1]
        centroids.append(np.column_stack([u[ii[order], jj[order]], v[ii[order], jj[order]], rel[order]]))
    return centroids, np.asarray(smoothed_maps)


def nearest_tx_grating_lobe_distances(dirs, beam_id, centroids):
    """Distance in direction-cosine space to nearest active-beam lobe centroid."""
    dirs = np.asarray(dirs, dtype=np.float64)
    beam_id = np.asarray(beam_id, dtype=np.int64)
    dist = np.full(len(dirs), np.nan, dtype=np.float64)
    for beam in range(len(centroids)):
        use = beam_id == beam
        if not np.any(use) or len(centroids[beam]) == 0:
            continue
        centers = centroids[beam][:, :2]
        uv = dirs[use, :2]
        dist[use] = np.sqrt(np.min(np.sum((uv[:, None, :] - centers[None, :, :]) ** 2, axis=2), axis=1))
    return dist


def plot_tx_grating_lobe_centroids(u, v, valid, smoothed_maps, centroids, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    beam_names = ["Zenith", "North", "East", "South", "West"]
    fig, axes = plt.subplots(3, 3, figsize=(10.5, 10.0), constrained_layout=True)
    panel_locs = {1: (0, 1), 4: (1, 0), 0: (1, 1), 2: (1, 2), 3: (2, 1)}
    for ax in axes.ravel():
        ax.axis("off")
    im = None
    for beam, (row, col) in panel_locs.items():
        ax = axes[row, col]
        ax.axis("on")
        panel = np.where(valid, smoothed_maps[beam], np.nan)
        im = ax.pcolormesh(u, v, panel, shading="auto", cmap="turbo", vmin=0.0, vmax=1.0)
        ax.contour(u, v, panel, levels=[0.12, 0.25, 0.5, 0.75], colors="white", linewidths=0.55, alpha=0.75)
        if len(centroids[beam]):
            ax.scatter(centroids[beam][:, 0], centroids[beam][:, 1], s=22, facecolors="none", edgecolors="black", linewidths=0.9)
        ax.add_patch(plt.Circle((0.0, 0.0), 1.0, color="black", fill=False, linewidth=0.8))
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(-1.0, 1.0)
        ax.set_aspect("equal")
        ax.set_title(f"{beam}: {beam_names[beam]} ({len(centroids[beam])})")
        ax.set_xlabel("u direction cosine")
        if col == 0:
            ax.set_ylabel("v direction cosine")
    fig.colorbar(im, ax=axes.ravel().tolist(), label="Gaussian-smoothed normalized TX power")
    fig.suptitle("Full-horizon TX grating-lobe centroids in direction-cosine space")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def steering_matrix(dmat: np.ndarray, u: np.ndarray, v: np.ndarray, w: np.ndarray, valid: np.ndarray):
    """Precompute interferometer steering vectors on the horizon grid."""
    import pansy_config as pc

    uvw = np.column_stack([u[valid], v[valid], w[valid]])
    k0 = 2.0 * np.pi / pc.wavelength
    phase = -1j * k0 * (dmat @ uvw.T)
    return np.exp(phase).astype(np.complex64), uvw


def coherence_map(xc, beam_id, phasecal, ch_pairs, steering, valid, shape):
    """Compute normalized interferometric coherence on the u/v grid."""
    z = np.exp(1j * (np.angle(xc) + phasecal[beam_id, ch_pairs[:, 0]] - phasecal[beam_id, ch_pairs[:, 1]]))
    coh_flat = np.abs(z @ steering) / len(ch_pairs)
    out = np.full(shape, np.nan, dtype=np.float32)
    out[valid] = coh_flat.astype(np.float32)
    return out


def local_coherence_peaks(coh: np.ndarray, threshold: float, max_peaks: int | None = None):
    """Find isolated high-coherence local maxima."""
    filled = np.nan_to_num(coh, nan=-np.inf)
    maxima = ndi.maximum_filter(filled, size=7, mode="constant", cval=-np.inf)
    mask = (filled == maxima) & (filled >= threshold)
    ii, jj = np.where(mask)
    if max_peaks is not None and len(ii) > max_peaks:
        order = np.argsort(filled[ii, jj])[::-1][: int(max_peaks)]
        ii = ii[order]
        jj = jj[order]
    return ii, jj


def selected_pulse_interferometer_response(candidates, obs, phasecal, ch_pairs, steering, valid, shape):
    """Normalized summed interferometer coherence over selected candidate pulses."""
    if obs is None or len(obs.get("xc", [])) == 0 or not candidates:
        return None
    pulses = sorted({int(c["pulse"]) for c in candidates if 0 <= int(c["pulse"]) < len(obs["xc"])})
    if not pulses:
        return None
    summed = np.zeros(shape, dtype=np.float64)
    for pulse in pulses:
        coh = coherence_map(obs["xc"][pulse], int(obs["beam_id"][pulse]), phasecal, ch_pairs, steering, valid, shape)
        summed += np.nan_to_num(coh, nan=0.0)
    maxval = float(np.nanmax(summed)) if np.size(summed) else np.nan
    if not np.isfinite(maxval) or maxval <= 0.0:
        return None
    out = np.full(shape, np.nan, dtype=np.float32)
    out[valid] = (summed[valid] / maxval).astype(np.float32)
    return out


def plot_single_echo(u, v, coh, peak_ij, obs, pulse_idx, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    im = ax.pcolormesh(u, v, coh, shading="auto", cmap="viridis", vmin=0.0, vmax=1.0)
    ax.contour(u, v, coh, levels=[0.9], colors="white", linewidths=0.9)
    if len(peak_ij[0]) > 0:
        ax.scatter(u[peak_ij], v[peak_ij], s=18, facecolors="none", edgecolors="red", linewidths=0.8)
    ax.set_xlabel("u/east direction cosine")
    ax.set_ylabel("v/north direction cosine")
    ax.set_title(
        "Single-pulse interferometer coherence\n"
        f"tx_idx={int(obs['tx_idx'][pulse_idx])}, SNR={10*np.log10(obs['snr'][pulse_idx]):.1f} dB"
    )
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.2)
    fig.colorbar(im, ax=ax, label="Normalized coherence")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_all_candidates(u, v, obs, candidates, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    circle = plt.Circle((0.0, 0.0), 1.0, color="black", fill=False, linewidth=1.0)
    ax.add_patch(circle)
    if candidates:
        t = np.asarray([c["t_rel"] for c in candidates])
        coh = np.asarray([c["coherence"] for c in candidates])
        sc = ax.scatter(
            [c["u"] for c in candidates],
            [c["v"] for c in candidates],
            c=t,
            s=10.0 + 50.0 * np.clip(coh - 0.9, 0.0, 0.1) / 0.1,
            cmap="plasma",
            alpha=0.75,
            edgecolors="none",
        )
        fig.colorbar(sc, ax=ax, label="Time since cut start (s)")
    ax.set_xlabel("u/east direction cosine")
    ax.set_ylabel("v/north direction cosine")
    ax.set_title(r"All local interferometer maxima with coherence $\geq 0.9$")
    ax.set_xlim(-1.03, 1.03)
    ax.set_ylim(-1.03, 1.03)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.2)
    ax.text(
        0.02,
        0.02,
        f"{len(candidates)} candidates from {len(obs['snr'])} pulses",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75},
    )
    fig.savefig(output, dpi=220)
    plt.close(fig)


def split_observations_by_range_time(
    obs,
    range_jump_km=1.0,
    range_rate_km_s=80.0,
    range_margin_km=0.4,
    min_points=3,
):
    """Split pulse-wise meteor observables into simple range-time components."""
    tx_idx = np.asarray(obs.get("tx_idx", []), dtype=np.float64)
    range_km = np.asarray(obs.get("range_km", []), dtype=np.float64)
    if len(tx_idx) != len(range_km) or len(tx_idx) == 0:
        return []

    good = np.isfinite(tx_idx) & np.isfinite(range_km)
    order = np.flatnonzero(good)
    if len(order) == 0:
        return []
    order = order[np.argsort(tx_idx[order], kind="mergesort")]

    segments = []
    start = 0
    for j in range(1, len(order)):
        prev_i = order[j - 1]
        cur_i = order[j]
        dt_s = max(0.0, float((tx_idx[cur_i] - tx_idx[prev_i]) / 1e6))
        allowed_jump = max(float(range_jump_km), float(range_rate_km_s) * dt_s + float(range_margin_km))
        if abs(float(range_km[cur_i] - range_km[prev_i])) > allowed_jump:
            segment = order[start:j]
            if len(segment) >= min_points:
                segments.append(segment)
            start = j

    segment = order[start:]
    if len(segment) >= min_points:
        segments.append(segment)
    return segments


def subset_pulse_observations(obs, idx):
    idx = np.asarray(idx, dtype=np.int64)
    n_pulse = len(obs.get("snr", []))
    return {
        key: val[idx] if isinstance(val, np.ndarray) and len(val) == n_pulse else val
        for key, val in obs.items()
    }


def fit_candidate_tracks(candidates, min_unique_pulses=None, tol=0.035, max_tracks=24, n_trials=30000):
    """Group high-coherence maxima into approximate u/v tracks.

    This is a deliberately simple multi-hypothesis finder for memo diagnostics.
    It fits straight lines in direction-cosine/time space and extracts the
    strongest remaining tracks one at a time.
    """
    if not candidates:
        return []

    rng = np.random.default_rng(20260615)
    t = np.asarray([c["t_rel"] for c in candidates], dtype=np.float64)
    u = np.asarray([c["u"] for c in candidates], dtype=np.float64)
    v = np.asarray([c["v"] for c in candidates], dtype=np.float64)
    pulse = np.asarray([c["pulse"] for c in candidates], dtype=np.int64)
    total_unique_pulses = len(np.unique(pulse))
    if min_unique_pulses is None:
        min_unique_pulses = int(min(total_unique_pulses, 10, max(5, np.ceil(0.75 * total_unique_pulses))))
    duration_s = float(np.nanmax(t) - np.nanmin(t)) if len(t) else 0.0
    min_seed_dt = min(0.20, max(0.005, 0.20 * duration_s))
    active = np.ones(len(candidates), dtype=bool)
    tracks = []

    for _track_i in range(max_tracks):
        idx_active = np.flatnonzero(active)
        if len(idx_active) < min_unique_pulses:
            break

        best = None
        for _ in range(n_trials):
            i, j = rng.choice(idx_active, size=2, replace=False)
            dt = t[j] - t[i]
            if abs(dt) < min_seed_dt:
                continue
            bu = (u[j] - u[i]) / dt
            bv = (v[j] - v[i]) / dt
            au = u[i] - bu * t[i]
            av = v[i] - bv * t[i]
            resid = np.hypot(u[idx_active] - (au + bu * t[idx_active]), v[idx_active] - (av + bv * t[idx_active]))
            inliers = idx_active[resid < tol]
            unique_pulses = len(np.unique(pulse[inliers]))
            if unique_pulses < min_unique_pulses:
                continue
            score = unique_pulses + 0.01 * len(inliers)
            if best is None or score > best["score"]:
                best = {"score": score, "inliers": inliers}

        if best is None:
            break

        inliers = best["inliers"]
        design = np.column_stack([np.ones(len(inliers)), t[inliers]])
        u_coeff = np.linalg.lstsq(design, u[inliers], rcond=None)[0]
        v_coeff = np.linalg.lstsq(design, v[inliers], rcond=None)[0]
        resid = np.hypot(u[inliers] - (u_coeff[0] + u_coeff[1] * t[inliers]), v[inliers] - (v_coeff[0] + v_coeff[1] * t[inliers]))
        keep = inliers[resid < tol]
        if len(np.unique(pulse[keep])) < min_unique_pulses:
            active[inliers] = False
            continue

        tracks.append(
            {
                "hypothesis_id": len(tracks) + 1,
                "idx": keep,
                "u_coeff": u_coeff,
                "v_coeff": v_coeff,
                "n": int(len(keep)),
                "unique_pulses": int(len(np.unique(pulse[keep]))),
                "min_unique_pulses": int(min_unique_pulses),
            }
        )
        active[keep] = False

    return tracks


def complete_tracks_to_common_pulses(tracks, candidates, max_completion_distance_dc=0.04, min_completion_fraction=0.65):
    """Assign one candidate from every available pulse to every hypothesis.

    The RANSAC-style finder is allowed to discover a hypothesis from a subset of
    pulses, but model ranking must compare hypotheses on the same measurements.
    This completion step evaluates the fitted direction-cosine line at each
    pulse and assigns the nearest high-coherence peak from that pulse.
    """
    if not tracks or not candidates:
        return tracks

    pulse_values = np.asarray([c["pulse"] for c in candidates], dtype=np.int64)
    t_values = np.asarray([c["t_rel"] for c in candidates], dtype=np.float64)
    u_values = np.asarray([c["u"] for c in candidates], dtype=np.float64)
    v_values = np.asarray([c["v"] for c in candidates], dtype=np.float64)
    common_pulses = np.unique(pulse_values)
    if len(common_pulses) == 0:
        return tracks

    for track in tracks:
        seed_idx = np.asarray(track.get("idx", []), dtype=np.int64)
        seed_idx = seed_idx[(seed_idx >= 0) & (seed_idx < len(candidates))]
        seed_pulses = set(int(pulse_values[i]) for i in seed_idx)
        completed_idx = []
        completed_distance = []
        missed_pulses = []
        for pulse in common_pulses:
            pulse_idx = np.flatnonzero(pulse_values == pulse)
            if len(pulse_idx) == 0:
                continue
            # All candidates from a pulse share the same time.
            t_pulse = float(t_values[pulse_idx[0]])
            pred_u = float(track["u_coeff"][0] + track["u_coeff"][1] * t_pulse)
            pred_v = float(track["v_coeff"][0] + track["v_coeff"][1] * t_pulse)
            dist = np.hypot(u_values[pulse_idx] - pred_u, v_values[pulse_idx] - pred_v)
            best_local = int(pulse_idx[int(np.argmin(dist))])
            best_distance = float(np.min(dist))
            if pulse not in seed_pulses and best_distance > max_completion_distance_dc:
                missed_pulses.append(int(pulse))
                continue
            completed_idx.append(best_local)
            completed_distance.append(best_distance)

        if completed_idx:
            track["seed_idx"] = np.asarray(track["idx"], dtype=np.int64)
            track["seed_unique_pulses"] = int(track.get("unique_pulses", 0))
            track["idx"] = np.asarray(completed_idx, dtype=np.int64)
            track["n"] = int(len(completed_idx))
            track["unique_pulses"] = int(len(np.unique(pulse_values[track["idx"]])))
            track["common_pulse_count"] = int(len(common_pulses))
            track["common_pulse_completion_distance_dc"] = np.asarray(completed_distance, dtype=np.float64)
            track["common_pulse_completion_missed_pulses"] = np.asarray(missed_pulses, dtype=np.int64)
            track["common_pulse_completion_fraction"] = float(track["unique_pulses"] / max(1, len(common_pulses)))
            track["common_pulse_completion_rms_dc"] = float(np.sqrt(np.mean(np.asarray(completed_distance) ** 2)))
            track["common_pulse_completion_max_dc"] = float(np.max(completed_distance))
            track["common_pulse_completion_reject"] = bool(
                track["common_pulse_completion_fraction"] < float(min_completion_fraction)
                or track["common_pulse_completion_rms_dc"] > max_completion_distance_dc
            )
    return tracks


def classify_track_visibility(track, candidates, t_start, t_end, min_elevation_deg=5.0):
    """Reject tracks whose fitted 3D path is below or too close to the horizon."""
    rows = [candidates[i] for i in track["idx"]]
    t = np.asarray([r["t_rel"] for r in rows], dtype=np.float64)
    u = np.asarray([r["u"] for r in rows], dtype=np.float64)
    v = np.asarray([r["v"] for r in rows], dtype=np.float64)
    rg = np.asarray([r["range_km"] for r in rows], dtype=np.float64)
    up = np.sqrt(np.maximum(0.0, 1.0 - u**2 - v**2))
    points = np.column_stack([rg * u, rg * v, rg * up])

    design = np.column_stack([np.ones(len(t)), t])
    coeff = np.linalg.lstsq(design, points, rcond=None)[0]
    td = np.linspace(t_start, t_end, 400)
    model = np.column_stack([
        coeff[0, 0] + coeff[1, 0] * td,
        coeff[0, 1] + coeff[1, 1] * td,
        coeff[0, 2] + coeff[1, 2] * td,
    ])
    radius = np.linalg.norm(model, axis=1)
    direction = model / np.maximum(radius[:, None], 1e-9)
    uv_radius = np.hypot(direction[:, 0], direction[:, 1])
    angle_step = np.arccos(np.clip(np.sum(direction[:-1] * direction[1:], axis=1), -1.0, 1.0))

    min_up_direction = float(np.nanmin(direction[:, 2])) if len(direction) else np.nan
    below_horizon = bool(np.any(model[:, 2] <= 0.0))
    near_horizon = bool(np.isfinite(min_up_direction) and min_up_direction < np.sin(np.deg2rad(min_elevation_deg)))
    wraps = bool(np.min(radius) < 20.0 or np.nanmax(angle_step) > np.deg2rad(20.0) or np.nanmax(uv_radius) > 1.0)
    reason = "kept"
    if bool(track.get("common_pulse_completion_reject", False)):
        reason = "completion mismatch"
    elif below_horizon:
        reason = "below horizon"
    elif near_horizon:
        reason = "near horizon"
    elif wraps:
        reason = "wraps"

    track.update(
        {
            "fit_points": points,
            "fit_t": t,
            "dense_t": td,
            "dense_model": model,
            "dense_direction": direction,
            "below_horizon": below_horizon,
            "near_horizon": near_horizon,
            "min_elevation_deg": float(np.rad2deg(np.arcsin(np.clip(min_up_direction, -1.0, 1.0)))) if np.isfinite(min_up_direction) else np.nan,
            "min_elevation_threshold_deg": float(min_elevation_deg),
            "wraps": wraps,
            "reason": reason,
            "min_up_km": float(np.min(model[:, 2])),
            "min_range_km": float(np.min(radius)),
            "max_angle_step_deg": float(np.rad2deg(np.nanmax(angle_step))) if len(angle_step) else 0.0,
        }
    )
    return track


def classify_track_linearity(
    track,
    candidates,
    angular_sigma_deg=0.25,
    range_sigma_km=0.15,
    p_threshold=0.01,
    rms_threshold_km=0.6,
    max_threshold_km=2.0,
):
    """Test whether candidate points are statistically consistent with one 3D line."""
    rows = [candidates[i] for i in track["idx"]]
    u = np.asarray([r["u"] for r in rows], dtype=np.float64)
    v = np.asarray([r["v"] for r in rows], dtype=np.float64)
    rg = np.asarray([r["range_km"] for r in rows], dtype=np.float64)
    snr = np.asarray([r["snr"] for r in rows], dtype=np.float64)
    up = np.sqrt(np.maximum(0.0, 1.0 - u**2 - v**2))
    points = np.column_stack([rg * u, rg * v, rg * up])

    center = np.mean(points, axis=0)
    _, _, vh = np.linalg.svd(points - center, full_matrices=False)
    direction = vh[0]
    along = (points - center) @ direction
    model = center + np.outer(along, direction)
    residual = points - model
    perp = np.linalg.norm(residual, axis=1)

    snr_db = 10.0 * np.log10(np.maximum(snr, 1e-12))
    empirical_angular_sigma_deg = np.clip(0.85 / np.sqrt(np.maximum(snr, 1.0)), 0.08, 0.45)
    angular_sigma_rad = np.deg2rad(np.maximum(angular_sigma_deg, empirical_angular_sigma_deg))
    sigma_km = np.sqrt((rg * angular_sigma_rad) ** 2 + range_sigma_km**2)
    chi2 = float(np.sum((perp / np.maximum(sigma_km, 1e-6)) ** 2))
    dof = max(1, 2 * len(points) - 4)
    p_value = float(st.chi2.sf(chi2, dof))
    reduced_chi2 = chi2 / dof
    line_rms_km = float(np.sqrt(np.mean(perp**2)))
    line_max_km = float(np.max(perp))
    path_length_km = float(np.nanmax(along) - np.nanmin(along)) if len(along) else np.nan
    line_rms_fraction = float(line_rms_km / max(path_length_km, 1e-6)) if np.isfinite(path_length_km) else np.nan
    line_length_adjusted_redchi = float(reduced_chi2 / max(path_length_km / 10.0, 1.0)) if np.isfinite(path_length_km) else float(reduced_chi2)
    statistical_reject = bool(p_value < p_threshold)
    rms_reject = bool(line_rms_km > rms_threshold_km)
    max_reject = bool(line_max_km > max_threshold_km)
    linearity_reject = bool(statistical_reject or rms_reject or max_reject)
    reasons = []
    if statistical_reject:
        reasons.append("chi2")
    if rms_reject:
        reasons.append("rms")
    if max_reject:
        reasons.append("max")

    track.update(
        {
            "line_center": center,
            "line_direction": direction,
            "line_model_points": model,
            "line_perp_km": perp,
            "line_sigma_km": sigma_km,
            "line_empirical_angular_sigma_deg": empirical_angular_sigma_deg,
            "line_snr_db": snr_db,
            "line_path_length_km": path_length_km,
            "line_rms_fraction": line_rms_fraction,
            "line_length_adjusted_reduced_chi2": line_length_adjusted_redchi,
            "detection_min_altitude_km": float(np.min(points[:, 2])),
            "detection_max_altitude_km": float(np.max(points[:, 2])),
            "low_detection_altitude_reject": bool(np.min(points[:, 2]) < 25.0),
            "line_chi2": chi2,
            "line_dof": dof,
            "line_p_value": p_value,
            "line_reduced_chi2": reduced_chi2,
            "linearity_reject": linearity_reject,
            "linearity_reason": ",".join(reasons) if reasons else "linear",
            "line_rms_km": line_rms_km,
            "line_max_km": line_max_km,
            "line_rms_threshold_km": rms_threshold_km,
            "line_max_threshold_km": max_threshold_km,
            "line_median_sigma_km": float(np.median(sigma_km)),
        }
    )
    return track


def classify_track_descent(track):
    """Meteors should move downward in local altitude during the observed path."""
    t = np.asarray(track["fit_t"], dtype=np.float64)
    up = np.asarray(track["fit_points"][:, 2], dtype=np.float64)
    if len(t) < 2:
        slope = np.nan
    else:
        design = np.column_stack([np.ones(len(t)), t - np.mean(t)])
        slope = float(np.linalg.lstsq(design, up, rcond=None)[0][1])
    track["descent_rate_km_s"] = slope
    track["descent_reject"] = bool(np.isfinite(slope) and slope >= 0.0)
    return track


def robust_std(x, floor=1e-3):
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return floor
    med = np.median(x)
    sig = 1.4826 * np.median(np.abs(x - med))
    if not np.isfinite(sig) or sig < floor:
        sig = np.std(x)
    return float(max(sig, floor))


def propagate_drag_model(params, t, rho_of_alt_m):
    """Propagate an ENU state with MSIS drag and return km / km/s arrays."""
    t = np.asarray(t, dtype=np.float64)
    t_ref = float(np.min(t))
    tau = t - t_ref
    order = np.argsort(tau)
    pos_m, vel_mps, beta = pbal.propagate(params, tau[order], rho_of_alt_m)
    pos = np.empty_like(pos_m)
    vel = np.empty_like(vel_mps)
    pos[order] = pos_m
    vel[order] = vel_mps
    return pos / 1e3, vel / 1e3, beta, t_ref


def ballistic_residuals(params, t, points, doppler_km_s, rho_of_alt_m, sigma_pos_km, sigma_dop_km_s):
    model, vel, _beta, _t_ref = propagate_drag_model(params, t, rho_of_alt_m)
    rng = np.linalg.norm(model, axis=1)
    pred_dop = np.sum(model * vel, axis=1) / np.maximum(rng, 1e-9)
    pos_res = ((points - model) / sigma_pos_km).ravel()
    dop_res = (doppler_km_s - pred_dop) / sigma_dop_km_s
    return np.concatenate([pos_res, dop_res])


def model_pulse_pair_phase(params, t_now, t_prev, rho_of_alt_m):
    t_eval = np.concatenate([t_prev, t_now])
    model, _vel, _beta, _t_ref = propagate_drag_model(params, t_eval, rho_of_alt_m)
    n_pair = len(t_now)
    prev = model[:n_pair]
    cur = model[n_pair:]
    range_delta_m = (np.linalg.norm(prev, axis=1) - np.linalg.norm(cur, axis=1)) * 1e3
    return np.exp(-1j * 4.0 * np.pi * range_delta_m / pc.wavelength)


def ballistic_residuals_with_phase(
    params,
    t,
    points,
    doppler_km_s,
    zpp,
    zpp_prev_t,
    rho_of_alt_m,
    sigma_pos_km,
    sigma_dop_km_s,
    sigma_phase_rad,
    keep=None,
):
    if keep is None:
        keep = np.ones(len(t), dtype=bool)
    base = ballistic_residuals(
        params,
        t[keep],
        points[keep],
        doppler_km_s[keep],
        rho_of_alt_m,
        sigma_pos_km,
        sigma_dop_km_s,
    )
    phase_keep = (
        keep
        & np.isfinite(zpp.real)
        & np.isfinite(zpp.imag)
        & np.isfinite(zpp_prev_t)
    )
    if np.count_nonzero(phase_keep) == 0:
        return base
    model_zpp = model_pulse_pair_phase(params, t[phase_keep], zpp_prev_t[phase_keep], rho_of_alt_m)
    phase_res = circular_phase_residual(np.angle(zpp[phase_keep]), np.angle(model_zpp)) / sigma_phase_rad
    return np.concatenate([base, phase_res])


def apply_final_phase_refit(
    track,
    rows,
    t,
    points,
    doppler_km_s,
    rho_of_alt_m,
    bounds,
    x_scale,
    sigma_pos_km,
    sigma_dop_km_s,
):
    zpp = np.asarray([r.get("zpp", np.nan + 1j * np.nan) for r in rows], dtype=np.complex64)
    zpp_prev_t = np.asarray([r.get("zpp_prev_t_rel", np.nan) for r in rows], dtype=np.float64)
    zpp_doppler_mps = np.asarray([r.get("zpp_doppler_mps", np.nan) for r in rows], dtype=np.float64)
    zpp_coarse_doppler_mps = np.asarray([r.get("zpp_coarse_doppler_mps", np.nan) for r in rows], dtype=np.float64)
    phase_valid = (
        track["ballistic_keep"]
        & np.isfinite(zpp.real)
        & np.isfinite(zpp.imag)
        & np.isfinite(zpp_prev_t)
    )
    if np.count_nonzero(phase_valid) < 3:
        track["ballistic_phase_fit_available"] = False
        track["ballistic_phase_n"] = int(np.count_nonzero(phase_valid))
        return track

    initial_model_zpp = model_pulse_pair_phase(track["ballistic_params"], t[phase_valid], zpp_prev_t[phase_valid], rho_of_alt_m)
    initial_phase_res = circular_phase_residual(np.angle(zpp[phase_valid]), np.angle(initial_model_zpp))
    sigma_phase_rad = robust_std(initial_phase_res, floor=0.10)
    sigma_phase_rad = float(np.clip(sigma_phase_rad, 0.10, 1.0))

    result = opt.least_squares(
        lambda p: ballistic_residuals_with_phase(
            p,
            t,
            points,
            doppler_km_s,
            zpp,
            zpp_prev_t,
            rho_of_alt_m,
            sigma_pos_km,
            sigma_dop_km_s,
            sigma_phase_rad,
            keep=track["ballistic_keep"],
        ),
        track["ballistic_params"],
        bounds=bounds,
        x_scale=x_scale,
        loss="soft_l1",
        f_scale=1.0,
        max_nfev=350,
    )

    model, vel, beta, t_ref = propagate_drag_model(result.x, t, rho_of_alt_m)
    rng = np.linalg.norm(model, axis=1)
    pred_dop = np.sum(model * vel, axis=1) / np.maximum(rng, 1e-9)
    pos_res = points - model
    dop_res = doppler_km_s - pred_dop
    final_model_zpp = model_pulse_pair_phase(result.x, t[phase_valid], zpp_prev_t[phase_valid], rho_of_alt_m)
    final_phase_res = circular_phase_residual(np.angle(zpp[phase_valid]), np.angle(final_model_zpp))
    phase_dt_s = t[phase_valid] - zpp_prev_t[phase_valid]
    model_phase_rad = np.angle(final_model_zpp)
    model_zpp_doppler_mps = -model_phase_rad * pc.wavelength / (4.0 * np.pi * phase_dt_s)
    phase_order = np.argsort(t[phase_valid])
    model_phase_unwrapped = np.full(np.count_nonzero(phase_valid), np.nan, dtype=np.float64)
    observed_phase_model_branch = np.full(np.count_nonzero(phase_valid), np.nan, dtype=np.float64)
    model_unwrapped_ordered = np.unwrap(np.angle(final_model_zpp)[phase_order])
    observed_ordered = np.angle(zpp[phase_valid])[phase_order]
    observed_branch_ordered = observed_ordered + 2.0 * np.pi * np.round((model_unwrapped_ordered - observed_ordered) / (2.0 * np.pi))
    model_phase_unwrapped[phase_order] = model_unwrapped_ordered
    observed_phase_model_branch[phase_order] = observed_branch_ordered

    pre_params = track["ballistic_params"]
    track.update(
        {
            "ballistic_pre_phase_params": pre_params,
            "ballistic_phase_fit_available": True,
            "ballistic_phase_n": int(np.count_nonzero(phase_valid)),
            "ballistic_phase_sigma_rad": sigma_phase_rad,
            "ballistic_phase_t": t[phase_valid],
            "ballistic_phase_prev_t": zpp_prev_t[phase_valid],
            "ballistic_phase_observed_rad": np.angle(zpp[phase_valid]),
            "ballistic_phase_model_rad": model_phase_rad,
            "ballistic_phase_doppler_mps": zpp_doppler_mps[phase_valid],
            "ballistic_phase_coarse_doppler_mps": zpp_coarse_doppler_mps[phase_valid],
            "ballistic_phase_model_doppler_mps": model_zpp_doppler_mps,
            "ballistic_phase_model_unwrapped_rad": model_phase_unwrapped,
            "ballistic_phase_observed_model_branch_rad": observed_phase_model_branch,
            "ballistic_phase_initial_residual_rad": initial_phase_res,
            "ballistic_phase_residual_rad": final_phase_res,
            "ballistic_phase_rms_rad": float(np.sqrt(np.mean(final_phase_res**2))),
            "ballistic_params": result.x,
            "ballistic_t_ref_s": t_ref,
            "ballistic_coefficient_kg_m2": beta,
            "ballistic_model": model,
            "ballistic_velocity_km_s": vel,
            "ballistic_pred_doppler_km_s": pred_dop,
            "ballistic_pos_res_km": pos_res,
            "ballistic_dop_res_km_s": dop_res,
            "ballistic_pos_rms_km": float(np.sqrt(np.mean(np.sum(pos_res[track["ballistic_keep"]] ** 2, axis=1)))),
            "ballistic_dop_rms_km_s": float(np.sqrt(np.mean(dop_res[track["ballistic_keep"]] ** 2))),
        }
    )
    return track


def drag_initial_guess(t, points, log10_beta=0.0):
    t = np.asarray(t, dtype=np.float64)
    tau = t - np.min(t)
    design = np.column_stack([np.ones(len(tau)), tau, 0.5 * tau**2])
    coeff = np.linalg.lstsq(design, points * 1e3, rcond=None)[0]
    pos0 = coeff[0]
    vel0 = coeff[1]
    acc0 = coeff[2]
    # If the quadratic fit is poorly conditioned, fall back to the linear term.
    if not np.all(np.isfinite(pos0)) or not np.all(np.isfinite(vel0)) or not np.all(np.isfinite(acc0)):
        design = np.column_stack([np.ones(len(tau)), tau])
        coeff = np.linalg.lstsq(design, points * 1e3, rcond=None)[0]
        pos0 = coeff[0]
        vel0 = coeff[1]
    pos0[2] = np.clip(pos0[2], 30e3, 220e3)
    vel0 = np.clip(vel0, -90e3, 90e3)
    return np.concatenate([pos0, vel0, [float(log10_beta)]])


def fit_ballistic_track(track, candidates, rho_of_alt_m, sigma_pos_km=0.5, sigma_dop_km_s=1.0, clip_sigma=3.5):
    rows = [candidates[i] for i in track["idx"]]
    t = np.asarray([r["t_rel"] for r in rows], dtype=np.float64)
    doppler_km_s = np.asarray([r["doppler_mps"] for r in rows], dtype=np.float64) / 1e3
    snr_db = 10.0 * np.log10(np.maximum(np.asarray([r["snr"] for r in rows], dtype=np.float64), 1e-12))
    points = np.asarray(track["fit_points"], dtype=np.float64)
    p0 = drag_initial_guess(t, points, log10_beta=0.0)
    bounds = (
        np.array([-np.inf, -np.inf, 20e3, -90e3, -90e3, -90e3, -4.0]),
        np.array([np.inf, np.inf, 220e3, 90e3, 90e3, 90e3, 6.0]),
    )
    x_scale = np.array([1e5, 1e5, 1e5, 7e4, 7e4, 7e4, 1.0])

    def fun(p, keep=None):
        if keep is None:
            keep = np.ones(len(t), dtype=bool)
        return ballistic_residuals(
            p,
            t[keep],
            points[keep],
            doppler_km_s[keep],
            rho_of_alt_m,
            sigma_pos_km,
            sigma_dop_km_s,
        )

    keep = np.ones(len(t), dtype=bool)
    result = opt.least_squares(
        fun,
        p0,
        bounds=bounds,
        x_scale=x_scale,
        loss="soft_l1",
        f_scale=1.0,
        max_nfev=350,
    )

    for _ in range(2):
        tfit = t[keep]
        pfit = points[keep]
        dfit = doppler_km_s[keep]
        result = opt.least_squares(
            lambda p: ballistic_residuals(p, tfit, pfit, dfit, rho_of_alt_m, sigma_pos_km, sigma_dop_km_s),
            result.x,
            bounds=bounds,
            x_scale=x_scale,
            loss="soft_l1",
            f_scale=1.0,
            max_nfev=350,
        )
        model, vel, beta, t_ref = propagate_drag_model(result.x, t, rho_of_alt_m)
        rng = np.linalg.norm(model, axis=1)
        pred_dop = np.sum(model * vel, axis=1) / np.maximum(rng, 1e-9)
        pos_res = points - model
        dop_res = doppler_km_s - pred_dop
        point_norm = np.sqrt(np.sum((pos_res / sigma_pos_km) ** 2, axis=1) + (dop_res / sigma_dop_km_s) ** 2)
        new_keep = point_norm < clip_sigma
        if np.count_nonzero(new_keep) < 10 or np.array_equal(new_keep, keep):
            keep = new_keep if np.count_nonzero(new_keep) >= 10 else keep
            break
        keep = new_keep

    model, vel, beta, t_ref = propagate_drag_model(result.x, t, rho_of_alt_m)
    rng = np.linalg.norm(model, axis=1)
    pred_dop = np.sum(model * vel, axis=1) / np.maximum(rng, 1e-9)
    pos_res = points - model
    dop_res = doppler_km_s - pred_dop
    pos_res_keep = pos_res[keep]
    dop_res_keep = dop_res[keep]
    chi2 = float(np.sum((pos_res_keep / sigma_pos_km) ** 2) + np.sum((dop_res_keep / sigma_dop_km_s) ** 2))
    dof = max(1, 4 * np.count_nonzero(keep) - 7)
    p_value = float(st.chi2.sf(chi2, dof))
    param_cov, param_std, cov_ok, cov_dof, cov_res_var = pbal.covariance_from_lsq(result, 4 * np.count_nonzero(keep))

    speed_km_s = np.linalg.norm(vel, axis=1)
    order = np.argsort(t)
    speed_t = t[order][keep[order]]
    speed_fit = speed_km_s[order][keep[order]]
    if len(speed_fit) >= 3:
        speed_design = np.column_stack([np.ones(len(speed_t)), speed_t - np.mean(speed_t)])
        speed_slope_km_s2 = float(np.linalg.lstsq(speed_design, speed_fit, rcond=None)[0][1])
        speed_delta_km_s = float(speed_fit[-1] - speed_fit[0])
    elif len(speed_fit) >= 2:
        dt = max(float(speed_t[-1] - speed_t[0]), 1e-9)
        speed_delta_km_s = float(speed_fit[-1] - speed_fit[0])
        speed_slope_km_s2 = speed_delta_km_s / dt
    else:
        speed_delta_km_s = np.nan
        speed_slope_km_s2 = np.nan
    speedup_flag = bool(np.isfinite(speed_slope_km_s2) and speed_slope_km_s2 > 0.0)
    start_altitude_km = np.nan
    if np.any(keep):
        first_keep = np.flatnonzero(keep)[np.argmin(t[keep])]
        start_altitude_km = float(model[first_keep, 2])
    low_start_altitude_reject = bool(np.isfinite(start_altitude_km) and start_altitude_km < 50.0)

    out = {
        "ballistic_params": result.x,
        "ballistic_param_names": np.asarray(["east_m", "north_m", "up_m", "ve_mps", "vn_mps", "vu_mps", "log10_beta_kg_m2"]),
        "ballistic_parameter_covariance": param_cov,
        "ballistic_parameter_std": param_std,
        "ballistic_covariance_available": bool(cov_ok),
        "ballistic_covariance_dof": int(cov_dof),
        "ballistic_covariance_residual_variance": float(cov_res_var),
        "ballistic_t_ref_s": t_ref,
        "ballistic_initial_guess": "quadratic_polynomial_position_time",
        "ballistic_coefficient_kg_m2": beta,
        "ballistic_keep": keep,
        "ballistic_model": model,
        "ballistic_velocity_km_s": vel,
        "ballistic_start_altitude_km": start_altitude_km,
        "ballistic_low_start_altitude_reject": low_start_altitude_reject,
        "ballistic_speed_km_s": speed_km_s,
        "ballistic_speed_slope_km_s2": speed_slope_km_s2,
        "ballistic_speed_delta_km_s": speed_delta_km_s,
        "ballistic_speedup_flag": speedup_flag,
        "ballistic_pred_doppler_km_s": pred_dop,
        "ballistic_pos_res_km": pos_res,
        "ballistic_dop_res_km_s": dop_res,
        "ballistic_pos_rms_km": float(np.sqrt(np.mean(np.sum(pos_res_keep**2, axis=1)))),
        "ballistic_dop_rms_km_s": float(np.sqrt(np.mean(dop_res_keep**2))),
        "ballistic_chi2": chi2,
        "ballistic_dof": dof,
        "ballistic_reduced_chi2": chi2 / dof,
        "ballistic_p_value": p_value,
        "ballistic_n": int(np.count_nonzero(keep)),
        "ballistic_n_outliers": int(len(keep) - np.count_nonzero(keep)),
        "ballistic_t": t,
        "ballistic_doppler_km_s": doppler_km_s,
        "ballistic_snr_db": snr_db,
        "ballistic_sigma_pos_km": sigma_pos_km,
        "ballistic_sigma_dop_km_s": sigma_dop_km_s,
    }
    zpp = np.asarray([r.get("zpp", np.nan + 1j * np.nan) for r in rows], dtype=np.complex64)
    zpp_prev_t = np.asarray([r.get("zpp_prev_t_rel", np.nan) for r in rows], dtype=np.float64)
    phase_valid = (
        out["ballistic_keep"]
        & np.isfinite(zpp.real)
        & np.isfinite(zpp.imag)
        & np.isfinite(zpp_prev_t)
    )
    out["ballistic_phase_fit_available"] = False
    out["ballistic_phase_n"] = int(np.count_nonzero(phase_valid))
    if np.count_nonzero(phase_valid):
        model_zpp = model_pulse_pair_phase(out["ballistic_params"], t[phase_valid], zpp_prev_t[phase_valid], rho_of_alt_m)
        phase_res = circular_phase_residual(np.angle(zpp[phase_valid]), np.angle(model_zpp))
        phase_dt_s = t[phase_valid] - zpp_prev_t[phase_valid]
        ipp_s = pmm.get_m_mode()["ipp_us"] * 1e-6
        phase_ipp = phase_dt_s / ipp_s
        out.update(
            {
                "ballistic_phase_t": t[phase_valid],
                "ballistic_phase_prev_t": zpp_prev_t[phase_valid],
                "ballistic_phase_dt_s": phase_dt_s,
                "ballistic_phase_ipp": phase_ipp,
                "ballistic_phase_observed_rad": np.angle(zpp[phase_valid]),
                "ballistic_phase_model_rad": np.angle(model_zpp),
                "ballistic_phase_residual_rad": phase_res,
                "ballistic_phase_rms_rad": float(np.sqrt(np.mean(phase_res**2))),
            }
        )
    return out


def fit_ballistic_survivors(tracks, candidates, event_epoch_unix, p_threshold=0.01):
    survivors = [
        t
        for t in tracks
        if t["reason"] == "kept"
        and "fit_points" in t
    ]
    if not survivors:
        return np.nan, np.nan
    fallback_level = "horizon_clear"

    required_pulses = frozenset(int(c["pulse"]) for c in candidates)
    required_unique_pulses = len(required_pulses)
    if required_unique_pulses == 0:
        required_unique_pulses = max(int(track.get("unique_pulses", 0)) for track in survivors)
    rho_of_alt_m, msis_meta = pbal.density_interpolator(event_epoch_unix)
    preliminary = [
        fit_ballistic_track(track, candidates, rho_of_alt_m, sigma_pos_km=0.5, sigma_dop_km_s=1.0)
        for track in survivors
    ]
    best_i = int(np.argmin([p["ballistic_pos_rms_km"] / 0.5 + p["ballistic_dop_rms_km_s"] / 1.0 for p in preliminary]))
    best = preliminary[best_i]
    best_keep = best["ballistic_keep"]
    sigma_pos = float(np.sqrt(np.mean(best["ballistic_pos_res_km"][best_keep] ** 2)))
    sigma_dop = float(np.sqrt(np.mean(best["ballistic_dop_res_km_s"][best_keep] ** 2)))
    sigma_pos = max(sigma_pos, 0.05)
    sigma_dop = max(sigma_dop, 0.05)

    for track in survivors:
        fit = fit_ballistic_track(track, candidates, rho_of_alt_m, sigma_pos_km=sigma_pos, sigma_dop_km_s=sigma_dop)
        track_pulses = frozenset(int(candidates[i]["pulse"]) for i in track["idx"])
        unique_pulses = int(track.get("unique_pulses", 0))
        fit["ballistic_required_unique_pulses"] = int(required_unique_pulses)
        fit["ballistic_pulse_coverage_reject"] = bool(
            unique_pulses != required_unique_pulses
            or (required_pulses and track_pulses != required_pulses)
        )
        fit["ballistic_reject"] = bool(
            fit["ballistic_p_value"] < p_threshold
            or fit["ballistic_low_start_altitude_reject"]
            or fit["ballistic_pulse_coverage_reject"]
        )
        fit["ballistic_model_type"] = "msis_drag"
        fit["ballistic_selection_level"] = fallback_level
        fit["msis_alt_grid_km"] = msis_meta["msis_alt_grid_km"]
        fit["msis_density_kg_m3"] = msis_meta["msis_density_kg_m3"]
        track.update(fit)

    survivors.sort(key=lambda t: (t["ballistic_reduced_chi2"], t["ballistic_dop_rms_km_s"], t["ballistic_pos_rms_km"]))
    for rank, track in enumerate(survivors):
        track["ballistic_rank"] = rank
    return sigma_pos, sigma_dop


def fit_fixed_velocity_survivors(tracks, candidates, rho_of_alt_m, sigma_pos_km, sigma_dop_km_s, fit_fixed_am=False):
    """Fit a constant-velocity range/Doppler trajectory for alias scoring."""
    try:
        import fit_best_alias_physics_models as physics
    except Exception as exc:
        for track in tracks:
            if "ballistic_t" in track:
                track["fixed_velocity_fit_error"] = f"{type(exc).__name__}: {exc}"
        return

    for track in tracks:
        if "fit_points" not in track:
            continue
        try:
            rows = [candidates[i] for i in track["idx"]]
            t_s = np.asarray(track.get("ballistic_t", track.get("fit_t")), dtype=np.float64)
            points_km = np.asarray(track["fit_points"], dtype=np.float64)
            doppler_km_s = np.asarray(
                track.get("ballistic_doppler_km_s", [r["doppler_mps"] / 1e3 for r in rows]),
                dtype=np.float64,
            )
            snr_db = np.asarray(
                track.get("ballistic_snr_db", [10.0 * np.log10(max(r["snr"], 1e-6)) for r in rows]),
                dtype=np.float64,
            )
            track["ballistic_t"] = t_s
            track["ballistic_doppler_km_s"] = doppler_km_s
            track["ballistic_snr_db"] = snr_db
            track["ballistic_required_unique_pulses"] = int(track.get("common_pulse_count", track.get("unique_pulses", 0)))
            track["ballistic_pulse_coverage_reject"] = False
            track["ballistic_low_start_altitude_reject"] = False
            track["ballistic_reject"] = False
            p0 = physics.linear_initial_guess(t_s, points_km)
            fit = physics.fit_model(
                "fixed_velocity",
                t_s,
                points_km,
                doppler_km_s,
                rho_of_alt_m,
                sigma_pos_km,
                sigma_dop_km_s,
                [p0],
            )
        except Exception as exc:
            track["fixed_velocity_fit_error"] = f"{type(exc).__name__}: {exc}"
            continue
        keep = np.asarray(fit["keep_mask"], dtype=bool)
        track.update(
            {
                "fixed_velocity_params": np.asarray(fit["params"], dtype=np.float64),
                "fixed_velocity_parameter_covariance": np.asarray(fit["parameter_covariance"], dtype=np.float64),
                "fixed_velocity_parameter_std": np.asarray(fit["parameter_std"], dtype=np.float64),
                "fixed_velocity_covariance_available": bool(fit["covariance_available"]),
                "fixed_velocity_keep": keep,
                "fixed_velocity_model": np.asarray(fit["model_enu_km"], dtype=np.float64),
                "fixed_velocity_velocity_km_s": np.asarray(fit["velocity_km_s"], dtype=np.float64),
                "fixed_velocity_pred_doppler_km_s": np.asarray(fit["pred_doppler_km_s"], dtype=np.float64),
                "fixed_velocity_pos_res_km": np.asarray(fit["pos_res_km"], dtype=np.float64),
                "fixed_velocity_dop_res_km_s": np.asarray(fit["dop_res_km_s"], dtype=np.float64),
                "fixed_velocity_chi2": float(fit["chi2"]),
                "fixed_velocity_dof": int(fit["dof"]),
                "fixed_velocity_reduced_chi2": float(fit["reduced_chi2"]),
                "fixed_velocity_bic": float(fit["bic"]),
                "fixed_velocity_pos_rms_km": float(fit["pos_rms_km"]),
                "fixed_velocity_dop_rms_km_s": float(fit["dop_rms_km_s"]),
                "fixed_velocity_n": int(fit["n_points"]),
                "fixed_velocity_n_outliers": int(len(keep) - np.count_nonzero(keep)),
                "fixed_velocity_t_ref_s": float(np.min(t_s)) if len(t_s) else 0.0,
            }
        )
        if not fit_fixed_am:
            plausibility_redchi = float(track.get("fixed_velocity_reduced_chi2", np.inf))
            track["interstellar_alias_plausible"] = bool(np.isfinite(plausibility_redchi) and plausibility_redchi < 1.0)
            track["interstellar_alias_plausibility_redchi_threshold"] = 1.0
            track["interstellar_alias_plausibility_redchi"] = plausibility_redchi
            track["interstellar_alias_plausibility_model"] = "fixed_velocity" if track["interstellar_alias_plausible"] else ""
            continue
        try:
            fixed_am_starts = [
                np.concatenate([fit["params"][:6], [log_b]])
                for log_b in (-4.0, -2.0, 0.0, 1.0, 2.0)
            ]
            fit_am = physics.fit_model(
                "fixed_am",
                t_s,
                points_km,
                doppler_km_s,
                rho_of_alt_m,
                sigma_pos_km,
                sigma_dop_km_s,
                fixed_am_starts,
            )
        except Exception as exc:
            track["fixed_am_fit_error"] = f"{type(exc).__name__}: {exc}"
            fit_am = None
        if fit_am is not None:
            keep_am = np.asarray(fit_am["keep_mask"], dtype=bool)
            track.update(
                {
                    "fixed_am_params": np.asarray(fit_am["params"], dtype=np.float64),
                    "fixed_am_parameter_covariance": np.asarray(fit_am["parameter_covariance"], dtype=np.float64),
                    "fixed_am_parameter_std": np.asarray(fit_am["parameter_std"], dtype=np.float64),
                    "fixed_am_covariance_available": bool(fit_am["covariance_available"]),
                    "fixed_am_keep": keep_am,
                    "fixed_am_model": np.asarray(fit_am["model_enu_km"], dtype=np.float64),
                    "fixed_am_velocity_km_s": np.asarray(fit_am["velocity_km_s"], dtype=np.float64),
                    "fixed_am_pred_doppler_km_s": np.asarray(fit_am["pred_doppler_km_s"], dtype=np.float64),
                    "fixed_am_pos_res_km": np.asarray(fit_am["pos_res_km"], dtype=np.float64),
                    "fixed_am_dop_res_km_s": np.asarray(fit_am["dop_res_km_s"], dtype=np.float64),
                    "fixed_am_chi2": float(fit_am["chi2"]),
                    "fixed_am_dof": int(fit_am["dof"]),
                    "fixed_am_reduced_chi2": float(fit_am["reduced_chi2"]),
                    "fixed_am_bic": float(fit_am["bic"]),
                    "fixed_am_pos_rms_km": float(fit_am["pos_rms_km"]),
                    "fixed_am_dop_rms_km_s": float(fit_am["dop_rms_km_s"]),
                    "fixed_am_n": int(fit_am["n_points"]),
                    "fixed_am_n_outliers": int(len(keep_am) - np.count_nonzero(keep_am)),
                    "fixed_am_cd_a_over_m_m2_kg": float(fit_am.get("cd_a_over_m_m2_kg", np.nan)),
                }
            )
        plausibility_terms = [
            float(track.get("fixed_velocity_reduced_chi2", np.inf)),
            float(track.get("fixed_am_reduced_chi2", np.inf)),
        ]
        plausibility_redchi = float(np.nanmin(plausibility_terms))
        plausible = bool(np.isfinite(plausibility_redchi) and plausibility_redchi < 1.0)
        track["interstellar_alias_plausible"] = plausible
        track["interstellar_alias_plausibility_redchi_threshold"] = 1.0
        track["interstellar_alias_plausibility_redchi"] = plausibility_redchi
        if plausible:
            track["interstellar_alias_plausibility_model"] = (
                "fixed_velocity"
                if plausibility_terms[0] <= plausibility_terms[1]
                else "fixed_am"
            )


def attach_selection_trajectory(track):
    """Choose the trajectory fit used for alias hypothesis testing."""
    candidates = []
    if np.isfinite(track.get("ballistic_reduced_chi2", np.nan)):
        candidates.append(("msis_drag", float(track["ballistic_reduced_chi2"])))
    if np.isfinite(track.get("fixed_velocity_reduced_chi2", np.nan)):
        candidates.append(("fixed_velocity", float(track["fixed_velocity_reduced_chi2"])))
    if not candidates:
        return
    model_type, redchi = min(candidates, key=lambda item: item[1])
    prefix = "fixed_velocity" if model_type == "fixed_velocity" else "ballistic"
    label = "fixed velocity" if model_type == "fixed_velocity" else "MSIS drag"
    track["selection_model_type"] = model_type
    track["selection_model_label"] = label
    track["selection_reduced_chi2"] = float(redchi)
    track["selection_chi2"] = float(track.get(f"{prefix}_chi2", np.nan))
    track["selection_dof"] = int(track.get(f"{prefix}_dof", 0))
    track["selection_n"] = int(track.get(f"{prefix}_n", 0))
    track["selection_n_outliers"] = int(track.get(f"{prefix}_n_outliers", 0))
    track["selection_pos_rms_km"] = float(track.get(f"{prefix}_pos_rms_km", np.nan))
    track["selection_dop_rms_km_s"] = float(track.get(f"{prefix}_dop_rms_km_s", np.nan))
    track["selection_model"] = np.asarray(track.get(f"{prefix}_model", np.empty((0, 3))), dtype=np.float64)
    track["selection_velocity_km_s"] = np.asarray(track.get(f"{prefix}_velocity_km_s", np.empty((0, 3))), dtype=np.float64)
    track["selection_pred_doppler_km_s"] = np.asarray(track.get(f"{prefix}_pred_doppler_km_s", np.empty(0)), dtype=np.float64)
    track["selection_pos_res_km"] = np.asarray(track.get(f"{prefix}_pos_res_km", np.empty((0, 3))), dtype=np.float64)
    track["selection_dop_res_km_s"] = np.asarray(track.get(f"{prefix}_dop_res_km_s", np.empty(0)), dtype=np.float64)
    track["selection_keep"] = np.asarray(track.get(f"{prefix}_keep", np.ones(len(track.get("idx", [])), dtype=bool)), dtype=bool)
    if model_type == "fixed_velocity":
        track["selection_params"] = np.asarray(track.get("fixed_velocity_params", np.full(6, np.nan)), dtype=np.float64)
        track["selection_parameter_covariance"] = np.asarray(
            track.get("fixed_velocity_parameter_covariance", np.full((6, 6), np.nan)),
            dtype=np.float64,
        )
        track["selection_parameter_std"] = np.asarray(track.get("fixed_velocity_parameter_std", np.full(6, np.nan)), dtype=np.float64)
        track["selection_t_ref_s"] = float(track.get("fixed_velocity_t_ref_s", 0.0))
    else:
        track["selection_params"] = np.asarray(track.get("ballistic_params", np.full(7, np.nan)), dtype=np.float64)
        track["selection_parameter_covariance"] = np.asarray(
            track.get("ballistic_parameter_covariance", np.full((7, 7), np.nan)),
            dtype=np.float64,
        )
        track["selection_parameter_std"] = np.asarray(track.get("ballistic_parameter_std", np.full(7, np.nan)), dtype=np.float64)
        track["selection_t_ref_s"] = float(track.get("ballistic_t_ref_s", 0.0))


def score_tx_beam_consistency(tracks, candidates):
    """Score how well each candidate lies near the active transmit beam."""
    beam_vecs = tx_beam_unit_vectors()
    for track in tracks:
        if "fit_points" not in track:
            continue
        rows = [candidates[i] for i in track["idx"]]
        dirs = np.asarray([[r["u"], r["v"], r["w"]] for r in rows], dtype=np.float64)
        beam_id = np.asarray([r["beam_id"] for r in rows], dtype=np.int64)
        snr = np.asarray([r["snr"] for r in rows], dtype=np.float64)
        t = np.asarray([r["t_rel"] for r in rows], dtype=np.float64)
        active = beam_vecs[beam_id]
        dot = np.sum(dirs * active, axis=1)
        angle_deg = np.rad2deg(np.arccos(np.clip(dot, -1.0, 1.0)))
        beam_center_dcos = np.linalg.norm(dirs[:, :2] - active[:, :2], axis=1)
        weight = np.maximum(snr, 1e-6)
        track["tx_beam_angle_deg"] = angle_deg
        track["tx_beam_center_distance_dc"] = beam_center_dcos
        track["tx_beam_time"] = t
        track["tx_beam_id"] = beam_id
        track["tx_beam_snr"] = snr
        track["tx_beam_snr_weighted_mean_dc"] = float(np.sum(weight * beam_center_dcos) / np.sum(weight))
        track["tx_beam_snr_weighted_rms_dc"] = float(np.sqrt(np.sum(weight * beam_center_dcos**2) / np.sum(weight)))
        track["tx_beam_weighted_rms_deg"] = float(np.sqrt(np.sum(weight * angle_deg**2) / np.sum(weight)))
        track["tx_beam_weighted_mean_deg"] = float(np.sum(weight * angle_deg) / np.sum(weight))
        track["tx_beam_median_deg"] = float(np.median(angle_deg))
        track["tx_beam_p90_deg"] = float(np.percentile(angle_deg, 90.0))

    scored = [t for t in tracks if "tx_beam_weighted_rms_deg" in t]
    scored.sort(key=lambda t: t["tx_beam_weighted_rms_deg"])
    for rank, track in enumerate(scored):
        track["tx_beam_rank"] = rank


def score_candidate_diagnostics(tracks, candidates, tx_gain_maps=None, tx_lobe_centroids=None):
    """Attach non-decision diagnostic metrics to ranked candidate tracks."""
    beam_vecs = tx_beam_unit_vectors()
    for track in tracks:
        if "fit_points" not in track:
            continue
        rows = [candidates[i] for i in track["idx"]]
        coherence = np.asarray([r["coherence"] for r in rows], dtype=np.float64)
        snr = np.asarray([r["snr"] for r in rows], dtype=np.float64)
        weight = np.maximum(snr, 1e-6)
        track["coherence_weighted_mean"] = float(np.sum(weight * coherence) / np.sum(weight))
        track["coherence_p10"] = float(np.percentile(coherence, 10.0))
        track["coherence_min"] = float(np.min(coherence))
        track["track_unique_pulse_fraction"] = float(track["unique_pulses"] / max(1, len({r["pulse"] for r in candidates})))

        dirs = np.asarray([[r["u"], r["v"], r["w"]] for r in rows], dtype=np.float64)
        beam_id = np.asarray([r["beam_id"] for r in rows], dtype=np.int64)
        ranges = np.asarray([r["range_km"] for r in rows], dtype=np.float64)
        doppler_km_s = np.asarray([r["doppler_mps"] / 1e3 for r in rows], dtype=np.float64)
        t = np.asarray([r["t_rel"] for r in rows], dtype=np.float64)
        track["measurement_range_span_km"] = float(np.nanmax(ranges) - np.nanmin(ranges)) if len(ranges) else np.nan
        track["measurement_doppler_span_km_s"] = (
            float(np.nanmax(doppler_km_s) - np.nanmin(doppler_km_s)) if len(doppler_km_s) else np.nan
        )
        if tx_lobe_centroids is not None:
            lobe_dist = nearest_tx_grating_lobe_distances(dirs, beam_id, tx_lobe_centroids)
            good_lobe = np.isfinite(lobe_dist)
            track["tx_lobe_distance_dc"] = lobe_dist
            if np.any(good_lobe):
                lobe_weight = weight[good_lobe]
                track["tx_lobe_snr_weighted_mean_dc"] = float(np.sum(lobe_weight * lobe_dist[good_lobe]) / np.sum(lobe_weight))
                track["tx_lobe_snr_weighted_rms_dc"] = float(np.sqrt(np.sum(lobe_weight * lobe_dist[good_lobe] ** 2) / np.sum(lobe_weight)))
                track["tx_lobe_p90_dc"] = float(np.percentile(lobe_dist[good_lobe], 90.0))
            else:
                track["tx_lobe_snr_weighted_mean_dc"] = np.nan
                track["tx_lobe_snr_weighted_rms_dc"] = np.nan
                track["tx_lobe_p90_dc"] = np.nan
            per_beam_lobe = np.full(5, np.nan, dtype=np.float64)
            per_beam_lobe_n = np.zeros(5, dtype=np.int64)
            for beam in range(5):
                bgood = good_lobe & (beam_id == beam)
                per_beam_lobe_n[beam] = int(np.count_nonzero(bgood))
                if per_beam_lobe_n[beam] > 0:
                    bw = weight[bgood]
                    per_beam_lobe[beam] = float(np.sum(bw * lobe_dist[bgood]) / np.sum(bw))
            track["tx_lobe_snr_weighted_mean_dc_by_beam"] = per_beam_lobe
            track["tx_lobe_n_by_beam"] = per_beam_lobe_n
        if tx_gain_maps is not None and all("grid_row" in r and "grid_col" in r for r in rows):
            grid_row = np.asarray([r["grid_row"] for r in rows], dtype=np.int64)
            grid_col = np.asarray([r["grid_col"] for r in rows], dtype=np.int64)
            gain_db = tx_gain_maps[beam_id, grid_row, grid_col].astype(np.float64)
        else:
            gain_db = tx_array_gain_db(dirs, beam_id, beam_vecs=beam_vecs)
        snr_db = 10.0 * np.log10(np.maximum(snr, 1e-6))
        range_corr_snr_db = snr_db + 40.0 * np.log10(np.maximum(ranges, 1e-6) / np.nanmedian(ranges))
        good = np.isfinite(gain_db) & np.isfinite(range_corr_snr_db)
        track["tx_array_time"] = t
        track["tx_array_gain_db"] = gain_db
        track["tx_array_range_corrected_snr_db"] = range_corr_snr_db
        per_beam_rms = np.full(5, np.nan, dtype=np.float64)
        per_beam_mad = np.full(5, np.nan, dtype=np.float64)
        per_beam_corr = np.full(5, np.nan, dtype=np.float64)
        per_beam_n = np.zeros(5, dtype=np.int64)
        for beam in range(5):
            bgood = good & (beam_id == beam)
            per_beam_n[beam] = int(np.count_nonzero(bgood))
            if per_beam_n[beam] >= 3:
                gain_rel_b = gain_db[bgood] - np.nanmedian(gain_db[bgood])
                snr_rel_b = range_corr_snr_db[bgood] - np.nanmedian(range_corr_snr_db[bgood])
                resid_b = snr_rel_b - gain_rel_b
                per_beam_rms[beam] = float(np.sqrt(np.mean(resid_b**2)))
                per_beam_mad[beam] = float(1.4826 * np.median(np.abs(resid_b - np.median(resid_b))))
                per_beam_corr[beam] = (
                    float(np.corrcoef(snr_rel_b, gain_rel_b)[0, 1])
                    if np.std(gain_rel_b) > 1e-6 and np.std(snr_rel_b) > 1e-6
                    else np.nan
                )
        track["tx_array_snr_rms_db_by_beam"] = per_beam_rms
        track["tx_array_snr_mad_db_by_beam"] = per_beam_mad
        track["tx_array_snr_corr_by_beam"] = per_beam_corr
        track["tx_array_snr_n_by_beam"] = per_beam_n
        if np.count_nonzero(good) >= 8:
            gain_rel = gain_db[good] - np.nanmedian(gain_db[good])
            snr_rel = range_corr_snr_db[good] - np.nanmedian(range_corr_snr_db[good])
            resid = snr_rel - gain_rel
            track["tx_array_snr_rms_db"] = float(np.sqrt(np.mean(resid**2)))
            track["tx_array_snr_mad_db"] = float(1.4826 * np.median(np.abs(resid - np.median(resid))))
            track["tx_array_snr_corr"] = (
                float(np.corrcoef(snr_rel, gain_rel)[0, 1])
                if np.std(gain_rel) > 1e-6 and np.std(snr_rel) > 1e-6
                else np.nan
            )
            track["tx_array_gain_span_db"] = float(np.nanpercentile(gain_db[good], 95) - np.nanpercentile(gain_db[good], 5))
        else:
            track["tx_array_snr_rms_db"] = np.nan
            track["tx_array_snr_mad_db"] = np.nan
            track["tx_array_snr_corr"] = np.nan
            track["tx_array_gain_span_db"] = np.nan


def score_combined_hypotheses(tracks):
    """Rank hypotheses with the best supported trajectory model and TX diagnostics."""
    for track in tracks:
        attach_selection_trajectory(track)
    scored = [t for t in tracks if "selection_reduced_chi2" in t]
    tx_beam_scale_dc = 0.035
    tx_beam_weight = 1.0
    line_weight = 0.10
    good_fit_redchi_threshold = 1.5
    head_echo_min_speed_km_s = 5.0
    if not scored:
        fallback = [t for t in tracks if t["reason"] == "kept" and "line_reduced_chi2" in t]
        if not fallback:
            return np.nan
        for track in fallback:
            track["combined_tx_sigma_deg"] = np.nan
            track["combined_tx_term"] = np.nan
            track["combined_score"] = float(track["line_reduced_chi2"])
            track["combined_score_source"] = "line_reduced_chi2"
        fallback.sort(
            key=lambda t: (
                bool(t.get("low_detection_altitude_reject", False)),
                bool(t.get("linearity_reject", False)),
                float(t.get("line_reduced_chi2", np.inf)),
            )
        )
        for rank, track in enumerate(fallback):
            track["combined_rank"] = rank
            track["combined_good_fit"] = False
            track["combined_reject"] = True
        delta = fallback[1]["combined_score"] - fallback[0]["combined_score"] if len(fallback) > 1 else np.inf
        odds = float(np.exp(0.5 * delta)) if np.isfinite(delta) and delta < 1400 else np.inf
        for track in fallback:
            track["combined_delta_to_next"] = float(delta)
            track["combined_odds_best_vs_second"] = float(odds)
        return np.nan
    for track in scored:
        track["combined_tx_sigma_deg"] = np.nan
        tx_beam = float(track.get("tx_beam_snr_weighted_mean_dc", np.nan))
        tx_term = 0.0 if not np.isfinite(tx_beam) else (tx_beam / tx_beam_scale_dc) ** 2
        line_term = float(track.get("line_length_adjusted_reduced_chi2", track.get("line_reduced_chi2", 0.0)))
        line_term = min(line_term, 100.0) / 25.0 if np.isfinite(line_term) else 0.0
        redchi = float(track.get("selection_reduced_chi2", np.inf))
        full_coverage = not bool(track.get("ballistic_pulse_coverage_reject", False))
        low_start_reject = bool(track.get("ballistic_low_start_altitude_reject", False))
        if track.get("selection_model_type") == "fixed_velocity":
            low_start_reject = False
        selection_params = np.asarray(track.get("selection_params", []), dtype=np.float64)
        selection_speed_km_s = (
            float(np.linalg.norm(selection_params[3:6]) / 1e3)
            if selection_params.size >= 6 and np.all(np.isfinite(selection_params[3:6]))
            else np.nan
        )
        head_echo_speed_reject = bool(
            np.isfinite(selection_speed_km_s) and selection_speed_km_s < head_echo_min_speed_km_s
        )
        range_span_km = float(track.get("measurement_range_span_km", np.nan))
        short_static_measurement_reject = bool(
            int(track.get("unique_pulses", 0)) <= 4
            and np.isfinite(range_span_km)
            and range_span_km < 0.2
        )
        good_fit = bool(
            full_coverage
            and not bool(track.get("common_pulse_completion_reject", False))
            and not low_start_reject
            and not head_echo_speed_reject
            and not short_static_measurement_reject
            and np.isfinite(redchi)
            and redchi <= good_fit_redchi_threshold
        )
        track["combined_tx_term"] = float(tx_term)
        track["combined_line_term"] = float(line_term)
        track["combined_tx_beam_scale_dc"] = float(tx_beam_scale_dc)
        track["combined_tx_beam_weight"] = float(tx_beam_weight)
        track["combined_line_weight"] = float(line_weight)
        track["combined_good_fit_redchi_threshold"] = float(good_fit_redchi_threshold)
        track["head_echo_min_speed_km_s"] = float(head_echo_min_speed_km_s)
        track["selection_speed_km_s"] = float(selection_speed_km_s)
        track["head_echo_speed_reject"] = head_echo_speed_reject
        track["short_static_measurement_reject"] = short_static_measurement_reject
        track["combined_good_fit"] = good_fit
        track["combined_tx_distance_dc"] = tx_beam
        if good_fit:
            if np.isfinite(tx_beam):
                track["combined_score"] = float(tx_beam)
                track["combined_score_source"] = f"tx_beam_center_distance_among_{track.get('selection_model_type', 'trajectory')}_redchi_le_1p5"
            else:
                track["combined_score"] = float(redchi)
                track["combined_score_source"] = f"{track.get('selection_model_type', 'trajectory')}_redchi_no_tx_beam_center_metric"
        else:
            tx_penalty = tx_term if np.isfinite(tx_term) else 1e6
            track["combined_score"] = float(redchi + 0.01 * tx_penalty + line_weight * line_term)
            track["combined_score_source"] = f"fallback_{track.get('selection_model_type', 'trajectory')}_redchi_plus_weak_tx_beam_center_proximity"
    scored.sort(
        key=lambda t: (
            bool(t.get("ballistic_pulse_coverage_reject", False)),
            bool(t.get("common_pulse_completion_reject", False)),
            bool(t.get("ballistic_low_start_altitude_reject", False)) and t.get("selection_model_type") != "fixed_velocity",
            bool(t.get("head_echo_speed_reject", False)),
            not bool(t.get("combined_good_fit", False)),
            (
                float(t.get("combined_tx_distance_dc", np.inf))
                if t.get("combined_good_fit", False) and np.isfinite(t.get("combined_tx_distance_dc", np.nan))
                else float(t.get("selection_reduced_chi2", np.inf))
            ),
            t["combined_score"],
        )
    )
    for rank, track in enumerate(scored):
        track["combined_rank"] = rank
        ballistic_reject_blocks = bool(track.get("ballistic_reject", False)) and not bool(track.get("combined_good_fit", False))
        if track.get("selection_model_type") == "fixed_velocity":
            ballistic_reject_blocks = False
        track["combined_reject"] = (
            rank != 0
            or not bool(track.get("combined_good_fit", False))
            or bool(track.get("common_pulse_completion_reject", False))
            or bool(track.get("ballistic_pulse_coverage_reject", False))
            or ballistic_reject_blocks
            or (bool(track.get("ballistic_low_start_altitude_reject", False)) and track.get("selection_model_type") != "fixed_velocity")
            or bool(track.get("head_echo_speed_reject", False))
            or bool(track.get("short_static_measurement_reject", False))
        )
    if len(scored) > 1:
        delta = scored[1]["combined_score"] - scored[0]["combined_score"]
        odds = float(np.exp(0.5 * delta))
    else:
        delta = np.inf
        odds = np.inf
    for track in scored:
        track["combined_delta_to_next"] = float(delta)
        track["combined_odds_best_vs_second"] = float(odds)
    return np.nan


def fit_winning_physics_trajectory_model(tracks, event_epoch_unix, sigma_pos_km, sigma_dop_km_s):
    """Run the three physics trajectory models only for the final winning alias."""
    ranked = sorted(
        [t for t in tracks if "combined_rank" in t and "ballistic_t" in t],
        key=lambda t: t["combined_rank"],
    )
    if not ranked:
        return None
    winner = ranked[0]
    try:
        import fit_best_alias_physics_models as physics

        t_s = np.asarray(winner["ballistic_t"], dtype=np.float64)
        points_km = np.asarray(winner["fit_points"], dtype=np.float64)
        doppler_km_s = np.asarray(winner["ballistic_doppler_km_s"], dtype=np.float64)
        rho_of_alt_m, _msis_meta = pbal.density_interpolator(event_epoch_unix)
        sigma_pos = float(sigma_pos_km) if np.isfinite(sigma_pos_km) and sigma_pos_km > 0 else physics.DEFAULT_SIGMA_POS_KM
        sigma_dop = float(sigma_dop_km_s) if np.isfinite(sigma_dop_km_s) and sigma_dop_km_s > 0 else physics.DEFAULT_SIGMA_DOP_KM_S
        p0_linear = physics.linear_initial_guess(t_s, points_km)
        fixed_velocity = physics.fit_model(
            "fixed_velocity",
            t_s,
            points_km,
            doppler_km_s,
            rho_of_alt_m,
            sigma_pos,
            sigma_dop,
            [p0_linear],
        )
        fixed_am_starts = [np.concatenate([fixed_velocity["params"][:6], [log_b]]) for log_b in (-4.0, -2.0, 0.0, 1.0, 2.0)]
        fixed_am = physics.fit_model(
            "fixed_am",
            t_s,
            points_km,
            doppler_km_s,
            rho_of_alt_m,
            sigma_pos,
            sigma_dop,
            fixed_am_starts,
        )
        radius_guess = 3.0 / (
            4.0
            * physics.SHRINKING_METEOROID_DENSITY_KG_M3
            * max(float(fixed_am["cd_a_over_m_m2_kg"]), 1e-30)
        )
        radius_starts = sorted(
            set(
                float(np.clip(r, physics.SHRINKING_RADIUS_MIN_M, physics.SHRINKING_RADIUS_MAX_M))
                for r in [radius_guess, *physics.SHRINKING_RADIUS_START_GRID_M]
            )
        )
        shrinking_starts = [np.concatenate([fixed_am["params"][:6], [np.log10(radius)]]) for radius in radius_starts]
        shrinking_radius = physics.fit_model(
            "shrinking_radius",
            t_s,
            points_km,
            doppler_km_s,
            rho_of_alt_m,
            sigma_pos,
            sigma_dop,
            shrinking_starts,
        )
        fits = {
            "fixed_velocity": fixed_velocity,
            "fixed_am": fixed_am,
            "shrinking_radius": shrinking_radius,
        }
        best_name = physics.best_model_name(fits)
        winner.update(
            {
                "physics_best_model": best_name,
                "physics_best_model_label": physics.model_label(best_name),
                "physics_bic_fixed_velocity": float(fixed_velocity["bic"]),
                "physics_bic_fixed_am": float(fixed_am["bic"]),
                "physics_bic_shrinking_radius": float(shrinking_radius["bic"]),
                "physics_reduced_chi2_fixed_velocity": float(fixed_velocity["reduced_chi2"]),
                "physics_reduced_chi2_fixed_am": float(fixed_am["reduced_chi2"]),
                "physics_reduced_chi2_shrinking_radius": float(shrinking_radius["reduced_chi2"]),
                "physics_fixed_am_cd_a_over_m_m2_kg": float(fixed_am["cd_a_over_m_m2_kg"]),
                "physics_shrinking_radius0_m": float(shrinking_radius["radius_m"][0]),
                "physics_shrinking_mass0_kg": float(shrinking_radius["mass_kg"][0]),
            }
        )
    except Exception as exc:
        winner["physics_model_error"] = f"{type(exc).__name__}: {exc}"
    return winner


def candidate_number(track):
    """Human-facing one-based candidate number derived from the final rank."""
    rank = track.get("combined_rank", track.get("ballistic_rank", None))
    if rank is None:
        return "?"
    return int(rank) + 1


def hypothesis_label(track):
    return f"H{int(track.get('hypothesis_id', 0)):02d}"


def track_cluster_indices(track):
    """Candidate indices belonging to the discovered visual cluster."""
    return np.asarray(track.get("seed_idx", track.get("idx", [])), dtype=np.int64)


def track_cluster_direction(track, candidates, n_points=120):
    """Direction-cosine line over the time span of the discovered cluster."""
    idx = track_cluster_indices(track)
    if len(idx) == 0:
        return np.empty((0, 3), dtype=np.float64)
    t = np.asarray([candidates[i]["t_rel"] for i in idx], dtype=np.float64)
    if not np.any(np.isfinite(t)):
        return np.empty((0, 3), dtype=np.float64)
    t0 = float(np.nanmin(t))
    t1 = float(np.nanmax(t))
    if not np.isfinite(t0) or not np.isfinite(t1):
        return np.empty((0, 3), dtype=np.float64)
    if t1 <= t0:
        td = np.asarray([t0], dtype=np.float64)
    else:
        td = np.linspace(t0, t1, int(n_points))
    u = track["u_coeff"][0] + track["u_coeff"][1] * td
    v = track["v_coeff"][0] + track["v_coeff"][1] * td
    w = np.sqrt(np.maximum(0.0, 1.0 - u**2 - v**2))
    return np.column_stack([u, v, w])


def snr_weighted_track_coherence(track, candidates):
    rows = [candidates[j] for j in track["idx"]]
    coherence = np.asarray([r["coherence"] for r in rows], dtype=np.float64)
    snr = np.asarray([r["snr"] for r in rows], dtype=np.float64)
    good = np.isfinite(coherence) & np.isfinite(snr) & (snr > 0.0)
    if not np.any(good):
        return np.nan
    return float(np.sum(snr[good] * coherence[good]) / np.sum(snr[good]))


def plot_visibility_rejections(candidates, tracks, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.6), constrained_layout=True)

    ax = axes[0]
    circle = plt.Circle((0.0, 0.0), 1.0, color="black", fill=False, linewidth=1.0)
    ax.add_patch(circle)
    colors = {"kept": "tab:green", "below horizon": "tab:red", "wraps": "tab:orange"}
    labels_done = set()
    for track in tracks:
        rows = [candidates[j] for j in track_cluster_indices(track)]
        uu = np.asarray([r["u"] for r in rows])
        vv = np.asarray([r["v"] for r in rows])
        label = track["reason"] if track["reason"] not in labels_done else None
        labels_done.add(track["reason"])
        color = colors.get(track["reason"], "0.4")
        ax.plot(uu, vv, ".", ms=3.0, color=color, alpha=0.65, label=label)
        dd = track_cluster_direction(track, candidates)
        if len(dd):
            ax.plot(dd[:, 0], dd[:, 1], "-", lw=1.1, color=color, alpha=0.9)
            mid = len(dd) // 2
            ax.text(dd[mid, 0], dd[mid, 1], hypothesis_label(track), fontsize=7)
    ax.set_xlabel("u/east direction cosine")
    ax.set_ylabel("v/north direction cosine")
    ax.set_title("Path hypotheses in direction-cosine space")
    ax.set_xlim(-1.03, 1.03)
    ax.set_ylim(-1.03, 1.03)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.2)
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[1]
    for track in tracks:
        z = track["dense_model"][:, 2]
        ax.plot(track["dense_t"], z, color=colors[track["reason"]], lw=1.2, alpha=0.9)
        ax.text(track["dense_t"][-1], z[-1], hypothesis_label(track), fontsize=7, color=colors[track["reason"]])
    ax.axhline(0.0, color="black", lw=1.0)
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Fitted local up coordinate (km)")
    ax.set_title("Full-path horizon visibility test")
    ax.grid(True, alpha=0.25)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_linearity_rejections(candidates, tracks, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    visible = [t for t in tracks if t["reason"] == "kept" and "linearity_reject" in t]
    n_rejected = sum(t.get("linearity_reject", False) for t in visible)
    n_kept = len(visible) - n_rejected
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.3), constrained_layout=True)

    ax = axes[0]
    labels_done = set()
    for track in visible:
        color = "tab:red" if track.get("linearity_reject", False) else "tab:green"
        label = "rejected: nonlinear" if track.get("linearity_reject", False) else "kept: linear"
        if label in labels_done:
            label = None
        else:
            labels_done.add(label)
        pts = track["fit_points"]
        model = track["line_model_points"]
        order = np.argsort((model - track["line_center"]) @ track["line_direction"])
        ax.plot(pts[:, 0], pts[:, 1], ".", ms=3.0, color=color, alpha=0.65, label=label)
        ax.plot(model[order, 0], model[order, 1], "-", lw=1.0, color=color, alpha=0.85)
        mid = len(order) // 2
        ax.text(model[order[mid], 0], model[order[mid], 1], hypothesis_label(track), fontsize=7)
    ax.set_xlabel("East (km)")
    ax.set_ylabel("North (km)")
    ax.set_title("3D candidate positions: horizontal projection")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)

    ax = axes[1]
    for track in visible:
        color = "tab:red" if track.get("linearity_reject", False) else "tab:green"
        pts = track["fit_points"]
        model = track["line_model_points"]
        along = (model - track["line_center"]) @ track["line_direction"]
        order = np.argsort(along)
        ax.plot(along, pts[:, 2], ".", ms=3.0, color=color, alpha=0.65)
        ax.plot(along[order], model[order, 2], "-", lw=1.0, color=color, alpha=0.85)
        mid = len(order) // 2
        ax.text(along[order[mid]], model[order[mid], 2], hypothesis_label(track), fontsize=7)
    ax.set_xlabel("Along best-fit 3D line (km)")
    ax.set_ylabel("Up (km)")
    ax.set_title("3D candidate positions: vertical projection")
    ax.grid(True, alpha=0.25)

    ax = axes[2]
    if visible:
        labels = np.arange(len(visible))
        hypothesis_labels = [hypothesis_label(t) for t in visible]
        rms = [t["line_rms_km"] for t in visible]
        colors = ["tab:red" if t.get("linearity_reject", False) else "tab:green" for t in visible]
        ax.scatter(labels, rms, c=colors, s=70, zorder=3)
        for lab, y, hlabel in zip(labels, rms, hypothesis_labels):
            ax.text(lab, y, hlabel, fontsize=8, ha="center", va="bottom")
        ax.set_xticks(labels)
        ax.set_xticklabels(hypothesis_labels, rotation=45, ha="right")
        ax.axhline(visible[0]["line_rms_threshold_km"], color="black", lw=1.0, ls="--", label="RMS threshold")
    ax.set_xlabel("Visible path hypothesis")
    ax.set_ylabel("3D line RMS residual (km)")
    ax.set_title(f"3D straight-line rejection: {n_rejected} rejected, {n_kept} kept")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_descent_rejections(tracks, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    scored = [
        t
        for t in tracks
        if t["reason"] == "kept" and not t.get("linearity_reject", False) and "descent_rate_km_s" in t
    ]
    scored.sort(key=lambda t: t.get("hypothesis_id", 999))
    if not scored:
        return

    n_rejected = sum(t.get("descent_reject", False) for t in scored)
    n_kept = len(scored) - n_rejected
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.2), constrained_layout=True)

    ax = axes[0]
    for track in scored:
        color = "tab:red" if track.get("descent_reject", False) else "tab:green"
        t = np.asarray(track["fit_t"], dtype=np.float64)
        z = np.asarray(track["fit_points"][:, 2], dtype=np.float64)
        order = np.argsort(t)
        design = np.column_stack([np.ones(len(t)), t - np.mean(t)])
        coeff = np.linalg.lstsq(design, z, rcond=None)[0]
        zfit = coeff[0] + coeff[1] * (t - np.mean(t))
        ax.plot(t[order], zfit[order], "-", color=color, lw=1.3, alpha=0.9)
        ax.plot(t[order], z[order], ".", color=color, ms=2.8, alpha=0.55)
        ax.text(t[order][-1], zfit[order][-1], hypothesis_label(track), fontsize=8, color=color)
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Local up coordinate (km)")
    ax.set_title("Altitude trend for line-valid hypotheses")
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    labels = [hypothesis_label(t) for t in scored]
    rates = np.asarray([t["descent_rate_km_s"] for t in scored], dtype=np.float64)
    colors = ["tab:red" if t.get("descent_reject", False) else "tab:green" for t in scored]
    x = np.arange(len(scored))
    ax.bar(x, rates, color=colors, alpha=0.85)
    ax.axhline(0.0, color="black", lw=1.0, ls="--", label="upward rejection threshold")
    for xx, yy in zip(x, rates):
        va = "bottom" if yy >= 0 else "top"
        offset = 0.25 if yy >= 0 else -0.25
        ax.text(xx, yy + offset, f"{yy:.2f}", ha="center", va=va, fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_xlabel("Line-valid path hypothesis")
    ax.set_ylabel("Fitted vertical velocity (km/s)")
    ax.set_title(f"Downward-motion test: {n_rejected} rejected, {n_kept} kept")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=8)

    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_ballistic_ranking(tracks, candidates, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    survivors = final_trajectory_tracks_for_view(tracks)
    survivors.sort(key=lambda t: t.get("combined_rank", t.get("ballistic_rank", 999)))
    n_rej = sum(t.get("combined_reject", False) for t in survivors)
    n_keep = len(survivors) - n_rej
    cmap = plt.get_cmap("plasma")
    all_snr_parts = [np.asarray(t.get("ballistic_snr_db", []), dtype=np.float64) for t in survivors]
    all_snr = np.concatenate(all_snr_parts) if all_snr_parts else np.array([], dtype=np.float64)
    all_snr = all_snr[np.isfinite(all_snr)]
    if len(all_snr):
        snr_norm = mpl.colors.Normalize(vmin=float(np.nanmin(all_snr)), vmax=float(np.nanmax(all_snr)))
    else:
        snr_norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
    snr_cmap = plt.get_cmap("viridis")
    chi_values = np.asarray([t.get("selection_reduced_chi2", t.get("ballistic_reduced_chi2", np.nan)) for t in survivors], dtype=np.float64)
    finite_chi = chi_values[np.isfinite(chi_values) & (chi_values > 0.0)]
    if len(finite_chi):
        norm = mpl.colors.LogNorm(vmin=max(np.nanmin(finite_chi), 1e-2), vmax=max(np.nanmax(finite_chi), 1e-1))
    else:
        norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)

    def track_color(track):
        return cmap(norm(max(float(track["ballistic_reduced_chi2"]), 1e-6)))

    def track_linestyle(track):
        return "--" if track.get("combined_reject", False) else "-"

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.2), constrained_layout=True)
    all_snr_parts = [np.asarray(t.get("ballistic_snr_db", []), dtype=np.float64) for t in survivors]
    all_snr = np.concatenate(all_snr_parts) if all_snr_parts else np.array([], dtype=np.float64)
    all_snr = all_snr[np.isfinite(all_snr)]
    if len(all_snr):
        snr_norm = mpl.colors.Normalize(vmin=float(np.nanmin(all_snr)), vmax=float(np.nanmax(all_snr)))
    else:
        snr_norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
    snr_cmap = plt.get_cmap("viridis")

    ax = axes[0]
    snr_scatter = None
    for track in survivors:
        label = hypothesis_label(track)
        color = track_color(track)
        t = track["ballistic_t"]
        dop = track["ballistic_doppler_km_s"]
        pred = track["ballistic_pred_doppler_km_s"]
        snr = np.asarray(track.get("ballistic_snr_db", np.full_like(dop, np.nan)), dtype=np.float64)
        order = np.argsort(t)
        ax.plot(t[order], pred[order], color=color, lw=1.4, ls=track_linestyle(track), alpha=0.92)
        snr_scatter = ax.scatter(
            t[order],
            dop[order],
            c=snr[order],
            cmap=snr_cmap,
            norm=snr_norm,
            s=13,
            marker="o",
            edgecolors="none",
            alpha=0.75,
        )
        ax.text(t[order][-1], pred[order][-1], label, fontsize=8, color=color)
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Doppler/range rate (km/s)")
    ax.set_title("Doppler fits; points colored by SNR")
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    for track in survivors:
        label = hypothesis_label(track)
        color = track_color(track)
        pts = track["fit_points"]
        model = track["ballistic_model"]
        keep = track["ballistic_keep"]
        order = np.argsort(track["ballistic_t"])
        ax.plot(model[order, 0], model[order, 1], color=color, lw=1.6, ls=track_linestyle(track), alpha=0.94)
        rows = [candidates[i] for i in track["idx"]]
        beam_id = np.asarray([r["beam_id"] for r in rows], dtype=np.int64)
        scatter_points_by_tx_beam(ax, pts, beam_id, mask=keep, size=2.8, alpha=0.65)
        for beam in tx_beam_center_projection_points(track, candidates):
            draw_tx_beam_projection_marker(ax, beam, color=color, fontsize=7)
        mid = order[len(order) // 2]
        ax.text(
            model[mid, 0],
            model[mid, 1],
            label,
            fontsize=8,
            color="black",
            bbox={"facecolor": "white", "edgecolor": color, "alpha": 0.85},
        )
    ax.set_xlabel("East (km)")
    ax.set_ylabel("North (km)")
    ax.set_title(f"Combined-ranked candidate paths: {n_rej} rejected, {n_keep} kept")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    ax.plot([], [], "-", color="black", label="accepted")
    ax.plot([], [], "--", color="black", label="rejected")
    for bid, name in enumerate(TX_BEAM_SHORT_NAMES):
        ax.plot([], [], ".", color=f"C{bid}", label=f"beam {name}")

    ax = axes[2]
    coh = [t.get("coherence_weighted_mean", np.nan) for t in survivors]
    bal_term = [t["ballistic_reduced_chi2"] for t in survivors]
    colors = [track_color(t) for t in survivors]
    markers = ["o" if not t.get("combined_reject", False) else "x" for t in survivors]
    for track, x, y, color, marker in zip(survivors, coh, bal_term, colors, markers):
        ax.scatter([x], [y], c=[color], s=95, marker=marker, linewidths=1.6)
        ax.text(x, y, hypothesis_label(track), fontsize=8)
    ax.set_xlabel("SNR-weighted mean coherence")
    ax.set_ylabel("Ballistic reduced $\\chi^2$")
    ax.set_title("Ballistic score and coherence diagnostic")
    ax.grid(True, alpha=0.25)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    fig.colorbar(sm, ax=axes.ravel().tolist(), label="Ballistic reduced $\\chi^2$")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def rotation_matrix_313(node_deg: float, inc_deg: float, argp_deg: float) -> np.ndarray:
    node = np.deg2rad(node_deg)
    inc = np.deg2rad(inc_deg)
    argp = np.deg2rad(argp_deg)
    cn, sn = np.cos(node), np.sin(node)
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(argp), np.sin(argp)
    rz_node = np.array([[cn, -sn, 0.0], [sn, cn, 0.0], [0.0, 0.0, 1.0]])
    rx_inc = np.array([[1.0, 0.0, 0.0], [0.0, ci, -si], [0.0, si, ci]])
    rz_argp = np.array([[cw, -sw, 0.0], [sw, cw, 0.0], [0.0, 0.0, 1.0]])
    return rz_node @ rx_inc @ rz_argp


def orbit_xy_from_elements(elements: np.ndarray, n_points: int = 360) -> tuple[np.ndarray, np.ndarray]:
    a_au, e, inc_deg, node_deg, argp_deg, _nu_deg, q_au = [float(x) for x in elements]
    if not np.all(np.isfinite(elements[:6])) or e < 0:
        return np.array([]), np.array([])
    if e < 1.0:
        f = np.linspace(0.0, 2.0 * np.pi, n_points)
        p = max(abs(a_au) * max(1.0 - e * e, 1e-12), 1e-9)
    else:
        fmax = np.arccos(np.clip(-1.0 / max(e, 1.0 + 1e-9), -1.0, 1.0)) - 1e-3
        f = np.linspace(-fmax, fmax, n_points)
        p = max(abs(q_au) * (1.0 + e), 1e-9)
    denom = 1.0 + e * np.cos(f)
    good = denom > 1e-6
    r = p / denom[good]
    f = f[good]
    perifocal = np.vstack([r * np.cos(f), r * np.sin(f), np.zeros_like(r)])
    xyz = rotation_matrix_313(node_deg, inc_deg, argp_deg) @ perifocal
    return xyz[0], xyz[1]


def plot_summary_orbit_panel(ax, orbit_h5_path: Path | None, hypothesis: str = "H03"):
    planets = [
        ("Mercury", 0.387, "0.62"),
        ("Venus", 0.723, "#b08a5b"),
        ("Earth", 1.000, "#4f8fd6"),
        ("Mars", 1.524, "#c46852"),
        ("Jupiter", 5.203, "#8d7658"),
    ]
    theta = np.linspace(0.0, 2.0 * np.pi, 361)
    for name, radius, color in planets:
        ax.plot(radius * np.cos(theta), radius * np.sin(theta), color=color, lw=0.8, alpha=0.75)
        if name in {"Earth", "Jupiter"}:
            ax.text(radius, 0.0, name, fontsize=6, color=color, ha="left", va="bottom")
    ax.scatter([0.0], [0.0], s=45, color="#ffd21f", edgecolor="black", zorder=5)
    result = read_dasst_orbit_group(orbit_h5_path, hypothesis)
    if result is None:
        ax.text(
            0.5,
            0.5,
            f"DASST orbit not available\nfor {hypothesis}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=10,
        )
    else:
        try:
            kep = np.asarray(result["kepler"], dtype=np.float64)
            std = np.asarray(result.get("kepler_std", np.full(7, np.nan)), dtype=np.float64)
            samples = np.asarray(result.get("kepler_samples", np.empty((0, 7))), dtype=np.float64)
            frac_e_gt_1 = float(result.get("frac_e_gt_1", np.nan))
            finite = samples[np.all(np.isfinite(samples), axis=1)][:120]
            for sample in finite:
                x, y = orbit_xy_from_elements(sample)
                if len(x):
                    ax.plot(x, y, color="0.25", alpha=0.045, lw=0.6)
            x, y = orbit_xy_from_elements(kep)
            if len(x):
                ax.plot(x, y, color="black", alpha=0.88, lw=1.8)
            ax.text(
                0.02,
                0.03,
                (
                    f"{hypothesis} DASST orbit\n"
                    f"a={kep[0]:.2g}+/-{std[0]:.1g} AU\n"
                    f"e={kep[1]:.4f}+/-{std[1]:.2g}\n"
                    f"i={kep[2]:.1f}+/-{std[2]:.1g} deg\n"
                    f"Omega={kep[3]:.1f}+/-{std[3]:.1g} deg\n"
                    f"omega={kep[4]:.1f}+/-{std[4]:.1g} deg\n"
                    f"nu={kep[5]:.1f}+/-{std[5]:.1g} deg\n"
                    f"q={kep[6]:.3f}+/-{std[6]:.2g} AU\n"
                    f"P(e>1)={frac_e_gt_1:.2f}"
                ),
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=7,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78},
            )
        except Exception as exc:
            ax.text(0.5, 0.5, f"Orbit plot failed:\n{exc}", transform=ax.transAxes, ha="center", va="center", fontsize=7)
    ax.set_title("3c. DASST-corrected orbit samples")
    ax.set_xlabel("Heliocentric ecliptic x (AU)")
    ax.set_ylabel("Heliocentric ecliptic y (AU)")
    ax.set_xlim(-8.0, 8.0)
    ax.set_ylim(-8.0, 8.0)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.90", lw=0.8)


def plot_disambiguation_summary(
    candidates,
    tracks,
    output: Path,
    orbit_h5_path: Path | None = None,
    obs=None,
    cut=None,
    obs_rti=None,
    event_epoch_unix=None,
    response_u=None,
    response_v=None,
    selected_response=None,
):
    """Combine the main disambiguation tests into one 3x3 summary figure."""
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 3, figsize=(16.0, 14.0), constrained_layout=False)
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.055, top=0.955, wspace=0.24, hspace=0.35)
    range_fit_ax = axes[0, 1]
    fit_summary_ax = axes[1, 1]

    visibility_colors = {"kept": "tab:green", "below horizon": "tab:red", "wraps": "tab:orange"}

    def annotate_panel(ax, text):
        ax.text(
            0.5,
            0.5,
            text,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9,
            color="0.25",
            bbox={"facecolor": "white", "edgecolor": "0.7", "alpha": 0.88, "pad": 4.0},
        )

    def legend_if_any(ax, **kwargs):
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(**kwargs)

    ax = axes[0, 0]
    if response_u is not None and response_v is not None and selected_response is not None:
        ax.pcolormesh(
            response_u,
            response_v,
            selected_response,
            shading="auto",
            cmap="Blues",
            vmin=0.45,
            vmax=1.0,
            alpha=0.34,
            zorder=0,
        )
    ax.add_patch(plt.Circle((0.0, 0.0), 1.0, color="black", fill=False, linewidth=1.0))
    if candidates:
        cand_u = np.asarray([c["u"] for c in candidates], dtype=np.float64)
        cand_v = np.asarray([c["v"] for c in candidates], dtype=np.float64)
        ax.plot(cand_u, cand_v, ".", ms=1.4, color="0.65", alpha=0.25, label="high-coherence peaks", zorder=1)
    beam_vecs = tx_beam_unit_vectors()
    for bid, (u_b, v_b, _w_b) in enumerate(beam_vecs):
        ax.scatter([u_b], [v_b], marker="x", s=48, color=f"C{bid}", linewidths=1.4, zorder=6)
        ax.text(
            u_b,
            v_b,
            f" {TX_BEAM_SHORT_NAMES[bid] if bid < len(TX_BEAM_SHORT_NAMES) else bid}",
            fontsize=7,
            color=f"C{bid}",
            va="center",
            ha="left",
            zorder=7,
        )
    labels_done = set()
    for track in sorted(tracks, key=lambda t: t.get("hypothesis_id", 999)):
        actual_idx = np.asarray(track.get("idx", []), dtype=np.int64)
        actual_idx = actual_idx[(actual_idx >= 0) & (actual_idx < len(candidates))]
        rows = [candidates[j] for j in actual_idx]
        if len(rows) == 0:
            rows = [candidates[j] for j in track_cluster_indices(track)]
        uu = np.asarray([r["u"] for r in rows])
        vv = np.asarray([r["v"] for r in rows])
        color = visibility_colors.get(track["reason"], "0.4")
        label = track["reason"] if track["reason"] not in labels_done else None
        labels_done.add(track["reason"])
        ax.scatter(
            uu,
            vv,
            s=1.0,
            marker="o",
            color=color,
            alpha=0.72,
            linewidths=0.0,
            label=label,
            zorder=5,
        )
        if len(uu):
            finite_uv = np.isfinite(uu) & np.isfinite(vv)
            if np.any(finite_uv):
                u0 = float(np.nanmin(uu[finite_uv]))
                u1 = float(np.nanmax(uu[finite_uv]))
                v0 = float(np.nanmin(vv[finite_uv]))
                v1 = float(np.nanmax(vv[finite_uv]))
                pad = 0.018
                width = max(u1 - u0, 0.018)
                height = max(v1 - v0, 0.018)
                cx = 0.5 * (u0 + u1)
                cy = 0.5 * (v0 + v1)
                ax.add_patch(
                    plt.Rectangle(
                        (cx - 0.5 * width - pad, cy - 0.5 * height - pad),
                        width + 2.0 * pad,
                        height + 2.0 * pad,
                        fill=False,
                        edgecolor=color,
                        linewidth=0.75,
                        alpha=0.72,
                        zorder=4,
                    )
                )
        dd = track_cluster_direction(track, candidates)
        if len(dd):
            mid = len(dd) // 2
            label_xy = dd[mid, :2].astype(np.float64)
            direction = dd[min(len(dd) - 1, mid + 1), :2] - dd[max(0, mid - 1), :2]
            normal = np.asarray([-direction[1], direction[0]], dtype=np.float64)
            norm = float(np.linalg.norm(normal))
            if norm < 1e-9:
                normal = label_xy.copy()
                norm = float(np.linalg.norm(normal))
            if norm < 1e-9:
                normal = np.asarray([0.0, 1.0], dtype=np.float64)
                norm = 1.0
            normal /= norm
            if len(uu):
                centroid = np.asarray([float(np.nanmean(uu)), float(np.nanmean(vv))])
                if np.dot(label_xy + 0.055 * normal - centroid, label_xy - centroid) < 0:
                    normal *= -1.0
            radial = label_xy.copy()
            radial_norm = float(np.linalg.norm(radial))
            if radial_norm < 1e-9:
                radial = normal.copy()
            else:
                radial /= radial_norm
            offset = 0.085
            directions = [
                normal,
                -normal,
                radial,
                -radial,
                normal + radial,
                normal - radial,
                -normal + radial,
                -normal - radial,
            ]
            point_xy = np.column_stack([uu, vv]) if len(uu) else np.empty((0, 2), dtype=np.float64)
            best_xy = label_xy + offset * normal
            best_score = -np.inf
            for direction_vec in directions:
                direction_vec = np.asarray(direction_vec, dtype=np.float64)
                direction_norm = float(np.linalg.norm(direction_vec))
                if direction_norm < 1e-9:
                    continue
                candidate_xy = label_xy + offset * direction_vec / direction_norm
                if float(np.linalg.norm(candidate_xy)) > 0.99:
                    continue
                if len(point_xy):
                    distance = np.sqrt(np.sum((point_xy - candidate_xy[None, :]) ** 2, axis=1))
                    score = float(np.nanmin(distance))
                else:
                    score = 0.0
                if score > best_score:
                    best_score = score
                    best_xy = candidate_xy
            label_xy = best_xy
            ax.text(
                label_xy[0],
                label_xy[1],
                hypothesis_label(track),
                fontsize=7,
                ha="center",
                va="center",
                zorder=4,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.30, "pad": 1.2},
            )
    if not tracks:
        n_pulses = len({c["pulse"] for c in candidates}) if candidates else 0
        annotate_panel(ax, f"No path hypotheses formed\n{len(candidates)} coherence peaks across {n_pulses} pulses")
    ax.set_xlabel("u/east direction cosine")
    ax.set_ylabel("v/north direction cosine")
    ax.set_title(f"1a. Provisional path hypotheses (N={len(tracks)})")
    ax.set_xlim(-1.03, 1.03)
    ax.set_ylim(-1.03, 1.03)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.2)
    legend_if_any(ax, loc="lower right", fontsize=7)

    ax = axes[0, 2]
    if obs is not None and cut is not None:
        rti_obs = obs if obs_rti is None else obs_rti
        fit_t0_s = float(obs["tx_idx"][0] / 1e6) if len(obs["tx_idx"]) else None
        plot_matched_filter_rti_panel(ax, rti_obs, cut, tracks, fit_t0_s=fit_t0_s, obs_overlay=obs)
    else:
        annotate_panel(ax, "No cut data available for RTI")
        ax.set_title("1c. Range-Doppler matched-filter RTI")

    horizon_visible = [t for t in tracks if t["reason"] == "kept"]
    visible = [t for t in horizon_visible if not t.get("descent_reject", False)]

    def path_stage_rejected(track):
        return bool(track.get("linearity_reject", False) or track.get("descent_reject", False))

    def path_stage_label(track):
        if track.get("linearity_reject", False) and track.get("descent_reject", False):
            return "rejected: nonlinear/upward"
        if track.get("linearity_reject", False):
            return "rejected: nonlinear"
        if track.get("descent_reject", False):
            return "rejected: upward"
        return "kept"

    line_survivors = [t for t in visible if not path_stage_rejected(t)]

    ax = axes[1, 0]
    tested = [t for t in tracks if "fit_points" in t and ("selection_reduced_chi2" in t or "fixed_velocity_reduced_chi2" in t)]
    if tested:
        x = np.asarray([straight_line_sumsq_km2(t) for t in tested], dtype=np.float64)
        y = np.asarray([t.get("selection_reduced_chi2", t.get("fixed_velocity_reduced_chi2", np.nan)) for t in tested], dtype=np.float64)
        colors = ["tab:green" if not t.get("combined_reject", True) else "tab:red" for t in tested]
        markers = ["o" if not t.get("combined_reject", True) else "x" for t in tested]
        for track, xx, yy, color, marker in zip(tested, x, y, colors, markers):
            if not (np.isfinite(xx) and np.isfinite(yy)):
                continue
            ax.scatter([max(xx, 1e-8)], [max(yy, 1e-8)], c=[color], marker=marker, s=58, linewidths=1.3, zorder=3)
            ax.text(max(xx, 1e-8), max(yy, 1e-8), hypothesis_label(track), fontsize=7, ha="left", va="bottom")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axhline(1.0, color="black", lw=0.9, ls="--", alpha=0.8, label=r"$\chi^2_\nu=1$")
        ax.plot([], [], "o", color="tab:green", label="accepted")
        ax.plot([], [], "x", color="tab:red", label="rejected")
    else:
        annotate_panel(ax, "No trajectory-fit hypotheses")
    ax.set_xlabel("Straight-line residual sum of squares (km$^2$)")
    ax.set_ylabel("Constant-velocity reduced $\\chi^2$")
    ax.set_title("2a. Trajectory fit vs. straight-line norm")
    ax.grid(True, alpha=0.25, which="both")
    legend_if_any(ax, loc="best", fontsize=7)

    ax = axes[1, 2]
    ranking_visible = sorted(
        [
            t
            for t in tracks
            if "combined_rank" in t
            and "selection_reduced_chi2" in t
        ],
        key=lambda t: t.get("combined_rank", 999),
    )
    if ranking_visible:
        labels = [hypothesis_label(t) for t in ranking_visible]
        redchi = np.asarray([t["selection_reduced_chi2"] for t in ranking_visible], dtype=np.float64)
        log_redchi = np.log10(np.maximum(redchi, 1e-12))
        tx_dcos = np.asarray(
            [t.get("tx_beam_snr_weighted_mean_dc", np.nan) for t in ranking_visible],
            dtype=np.float64,
        )
        valid_tx = np.isfinite(tx_dcos) & (tx_dcos > 0.0)
        if not np.any(valid_tx):
            tx_dcos = np.asarray(
                [t.get("tx_lobe_snr_weighted_mean_dc", np.nan) for t in ranking_visible],
                dtype=np.float64,
            )
        valid_tx = np.isfinite(tx_dcos) & (tx_dcos > 0.0)
        if not np.any(valid_tx):
            tx_dcos = np.asarray([t.get("tx_beam_weighted_rms_deg", np.nan) for t in ranking_visible], dtype=np.float64)
            # Keep the panel defined even for legacy products missing the dcos metric.
            tx_dcos = np.deg2rad(np.maximum(tx_dcos, 1e-6))
        log_tx_dcos = np.log10(np.maximum(tx_dcos, 1e-6))
        colors = ["tab:red" if t.get("combined_reject", False) else "tab:green" for t in ranking_visible]
        ax.scatter(log_tx_dcos, log_redchi, c=colors, s=60, zorder=3)
        for xx_i, yy, label, chi, track in zip(log_tx_dcos, log_redchi, labels, redchi, ranking_visible):
            ax.text(xx_i, yy, label, fontsize=7, ha="center", va="bottom", clip_on=True)
            ax.text(
                xx_i,
                yy - 0.08,
                f"{chi:.1f} {track.get('selection_model_type', '')[:3]}",
                fontsize=6,
                ha="center",
                va="top",
                color="0.25",
                clip_on=True,
            )
        ax.axhline(0.0, color="black", lw=1.0, ls="--", label=r"$\chi^2_\nu=1$")
    ax.set_xlabel(r"$\log_{10}$ SNR-weighted TX beam-center distance (dcos)")
    ax.set_ylabel(r"$\log_{10}$ selected-model reduced $\chi^2$")
    ax.set_title("2c. TX beam-center distance vs. selected trajectory fit")
    ax.grid(True, alpha=0.25)
    if not ranking_visible:
        annotate_panel(ax, "No trajectory-ranked candidates")
    legend_if_any(ax, fontsize=7)

    survivors = final_trajectory_tracks_for_view(tracks)
    survivors.sort(key=lambda t: t.get("combined_rank", t.get("ballistic_rank", 999)))
    chi_values = np.asarray([t.get("selection_reduced_chi2", t.get("ballistic_reduced_chi2", np.nan)) for t in survivors], dtype=np.float64)
    finite_chi = chi_values[np.isfinite(chi_values) & (chi_values > 0.0)]
    if len(finite_chi):
        norm = mpl.colors.LogNorm(vmin=max(np.nanmin(finite_chi), 1e-2), vmax=max(np.nanmax(finite_chi), 1e-1))
    else:
        norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
    cmap = plt.get_cmap("plasma")
    all_snr_parts = [np.asarray(t.get("ballistic_snr_db", []), dtype=np.float64) for t in survivors]
    all_snr = np.concatenate(all_snr_parts) if all_snr_parts else np.array([], dtype=np.float64)
    all_snr = all_snr[np.isfinite(all_snr)]
    if len(all_snr):
        snr_norm = mpl.colors.Normalize(vmin=float(np.nanmin(all_snr)), vmax=float(np.nanmax(all_snr)))
    else:
        snr_norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
    snr_cmap = plt.get_cmap("viridis")

    winner = winning_track_for_view(survivors)

    def is_winner(track):
        return winner is not None and track.get("hypothesis_id") == winner.get("hypothesis_id")

    def bcolor(track):
        return "black" if is_winner(track) else "0.65"

    def bls(track):
        return "--" if track.get("combined_reject", False) else "-"

    if event_epoch_unix is None:
        event_epoch_unix = np.nan
    plot_best_range_fit_panel(range_fit_ax, winner, survivors, event_epoch_unix, annotate_panel)
    plot_best_fit_summary_panel(fit_summary_ax, winner, orbit_h5_path, annotate_panel)

    ax = axes[2, 0]
    snr_scatter = None
    for track in survivors:
        order = np.argsort(track["ballistic_t"])
        color = bcolor(track)
        snr = np.asarray(track.get("ballistic_snr_db", np.full_like(track["ballistic_doppler_km_s"], np.nan)), dtype=np.float64)
        pred_dop = np.asarray(
            track["selection_pred_doppler_km_s"]
            if "selection_pred_doppler_km_s" in track
            else track["ballistic_pred_doppler_km_s"],
            dtype=np.float64,
        )
        ax.plot(
            track["ballistic_t"][order],
            pred_dop[order],
            color=color,
            lw=1.7 if is_winner(track) else 1.1,
            ls=bls(track),
            alpha=0.96 if is_winner(track) else 0.75,
        )
        if is_winner(track):
            band_95 = doppler_uncertainty_from_ballistic_covariance(track, np.asarray(track["ballistic_t"], dtype=np.float64), event_epoch_unix)
            if band_95 is not None and np.any(np.isfinite(band_95)):
                ax.fill_between(
                    track["ballistic_t"][order],
                    (pred_dop - band_95)[order],
                    (pred_dop + band_95)[order],
                    color="black",
                    alpha=0.13,
                    linewidth=0.0,
                    label="95% fit region",
                )
        snr_scatter = ax.scatter(
            track["ballistic_t"][order],
            track["ballistic_doppler_km_s"][order],
            c=snr[order],
            cmap=snr_cmap,
            norm=snr_norm,
            s=24,
            marker="o",
            edgecolors="none",
            alpha=0.75,
            zorder=100,
        )
        ax.text(track["ballistic_t"][order][-1], pred_dop[order][-1], hypothesis_label(track), fontsize=7, color=color)
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Doppler/range rate (km/s)")
    ax.set_title("3a. Selected-model Doppler fit")
    ax.grid(True, alpha=0.25)
    if not survivors:
        annotate_panel(ax, "No final trajectory-fit candidates")
    ax = axes[2, 1]
    summary_tracks = survivors
    set_winning_path_view(ax, survivors, candidates)

    def beam_center_inside_current_axes(beam):
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        x = float(beam["east_km"])
        y = float(beam["north_km"])
        return min(x0, x1) <= x <= max(x0, x1) and min(y0, y1) <= y <= max(y0, y1)

    for track in summary_tracks:
        order = np.argsort(track["ballistic_t"])
        color = bcolor(track)
        pts = track["fit_points"]
        model = np.asarray(track.get("selection_model", track.get("ballistic_model")), dtype=np.float64)
        keep = np.asarray(track.get("selection_keep", track.get("ballistic_keep")), dtype=bool)
        ax.plot(
            model[order, 0],
            model[order, 1],
            color=color,
            lw=1.8 if is_winner(track) else 1.1,
            ls=bls(track),
            alpha=0.97 if is_winner(track) else 0.7,
        )
        rows = [candidates[i] for i in track["idx"]]
        beam_id = np.asarray([r["beam_id"] for r in rows], dtype=np.int64)
        scatter_points_by_tx_beam(ax, pts, beam_id, mask=keep, size=5.2, alpha=0.75)
        for beam in tx_beam_center_projection_points(track, candidates):
            if beam_center_inside_current_axes(beam):
                draw_tx_beam_projection_marker(ax, beam, color=color, fontsize=6.5)
        mid = order[len(order) // 2]
        ax.text(
            model[mid, 0],
            model[mid, 1],
            hypothesis_label(track),
            fontsize=6.4,
            bbox={"facecolor": "white", "edgecolor": color, "alpha": 0.85, "pad": 1.2},
            clip_on=True,
        )
    ax.set_xlabel("East (km)")
    ax.set_ylabel("North (km)")
    ax.set_title("3b. Final ballistic-ranked paths")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    for bid, name in enumerate(TX_BEAM_SHORT_NAMES):
        ax.plot([], [], ".", color=f"C{bid}", label=f"beam {name}")
    if not summary_tracks:
        annotate_panel(ax, "No final ballistic ranking\nSee rejection stages above")
    ax = axes[2, 2]
    orbit_hypothesis = hypothesis_label(winner) if winner is not None else "H03"
    plot_summary_orbit_panel(ax, orbit_h5_path, hypothesis=orbit_hypothesis)

    fig.savefig(output, dpi=240)
    plt.close(fig)


def plot_rejected_event_summary(output: Path, sample_idx: int, reason: str, details: str):
    """Write the standard overview filename for events rejected before AOI search."""
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 3, figsize=(16.0, 14.0), constrained_layout=True)
    message = (
        f"Event {sample_idx}\n"
        "Rejected before interferometric disambiguation\n"
        f"{reason}\n"
        f"{details}"
    )
    for ax in axes.flat:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(
            0.5,
            0.5,
            message,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=10,
            color="0.25",
            bbox={"facecolor": "white", "edgecolor": "0.7", "alpha": 0.9, "pad": 5.0},
        )
        ax.set_frame_on(True)
    axes[0, 0].set_title("1a. Provisional path hypotheses")
    axes[0, 1].set_title("1b. Full-path horizon test")
    axes[0, 2].set_title("1c. SNR-weighted coherence averages")
    axes[1, 0].set_title("2a. 3D line horizontal projection")
    axes[1, 1].set_title("2b. 3D line vertical projection")
    axes[1, 2].set_title("2c. Straight-line chi-square metric")
    axes[2, 0].set_title("3a. Ballistic Doppler fit")
    axes[2, 1].set_title("3b. Final ballistic-ranked paths")
    axes[2, 2].set_title("3c. Orbit samples")
    fig.savefig(output, dpi=240)
    plt.close(fig)


def plot_empty_event_summary(output: Path, sample_idx: int, n_raw: int, snr_threshold: float):
    """Write the standard overview filename for events with no usable detections."""
    plot_rejected_event_summary(
        output,
        sample_idx,
        f"0 measurement points above SNR threshold {snr_threshold:g}",
        f"{n_raw} raw cut samples before thresholding",
    )


def _h5_set_scalar_attrs(group, attrs):
    """Set HDF5 scalar attributes while skipping values that are awkward to serialize."""
    for key, value in attrs.items():
        if value is None:
            continue
        if isinstance(value, (str, bytes, bool, int, float, np.integer, np.floating, np.bool_)):
            group.attrs[key] = value


def _h5_create_array(group, name, value, dtype=np.float64):
    arr = np.asarray(value, dtype=dtype)
    group.create_dataset(name, data=arr)


def tx_phase_quality_attrs(quality: dict | None) -> dict:
    if not quality:
        return {}
    out = {}
    for key, value in quality.items():
        if isinstance(value, (str, bytes, bool, int, float, np.integer, np.floating, np.bool_)):
            out[f"tx_phase_quality_{key}"] = value
    return out


def cut_tx_waveform_quality(cut: dict, min_corr: float = 0.8, max_phase_std_deg: float = 20.0) -> dict:
    z_tx = np.asarray(cut["ztx_pulses_re"], dtype=np.float32) + 1j * np.asarray(cut["ztx_pulses_im"], dtype=np.float32)
    if z_tx.ndim != 2 or z_tx.shape[0] < 2 or z_tx.shape[1] < 4:
        return {
            "available": False,
            "good": False,
            "reason": "cut_tx_waveform_unavailable",
        }
    pulse_power = np.sqrt(np.sum(np.abs(z_tx) ** 2, axis=1))
    good_pulses = np.isfinite(pulse_power) & (pulse_power > 0)
    if np.count_nonzero(good_pulses) < 2:
        return {
            "available": False,
            "good": False,
            "reason": "too_few_valid_cut_tx_pulses",
        }
    z_norm = z_tx[good_pulses] / pulse_power[good_pulses, None]
    template = np.mean(z_norm, axis=0)
    template_norm = np.sqrt(np.sum(np.abs(template) ** 2))
    if not np.isfinite(template_norm) or template_norm <= 0:
        return {
            "available": False,
            "good": False,
            "reason": "cut_tx_waveform_template_failed",
        }
    template = template / template_norm
    corr = np.sum(np.conj(template)[None, :] * z_norm, axis=1)
    corr_abs = np.abs(corr)
    phase = np.angle(corr)
    mean_phase = np.angle(np.mean(np.exp(1j * phase)))
    phase_resid = np.angle(np.exp(1j * (phase - mean_phase)))
    phase_std_deg = float(np.rad2deg(np.sqrt(np.mean(phase_resid**2))))
    median_corr = float(np.nanmedian(corr_abs))
    good = bool(median_corr >= float(min_corr) and phase_std_deg <= float(max_phase_std_deg))
    return {
        "available": True,
        "good": good,
        "raw_good": good,
        "reason": "ok_cut_tx_waveform" if good else "cut_tx_waveform_unstable",
        "computed_available": True,
        "computed_from_cut_tx_waveform": True,
        "valid_channel_count": int(np.count_nonzero(good_pulses)),
        "mean_amplitude": median_corr,
        "max_abs_diff_deg": phase_std_deg,
        "rms_diff_deg": phase_std_deg,
        "threshold_deg": float(max_phase_std_deg),
    }


def write_rejected_disambiguation_diagnostics_h5(
    output: Path,
    sample_idx: int,
    reason: str,
    details: str,
    extra_attrs: dict | None = None,
):
    """Write a minimal per-event diagnostic file for events rejected before path fitting."""
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as h:
        h.attrs["schema_version"] = "pansy_interferometric_disambiguation_diagnostics_v1"
        h.attrs["source_program"] = "plot_interferometric_disambiguation.py"
        h.attrs["sample_idx"] = int(sample_idx)
        h.attrs["sample_epoch_unix"] = float(sample_idx / 1e6)
        h.attrs["event_status"] = "rejected_before_interferometric_disambiguation"
        h.attrs["rejection_reason"] = str(reason)
        h.attrs["rejection_details"] = str(details)
        h.attrs["position_frame"] = "PANSY local ENU"
        h.attrs["position_units"] = "km"
        h.attrs["direction_cosine_units"] = "dimensionless"
        _h5_set_scalar_attrs(h, extra_attrs or {})
        h.create_group("candidates")
        h.create_group("hypotheses")


def write_disambiguation_diagnostics_h5(
    candidates,
    tracks,
    output: Path,
    sample_idx: int,
    state_h5_path: Path | None = None,
    dasst_orbit_h5_path: Path | None = None,
    sigma_pos_km=np.nan,
    sigma_dop_km_s=np.nan,
    range_time_component_lengths=None,
    selected_range_time_component=-1,
    tx_phase_quality=None,
):
    """Persist compact per-event products needed for later statistical summaries.

    The key product for position histograms is the winning
    hypotheses/Hxx/position_enu_km dataset. Non-winning hypotheses retain
    scalar scores and flags, but not 3D position arrays.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as h:
        _h5_set_scalar_attrs(
            h,
            {
                "schema_version": "pansy_interferometric_disambiguation_diagnostics_v1",
                "source_program": "plot_interferometric_disambiguation.py",
                "sample_idx": int(sample_idx),
                "sample_epoch_unix": float(sample_idx / 1e6),
                "event_status": "processed",
                "position_frame": "PANSY local ENU",
                "position_units": "km",
                "direction_cosine_convention": "u=east, v=north, up=sqrt(1-u^2-v^2)",
                "direction_cosine_units": "dimensionless",
                "range_units": "km",
                "doppler_units": "m/s",
                "sigma_pos_km": float(sigma_pos_km),
                "sigma_dop_km_s": float(sigma_dop_km_s),
                "orbit_state_h5": str(state_h5_path) if state_h5_path is not None else "",
                "dasst_orbit_h5": str(dasst_orbit_h5_path) if dasst_orbit_h5_path is not None else "",
                "selection_metric": "combined_score",
                "range_time_component_count": int(len(range_time_component_lengths or [])),
                "selected_range_time_component": int(selected_range_time_component),
                **tx_phase_quality_attrs(tx_phase_quality),
            },
        )
        if range_time_component_lengths is not None:
            h.create_dataset(
                "range_time_component_lengths",
                data=np.asarray(range_time_component_lengths, dtype=np.int64),
            )

        cand_grp = h.create_group("candidates")
        if candidates:
            fields = {
                "u": np.float64,
                "v": np.float64,
                "w": np.float64,
                "grid_row": np.int64,
                "grid_col": np.int64,
                "coherence": np.float64,
                "t_rel_s": np.float64,
                "pulse": np.int64,
                "range_km": np.float64,
                "doppler_mps": np.float64,
                "snr": np.float64,
                "beam_id": np.int64,
                "zpp": np.complex64,
                "zpp_phase_rad": np.float64,
                "zpp_phase_dt_s": np.float64,
                "zpp_doppler_mps": np.float64,
                "zpp_coarse_doppler_mps": np.float64,
                "zpp_coarse_phase_rad": np.float64,
                "zpp_prev_t_rel": np.float64,
                "zpp_prev_pulse": np.int64,
            }
            for field, dtype in fields.items():
                key = "t_rel" if field == "t_rel_s" else field
                cand_grp.create_dataset(field, data=np.asarray([c[key] for c in candidates], dtype=dtype))

        scored = [t for t in tracks if "combined_score" in t]
        winners = sorted(
            [t for t in scored if not t.get("combined_reject", False)],
            key=lambda t: t.get("combined_rank", 999999),
        )
        selected = winners[0] if winners else None
        if selected is not None:
            _h5_set_scalar_attrs(
                h,
                {
                    "selected_hypothesis": hypothesis_label(selected),
                    "selected_hypothesis_id": int(selected.get("hypothesis_id", -1)),
                    "selected_candidate_number": int(candidate_number(selected)),
                    "selected_combined_rank": int(selected.get("combined_rank", -1)),
                    "selected_combined_score": float(selected.get("combined_score", np.nan)),
                },
            )
        else:
            _h5_set_scalar_attrs(
                h,
                {
                    "selected_hypothesis": "",
                    "selected_hypothesis_id": -1,
                    "selected_candidate_number": -1,
                    "selected_combined_rank": -1,
                    "selected_combined_score": np.nan,
                },
            )

        hyp_grp = h.create_group("hypotheses")
        for track in sorted(tracks, key=lambda t: t.get("hypothesis_id", 999999)):
            label = hypothesis_label(track)
            grp = hyp_grp.create_group(label)
            _h5_set_scalar_attrs(
                grp,
                {
                    "hypothesis_label": label,
                    "hypothesis_id": int(track.get("hypothesis_id", -1)),
                    "candidate_number": int(candidate_number(track)) if "combined_rank" in track else -1,
                    "reason": str(track.get("reason", "unknown")),
                    "unique_pulses": int(track.get("unique_pulses", 0)),
                    "seed_unique_pulses": int(track.get("seed_unique_pulses", track.get("unique_pulses", 0))),
                    "common_pulse_count": int(track.get("common_pulse_count", track.get("unique_pulses", 0))),
                    "common_pulse_completion_rms_dc": float(track.get("common_pulse_completion_rms_dc", np.nan)),
                    "common_pulse_completion_max_dc": float(track.get("common_pulse_completion_max_dc", np.nan)),
                    "common_pulse_completion_reject": bool(track.get("common_pulse_completion_reject", False)),
                    "below_horizon": bool(track.get("below_horizon", False)),
                    "near_horizon": bool(track.get("near_horizon", False)),
                    "min_elevation_deg": float(track.get("min_elevation_deg", np.nan)),
                    "min_elevation_threshold_deg": float(track.get("min_elevation_threshold_deg", np.nan)),
                    "wraps": bool(track.get("wraps", False)),
                    "linearity_reject": bool(track.get("linearity_reject", False)),
                    "descent_reject": bool(track.get("descent_reject", False)),
                    "low_detection_altitude_reject": bool(track.get("low_detection_altitude_reject", False)),
                    "ballistic_reject": bool(track.get("ballistic_reject", False)),
                    "ballistic_pulse_coverage_reject": bool(track.get("ballistic_pulse_coverage_reject", False)),
                    "ballistic_required_unique_pulses": int(track.get("ballistic_required_unique_pulses", 0)),
                    "ballistic_low_start_altitude_reject": bool(track.get("ballistic_low_start_altitude_reject", False)),
                    "combined_reject": bool(track.get("combined_reject", False)),
                    "line_reduced_chi2": float(track.get("line_reduced_chi2", np.nan)),
                    "line_rms_km": float(track.get("line_rms_km", np.nan)),
                    "line_max_km": float(track.get("line_max_km", np.nan)),
                    "descent_rate_km_s": float(track.get("descent_rate_km_s", np.nan)),
                    "ballistic_rank": int(track.get("ballistic_rank", -1)) if "ballistic_rank" in track else -1,
                    "ballistic_reduced_chi2": float(track.get("ballistic_reduced_chi2", np.nan)),
                    "ballistic_pos_rms_km": float(track.get("ballistic_pos_rms_km", np.nan)),
                    "ballistic_dop_rms_km_s": float(track.get("ballistic_dop_rms_km_s", np.nan)),
                    "ballistic_phase_fit_available": bool(track.get("ballistic_phase_fit_available", False)),
                    "ballistic_phase_n": int(track.get("ballistic_phase_n", 0)),
                    "ballistic_phase_sigma_rad": float(track.get("ballistic_phase_sigma_rad", np.nan)),
                    "ballistic_phase_rms_rad": float(track.get("ballistic_phase_rms_rad", np.nan)),
                    "ballistic_start_altitude_km": float(track.get("ballistic_start_altitude_km", np.nan)),
                    "fixed_velocity_reduced_chi2": float(track.get("fixed_velocity_reduced_chi2", np.nan)),
                    "fixed_velocity_bic": float(track.get("fixed_velocity_bic", np.nan)),
                    "fixed_velocity_pos_rms_km": float(track.get("fixed_velocity_pos_rms_km", np.nan)),
                    "fixed_velocity_dop_rms_km_s": float(track.get("fixed_velocity_dop_rms_km_s", np.nan)),
                    "fixed_am_reduced_chi2": float(track.get("fixed_am_reduced_chi2", np.nan)),
                    "fixed_am_bic": float(track.get("fixed_am_bic", np.nan)),
                    "fixed_am_pos_rms_km": float(track.get("fixed_am_pos_rms_km", np.nan)),
                    "fixed_am_dop_rms_km_s": float(track.get("fixed_am_dop_rms_km_s", np.nan)),
                    "fixed_am_cd_a_over_m_m2_kg": float(track.get("fixed_am_cd_a_over_m_m2_kg", np.nan)),
                    "interstellar_alias_plausible": bool(track.get("interstellar_alias_plausible", False)),
                    "interstellar_alias_plausibility_redchi_threshold": float(track.get("interstellar_alias_plausibility_redchi_threshold", np.nan)),
                    "interstellar_alias_plausibility_redchi": float(track.get("interstellar_alias_plausibility_redchi", np.nan)),
                    "interstellar_alias_plausibility_model": str(track.get("interstellar_alias_plausibility_model", "")),
                    "selection_model_type": str(track.get("selection_model_type", "")),
                    "selection_model_label": str(track.get("selection_model_label", "")),
                    "selection_reduced_chi2": float(track.get("selection_reduced_chi2", np.nan)),
                    "selection_chi2": float(track.get("selection_chi2", np.nan)),
                    "selection_dof": int(track.get("selection_dof", 0)),
                    "selection_pos_rms_km": float(track.get("selection_pos_rms_km", np.nan)),
                    "selection_dop_rms_km_s": float(track.get("selection_dop_rms_km_s", np.nan)),
                    "selection_speed_km_s": float(track.get("selection_speed_km_s", np.nan)),
                    "head_echo_min_speed_km_s": float(track.get("head_echo_min_speed_km_s", np.nan)),
                    "head_echo_speed_reject": bool(track.get("head_echo_speed_reject", False)),
                    "measurement_range_span_km": float(track.get("measurement_range_span_km", np.nan)),
                    "measurement_doppler_span_km_s": float(track.get("measurement_doppler_span_km_s", np.nan)),
                    "short_static_measurement_reject": bool(track.get("short_static_measurement_reject", False)),
                    "combined_rank": int(track.get("combined_rank", -1)) if "combined_rank" in track else -1,
                    "combined_score": float(track.get("combined_score", np.nan)),
                    "combined_score_source": str(track.get("combined_score_source", "")),
                    "combined_good_fit": bool(track.get("combined_good_fit", False)),
                    "combined_good_fit_redchi_threshold": float(track.get("combined_good_fit_redchi_threshold", np.nan)),
                    "combined_tx_distance_dc": float(track.get("combined_tx_distance_dc", np.nan)),
                    "combined_tx_metric_dc": float(track.get("combined_tx_metric_dc", np.nan)),
                    "combined_tx_metric_source": str(track.get("combined_tx_metric_source", "")),
                    "combined_completion_term": float(track.get("combined_completion_term", np.nan)),
                    "combined_tx_term": float(track.get("combined_tx_term", np.nan)),
                    "combined_line_term": float(track.get("combined_line_term", np.nan)),
                    "combined_tx_beam_scale_dc": float(track.get("combined_tx_beam_scale_dc", np.nan)),
                    "coherence_weighted_mean": float(track.get("coherence_weighted_mean", np.nan)),
                    "tx_beam_rank": int(track.get("tx_beam_rank", -1)) if "tx_beam_rank" in track else -1,
                    "tx_beam_snr_weighted_mean_dc": float(track.get("tx_beam_snr_weighted_mean_dc", np.nan)),
                    "tx_beam_snr_weighted_rms_dc": float(track.get("tx_beam_snr_weighted_rms_dc", np.nan)),
                    "tx_lobe_snr_weighted_mean_dc": float(track.get("tx_lobe_snr_weighted_mean_dc", np.nan)),
                    "tx_lobe_snr_weighted_rms_dc": float(track.get("tx_lobe_snr_weighted_rms_dc", np.nan)),
                },
            )
            idx = np.asarray(track.get("idx", []), dtype=np.int64)
            grp.create_dataset("candidate_indices", data=idx)
            if "seed_idx" in track:
                grp.create_dataset("seed_candidate_indices", data=np.asarray(track["seed_idx"], dtype=np.int64))
            if "common_pulse_completion_distance_dc" in track:
                grp.create_dataset(
                    "common_pulse_completion_distance_dc",
                    data=np.asarray(track["common_pulse_completion_distance_dc"], dtype=np.float64),
                )
            if len(idx) and candidates:
                rows = [candidates[int(i)] for i in idx]
                for name, key, dtype in [
                    ("pulse", "pulse", np.int64),
                    ("t_rel_s", "t_rel", np.float64),
                    ("range_km", "range_km", np.float64),
                    ("doppler_mps", "doppler_mps", np.float64),
                    ("snr", "snr", np.float64),
                    ("beam_id", "beam_id", np.int64),
                    ("zpp", "zpp", np.complex64),
                    ("zpp_phase_rad", "zpp_phase_rad", np.float64),
                    ("zpp_phase_dt_s", "zpp_phase_dt_s", np.float64),
                    ("zpp_doppler_mps", "zpp_doppler_mps", np.float64),
                    ("zpp_coarse_doppler_mps", "zpp_coarse_doppler_mps", np.float64),
                    ("zpp_coarse_phase_rad", "zpp_coarse_phase_rad", np.float64),
                    ("zpp_prev_t_rel", "zpp_prev_t_rel", np.float64),
                    ("zpp_prev_pulse", "zpp_prev_pulse", np.int64),
                    ("coherence", "coherence", np.float64),
                ]:
                    grp.create_dataset(name, data=np.asarray([r.get(key, np.nan) for r in rows], dtype=dtype))
                grp.create_dataset(
                    "direction_cosines_uvw",
                    data=np.asarray([[r["u"], r["v"], r["w"]] for r in rows], dtype=np.float64),
                )

            if selected is not None and track is selected and "fit_points" in track:
                grp.create_dataset("position_enu_km", data=np.asarray(track["fit_points"], dtype=np.float64))

            for name, value in [
                ("fit_t_s", track.get("fit_t")),
                ("dense_t_s", track.get("dense_t")),
                ("dense_direction_uvw", track.get("dense_direction")),
                ("line_perpendicular_residual_km", track.get("line_perp_km")),
                ("line_sigma_km", track.get("line_sigma_km")),
                ("ballistic_keep", track.get("ballistic_keep")),
                ("ballistic_params", track.get("ballistic_params")),
                ("ballistic_parameter_covariance", track.get("ballistic_parameter_covariance")),
                ("ballistic_parameter_std", track.get("ballistic_parameter_std")),
                ("ballistic_velocity_km_s", track.get("ballistic_velocity_km_s")),
                ("ballistic_speed_km_s", track.get("ballistic_speed_km_s")),
                ("ballistic_t_s", track.get("ballistic_t")),
                ("ballistic_doppler_km_s", track.get("ballistic_doppler_km_s")),
                ("ballistic_pos_res_km", track.get("ballistic_pos_res_km")),
                ("ballistic_dop_res_km_s", track.get("ballistic_dop_res_km_s")),
                ("ballistic_pred_doppler_km_s", track.get("ballistic_pred_doppler_km_s")),
                ("ballistic_pre_phase_params", track.get("ballistic_pre_phase_params")),
                ("ballistic_phase_t_s", track.get("ballistic_phase_t")),
                ("ballistic_phase_prev_t_s", track.get("ballistic_phase_prev_t")),
                ("ballistic_phase_observed_rad", track.get("ballistic_phase_observed_rad")),
                ("ballistic_phase_model_rad", track.get("ballistic_phase_model_rad")),
                ("ballistic_phase_doppler_mps", track.get("ballistic_phase_doppler_mps")),
                ("ballistic_phase_coarse_doppler_mps", track.get("ballistic_phase_coarse_doppler_mps")),
                ("ballistic_phase_model_doppler_mps", track.get("ballistic_phase_model_doppler_mps")),
                ("ballistic_phase_model_unwrapped_rad", track.get("ballistic_phase_model_unwrapped_rad")),
                ("ballistic_phase_observed_model_branch_rad", track.get("ballistic_phase_observed_model_branch_rad")),
                ("ballistic_phase_initial_residual_rad", track.get("ballistic_phase_initial_residual_rad")),
                ("ballistic_phase_residual_rad", track.get("ballistic_phase_residual_rad")),
                ("fixed_velocity_params", track.get("fixed_velocity_params")),
                ("fixed_velocity_parameter_covariance", track.get("fixed_velocity_parameter_covariance")),
                ("fixed_velocity_parameter_std", track.get("fixed_velocity_parameter_std")),
                ("fixed_velocity_keep", track.get("fixed_velocity_keep")),
                ("fixed_velocity_model", track.get("fixed_velocity_model")),
                ("fixed_velocity_velocity_km_s", track.get("fixed_velocity_velocity_km_s")),
                ("fixed_velocity_pred_doppler_km_s", track.get("fixed_velocity_pred_doppler_km_s")),
                ("fixed_velocity_pos_res_km", track.get("fixed_velocity_pos_res_km")),
                ("fixed_velocity_dop_res_km_s", track.get("fixed_velocity_dop_res_km_s")),
                ("fixed_am_params", track.get("fixed_am_params")),
                ("fixed_am_parameter_covariance", track.get("fixed_am_parameter_covariance")),
                ("fixed_am_parameter_std", track.get("fixed_am_parameter_std")),
                ("fixed_am_keep", track.get("fixed_am_keep")),
                ("fixed_am_model", track.get("fixed_am_model")),
                ("fixed_am_velocity_km_s", track.get("fixed_am_velocity_km_s")),
                ("fixed_am_pred_doppler_km_s", track.get("fixed_am_pred_doppler_km_s")),
                ("fixed_am_pos_res_km", track.get("fixed_am_pos_res_km")),
                ("fixed_am_dop_res_km_s", track.get("fixed_am_dop_res_km_s")),
                ("selection_params", track.get("selection_params")),
                ("selection_parameter_covariance", track.get("selection_parameter_covariance")),
                ("selection_parameter_std", track.get("selection_parameter_std")),
                ("selection_keep", track.get("selection_keep")),
                ("selection_model", track.get("selection_model")),
                ("selection_velocity_km_s", track.get("selection_velocity_km_s")),
                ("selection_pred_doppler_km_s", track.get("selection_pred_doppler_km_s")),
                ("selection_pos_res_km", track.get("selection_pos_res_km")),
                ("selection_dop_res_km_s", track.get("selection_dop_res_km_s")),
                ("tx_beam_center_distance_dc", track.get("tx_beam_center_distance_dc")),
                ("tx_lobe_distance_dc", track.get("tx_lobe_distance_dc")),
            ]:
                if value is not None:
                    grp.create_dataset(name, data=np.asarray(value))


def plot_ballistic_horizontal_ranking(tracks, candidates, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    survivors = [
        t
        for t in tracks
        if t["reason"] == "kept" and not t.get("linearity_reject", False) and not t.get("descent_reject", False)
    ]
    survivors.sort(key=lambda t: t.get("ballistic_rank", 999))
    fig, ax = plt.subplots(figsize=(7.6, 6.8), constrained_layout=True)

    for track in survivors:
        label = hypothesis_label(track)
        color = "tab:red" if track.get("ballistic_reject", False) else "tab:green"
        pts = track["fit_points"]
        model = track["ballistic_model"]
        keep = track["ballistic_keep"]
        order = np.argsort(track["ballistic_t"])
        ax.plot(model[order, 0], model[order, 1], "-", color=color, lw=1.5, alpha=0.9)
        ax.plot(pts[keep, 0], pts[keep, 1], ".", color=color, ms=3.0, alpha=0.65)
        if np.any(~keep):
            ax.plot(pts[~keep, 0], pts[~keep, 1], "x", color=color, ms=4.0, alpha=0.45)
        for beam in tx_beam_center_projection_points(track, candidates):
            draw_tx_beam_projection_marker(ax, beam, color=color, fontsize=7)
        mid = order[len(order) // 2]
        ax.text(
            model[mid, 0],
            model[mid, 1],
            label,
            fontsize=8,
            color="black",
            bbox={"facecolor": "white", "edgecolor": color, "alpha": 0.85},
        )

    ax.set_xlabel("East (km)")
    ax.set_ylabel("North (km)")
    ax.set_title("Candidate paths ranked by robust ballistic statistical score")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    ax.plot([], [], "-", color="tab:green", label="survives ballistic fit")
    ax.plot([], [], "-", color="tab:red", label="rejected by ballistic fit")
    ax.plot([], [], "x", color="0.4", label="robust-fit outlier")
    ax.legend(loc="best", fontsize=8)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_tx_beam_consistency(tracks, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    scored = [t for t in tracks if "tx_beam_weighted_rms_deg" in t]
    scored.sort(key=lambda t: t.get("combined_rank", t.get("ballistic_rank", 999)))

    markers = ["o", "s", "^", "D", "v", "P", "X"]

    def track_style(track, idx):
        if not track.get("combined_reject", False):
            return "#009E73", markers[idx % len(markers)], "accepted"
        return "#CC3311", markers[idx % len(markers)], "rejected"

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0), constrained_layout=True)

    ax = axes[0, 0]
    for idx, track in enumerate(scored):
        color, marker, status = track_style(track, idx)
        order = np.argsort(track["tx_beam_time"])
        ax.plot(
            track["tx_beam_time"][order],
            track["tx_beam_angle_deg"][order],
            linestyle="None",
            marker=marker,
            ms=3.8,
            color=color,
            alpha=0.70,
            label=f"{hypothesis_label(track)} ({status})",
        )
        ax.text(
            track["tx_beam_time"][order][-1],
            track["tx_beam_angle_deg"][order][-1],
            hypothesis_label(track),
            color=color,
            fontsize=8,
        )
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Angular separation from active TX beam (deg)")
    ax.set_title("Per-pulse TX beam pointing consistency")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)

    ax = axes[0, 1]
    lobe_scored = [t for t in scored if "tx_lobe_distance_dc" in t]
    for idx, track in enumerate(lobe_scored):
        color, marker, status = track_style(track, idx)
        t = np.asarray(track["tx_beam_time"], dtype=np.float64)
        dist = np.asarray(track["tx_lobe_distance_dc"], dtype=np.float64)
        good = np.isfinite(t) & np.isfinite(dist)
        if np.count_nonzero(good) < 2:
            continue
        order = np.argsort(t[good])
        tt = t[good][order]
        dd = dist[good][order]
        ax.plot(
            tt,
            dd,
            linestyle="None",
            marker=marker,
            ms=3.8,
            color=color,
            alpha=0.70,
            label=f"{hypothesis_label(track)} ({status})",
        )
        ax.text(tt[-1], dd[-1], hypothesis_label(track), color=color, fontsize=8)
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Nearest TX lobe centroid distance")
    ax.set_title("Per-pulse nearest grating-lobe centroid distance")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)

    ax = axes[1, 0]
    x = np.arange(len(scored))
    labels = [hypothesis_label(t) for t in scored]
    rms = [t["tx_beam_weighted_rms_deg"] for t in scored]
    colors = [track_style(t, i)[0] for i, t in enumerate(scored)]
    bars = ax.bar(x, rms, color=colors, alpha=0.88, edgecolor="black", linewidth=0.7)
    for bar, track in zip(bars, scored):
        if track.get("combined_reject", False):
            bar.set_hatch("//")
    for xx, y in zip(x, rms):
        ax.text(xx, y, f"{y:.1f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_xlabel("Path hypothesis")
    ax.set_ylabel("SNR-weighted RMS beam angle (deg)")
    ax.set_title("TX beam pointing score")
    ax.grid(True, axis="y", alpha=0.25)

    ax = axes[1, 1]
    lobe_labels = [hypothesis_label(t) for t in lobe_scored]
    lobe_x = np.arange(len(lobe_scored))
    lobe_mean = [t.get("tx_lobe_snr_weighted_mean_dc", np.nan) for t in lobe_scored]
    lobe_rms = [t.get("tx_lobe_snr_weighted_rms_dc", np.nan) for t in lobe_scored]
    lobe_colors = [track_style(t, i)[0] for i, t in enumerate(lobe_scored)]
    bars_mean = ax.bar(lobe_x - 0.18, lobe_mean, width=0.36, color=lobe_colors, alpha=0.88, edgecolor="black", label="mean")
    bars_rms = ax.bar(lobe_x + 0.18, lobe_rms, width=0.36, color=lobe_colors, alpha=0.35, edgecolor="black", label="RMS")
    for bars in (bars_mean, bars_rms):
        for bar, track in zip(bars, lobe_scored):
            if track.get("combined_reject", False):
                bar.set_hatch("//")
    for x, y in zip(lobe_x, lobe_mean):
        if np.isfinite(y):
            ax.text(x - 0.18, y, f"{y:.3f}", ha="center", va="bottom", fontsize=8, rotation=90)
    ax.set_xticks(lobe_x)
    ax.set_xticklabels(lobe_labels, rotation=45, ha="right")
    ax.set_xlabel("Path hypothesis")
    ax.set_ylabel("Direction-cosine distance")
    ax.set_title("SNR-weighted nearest-lobe centroid distance")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_tx_array_snr_consistency(tracks, output: Path):
    """Plot full TX array-factor gain against range-corrected SNR for aliases."""
    output.parent.mkdir(parents=True, exist_ok=True)
    scored = [t for t in tracks if "tx_array_snr_rms_db" in t]
    scored.sort(key=lambda t: t.get("combined_rank", t.get("ballistic_rank", 999)))
    if not scored:
        return

    cmap = plt.get_cmap("viridis")

    def track_color(track):
        if track.get("combined_reject", False):
            return "tab:red"
        denom = max(1, len(scored) - 1)
        return cmap(track.get("combined_rank", track.get("ballistic_rank", 0)) / denom)

    fig, axes = plt.subplots(1, 3, figsize=(15.2, 5.2), constrained_layout=True)

    ax = axes[0]
    for track in scored:
        color = track_color(track)
        good = np.isfinite(track["tx_array_gain_db"]) & np.isfinite(track["tx_array_range_corrected_snr_db"])
        if np.count_nonzero(good) < 2:
            continue
        gain_rel = track["tx_array_gain_db"][good] - np.nanmedian(track["tx_array_gain_db"][good])
        snr_rel = track["tx_array_range_corrected_snr_db"][good] - np.nanmedian(
            track["tx_array_range_corrected_snr_db"][good]
        )
        ax.scatter(gain_rel, snr_rel, s=8, color=color, alpha=0.45)
        ax.text(
            np.nanmedian(gain_rel),
            np.nanmedian(snr_rel),
            hypothesis_label(track),
            color=color,
            fontsize=9,
        )
    lim = np.nanmax(np.abs(ax.axis())) if np.all(np.isfinite(ax.axis())) else 20.0
    lim = max(5.0, min(40.0, lim))
    ax.plot([-lim, lim], [-lim, lim], color="black", lw=1.0, ls="--")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Relative TX array-factor gain (dB)")
    ax.set_ylabel("Relative range-corrected SNR (dB)")
    ax.set_title("Amplitude consistency scatter")
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    ranks = [t.get("combined_rank", t.get("ballistic_rank", 999)) + 1 for t in scored]
    rms = [t["tx_array_snr_rms_db"] for t in scored]
    colors = [track_color(t) for t in scored]
    ax.bar(ranks, rms, color=colors, alpha=0.85)
    for x, y in zip(ranks, rms):
        if np.isfinite(y):
            ax.text(x, y, f"{y:.1f}", ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("Combined statistical candidate number")
    ax.set_ylabel("Array-factor/SNR RMS mismatch (dB)")
    ax.set_title("Lower is more amplitude-consistent")
    ax.grid(True, axis="y", alpha=0.25)

    ax = axes[2]
    corr = [t["tx_array_snr_corr"] for t in scored]
    gain_span = [t["tx_array_gain_span_db"] for t in scored]
    sc = ax.scatter(rms, corr, c=gain_span, cmap="plasma", s=90)
    for track, x, y in zip(scored, rms, corr):
        ax.text(x, y, hypothesis_label(track), fontsize=8)
    ax.axhline(0.0, color="black", lw=1.0, ls="--")
    ax.set_xlabel("RMS mismatch (dB)")
    ax.set_ylabel("SNR/gain correlation")
    ax.set_title("Diagnostic strength and consistency")
    ax.grid(True, alpha=0.25)
    fig.colorbar(sc, ax=ax, label="Predicted gain span (dB)")

    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_tx_array_snr_overlay(tracks, output: Path, n_best=3):
    """Line-overlay TX array factor and range-corrected SNR for best fits."""
    output.parent.mkdir(parents=True, exist_ok=True)
    ranked = [t for t in tracks if "tx_array_snr_rms_db" in t and "ballistic_rank" in t]
    ranked.sort(key=lambda t: t["ballistic_rank"])
    ranked = ranked[:n_best]
    if not ranked:
        return

    fig, axes = plt.subplots(len(ranked), 1, figsize=(10.5, 3.0 * len(ranked)), sharex=True, constrained_layout=True)
    if len(ranked) == 1:
        axes = [axes]
    colors = plt.get_cmap("tab10").colors

    for ax, track in zip(axes, ranked):
        t = np.asarray(track["tx_array_time"], dtype=np.float64)
        gain_db = np.asarray(track["tx_array_gain_db"], dtype=np.float64)
        snr_db = np.asarray(track["tx_array_range_corrected_snr_db"], dtype=np.float64)
        good = np.isfinite(t) & np.isfinite(gain_db) & np.isfinite(snr_db)
        if np.count_nonzero(good) < 2:
            continue

        t = t[good]
        gain_rel = gain_db[good] - np.nanmedian(gain_db[good])
        snr_rel = snr_db[good] - np.nanmedian(snr_db[good])
        order = np.argsort(t)
        color = colors[int(track["ballistic_rank"]) % len(colors)]

        ax.plot(t[order], gain_rel[order], "-", color=color, lw=1.6, label="TX array-factor gain")
        ax.plot(t[order], snr_rel[order], "-", color="black", lw=1.2, alpha=0.85, label="range-corrected SNR")
        ax.axhline(0.0, color="0.5", lw=0.8, ls="--")
        ax.set_ylabel("Relative level (dB)")
        ax.set_title(
            f"{hypothesis_label(track)}: ballistic rank {track['ballistic_rank']}, "
            f"candidate {track.get('combined_rank', track['ballistic_rank']) + 1}, "
            f"$\\chi^2_\\nu$={track['ballistic_reduced_chi2']:.3f}, "
            f"RMS={track['tx_array_snr_rms_db']:.2f} dB, "
            f"corr={track['tx_array_snr_corr']:.2f}"
        )
        ax.grid(True, alpha=0.25)
        ax.legend(loc="upper right", fontsize=8)

    axes[-1].set_xlabel("Time since cut start (s)")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_tx_array_snr_overlay_by_beam(tracks, output: Path, n_best=3):
    """Overlay TX gain and SNR separately for each active transmit beam."""
    output.parent.mkdir(parents=True, exist_ok=True)
    ranked = [t for t in tracks if "tx_array_snr_rms_db" in t and "ballistic_rank" in t]
    ranked.sort(key=lambda t: t["ballistic_rank"])
    ranked = ranked[:n_best]
    if not ranked:
        return

    beam_names = ["Zenith", "North", "East", "South", "West"]
    fig, axes = plt.subplots(
        len(ranked),
        len(beam_names),
        figsize=(18.0, 3.0 * len(ranked)),
        sharex=False,
        sharey=True,
        constrained_layout=True,
    )
    if len(ranked) == 1:
        axes = np.asarray([axes])
    colors = plt.get_cmap("tab10").colors

    for row, track in enumerate(ranked):
        t = np.asarray(track["tx_array_time"], dtype=np.float64)
        gain_db = np.asarray(track["tx_array_gain_db"], dtype=np.float64)
        snr_db = np.asarray(track["tx_array_range_corrected_snr_db"], dtype=np.float64)
        beam_id = np.asarray(track["tx_beam_id"], dtype=np.int64)
        good_all = np.isfinite(t) & np.isfinite(gain_db) & np.isfinite(snr_db)
        color = colors[int(track["ballistic_rank"]) % len(colors)]

        for beam in range(len(beam_names)):
            ax = axes[row, beam]
            good = good_all & (beam_id == beam)
            if np.count_nonzero(good) >= 2:
                tt = t[good]
                gain_rel = gain_db[good] - np.nanmedian(gain_db[good])
                snr_rel = snr_db[good] - np.nanmedian(snr_db[good])
                order = np.argsort(tt)
                ax.plot(tt[order], gain_rel[order], "-", color=color, lw=1.5, label="TX gain")
                ax.plot(tt[order], snr_rel[order], "-", color="black", lw=1.2, alpha=0.9, label="SNR")
            ax.axhline(0.0, color="0.55", lw=0.8, ls="--")
            ax.grid(True, alpha=0.25)
            if row == 0:
                ax.set_title(f"{beam}: {beam_names[beam]}")
            if beam == 0:
                ax.set_ylabel(
                    f"{hypothesis_label(track)}\n"
                    f"ballistic rank {track['ballistic_rank']}\n"
                    f"$\\chi^2_\\nu$={track['ballistic_reduced_chi2']:.2f}\n"
                    "Relative dB"
                )
            if row == len(ranked) - 1:
                ax.set_xlabel("Time since cut start (s)")
            if row == 0 and beam == len(beam_names) - 1:
                ax.legend(loc="upper right", fontsize=8)

    fig.suptitle("TX array-factor gain and range-corrected SNR by active transmit beam")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_tx_array_snr_overlay_separate_by_beam(tracks, output_dir: Path, sample_idx: int, n_best=3):
    """Write one TX-gain/SNR overlay file per hypothesis and transmit beam."""
    output_dir.mkdir(parents=True, exist_ok=True)
    ranked = [t for t in tracks if "tx_array_snr_rms_db" in t and "ballistic_rank" in t]
    ranked.sort(key=lambda t: t["ballistic_rank"])
    ranked = ranked[:n_best]
    beam_names = ["Zenith", "North", "East", "South", "West"]
    colors = plt.get_cmap("tab10").colors
    paths = []

    for track in ranked:
        t = np.asarray(track["tx_array_time"], dtype=np.float64)
        gain_db = np.asarray(track["tx_array_gain_db"], dtype=np.float64)
        snr_db = np.asarray(track["tx_array_range_corrected_snr_db"], dtype=np.float64)
        beam_id = np.asarray(track["tx_beam_id"], dtype=np.int64)
        good_all = np.isfinite(t) & np.isfinite(gain_db) & np.isfinite(snr_db)
        color = colors[int(track["ballistic_rank"]) % len(colors)]

        for beam, beam_name in enumerate(beam_names):
            good = good_all & (beam_id == beam)
            if np.count_nonzero(good) < 2:
                continue
            tt = t[good]
            gain_rel = gain_db[good] - np.nanmedian(gain_db[good])
            snr_rel = snr_db[good] - np.nanmedian(snr_db[good])
            order = np.argsort(tt)

            fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
            ax.plot(tt[order], gain_rel[order], "-", color=color, lw=1.8, label="TX array-factor gain")
            ax.plot(tt[order], snr_rel[order], "-", color="black", lw=1.5, alpha=0.9, label="range-corrected SNR")
            ax.axhline(0.0, color="0.55", lw=0.9, ls="--")
            ax.set_xlabel("Time since cut start (s)")
            ax.set_ylabel("Relative level within beam (dB)")
            ax.set_title(
                f"Sample {sample_idx}, {hypothesis_label(track)}, "
                f"ballistic rank {track['ballistic_rank']}, "
                f"TX beam {beam}: {beam_name}\n"
                f"$\\chi^2_\\nu$={track['ballistic_reduced_chi2']:.3f}, "
                f"beam RMS={track['tx_array_snr_rms_db_by_beam'][beam]:.2f} dB, "
                f"beam corr={track['tx_array_snr_corr_by_beam'][beam]:.2f}, "
                f"n={track['tx_array_snr_n_by_beam'][beam]}"
            )
            ax.grid(True, alpha=0.25)
            ax.legend(loc="best", fontsize=9)
            path = output_dir / (
                f"pansy_interferometer_tx_array_snr_overlay_{sample_idx}"
                f"_{hypothesis_label(track).lower()}_rank{track['ballistic_rank']}_beam{beam}_{beam_name.lower()}.png"
            )
            fig.savefig(path, dpi=220)
            plt.close(fig)
            paths.append(path)
    return paths


def state_at_first_detection_from_ballistic_track(track):
    """Return fitted ENU state at the first-detection fit epoch."""
    params = np.asarray(
        track["selection_params"] if "selection_params" in track else track["ballistic_params"],
        dtype=np.float64,
    )
    t_ref = float(track.get("selection_t_ref_s", track.get("ballistic_t_ref_s", 0.0)))
    state_enu_m_mps = np.asarray(params[:6], dtype=np.float64)
    return state_enu_m_mps, t_ref, float(state_enu_m_mps[2] / 1e3)


def orbit_for_candidate_track(track, sample_epoch_unix):
    state_enu, epoch_offset_s, alt_km = state_at_first_detection_from_ballistic_track(track)
    epoch_unix = float(sample_epoch_unix) + epoch_offset_s
    state_gcrs = pbal.enu_state_to_gcrs(state_enu, epoch_unix)
    r_km, v_km_s = porb.heliocentric_state_from_gcrs(state_gcrs, epoch_unix)
    kepler = porb.kepler_from_state(r_km, v_km_s)
    return {
        "epoch_unix": epoch_unix,
        "alt_km": alt_km,
        "state_epoch": "first_detection",
        "state_enu_m_mps": state_enu,
        "state_gcrs_m_mps": state_gcrs,
        "heliocentric_state_km_kms": np.concatenate([r_km, v_km_s]),
        "kepler": kepler,
    }


def first_detection_gcrs_state_samples(track, sample_epoch_unix, n_samples=300, seed=20260615):
    cov = np.asarray(track.get("selection_parameter_covariance", track.get("ballistic_parameter_covariance", np.nan)), dtype=np.float64)
    params = np.asarray(track.get("selection_params", track.get("ballistic_params", np.nan)), dtype=np.float64)
    if params.shape not in ((6,), (7,)) or cov.shape != (len(params), len(params)) or not np.all(np.isfinite(cov)) or not np.all(np.isfinite(params)):
        return np.empty((0, 6), dtype=np.float64)
    cov = 0.5 * (cov + cov.T)
    rng = np.random.default_rng(seed + int(track.get("hypothesis_id", 0)))
    try:
        draws = rng.multivariate_normal(params, cov, size=n_samples, check_valid="ignore")
    except np.linalg.LinAlgError:
        return np.empty((0, 6), dtype=np.float64)
    t_ref = float(track.get("selection_t_ref_s", track.get("ballistic_t_ref_s", 0.0)))
    epoch_unix = float(sample_epoch_unix) + t_ref
    states = []
    for sample_params in draws:
        if not np.all(np.isfinite(sample_params)):
            continue
        sample_params[2] = np.clip(sample_params[2], 20e3, 220e3)
        if sample_params.shape == (7,):
            sample_params[6] = np.clip(sample_params[6], -4.0, 6.0)
        try:
            state_gcrs = pbal.enu_state_to_gcrs(sample_params[:6], epoch_unix)
        except Exception:
            continue
        if np.all(np.isfinite(state_gcrs)):
            states.append(state_gcrs)
    return np.asarray(states, dtype=np.float64)


def orbit_uncertainty_for_candidate_track(track, sample_epoch_unix, nominal_kepler, n_samples=300, seed=20260615):
    t_ref = float(track.get("selection_t_ref_s", track.get("ballistic_t_ref_s", 0.0)))
    epoch_unix = float(sample_epoch_unix) + t_ref
    state_samples = first_detection_gcrs_state_samples(track, sample_epoch_unix, n_samples=n_samples, seed=seed)
    samples = []
    for state_gcrs in state_samples:
        try:
            r_km, v_km_s = porb.heliocentric_state_from_gcrs(state_gcrs, epoch_unix)
            kep = porb.kepler_from_state(r_km, v_km_s)
        except Exception:
            continue
        if np.all(np.isfinite(kep)):
            samples.append(kep)
    if not samples:
        return np.full(7, np.nan), 0
    std, _cov = porb.summarize_samples(np.asarray(samples), np.asarray(nominal_kepler, dtype=np.float64))
    return std, len(samples)


def write_candidate_orbit_state_h5(orbits, output_path, sample_epoch_unix, n_samples=300):
    if output_path is None or not orbits:
        return
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_path, "w") as h:
        h.attrs["source_program"] = "plot_interferometric_disambiguation.py"
        h.attrs["sample_epoch_unix"] = float(sample_epoch_unix)
        h.attrs["reference_frame"] = "GCRS"
        h.attrs["state_epoch"] = "first_detection"
        h.attrs["state_units"] = "m,m,m,m/s,m/s,m/s"
        h.attrs["ballistic_coefficient_parameter"] = "log10_beta_kg_m2"
        h.attrs["state_source_parameter_preference"] = "selection_params, falling back to ballistic_params"
        h.attrs["orbit_note"] = "Input states for DASST zenithal-attraction removal; no above-atmosphere back-propagation."
        for track, orbit in orbits:
            label = hypothesis_label(track)
            grp = h.create_group(label)
            params = np.asarray(track.get("selection_params", track.get("ballistic_params", np.full(7, np.nan))), dtype=np.float64)
            param_std = np.asarray(track.get("selection_parameter_std", track.get("ballistic_parameter_std", np.full(len(params), np.nan))), dtype=np.float64)
            nsamp = int(n_samples) if int(track.get("combined_rank", -1)) == 0 else 0
            state_samples = first_detection_gcrs_state_samples(track, sample_epoch_unix, n_samples=nsamp)
            grp.attrs["hypothesis_label"] = label
            grp.attrs["hypothesis_id"] = int(track.get("hypothesis_id", -1))
            grp.attrs["candidate_number"] = int(candidate_number(track))
            grp.attrs["combined_rank"] = int(track.get("combined_rank", -1))
            grp.attrs["combined_score"] = float(track.get("combined_score", np.nan))
            grp.attrs["selection_model_type"] = str(track.get("selection_model_type", ""))
            grp.attrs["interstellar_alias_plausible"] = bool(track.get("interstellar_alias_plausible", False))
            grp.attrs["interstellar_alias_plausibility_model"] = str(track.get("interstellar_alias_plausibility_model", ""))
            grp.attrs["interstellar_alias_plausibility_redchi"] = float(track.get("interstellar_alias_plausibility_redchi", np.nan))
            grp.attrs["interstellar_alias_plausibility_redchi_threshold"] = float(track.get("interstellar_alias_plausibility_redchi_threshold", np.nan))
            grp.attrs["fixed_velocity_reduced_chi2"] = float(track.get("fixed_velocity_reduced_chi2", np.nan))
            grp.attrs["fixed_am_reduced_chi2"] = float(track.get("fixed_am_reduced_chi2", np.nan))
            grp.attrs["tx_beam_snr_weighted_mean_dc"] = float(track.get("tx_beam_snr_weighted_mean_dc", np.nan))
            grp.attrs["tx_beam_snr_weighted_rms_dc"] = float(track.get("tx_beam_snr_weighted_rms_dc", np.nan))
            grp.attrs["tx_beam_weighted_mean_deg"] = float(track.get("tx_beam_weighted_mean_deg", np.nan))
            grp.attrs["tx_beam_weighted_rms_deg"] = float(track.get("tx_beam_weighted_rms_deg", np.nan))
            grp.attrs["tx_lobe_snr_weighted_mean_dc"] = float(track.get("tx_lobe_snr_weighted_mean_dc", np.nan))
            grp.attrs["tx_lobe_snr_weighted_rms_dc"] = float(track.get("tx_lobe_snr_weighted_rms_dc", np.nan))
            grp.attrs["tx_lobe_p90_dc"] = float(track.get("tx_lobe_p90_dc", np.nan))
            grp.attrs["first_alt_km"] = float(orbit["alt_km"])
            grp.attrs["epoch_unix"] = float(orbit["epoch_unix"])
            grp.attrs["log10_beta_kg_m2"] = float(params[6]) if params.shape == (7,) else np.nan
            grp.attrs["sigma_log10_beta"] = float(param_std[6]) if param_std.shape == (7,) else np.nan
            grp.create_dataset("state_gcrs_m_mps", data=np.asarray(orbit["state_gcrs_m_mps"], dtype=np.float64))
            grp.create_dataset("state_gcrs_samples_m_mps", data=state_samples)
            grp.create_dataset("selection_params", data=params)
            grp.create_dataset("selection_parameter_covariance", data=np.asarray(track.get("selection_parameter_covariance", np.full((len(params), len(params)), np.nan)), dtype=np.float64))
            fit_points = np.asarray(track.get("fit_points", []), dtype=np.float64)
            if fit_points.ndim == 2 and fit_points.shape[1] == 3 and len(fit_points):
                path_range_km = np.linalg.norm(fit_points, axis=1)
                direction = np.full_like(fit_points, np.nan, dtype=np.float64)
                good_range = path_range_km > 0.0
                direction[good_range, 0] = fit_points[good_range, 0] / path_range_km[good_range]
                direction[good_range, 1] = fit_points[good_range, 1] / path_range_km[good_range]
                direction[good_range, 2] = -fit_points[good_range, 2] / path_range_km[good_range]
                path_payload = [
                    ("path_t_rel_s", track.get("fit_t")),
                    ("path_position_enu_km", fit_points),
                    ("path_direction_cosines_uvw", direction),
                    ("path_range_km", path_range_km),
                    ("path_snr", track.get("tx_beam_snr")),
                    ("path_beam_id", track.get("tx_beam_id")),
                    ("path_selection_keep", track.get("selection_keep")),
                ]
            else:
                path_payload = []
            for name, value in path_payload:
                if value is not None:
                    grp.create_dataset(name, data=np.asarray(value))
            if "ballistic_params" in track:
                grp.create_dataset("ballistic_params", data=np.asarray(track["ballistic_params"], dtype=np.float64))
                grp.create_dataset("ballistic_parameter_covariance", data=np.asarray(track.get("ballistic_parameter_covariance", np.full((7, 7), np.nan)), dtype=np.float64))
    print(f"candidate_orbit_state_h5 {output_path}")


def print_candidate_orbits(tracks, sample_epoch_unix, n_candidates=3, state_h5_path=None, n_samples=300):
    """Compute and print preliminary orbital elements for top-ranked aliases."""
    ranked = sorted(
        [t for t in tracks if "combined_score" in t and ("selection_params" in t or "ballistic_params" in t)],
        key=lambda t: t["combined_rank"],
    )
    top = ranked[:n_candidates]
    if not top:
        return []

    print(
        "candidate_orbit_columns "
        "hypothesis_id candidate_number combined_rank combined_score log10_beta_kg_m2 first_alt_km v_gcrs_km_s "
        "a_au e i_deg raan_deg argp_deg nu_deg q_au"
    )
    print(
        "candidate_orbit_uncertainty_columns "
        "hypothesis_id n_samples sigma_log10_beta sigma_a_au sigma_e sigma_i_deg sigma_raan_deg "
        "sigma_argp_deg sigma_nu_deg sigma_q_au"
    )
    orbits = []
    for track in top:
        orbit = orbit_for_candidate_track(track, sample_epoch_unix)
        track["candidate_orbit"] = orbit
        kep = orbit["kepler"]
        nsamp = int(n_samples) if int(track.get("combined_rank", -1)) == 0 else 0
        kep_std, n_unc = orbit_uncertainty_for_candidate_track(track, sample_epoch_unix, kep, n_samples=nsamp)
        orbit["kepler_std"] = kep_std
        orbit["n_uncertainty_samples"] = int(n_unc)
        gcrs_speed_km_s = np.linalg.norm(np.asarray(orbit["state_gcrs_m_mps"][3:], dtype=np.float64)) / 1e3
        params = np.asarray(track.get("selection_params", track.get("ballistic_params", np.full(7, np.nan))), dtype=np.float64)
        param_std = np.asarray(track.get("selection_parameter_std", track.get("ballistic_parameter_std", np.full(len(params), np.nan))), dtype=np.float64)
        log10_beta = float(params[6]) if params.shape == (7,) else float("nan")
        sigma_log10_beta = float(param_std[6]) if param_std.shape == (7,) else float("nan")
        print(
            "candidate_orbit",
            hypothesis_label(track),
            candidate_number(track),
            track["combined_rank"],
            f"{track['combined_score']:.3f}",
            f"{log10_beta:.3f}",
            f"{orbit['alt_km']:.3f}",
            f"{gcrs_speed_km_s:.3f}",
            f"{kep[0]:.6g}",
            f"{kep[1]:.6f}",
            f"{kep[2]:.3f}",
            f"{kep[3]:.3f}",
            f"{kep[4]:.3f}",
            f"{kep[5]:.3f}",
            f"{kep[6]:.6g}",
        )
        print(
            "candidate_orbit_uncertainty",
            hypothesis_label(track),
            int(n_unc),
            f"{sigma_log10_beta:.3f}",
            f"{kep_std[0]:.6g}",
            f"{kep_std[1]:.6f}",
            f"{kep_std[2]:.3f}",
            f"{kep_std[3]:.3f}",
            f"{kep_std[4]:.3f}",
            f"{kep_std[5]:.3f}",
            f"{kep_std[6]:.6g}",
        )
        orbits.append((track, orbit))
    write_candidate_orbit_state_h5(orbits, state_h5_path, sample_epoch_unix, n_samples=n_samples)
    return orbits


def main():
    parser = argparse.ArgumentParser(description="Plot PANSY interferometric alias disambiguation memo figures.")
    parser.add_argument("--sample-idx", type=int, default=1746494797212678)
    parser.add_argument("--cut-dir", type=Path, default=Path("data/metadata/cut"))
    parser.add_argument("--output-dir", type=Path, default=Path("../pansy_paper/memos/figures"))
    parser.add_argument("--grid-n", type=int, default=501)
    parser.add_argument("--coherence-threshold", type=float, default=0.80)
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    parser.add_argument("--linearity-angular-sigma-deg", type=float, default=0.25)
    parser.add_argument("--linearity-range-sigma-km", type=float, default=0.15)
    parser.add_argument("--linearity-p-threshold", type=float, default=0.01)
    parser.add_argument("--ballistic-p-threshold", type=float, default=0.01)
    parser.add_argument("--orbit-state-h5", type=Path, default=None)
    parser.add_argument("--dasst-orbit-h5", type=Path, default=None)
    parser.add_argument("--orbit-metadata-dir", type=Path, default=Path("data/metadata/orbit"))
    parser.add_argument("--run-dasst", action="store_true", help="Run DASST orbit determination after writing candidate state HDF5.")
    parser.add_argument("--orbit-samples", type=int, default=300, help="Number of covariance samples for DASST orbit uncertainty.")
    parser.add_argument("--overview-only", action="store_true", help="Only write the final disambiguation overview and orbit HDF5 products.")
    parser.add_argument("--path-max-tracks", type=int, default=24, help="Maximum direction-cosine path hypotheses to extract.")
    parser.add_argument("--path-trials", type=int, default=5000, help="Random two-point path trials per extracted hypothesis.")
    parser.add_argument("--max-peaks-per-pulse", type=int, default=32, help="Maximum interferometer local maxima retained per pulse.")
    parser.add_argument(
        "--alias-fit-models",
        choices=("fixed_velocity", "fixed_velocity_and_drag", "drag_and_fixed_am"),
        default="fixed_velocity",
        help="Models fit to every interferometric alias for ranking.",
    )
    parser.add_argument("--skip-secondary-models", action="store_true", help="Skip final winning-alias fixed-A/m and shrinking-radius refits.")
    parser.add_argument("--skip-orbit-products", action="store_true", help="Skip candidate orbit computation and orbit-state HDF5 output.")
    parser.add_argument("--recompute-cut-observables", action="store_true", help="Recompute full per-pulse range-Doppler observables instead of using cached cut detections.")
    parser.add_argument("--tx-phase-quality-h5", type=Path, default=None, help="HDF5 TX cross-phase quality table produced by tx_phase_quality.py.")
    parser.add_argument("--require-good-tx-phase", action="store_true", help="Reject events whose nearest TX cross-phase sample is missing or misaligned.")
    parser.add_argument("--tx-phase-max-age-s", type=float, default=None, help="Override maximum nearest TX phase sample age in seconds.")
    parser.add_argument("--compute-missing-tx-phase", action="store_true", help="Compute TX cross-phase from raw voltage when no nearby txphase metadata exists.")
    parser.add_argument("--tx-phase-raw-voltage-dir", type=Path, default=Path(pc.raw_voltage_dir), help="Digital RF raw voltage directory for missing TX phase fallback.")
    parser.add_argument("--tx-phase-tx-metadata-dir", type=Path, default=Path(pc.tx_metadata_dir), help="TX metadata directory for missing TX phase fallback.")
    parser.add_argument("--tx-phase-fallback-search-radius-s", type=float, default=txpq.DEFAULT_FALLBACK_SEARCH_RADIUS_S, help="Search radius for TX pulse metadata when computing missing TX phase.")
    parser.add_argument("--tx-phase-fallback-samples", type=int, default=txpq.DEFAULT_FALLBACK_TX_SAMPLES, help="Raw TX samples used for missing TX phase cross-correlation.")
    parser.add_argument("--stage-timing", action="store_true", help="Print wall-clock timing for major processing stages.")
    args = parser.parse_args()

    stage_t0 = time.perf_counter()
    stage_last = stage_t0

    def stage(label):
        nonlocal stage_last
        if not args.stage_timing:
            return
        now = time.perf_counter()
        print(f"stage_timing {label} dt_s {now - stage_last:.3f} elapsed_s {now - stage_t0:.3f}", flush=True)
        stage_last = now

    disambiguation_summary_out = args.output_dir / f"pansy_interferometer_disambiguation_summary_{args.sample_idx}.png"
    diagnostics_h5 = args.output_dir / f"pansy_disambiguation_diagnostics_{args.sample_idx}.h5"
    tx_quality = None
    cut = None
    if args.tx_phase_quality_h5 is not None:
        try:
            tx_quality = txpq.quality_for_sample(
                args.tx_phase_quality_h5,
                args.sample_idx,
                max_age_s=args.tx_phase_max_age_s,
                compute_if_missing=args.compute_missing_tx_phase,
                raw_voltage_dir=args.tx_phase_raw_voltage_dir,
                tx_metadata_dir=args.tx_phase_tx_metadata_dir,
                fallback_search_radius_s=args.tx_phase_fallback_search_radius_s,
                fallback_tx_samples=args.tx_phase_fallback_samples,
            )
        except Exception as exc:
            tx_quality = {
                "available": False,
                "good": False,
                "reason": f"tx_phase_quality_read_failed:{type(exc).__name__}",
            }
    elif args.require_good_tx_phase:
        tx_quality = {
            "available": False,
            "good": False,
            "reason": "tx_phase_quality_h5_not_provided",
        }
    if (
        args.require_good_tx_phase
        and args.compute_missing_tx_phase
        and tx_quality
        and not tx_quality.get("good", False)
        and tx_quality.get("reason") != "tx_phase_mismatch"
    ):
        try:
            cut = load_cut(args.cut_dir, args.sample_idx)
            cut_quality = cut_tx_waveform_quality(cut)
            if cut_quality.get("good", False):
                inherited = {
                    key: value
                    for key, value in tx_quality.items()
                    if key.startswith("reference_") or key in {"max_age_s"}
                }
                tx_quality = {**inherited, **cut_quality}
            else:
                tx_quality = {
                    **tx_quality,
                    "computed_available": bool(cut_quality.get("available", False)),
                    "computed_reason": str(cut_quality.get("reason", "cut_tx_waveform_failed")),
                    "cut_tx_waveform_good": bool(cut_quality.get("good", False)),
                    "cut_tx_waveform_phase_std_deg": float(cut_quality.get("rms_diff_deg", np.nan)),
                    "cut_tx_waveform_median_corr": float(cut_quality.get("mean_amplitude", np.nan)),
                }
        except Exception as exc:
            tx_quality = {
                **(tx_quality or {}),
                "computed_available": False,
                "computed_reason": f"cut_tx_waveform_read_failed:{type(exc).__name__}",
            }
    if args.require_good_tx_phase and not bool(tx_quality and tx_quality.get("good", False)):
        reason = str((tx_quality or {}).get("reason", "bad_tx_phase_alignment"))
        details = (
            f"TX cross-phase quality gate failed: {reason}; "
            f"nearest_sample_idx={(tx_quality or {}).get('nearest_sample_idx', 'n/a')}; "
            f"age_s={(tx_quality or {}).get('age_s', 'n/a')}; "
            f"max_abs_diff_deg={(tx_quality or {}).get('max_abs_diff_deg', 'n/a')}; "
            f"rms_diff_deg={(tx_quality or {}).get('rms_diff_deg', 'n/a')}"
        )
        plot_rejected_event_summary(
            disambiguation_summary_out,
            args.sample_idx,
            "bad or missing TX cross-phase alignment",
            details,
        )
        write_rejected_disambiguation_diagnostics_h5(
            diagnostics_h5,
            args.sample_idx,
            "bad_tx_phase_alignment",
            details,
            extra_attrs=tx_phase_quality_attrs(tx_quality),
        )
        print(disambiguation_summary_out)
        print(diagnostics_h5)
        print("single_echo_high_coherence_peaks 0")
        print("all_high_coherence_candidates 0")
        print("tx_lobe_centroid_counts 0 0 0 0 0")
        print("path_hypotheses 0")
        print("path_hypotheses_kept 0")
        print("path_hypotheses_below_horizon 0")
        print("path_hypotheses_wrapping 0")
        print("path_hypotheses_linearity_rejected 0")
        print("path_hypotheses_after_linearity 0")
        print("path_hypotheses_low_detection_altitude_rejected 0")
        print("path_hypotheses_after_detection_altitude 0")
        print("path_hypotheses_descent_rejected 0")
        print("path_hypotheses_after_descent 0")
        print("path_hypotheses_low_start_altitude_rejected 0")
        print("ballistic_sigma_pos_km nan")
        print("ballistic_sigma_dop_km_s nan")
        print("tx_pointing_score_used False")
        print("combined_hypotheses_rejected 0")
        print("combined_hypotheses_after 0")
        print("event_rejected_reason bad_tx_phase_alignment")
        return
    try:
        if cut is None:
            cut = load_cut(args.cut_dir, args.sample_idx)
        if args.recompute_cut_observables:
            obs_all = recompute_cut_observables(cut)
        else:
            obs_all = cached_cut_observables(cut)
    except Exception as exc:
        plot_rejected_event_summary(
            disambiguation_summary_out,
            args.sample_idx,
            "cut observable recomputation failed",
            f"{type(exc).__name__}: {exc}",
        )
        write_rejected_disambiguation_diagnostics_h5(
            diagnostics_h5,
            args.sample_idx,
            "cut_observable_recomputation_failed",
            f"{type(exc).__name__}: {exc}",
        )
        print(disambiguation_summary_out)
        print(diagnostics_h5)
        print("single_echo_high_coherence_peaks 0")
        print("all_high_coherence_candidates 0")
        print("tx_lobe_centroid_counts 0 0 0 0 0")
        print("path_hypotheses 0")
        print("path_hypotheses_kept 0")
        print("path_hypotheses_below_horizon 0")
        print("path_hypotheses_wrapping 0")
        print("path_hypotheses_linearity_rejected 0")
        print("path_hypotheses_after_linearity 0")
        print("path_hypotheses_low_detection_altitude_rejected 0")
        print("path_hypotheses_after_detection_altitude 0")
        print("path_hypotheses_descent_rejected 0")
        print("path_hypotheses_after_descent 0")
        print("path_hypotheses_low_start_altitude_rejected 0")
        print("ballistic_sigma_pos_km nan")
        print("ballistic_sigma_dop_km_s nan")
        print("tx_pointing_score_used False")
        print("combined_hypotheses_rejected 0")
        print("combined_hypotheses_after 0")
        print("event_rejected_reason cut_observable_recomputation_failed")
        return
    stage("load_and_recompute_observables")
    good = obs_all["snr"] > args.snr_threshold
    obs = {
        key: val[good] if isinstance(val, np.ndarray) and len(val) == len(good) else val
        for key, val in obs_all.items()
    }
    if len(obs["snr"]) == 0:
        plot_empty_event_summary(disambiguation_summary_out, args.sample_idx, int(len(obs_all["snr"])), args.snr_threshold)
        write_rejected_disambiguation_diagnostics_h5(
            diagnostics_h5,
            args.sample_idx,
            "no_measurements_above_snr_threshold",
            f"{int(len(obs_all['snr']))} raw cut samples before SNR threshold {args.snr_threshold:g}",
        )
        print(disambiguation_summary_out)
        print(diagnostics_h5)
        print("single_echo_high_coherence_peaks 0")
        print("all_high_coherence_candidates 0")
        print("tx_lobe_centroid_counts 0 0 0 0 0")
        print("path_hypotheses 0")
        print("path_hypotheses_kept 0")
        print("path_hypotheses_below_horizon 0")
        print("path_hypotheses_wrapping 0")
        print("path_hypotheses_linearity_rejected 0")
        print("path_hypotheses_after_linearity 0")
        print("path_hypotheses_low_detection_altitude_rejected 0")
        print("path_hypotheses_after_detection_altitude 0")
        print("path_hypotheses_descent_rejected 0")
        print("path_hypotheses_after_descent 0")
        print("path_hypotheses_low_start_altitude_rejected 0")
        print("ballistic_sigma_pos_km nan")
        print("ballistic_sigma_dop_km_s nan")
        print("tx_pointing_score_used False")
        print("combined_hypotheses_rejected 0")
        print("combined_hypotheses_after 0")
        print("event_rejected_reason no_measurements_above_snr_threshold")
        return

    range_time_segments = split_observations_by_range_time(obs, min_points=3)
    range_time_component_lengths = [int(len(segment)) for segment in range_time_segments]
    selected_range_time_component = -1
    if len(range_time_segments) >= 2 and max(range_time_component_lengths, default=0) >= 3:
        weights = [float(np.nansum(np.asarray(obs["snr"], dtype=np.float64)[segment])) for segment in range_time_segments]
        selected_range_time_component = int(np.argmax(weights))
        obs = subset_pulse_observations(obs, range_time_segments[selected_range_time_component])
    stage("snr_and_range_time_filter")

    antpos = pint.get_antpos()
    ch_pairs = np.asarray(list(itertools.combinations(np.arange(7), 2)))
    dmat = pint.pair_mat(ch_pairs, antpos)
    phasecal = pint.get_phasecal()
    u, v, w, valid = horizon_grid(args.grid_n)
    steering, _uvw = steering_matrix(dmat, u, v, w, valid)
    tx_gain_maps = precompute_tx_array_gain_maps(u, v, w, valid)
    tx_lobe_centroids, tx_lobe_smoothed_maps = tx_grating_lobe_centroids_from_gain_maps(tx_gain_maps, u, v, valid)
    stage("grid_and_tx_maps")

    strongest = int(np.argmax(obs["snr"]))
    single = coherence_map(obs["xc"][strongest], int(obs["beam_id"][strongest]), phasecal, ch_pairs, steering, valid, u.shape)
    peak_ij = local_coherence_peaks(single, args.coherence_threshold, max_peaks=args.max_peaks_per_pulse)

    candidates = []
    t_rel = obs["tx_idx"] / 1e6 - obs["tx_idx"][0] / 1e6
    for i in range(len(obs["snr"])):
        coh = coherence_map(obs["xc"][i], int(obs["beam_id"][i]), phasecal, ch_pairs, steering, valid, u.shape)
        ii, jj = local_coherence_peaks(coh, args.coherence_threshold, max_peaks=args.max_peaks_per_pulse)
        for row, col in zip(ii, jj):
            candidates.append(
                {
                    "u": float(u[row, col]),
                    "v": float(v[row, col]),
                    "w": float(w[row, col]),
                    "grid_row": int(row),
                    "grid_col": int(col),
                    "coherence": float(coh[row, col]),
                    "t_rel": float(t_rel[i]),
                    "pulse": int(i),
                    "range_km": float(obs["range_km"][i]),
                    "doppler_mps": float(obs["doppler_mps"][i]),
                    "snr": float(obs["snr"][i]),
                    "beam_id": int(obs["beam_id"][i]),
                    "zpp": complex(obs["zpp"][i]) if "zpp" in obs else complex(np.nan, np.nan),
                    "zpp_phase_rad": float(obs["zpp_phase_rad"][i]) if "zpp_phase_rad" in obs else np.nan,
                    "zpp_phase_dt_s": float(obs["zpp_phase_dt_s"][i]) if "zpp_phase_dt_s" in obs else np.nan,
                    "zpp_doppler_mps": float(obs["zpp_doppler_mps"][i]) if "zpp_doppler_mps" in obs else np.nan,
                    "zpp_coarse_doppler_mps": float(obs["zpp_coarse_doppler_mps"][i]) if "zpp_coarse_doppler_mps" in obs else np.nan,
                    "zpp_coarse_phase_rad": float(obs["zpp_coarse_phase_rad"][i]) if "zpp_coarse_phase_rad" in obs else np.nan,
                    "zpp_prev_t_rel": float(obs["zpp_prev_t_rel"][i]) if "zpp_prev_t_rel" in obs else np.nan,
                    "zpp_prev_pulse": int(obs["zpp_prev_pulse"][i]) if "zpp_prev_pulse" in obs else -1,
                }
            )
    stage("coherence_candidates")

    single_out = args.output_dir / f"pansy_interferometer_single_echo_{args.sample_idx}.png"
    all_out = args.output_dir / f"pansy_interferometer_alias_candidates_{args.sample_idx}.png"
    visibility_out = args.output_dir / f"pansy_interferometer_visibility_rejections_{args.sample_idx}.png"
    linearity_out = args.output_dir / f"pansy_interferometer_linearity_rejections_{args.sample_idx}.png"
    descent_out = args.output_dir / f"pansy_interferometer_descent_rejections_{args.sample_idx}.png"
    ballistic_out = args.output_dir / f"pansy_interferometer_ballistic_ranking_{args.sample_idx}.png"
    dasst_orbit_h5 = args.dasst_orbit_h5 or (args.output_dir / f"pansy_candidate_orbits_dasst_{args.sample_idx}.h5")
    ballistic_horizontal_out = args.output_dir / f"pansy_interferometer_ballistic_horizontal_ranking_{args.sample_idx}.png"
    tx_beam_out = args.output_dir / f"pansy_interferometer_tx_beam_consistency_{args.sample_idx}.png"
    tx_array_snr_out = args.output_dir / f"pansy_interferometer_tx_array_snr_consistency_{args.sample_idx}.png"
    tx_array_snr_overlay_out = args.output_dir / f"pansy_interferometer_tx_array_snr_overlay_{args.sample_idx}.png"
    tx_array_snr_overlay_by_beam_out = args.output_dir / f"pansy_interferometer_tx_array_snr_overlay_by_beam_{args.sample_idx}.png"
    tx_lobe_centroids_out = args.output_dir / "pansy_tx_grating_lobe_centroids.png"
    antennas_out = args.output_dir / "pansy_interferometer_antenna_positions.png"
    if not args.overview_only:
        plot_antenna_positions(antennas_out)
        plot_tx_grating_lobe_centroids(u, v, valid, tx_lobe_smoothed_maps, tx_lobe_centroids, tx_lobe_centroids_out)
        plot_single_echo(u, v, single, peak_ij, obs, strongest, single_out)
        plot_all_candidates(u, v, obs, candidates, all_out)
    tracks = complete_tracks_to_common_pulses(
        fit_candidate_tracks(
            candidates,
            max_tracks=args.path_max_tracks,
            n_trials=args.path_trials,
        ),
        candidates,
    )
    stage("path_hypotheses")
    t_start = float(np.min(t_rel))
    t_end = float(np.max(t_rel))
    tracks = [classify_track_visibility(track, candidates, t_start, t_end) for track in tracks]
    tracks = [
        classify_track_descent(track)
        if track["reason"] == "kept"
        else track
        for track in tracks
    ]
    tracks = [
        classify_track_linearity(
            track,
            candidates,
            angular_sigma_deg=args.linearity_angular_sigma_deg,
            range_sigma_km=args.linearity_range_sigma_km,
            p_threshold=args.linearity_p_threshold,
        )
        if track["reason"] == "kept"
        else track
        for track in tracks
    ]
    stage("path_diagnostics")
    rho_of_alt_m, _msis_meta = pbal.density_interpolator(args.sample_idx / 1e6)
    if args.alias_fit_models in {"fixed_velocity_and_drag", "drag_and_fixed_am"}:
        sigma_pos, sigma_dop = fit_ballistic_survivors(
            tracks,
            candidates,
            event_epoch_unix=args.sample_idx / 1e6,
            p_threshold=args.ballistic_p_threshold,
        )
        stage("msis_drag_fits")
        fit_fixed_velocity_survivors(
            tracks,
            candidates,
            rho_of_alt_m,
            sigma_pos,
            sigma_dop,
            fit_fixed_am=args.alias_fit_models == "drag_and_fixed_am",
        )
        stage("fixed_velocity_alias_fits")
    else:
        sigma_pos, sigma_dop = 0.5, 1.0
        fit_fixed_velocity_survivors(tracks, candidates, rho_of_alt_m, sigma_pos, sigma_dop, fit_fixed_am=False)
        stage("fixed_velocity_alias_fits")
    score_tx_beam_consistency(tracks, candidates)
    score_candidate_diagnostics(tracks, candidates, tx_gain_maps=tx_gain_maps, tx_lobe_centroids=tx_lobe_centroids)
    tx_sigma = score_combined_hypotheses(tracks)
    stage("scoring")
    if not args.skip_secondary_models:
        fit_winning_physics_trajectory_model(tracks, args.sample_idx / 1e6, sigma_pos, sigma_dop)
        stage("winning_secondary_models")
    best = sorted([t for t in tracks if "combined_score" in t], key=lambda t: t["combined_rank"])
    response_selected = selected_pulse_interferometer_response(
        candidates,
        obs,
        phasecal,
        ch_pairs,
        steering,
        valid,
        u.shape,
    )
    state_h5 = args.orbit_state_h5 or (args.output_dir / f"pansy_candidate_orbit_states_{args.sample_idx}.h5")
    plausible_orbit_aliases = [
        t
        for t in best
        if ("selection_params" in t or "ballistic_params" in t)
        and bool(t.get("interstellar_alias_plausible", False))
    ]
    selected_orbit_aliases = [t for t in best if "selection_params" in t or "ballistic_params" in t][:1]
    orbit_aliases = []
    seen_hypothesis_ids = set()
    for track in selected_orbit_aliases + plausible_orbit_aliases:
        hyp_id = int(track.get("hypothesis_id", -1))
        if hyp_id in seen_hypothesis_ids:
            continue
        orbit_aliases.append(track)
        seen_hypothesis_ids.add(hyp_id)
    if orbit_aliases and not args.skip_orbit_products:
        print_candidate_orbits(
            orbit_aliases,
            args.sample_idx / 1e6,
            n_candidates=len(orbit_aliases),
            state_h5_path=state_h5,
            n_samples=args.orbit_samples,
        )
        if args.run_dasst and state_h5.exists():
            cmd = [
                sys.executable,
                str(Path(__file__).with_name("dasst_orbits_from_candidate_states.py")),
                str(state_h5),
                "--output-h5",
                str(dasst_orbit_h5),
                "--metadata-dir",
                str(args.orbit_metadata_dir),
            ]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as exc:
                print(f"dasst_failed sample_idx={args.sample_idx} returncode={exc.returncode}")
        stage("orbit_products")
    if not args.overview_only:
        plot_visibility_rejections(candidates, tracks, visibility_out)
        plot_linearity_rejections(candidates, tracks, linearity_out)
        plot_descent_rejections(tracks, descent_out)
        plot_ballistic_ranking(tracks, candidates, ballistic_out)
        plot_ballistic_horizontal_ranking(tracks, candidates, ballistic_horizontal_out)
        plot_tx_beam_consistency(tracks, tx_beam_out)
        plot_tx_array_snr_consistency(tracks, tx_array_snr_out)
        plot_tx_array_snr_overlay(tracks, tx_array_snr_overlay_out, n_best=3)
        plot_tx_array_snr_overlay_by_beam(tracks, tx_array_snr_overlay_by_beam_out, n_best=3)
        separate_overlay_paths = plot_tx_array_snr_overlay_separate_by_beam(
            tracks,
            args.output_dir / "tx_array_snr_overlays",
            args.sample_idx,
            n_best=3,
        )
    else:
        separate_overlay_paths = []
    plot_disambiguation_summary(
        candidates,
        tracks,
        disambiguation_summary_out,
        orbit_h5_path=dasst_orbit_h5,
        obs=obs,
        cut=cut,
        obs_rti=obs_all,
        event_epoch_unix=args.sample_idx / 1e6,
        response_u=u,
        response_v=v,
        selected_response=response_selected,
    )
    stage("summary_plot")
    write_disambiguation_diagnostics_h5(
        candidates,
        tracks,
        diagnostics_h5,
        args.sample_idx,
        state_h5_path=state_h5,
        dasst_orbit_h5_path=dasst_orbit_h5,
        sigma_pos_km=sigma_pos,
        sigma_dop_km_s=sigma_dop,
        range_time_component_lengths=range_time_component_lengths,
        selected_range_time_component=selected_range_time_component,
        tx_phase_quality=tx_quality,
    )
    stage("diagnostics_h5")

    if not args.overview_only:
        print(antennas_out)
        print(tx_lobe_centroids_out)
        print(single_out)
        print(all_out)
        print(visibility_out)
        print(linearity_out)
        print(descent_out)
        print(ballistic_out)
    print(disambiguation_summary_out)
    if not args.overview_only:
        print(ballistic_horizontal_out)
        print(tx_beam_out)
        print(tx_array_snr_out)
        print(tx_array_snr_overlay_out)
        print(tx_array_snr_overlay_by_beam_out)
    for path in separate_overlay_paths:
            print(path)
    print(diagnostics_h5)
    print(f"single_echo_high_coherence_peaks {len(peak_ij[0])}")
    print(f"all_high_coherence_candidates {len(candidates)}")
    print(f"range_time_components {len(range_time_component_lengths)}")
    print("range_time_component_lengths " + " ".join(str(n) for n in range_time_component_lengths))
    print(f"selected_range_time_component {selected_range_time_component}")
    print("tx_lobe_centroid_counts " + " ".join(str(len(c)) for c in tx_lobe_centroids))
    print(f"path_hypotheses {len(tracks)}")
    print(f"path_hypotheses_kept {sum(t['reason'] == 'kept' for t in tracks)}")
    print(f"path_hypotheses_below_horizon {sum(t['reason'] == 'below horizon' for t in tracks)}")
    print(f"path_hypotheses_wrapping {sum(t['reason'] == 'wraps' for t in tracks)}")
    print(f"path_hypotheses_linearity_rejected {sum(t.get('linearity_reject', False) for t in tracks)}")
    print(f"path_hypotheses_after_linearity {sum(t['reason'] == 'kept' and not t.get('linearity_reject', False) for t in tracks)}")
    print(f"path_hypotheses_low_detection_altitude_rejected {sum(t.get('low_detection_altitude_reject', False) for t in tracks)}")
    print(f"path_hypotheses_after_detection_altitude {sum(t['reason'] == 'kept' and not t.get('linearity_reject', False) and not t.get('low_detection_altitude_reject', False) for t in tracks)}")
    print(f"path_hypotheses_descent_rejected {sum(t.get('descent_reject', False) for t in tracks)}")
    print(f"path_hypotheses_after_descent {sum(t['reason'] == 'kept' and not t.get('linearity_reject', False) and not t.get('low_detection_altitude_reject', False) and not t.get('descent_reject', False) for t in tracks)}")
    print(f"path_hypotheses_low_start_altitude_rejected {sum(t.get('ballistic_low_start_altitude_reject', False) for t in tracks)}")
    print(f"ballistic_sigma_pos_km {sigma_pos:.4f}")
    print(f"ballistic_sigma_dop_km_s {sigma_dop:.4f}")
    print("tx_pointing_score_used False")
    ranked_tracks = [t for t in tracks if "combined_rank" in t]
    accepted_tracks = [t for t in ranked_tracks if not bool(t.get("combined_reject", False))]
    print(f"combined_hypotheses_ranked {len(ranked_tracks)}")
    print(f"combined_hypotheses_rejected {sum(t.get('combined_reject', False) for t in ranked_tracks)}")
    print(f"combined_hypotheses_after {len(accepted_tracks)}")
    print(
        "hypothesis_status_columns "
        "hypothesis reason unique_pulses below_horizon wraps linearity_reject "
        "line_redchi line_rms_km line_max_km min_detect_alt_km low_detection_altitude_reject "
        "descent_rate_km_s descent_reject speedup_flag start_alt_km low_start_altitude_reject "
        "pulse_coverage_reject required_unique_pulses ballistic_rank combined_rank candidate "
        "selection_model selection_redchi fixed_velocity_redchi fixed_am_redchi "
        "short_static_measurement_reject range_span_km interstellar_plausible"
    )
    for track in sorted(tracks, key=lambda t: t.get("hypothesis_id", 999)):
        print(
            "hypothesis_status",
            hypothesis_label(track),
            track.get("reason", "unknown"),
            track.get("unique_pulses", 0),
            track.get("below_horizon", False),
            track.get("wraps", False),
            track.get("linearity_reject", False),
            f"{track.get('line_reduced_chi2', np.nan):.3f}" if "line_reduced_chi2" in track else "NA",
            f"{track.get('line_rms_km', np.nan):.3f}" if "line_rms_km" in track else "NA",
            f"{track.get('line_max_km', np.nan):.3f}" if "line_max_km" in track else "NA",
            f"{track.get('detection_min_altitude_km', np.nan):.3f}" if "detection_min_altitude_km" in track else "NA",
            track.get("low_detection_altitude_reject", False),
            f"{track.get('descent_rate_km_s', np.nan):.3f}" if "descent_rate_km_s" in track else "NA",
            track.get("descent_reject", False),
            track.get("ballistic_speedup_flag", False),
            f"{track.get('ballistic_start_altitude_km', np.nan):.3f}" if "ballistic_start_altitude_km" in track else "NA",
            track.get("ballistic_low_start_altitude_reject", False),
            track.get("ballistic_pulse_coverage_reject", False),
            track.get("ballistic_required_unique_pulses", "NA"),
            track.get("ballistic_rank", "NA"),
            track.get("combined_rank", "NA"),
            candidate_number(track) if "combined_rank" in track else "NA",
            track.get("selection_model_type", "NA"),
            f"{track.get('selection_reduced_chi2', np.nan):.3f}" if "selection_reduced_chi2" in track else "NA",
            f"{track.get('fixed_velocity_reduced_chi2', np.nan):.3f}" if "fixed_velocity_reduced_chi2" in track else "NA",
            f"{track.get('fixed_am_reduced_chi2', np.nan):.3f}" if "fixed_am_reduced_chi2" in track else "NA",
            track.get("short_static_measurement_reject", False),
            f"{track.get('measurement_range_span_km', np.nan):.3f}",
            track.get("interstellar_alias_plausible", False),
        )
    for track in sorted(
        [
            t
            for t in tracks
            if t["reason"] == "kept" and not t.get("linearity_reject", False) and not t.get("descent_reject", False)
            and "ballistic_keep" in t
        ],
        key=lambda t: t.get("ballistic_rank", 999),
    ):
        fit_keep = track["ballistic_keep"]
        chi_pos = float(np.sum((track["ballistic_pos_res_km"][fit_keep] / sigma_pos) ** 2))
        chi_dop = float(np.sum((track["ballistic_dop_res_km_s"][fit_keep] / sigma_dop) ** 2))
        print(
            "ballistic_rank",
            track["ballistic_rank"],
            "hypothesis",
            hypothesis_label(track),
            "candidate",
            candidate_number(track) if "combined_rank" in track else "unranked",
            "redchi",
            f"{track['ballistic_reduced_chi2']:.3f}",
            "chi_pos",
            f"{chi_pos:.1f}",
            "chi_dop",
            f"{chi_dop:.1f}",
            "pos_rms_km",
            f"{track['ballistic_pos_rms_km']:.3f}",
            "dop_rms_km_s",
            f"{track['ballistic_dop_rms_km_s']:.3f}",
            "start_alt_km",
            f"{track.get('ballistic_start_altitude_km', np.nan):.3f}",
            "low_start_alt_reject",
            track.get("ballistic_low_start_altitude_reject", False),
            "speed_slope_km_s2",
            f"{track.get('ballistic_speed_slope_km_s2', np.nan):.3f}",
            "speed_delta_km_s",
            f"{track.get('ballistic_speed_delta_km_s', np.nan):.3f}",
            "speedup_flag",
            track.get("ballistic_speedup_flag", False),
            "n_out",
            track["ballistic_n_outliers"],
            "pulse_coverage_reject",
            track.get("ballistic_pulse_coverage_reject", False),
            "required_unique_pulses",
            track.get("ballistic_required_unique_pulses", "NA"),
            "ballistic_reject",
            track["ballistic_reject"],
        )
    for track in sorted(
        [t for t in tracks if "combined_score" in t and "ballistic_reduced_chi2" in t],
        key=lambda t: t.get("combined_rank", 999),
    ):
        print(
            "combined_rank",
            track["combined_rank"],
            "hypothesis",
            hypothesis_label(track),
            "candidate",
            candidate_number(track),
            "ballistic_rank",
            track.get("ballistic_rank", "NA"),
            "tx_beam_rank",
            track.get("tx_beam_rank", "NA"),
            "combined_score",
            f"{track['combined_score']:.3f}",
            "good_fit",
            track.get("combined_good_fit", False),
            "good_fit_redchi_threshold",
            f"{track.get('combined_good_fit_redchi_threshold', np.nan):.3f}",
            "selection_model",
            track.get("selection_model_type", "NA"),
            "selection_term",
            f"{track.get('selection_reduced_chi2', np.nan):.3f}" if "selection_reduced_chi2" in track else "NA",
            "ballistic_term",
            f"{track['ballistic_reduced_chi2']:.3f}" if "ballistic_reduced_chi2" in track else "NA",
            "fixed_velocity_term",
            f"{track.get('fixed_velocity_reduced_chi2', np.nan):.3f}" if "fixed_velocity_reduced_chi2" in track else "NA",
            "fixed_am_term",
            f"{track.get('fixed_am_reduced_chi2', np.nan):.3f}" if "fixed_am_reduced_chi2" in track else "NA",
            "interstellar_plausible",
            track.get("interstellar_alias_plausible", False),
            "coh_mean",
            f"{track.get('coherence_weighted_mean', np.nan):.3f}",
            "tx_af_snr_rms_db",
            f"{track.get('tx_array_snr_rms_db', np.nan):.2f}",
            "tx_af_snr_corr",
            f"{track.get('tx_array_snr_corr', np.nan):.3f}",
            "tx_af_gain_span_db",
            f"{track.get('tx_array_gain_span_db', np.nan):.2f}",
            "tx_rms_deg_diagnostic",
            f"{track.get('tx_beam_weighted_rms_deg', np.nan):.2f}",
            "tx_beam_mean_dc",
            f"{track.get('tx_beam_snr_weighted_mean_dc', np.nan):.4f}",
            "tx_beam_rms_dc",
            f"{track.get('tx_beam_snr_weighted_rms_dc', np.nan):.4f}",
            "tx_lobe_mean_dc",
            f"{track.get('tx_lobe_snr_weighted_mean_dc', np.nan):.4f}",
            "tx_lobe_rms_dc",
            f"{track.get('tx_lobe_snr_weighted_rms_dc', np.nan):.4f}",
            "speedup_flag",
            track.get("ballistic_speedup_flag", False),
            "low_start_alt_reject",
            track.get("ballistic_low_start_altitude_reject", False),
            "pulse_coverage_reject",
            track.get("ballistic_pulse_coverage_reject", False),
            "short_static_measurement_reject",
            track.get("short_static_measurement_reject", False),
            "range_span_km",
            f"{track.get('measurement_range_span_km', np.nan):.3f}",
            "reject",
            track["combined_reject"],
        )
        if "tx_array_snr_corr_by_beam" in track:
            corr = track["tx_array_snr_corr_by_beam"]
            rms = track["tx_array_snr_rms_db_by_beam"]
            n_by_beam = track["tx_array_snr_n_by_beam"]
            print(
                "tx_af_snr_by_beam",
                "combined_rank",
                track["combined_rank"],
                "corr",
                ",".join(f"{x:.3f}" if np.isfinite(x) else "nan" for x in corr),
                "rms_db",
                ",".join(f"{x:.2f}" if np.isfinite(x) else "nan" for x in rms),
                "n",
                ",".join(str(int(x)) for x in n_by_beam),
            )
        if "tx_lobe_snr_weighted_mean_dc_by_beam" in track:
            lobe = track["tx_lobe_snr_weighted_mean_dc_by_beam"]
            n_by_beam = track["tx_lobe_n_by_beam"]
            print(
                "tx_lobe_distance_by_beam",
                "combined_rank",
                track["combined_rank"],
                "mean_dc",
                ",".join(f"{x:.4f}" if np.isfinite(x) else "nan" for x in lobe),
                "n",
                ",".join(str(int(x)) for x in n_by_beam),
            )
    if len(best) > 1:
        print(f"combined_delta_best_second {best[0]['combined_delta_to_next']:.3f}")
        print(f"combined_odds_best_vs_second {best[0]['combined_odds_best_vs_second']:.3f}")


if __name__ == "__main__":
    main()
