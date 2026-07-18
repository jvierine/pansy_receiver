#!/usr/bin/env python3
"""Compare shrinking-radius profiles with and without decoded beat-phase acceleration."""

from __future__ import annotations

import argparse
import datetime as dt
import itertools
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from scipy.sparse.csgraph import minimum_spanning_tree

import fit_best_alias_physics_models as physics
import pansy_config as pc
import pansy_interferometry as interferometry
from run_catalogue_mass_profiles import (
    load_selected_fit,
    log_grid_probability,
    profile_interval,
    radius_um_to_mass_kg,
    weighted_quantile,
)


def circular_residual(observed: np.ndarray, predicted: np.ndarray) -> np.ndarray:
    return np.angle(np.exp(1j * (observed - predicted)))


def load_nonoverlapping_phase_acceleration(path: Path) -> dict:
    with h5py.File(path, "r") as handle:
        group = handle["baud_averaged_beat_pairs"]
        pair = {name: np.asarray(values) for name, values in group.items()}
        measurement_tx_s = np.asarray(handle["measurement/tx_idx"], dtype=float) / 1e6
    previous = pair["previous_index"].astype(int)
    current = pair["current_index"].astype(int)
    beam = pair["beam_id"].astype(int)
    beat = pair["beat_real"] + 1j * pair["beat_imag"]
    rows = []
    for beam_id in np.unique(beam):
        indices = np.flatnonzero(beam == beam_id)
        indices = indices[np.argsort(pair["time_s"][indices])]
        candidates = []
        for first, second in zip(indices[:-1], indices[1:]):
            if previous[second] != current[first]:
                continue
            if abs(pair["delta_t_s"][first] - 0.008) > 5e-6:
                continue
            if abs(pair["delta_t_s"][second] - 0.008) > 5e-6:
                continue
            if abs(pair["delta_t_s"][second] - pair["delta_t_s"][first]) > 5e-6:
                continue
            observed_delta_phase = float(np.angle(beat[second] * np.conj(beat[first])))
            phase_std = float(np.hypot(pair["beat_phase_std_rad"][first], pair["beat_phase_std_rad"][second]))
            rows.append(
                (
                    beam_id,
                    first,
                    second,
                    pair["time_s"][first],
                    pair["time_s"][second],
                    pair["delta_t_s"][first],
                    pair["delta_t_s"][second],
                    observed_delta_phase,
                    phase_std,
                )
            )
            candidates.append(len(rows) - 1)
        # Each acceleration estimate uses three pulses. Retain every third
        # candidate so no raw pulse appears in more than one likelihood term.
        keep_candidates = set(candidates[::3])
        reject = set(candidates) - keep_candidates
        rows = [row for index, row in enumerate(rows) if index not in reject]
    dtype = [
        ("beam_id", "i4"),
        ("first_pair", "i4"),
        ("second_pair", "i4"),
        ("first_time_s", "f8"),
        ("second_time_s", "f8"),
        ("first_dt_s", "f8"),
        ("second_dt_s", "f8"),
        ("observed_delta_phase_rad", "f8"),
        ("formal_phase_std_rad", "f8"),
    ]
    return {"samples": np.asarray(rows, dtype=dtype), "measurement_tx_s": measurement_tx_s}


def predicted_delta_phase(predicted_doppler_mps: np.ndarray, t_s: np.ndarray, phase_data: dict) -> np.ndarray:
    samples = phase_data["samples"]
    absolute_zero = phase_data["measurement_tx_s"][0] - t_s[0]
    first_t = samples["first_time_s"] - absolute_zero
    second_t = samples["second_time_s"] - absolute_zero
    first_velocity = np.interp(first_t, t_s, predicted_doppler_mps)
    second_velocity = np.interp(second_t, t_s, predicted_doppler_mps)
    first_phase = 4.0 * np.pi * first_velocity * samples["first_dt_s"] / pc.wavelength
    second_phase = 4.0 * np.pi * second_velocity * samples["second_dt_s"] / pc.wavelength
    return np.angle(np.exp(1j * (second_phase - first_phase)))


def predicted_radial_acceleration(predicted_doppler_mps: np.ndarray, t_s: np.ndarray, phase_data: dict) -> np.ndarray:
    samples = phase_data["samples"]
    absolute_zero = phase_data["measurement_tx_s"][0] - t_s[0]
    first_t = samples["first_time_s"] - absolute_zero
    second_t = samples["second_time_s"] - absolute_zero
    first_velocity = np.interp(first_t, t_s, predicted_doppler_mps)
    second_velocity = np.interp(second_t, t_s, predicted_doppler_mps)
    return (second_velocity - first_velocity) / (second_t - first_t)


def load_nonoverlapping_cross_phase(path: Path) -> dict:
    with h5py.File(path, "r") as handle:
        measurement = handle["measurement"]
        xc_name = "xc_calibrated" if "xc_calibrated" in measurement else "xc"
        xc = np.asarray(measurement[xc_name], dtype=np.complex128)
        tx_s = np.asarray(measurement["tx_idx"], dtype=float) / 1e6
        beam_id = np.asarray(measurement["beam_id"], dtype=int)
        snr = np.asarray(measurement["snr"], dtype=float)
        range_km = np.asarray(measurement["range_km"], dtype=float)
    xc /= np.maximum(np.abs(xc), 1e-300)
    previous = []
    current = []
    for beam in np.unique(beam_id):
        indices = np.flatnonzero(beam_id == beam)
        candidates = []
        for prev, cur in zip(indices[:-1], indices[1:]):
            delta_t = tx_s[cur] - tx_s[prev]
            if np.isfinite(delta_t) and abs(delta_t - 0.008) <= 5e-6:
                candidates.append((prev, cur))
        for prev, cur in candidates[::2]:
            previous.append(prev)
            current.append(cur)
    previous = np.asarray(previous, dtype=int)
    current = np.asarray(current, dtype=int)
    observed = np.angle(xc[current] * np.conj(xc[previous]))
    pair_snr = np.minimum(snr[previous], snr[current])
    weight = np.clip(pair_snr, 0.0, 100.0)
    good_weight = np.isfinite(weight) & (weight > 0.0)
    if not np.any(good_weight):
        weight = np.ones_like(pair_snr)
    else:
        weight = np.where(good_weight, weight, np.nanmin(weight[good_weight]))
    weight /= np.mean(weight)
    n_channels = int((1.0 + np.sqrt(1.0 + 8.0 * xc.shape[1])) / 2.0)
    channel_pairs = np.asarray(list(itertools.combinations(range(n_channels), 2)), dtype=int)
    if len(channel_pairs) != xc.shape[1]:
        raise RuntimeError(f"cannot infer receiver channels from {xc.shape[1]} cross-products")
    antenna_position_m = np.asarray(interferometry.get_antpos(), dtype=float)
    baseline_m = antenna_position_m[channel_pairs[:, 1]] - antenna_position_m[channel_pairs[:, 0]]
    return {
        "previous": previous,
        "current": current,
        "time_s": 0.5 * (tx_s[previous] + tx_s[current]),
        "delta_t_s": tx_s[current] - tx_s[previous],
        "range_km": 0.5 * (range_km[previous] + range_km[current]),
        "observed_phase_rad": observed,
        "snr_linear": pair_snr,
        "weight": weight,
        "channel_pairs": channel_pairs,
        "baseline_m": baseline_m,
        "beam_id": beam_id[current],
    }


def predicted_cross_phase(position_km: np.ndarray, cross_data: dict) -> np.ndarray:
    direction = position_km / np.linalg.norm(position_km, axis=1)[:, None]
    delta_direction = direction[cross_data["current"]] - direction[cross_data["previous"]]
    return (2.0 * np.pi / pc.wavelength) * (delta_direction @ cross_data["baseline_m"].T)


def independent_baselines(channel_pairs: np.ndarray, sigma_rad: np.ndarray) -> np.ndarray:
    n_channels = int(np.max(channel_pairs)) + 1
    graph = np.zeros((n_channels, n_channels), dtype=float)
    for (left, right), sigma in zip(channel_pairs, sigma_rad):
        graph[left, right] = graph[right, left] = max(float(sigma), 1e-6)
    tree = minimum_spanning_tree(graph).toarray()
    selected = []
    for index, (left, right) in enumerate(channel_pairs):
        if tree[left, right] != 0.0 or tree[right, left] != 0.0:
            selected.append(index)
    return np.asarray(selected, dtype=int)


def fit_profile(diagnostics_h5: Path, baseline_h5: Path, beat_h5: Path, output_h5: Path, output_png: Path) -> dict:
    observations, stored = load_selected_fit(diagnostics_h5)
    t_s = observations["t_s"]
    points = observations["points_km"]
    doppler = observations["doppler_km_s"]
    finite = np.isfinite(t_s) & np.all(np.isfinite(points), axis=1) & np.isfinite(doppler)
    keep = stored["keep"] & finite
    if np.count_nonzero(keep) < 10:
        keep = finite
    stored_prediction = physics.predicted_doppler(stored["model_km"], stored["velocity_km_s"])
    sigma_position = float(np.sqrt(np.mean((points[keep] - stored["model_km"][keep]) ** 2)))
    sigma_doppler = float(np.sqrt(np.mean((doppler[keep] - stored_prediction[keep]) ** 2)))
    phase_data = load_nonoverlapping_phase_acceleration(beat_h5)
    stored_phase_prediction = predicted_delta_phase(stored_prediction * 1e3, t_s, phase_data)
    phase_residual = circular_residual(phase_data["samples"]["observed_delta_phase_rad"], stored_phase_prediction)
    phase_samples = phase_data["samples"]
    absolute_zero = phase_data["measurement_tx_s"][0] - t_s[0]
    phase_midpoint_t = 0.5 * (phase_samples["first_time_s"] + phase_samples["second_time_s"]) - absolute_zero
    phase_snr = np.interp(phase_midpoint_t, t_s, observations["snr"])
    phase_weight = np.clip(phase_snr, 0.0, 100.0)
    finite_positive_weight = np.isfinite(phase_weight) & (phase_weight > 0.0)
    if not np.any(finite_positive_weight):
        phase_weight = np.ones_like(phase_residual)
    else:
        floor = float(np.nanmin(phase_weight[finite_positive_weight]))
        phase_weight = np.where(finite_positive_weight, phase_weight, floor)
    phase_weight /= np.mean(phase_weight)
    sigma_phase = float(np.sqrt(np.sum(phase_weight * phase_residual**2) / np.sum(phase_weight)))
    if not np.isfinite(sigma_phase) or sigma_phase <= 0.0:
        raise RuntimeError("invalid beat-phase residual RMS")
    mean_dt = 0.5 * (phase_samples["first_dt_s"] + phase_samples["second_dt_s"])
    measured_acceleration_principal = (
        pc.wavelength
        * phase_samples["observed_delta_phase_rad"]
        / (4.0 * np.pi * mean_dt**2)
    )
    stored_acceleration = predicted_radial_acceleration(stored_prediction * 1e3, t_s, phase_data)
    acceleration_ambiguity = pc.wavelength / (2.0 * mean_dt**2)
    acceleration_residual = measured_acceleration_principal - stored_acceleration
    acceleration_residual -= np.rint(acceleration_residual / acceleration_ambiguity) * acceleration_ambiguity
    sigma_acceleration = float(np.sqrt(np.mean(acceleration_residual**2)))
    if not np.isfinite(sigma_acceleration) or sigma_acceleration <= 0.0:
        raise RuntimeError("invalid radial-acceleration residual RMS")
    cross_data = load_nonoverlapping_cross_phase(beat_h5)
    stored_cross_prediction = predicted_cross_phase(stored["model_km"], cross_data)
    stored_cross_residual = circular_residual(
        cross_data["observed_phase_rad"], stored_cross_prediction
    )
    cross_sigma = np.sqrt(
        np.sum(cross_data["weight"][:, None] * stored_cross_residual**2, axis=0)
        / np.sum(cross_data["weight"])
    )
    cross_sigma = np.maximum(cross_sigma, 0.10)
    selected_baselines = independent_baselines(cross_data["channel_pairs"], cross_sigma)
    if len(selected_baselines) != len(np.asarray(interferometry.get_antpos())) - 1:
        raise RuntimeError("cross-phase baseline selection did not produce a spanning tree")
    density, _metadata = physics.pbal.density_interpolator(observations["sample_epoch_unix"])

    with h5py.File(baseline_h5, "r") as baseline:
        radius_um = np.asarray(baseline["profile/radius_um"], dtype=float)
        baseline_delta = np.asarray(baseline["profile/delta_chi2"], dtype=float)
        baseline_params6 = np.asarray(baseline["profile/parameters6"], dtype=float)
        baseline_best = float(baseline["result/free_best_radius_um"][()])
        baseline_interval = (
            float(baseline["result/profile_ci95_lower_radius_um"][()]),
            float(baseline["result/profile_ci95_upper_radius_um"][()]),
        )
    log_radius_m = np.log10(radius_um * 1e-6)
    lower6 = np.asarray([-np.inf, -np.inf, 20e3, -90e3, -90e3, -90e3])
    upper6 = np.asarray([np.inf, np.inf, 220e3, 90e3, 90e3, 90e3])
    scale6 = np.asarray([1e5, 1e5, 1e5, 7e4, 7e4, 7e4])

    def residual6(parameters6, fixed_log_radius):
        position, velocity, _radius, _mass, _success, _message = physics.propagate_shrinking_radius_model(
            np.r_[parameters6, fixed_log_radius], t_s, density
        )
        prediction = physics.predicted_doppler(position, velocity)
        phase_prediction = predicted_delta_phase(prediction * 1e3, t_s, phase_data)
        phase_fit_residual = circular_residual(
            phase_samples["observed_delta_phase_rad"], phase_prediction
        )
        cross_prediction = predicted_cross_phase(position, cross_data)
        cross_fit_residual = circular_residual(
            cross_data["observed_phase_rad"], cross_prediction
        )[:, selected_baselines]
        return np.r_[
            ((points[keep] - position[keep]) / sigma_position).ravel(),
            (doppler[keep] - prediction[keep]) / sigma_doppler,
            np.sqrt(phase_weight) * phase_fit_residual / sigma_phase,
            (
                np.sqrt(cross_data["weight"][:, None])
                * cross_fit_residual
                / cross_sigma[selected_baselines][None, :]
            ).ravel(),
        ]

    chi2 = np.full(len(radius_um), np.nan)
    parameters6 = np.full((len(radius_um), 6), np.nan)
    success = np.zeros(len(radius_um), dtype=bool)
    for index in np.argsort(np.abs(np.log(radius_um) - np.log(baseline_best))):
        starts = [baseline_params6[index]]
        if index > 0 and np.all(np.isfinite(parameters6[index - 1])):
            starts.append(parameters6[index - 1])
        candidates = []
        for start in starts:
            if not np.all(np.isfinite(start)):
                continue
            try:
                result = opt.least_squares(
                    lambda values: residual6(values, log_radius_m[index]),
                    start,
                    bounds=(lower6, upper6),
                    x_scale=scale6,
                    loss="linear",
                    max_nfev=60,
                )
                candidates.append(result)
            except Exception:
                pass
        if candidates:
            best = min(candidates, key=lambda candidate: candidate.cost)
            parameters6[index] = best.x
            chi2[index] = float(np.sum(residual6(best.x, log_radius_m[index]) ** 2))
            success[index] = bool(best.success)
    minimum = float(np.nanmin(chi2))
    delta = chi2 - minimum
    probability, weights = log_grid_probability(radius_um, delta)
    quantiles = weighted_quantile(radius_um, weights, [0.025, 0.5, 0.975])
    threshold = 3.841458820694124
    lower, upper, lower_bounded, upper_bounded, lower_status, upper_status = profile_interval(radius_um, delta, threshold)
    best_radius = float(radius_um[np.nanargmin(chi2)])
    best_index = int(np.nanargmin(chi2))
    best_position, best_velocity, *_ = physics.propagate_shrinking_radius_model(
        np.r_[parameters6[best_index], log_radius_m[best_index]], t_s, density
    )
    best_prediction = physics.predicted_doppler(best_position, best_velocity)
    best_phase_prediction = predicted_delta_phase(best_prediction * 1e3, t_s, phase_data)
    best_acceleration = predicted_radial_acceleration(best_prediction * 1e3, t_s, phase_data)
    best_cross_prediction = predicted_cross_phase(best_position, cross_data)
    display_cross_phase = cross_data["observed_phase_rad"] + 2.0 * np.pi * np.rint(
        (best_cross_prediction - cross_data["observed_phase_rad"]) / (2.0 * np.pi)
    )
    selected_baseline_m = cross_data["baseline_m"][selected_baselines, :2]
    cross_design = (2.0 * np.pi / pc.wavelength) * selected_baseline_m
    cross_design_weighted = cross_design / cross_sigma[selected_baselines, None]
    measured_horizontal_velocity_km_s = np.full((len(display_cross_phase), 2), np.nan)
    for row in range(len(display_cross_phase)):
        delta_direction, *_ = np.linalg.lstsq(
            cross_design_weighted,
            display_cross_phase[row, selected_baselines] / cross_sigma[selected_baselines],
            rcond=None,
        )
        measured_horizontal_velocity_km_s[row] = (
            cross_data["range_km"][row]
            * delta_direction
            / cross_data["delta_t_s"][row]
        )
    best_direction = best_position / np.linalg.norm(best_position, axis=1)[:, None]
    best_horizontal_velocity_km_s = (
        cross_data["range_km"][:, None]
        * (best_direction[cross_data["current"], :2] - best_direction[cross_data["previous"], :2])
        / cross_data["delta_t_s"][:, None]
    )
    display_acceleration = measured_acceleration_principal + np.rint(
        (best_acceleration - measured_acceleration_principal) / acceleration_ambiguity
    ) * acceleration_ambiguity

    output_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_h5, "w") as handle:
        handle.attrs["sample_idx"] = int(observations["sample_idx"])
        handle.attrs["sigma_phase_rad"] = sigma_phase
        handle.attrs["phase_weighting"] = "linear SNR interpolated to each phase sample, capped at 100 and normalized to unit mean"
        handle.attrs["cross_phase_likelihood"] = "same-transmit-beam 8-ms non-overlapping echo pairs; six independent spanning-tree receiver baselines; modulo-2pi residual with fixed per-baseline RMS and SNR weights"
        handle.attrs["sigma_radial_acceleration_mps2"] = sigma_acceleration
        handle.attrs["n_nonoverlapping_phase_acceleration"] = len(phase_data["samples"])
        handle.attrs["phase_likelihood"] = "beat-phase residual wrapped to [-pi, pi) independently for each model"
        profile = handle.create_group("profile")
        profile["radius_um"] = radius_um
        profile["mass_kg"] = radius_um_to_mass_kg(radius_um)
        profile["chi2"] = chi2
        profile["delta_chi2"] = delta
        profile["probability_density_log_radius"] = probability
        profile["parameters6"] = parameters6
        profile["success"] = success
        result = handle.create_group("result")
        result["best_radius_um"] = best_radius
        result["ci95_lower_radius_um"] = lower
        result["ci95_upper_radius_um"] = upper
        result["ci95_lower_bounded"] = lower_bounded
        result["ci95_upper_bounded"] = upper_bounded
        result.attrs["ci95_lower_status"] = lower_status
        result.attrs["ci95_upper_status"] = upper_status
        result["marginal_radius_quantiles_um"] = quantiles
        phase = handle.create_group("phase_acceleration")
        phase["samples"] = phase_data["samples"]
        phase["snr_linear"] = phase_snr
        phase["normalized_snr_weight"] = phase_weight
        phase["stored_model_prediction_rad"] = stored_phase_prediction
        phase["stored_model_residual_rad"] = phase_residual
        phase["best_model_prediction_rad"] = best_phase_prediction
        phase["measured_radial_acceleration_principal_mps2"] = measured_acceleration_principal
        phase["measured_radial_acceleration_display_mps2"] = display_acceleration
        phase["acceleration_ambiguity_period_mps2"] = acceleration_ambiguity
        phase["best_model_radial_acceleration_mps2"] = best_acceleration
        phase["stored_model_radial_acceleration_mps2"] = stored_acceleration
        phase["stored_model_radial_acceleration_residual_mps2"] = acceleration_residual
        cross = handle.create_group("cross_phase_velocity")
        for name in ("previous", "current", "time_s", "delta_t_s", "range_km", "snr_linear", "weight", "channel_pairs", "baseline_m", "beam_id"):
            cross[name] = cross_data[name]
        cross["observed_phase_rad"] = cross_data["observed_phase_rad"]
        cross["stored_model_prediction_rad"] = stored_cross_prediction
        cross["stored_model_residual_rad"] = stored_cross_residual
        cross["phase_sigma_rad"] = cross_sigma
        cross["selected_baseline_indices"] = selected_baselines
        cross["best_model_prediction_rad"] = best_cross_prediction
        cross["display_phase_rad"] = display_cross_phase
        cross["measured_horizontal_velocity_km_s"] = measured_horizontal_velocity_km_s
        cross["best_model_horizontal_velocity_km_s"] = best_horizontal_velocity_km_s
        cross.attrs["pairing"] = "same transmit beam, 8 ms, non-overlapping"

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3), constrained_layout=True)
    ax = axes[0]
    ax.plot(radius_um, baseline_delta, color="0.5", lw=1.2, label="position + Doppler")
    ax.plot(radius_um, delta, color="C0", lw=1.4, label="+ radial deceleration")
    ax.axhline(threshold, color="0.2", ls="--", lw=0.8, label="95% profile threshold")
    ax.set_xscale("log")
    ax.set_ylim(0, min(30.0, max(8.0, np.nanpercentile(delta, 90))))
    ax.set_xlabel(r"Initial radius $r_0$ ($\mu$m)")
    ax.set_ylabel(r"$\Delta\chi^2$")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2, which="both")
    ax = axes[1]
    accepted = np.isfinite(delta)
    ax.plot(radius_um[accepted], probability[accepted] / np.nanmax(probability), color="C0")
    ax.fill_between(radius_um[accepted], 0, probability[accepted] / np.nanmax(probability), color="C0", alpha=0.2)
    ax.set_xscale("log")
    ax.set_xlabel(r"Initial radius $r_0$ ($\mu$m)")
    ax.set_ylabel("Relative probability")
    ax.grid(alpha=0.2, which="both")
    baseline_text = f"old 95%: {baseline_interval[0]:.0f}--{baseline_interval[1]:.0f} $\\mu$m"
    new_text = f"new 95%: {lower:.0f}--{upper:.0f} $\\mu$m"
    ax.text(
        0.03,
        0.96,
        baseline_text
        + "\n"
        + new_text
        + f"\nN accel={len(phase_data['samples'])}, RMS={sigma_acceleration / 1e3:.2f} km s$^{{-2}}$",
        transform=ax.transAxes,
        va="top",
    )
    fig.savefig(output_png, dpi=190)
    plt.close(fig)
    return {
        "sample_idx": int(observations["sample_idx"]),
        "baseline_best_radius_um": baseline_best,
        "baseline_interval_um": baseline_interval,
        "new_best_radius_um": best_radius,
        "new_interval_um": (lower, upper),
        "new_quantiles_um": quantiles,
        "sigma_phase_rad": sigma_phase,
        "sigma_acceleration_mps2": sigma_acceleration,
        "n_phase": len(phase_data["samples"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics", type=Path)
    parser.add_argument("--baseline-profile", type=Path)
    parser.add_argument("--beat-h5", type=Path)
    parser.add_argument("--output-h5", type=Path)
    parser.add_argument("--output-png", type=Path)
    parser.add_argument("--sample-idx", type=int)
    parser.add_argument("--base", type=Path, default=Path("/mnt/data/juha/pansy"))
    parser.add_argument("--baseline-profile-dir", type=Path)
    parser.add_argument("--beat-h5-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if args.sample_idx is not None:
        if args.baseline_profile_dir is None or args.beat_h5_dir is None or args.output_dir is None:
            parser.error("--sample-idx requires --baseline-profile-dir, --beat-h5-dir, and --output-dir")
        day = dt.datetime.fromtimestamp(args.sample_idx / 1e6, tz=dt.timezone.utc).strftime("%Y-%m-%d")
        diagnostics = args.base / "events" / day / f"pansy_disambiguation_diagnostics_{args.sample_idx}.h5"
        baseline_profile = args.baseline_profile_dir / f"mass_profile_{args.sample_idx}.h5"
        beat_h5 = args.beat_h5_dir / f"inter_pulse_phase_{args.sample_idx}.h5"
        output_h5 = args.output_dir / f"mass_profile_with_acceleration_{args.sample_idx}.h5"
        output_png = args.output_dir / f"mass_profile_with_acceleration_{args.sample_idx}.png"
    else:
        required = (args.diagnostics, args.baseline_profile, args.beat_h5, args.output_h5, args.output_png)
        if any(value is None for value in required):
            parser.error("provide explicit input/output paths or use --sample-idx batch mode")
        diagnostics, baseline_profile, beat_h5, output_h5, output_png = required
    summary = fit_profile(diagnostics, baseline_profile, beat_h5, output_h5, output_png)
    print(summary, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
