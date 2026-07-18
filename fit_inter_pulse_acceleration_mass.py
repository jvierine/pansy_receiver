#!/usr/bin/env python3
"""Compare shrinking-radius profiles with and without decoded beat-phase acceleration."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

import fit_best_alias_physics_models as physics
import pansy_config as pc
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
    sigma_phase = float(np.sqrt(np.mean(phase_residual**2)))
    if not np.isfinite(sigma_phase) or sigma_phase <= 0.0:
        raise RuntimeError("invalid beat-phase residual RMS")
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
        return np.r_[
            ((points[keep] - position[keep]) / sigma_position).ravel(),
            (doppler[keep] - prediction[keep]) / sigma_doppler,
            circular_residual(phase_data["samples"]["observed_delta_phase_rad"], phase_prediction) / sigma_phase,
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

    output_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_h5, "w") as handle:
        handle.attrs["sample_idx"] = int(observations["sample_idx"])
        handle.attrs["sigma_phase_rad"] = sigma_phase
        handle.attrs["n_nonoverlapping_phase_acceleration"] = len(phase_data["samples"])
        handle.attrs["phase_likelihood"] = "wrapped circular residual of consecutive SNR-weighted decoded beat phasors"
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
        phase["stored_model_prediction_rad"] = stored_phase_prediction
        phase["stored_model_residual_rad"] = phase_residual

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3), constrained_layout=True)
    ax = axes[0]
    ax.plot(radius_um, baseline_delta, color="0.5", lw=1.2, label="position + Doppler")
    ax.plot(radius_um, delta, color="C0", lw=1.4, label="+ decoded beat phase")
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
    ax.text(0.03, 0.96, baseline_text + "\n" + new_text + f"\nN phase={len(phase_data['samples'])}, RMS={sigma_phase:.2f} rad", transform=ax.transAxes, va="top")
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
        "n_phase": len(phase_data["samples"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--baseline-profile", type=Path, required=True)
    parser.add_argument("--beat-h5", type=Path, required=True)
    parser.add_argument("--output-h5", type=Path, required=True)
    parser.add_argument("--output-png", type=Path, required=True)
    args = parser.parse_args()
    summary = fit_profile(args.diagnostics, args.baseline_profile, args.beat_h5, args.output_h5, args.output_png)
    print(summary, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
