#!/usr/bin/env python3
"""Decompose a phase-aware mass profile into measurement contributions."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

import fit_best_alias_physics_models as physics
import pansy_config as pc
from run_catalogue_mass_profiles import load_selected_fit


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--baseline-profile", type=Path, required=True)
    parser.add_argument("--phase-profile", type=Path, required=True)
    parser.add_argument("--beat-h5", type=Path, required=True)
    parser.add_argument("--output-h5", type=Path, required=True)
    parser.add_argument("--output-plot", type=Path, required=True)
    return parser.parse_args()


def circular_residual(observed, predicted):
    return np.angle(np.exp(1j * (observed - predicted)))


def predicted_phase(predicted_doppler_mps, t_s, absolute_time_zero, samples):
    first_t = samples["first_time_s"] - absolute_time_zero
    second_t = samples["second_time_s"] - absolute_time_zero
    first_velocity = np.interp(first_t, t_s, predicted_doppler_mps)
    second_velocity = np.interp(second_t, t_s, predicted_doppler_mps)
    first_phase = 4.0 * np.pi * first_velocity * samples["first_dt_s"] / pc.wavelength
    second_phase = 4.0 * np.pi * second_velocity * samples["second_dt_s"] / pc.wavelength
    return np.angle(np.exp(1j * (second_phase - first_phase)))


def phase_chi2(group, prediction):
    selected = np.asarray(group["shared_inlier_mask"], dtype=bool)
    observed = np.asarray(group["samples"])["observed_delta_phase_rad"]
    residual = circular_residual(observed[selected], prediction[selected])
    covariance = np.asarray(group["phase_difference_covariance_rad2"], dtype=float)
    covariance = covariance[np.ix_(selected, selected)]
    scale = float(group.attrs["covariance_scale"])
    return float(residual @ np.linalg.solve(covariance, residual) / scale**2)


def residual_correlation_summary(residual):
    residual = np.asarray(residual, dtype=float)
    correlations = []
    for lag in range(1, min(31, len(residual) // 3)):
        correlations.append(float(np.corrcoef(residual[:-lag], residual[lag:])[0, 1]))
    positive = []
    for value in correlations:
        if not np.isfinite(value) or value <= 0.0:
            break
        positive.append(value)
    correlation_time = 1.0 + 2.0 * float(np.sum(positive))
    effective_samples = len(residual) / correlation_time
    return np.asarray(correlations), correlation_time, effective_samples


def main():
    args = parse_args()
    observations, stored = load_selected_fit(args.diagnostics)
    t_s = np.asarray(observations["t_s"], dtype=float)
    points = np.asarray(observations["points_km"], dtype=float)
    doppler = np.asarray(observations["doppler_km_s"], dtype=float)
    sample_epoch = float(observations["sample_epoch_unix"])
    finite = np.isfinite(t_s) & np.all(np.isfinite(points), axis=1) & np.isfinite(doppler)
    stored_keep = np.asarray(stored["keep"], dtype=bool) & finite
    if np.count_nonzero(stored_keep) < 10:
        stored_keep = finite
    stored_prediction = physics.predicted_doppler(stored["model_km"], stored["velocity_km_s"])
    sigma_position = float(np.sqrt(np.mean((points[stored_keep] - stored["model_km"][stored_keep]) ** 2)))
    sigma_doppler = float(np.sqrt(np.mean((doppler[stored_keep] - stored_prediction[stored_keep]) ** 2)))
    density, _ = physics.pbal.density_interpolator(sample_epoch)

    with h5py.File(args.beat_h5, "r") as beat:
        absolute_time_zero = float(np.asarray(beat["measurement/tx_idx"])[0]) / 1e6 - t_s[0]

    with h5py.File(args.baseline_profile, "r") as baseline, h5py.File(args.phase_profile, "r") as phase:
        radius_um = np.asarray(phase["profile/radius_um"], dtype=float)
        parameters6 = np.asarray(phase["profile/parameters6"], dtype=float)
        stored_chi2 = np.asarray(phase["profile/chi2"], dtype=float)
        echo_keep = np.asarray(phase["result/echo_shared_inlier_mask"], dtype=bool) & finite
        phase_groups = (phase["phase_acceleration"], phase["phase_acceleration_16ms"])

        names = ("position_radial", "position_transverse", "doppler", "phase_8ms", "phase_16ms")
        contribution = {name: np.full(len(radius_um), np.nan) for name in names}
        coordinate_position_chi2 = np.full((3, len(radius_um)), np.nan)
        for index, (radius, parameters) in enumerate(zip(radius_um, parameters6)):
            if not np.all(np.isfinite(parameters)):
                continue
            model, velocity, *_ = physics.propagate_shrinking_radius_model(
                np.r_[parameters, np.log10(radius * 1e-6)], t_s, density
            )
            prediction = physics.predicted_doppler(model, velocity)
            position_residual = points - model
            line_of_sight = model / np.maximum(np.linalg.norm(model, axis=1)[:, None], 1e-12)
            radial_residual = np.sum(position_residual * line_of_sight, axis=1)
            transverse_squared = np.sum(position_residual**2, axis=1) - radial_residual**2
            contribution["position_radial"][index] = np.sum(
                (radial_residual[echo_keep] / sigma_position) ** 2
            )
            contribution["position_transverse"][index] = np.sum(
                np.maximum(transverse_squared[echo_keep], 0.0) / sigma_position**2
            )
            coordinate_position_chi2[:, index] = np.sum(
                (position_residual[echo_keep] / sigma_position) ** 2, axis=0
            )
            contribution["doppler"][index] = np.sum(
                ((doppler[echo_keep] - prediction[echo_keep]) / sigma_doppler) ** 2
            )
            for name, group in zip(("phase_8ms", "phase_16ms"), phase_groups):
                samples = np.asarray(group["samples"])
                phase_prediction = predicted_phase(prediction * 1e3, t_s, absolute_time_zero, samples)
                contribution[name][index] = phase_chi2(group, phase_prediction)

        reconstructed = sum(contribution.values())
        mismatch = reconstructed - stored_chi2
        finite_mismatch = np.isfinite(mismatch)
        maximum_mismatch = float(np.max(np.abs(mismatch[finite_mismatch])))
        if maximum_mismatch > 1e-3:
            raise RuntimeError(f"profile reconstruction mismatch: {maximum_mismatch:g}")

    minimum_index = int(np.nanargmin(stored_chi2))
    delta = {name: values - values[minimum_index] for name, values in contribution.items()}
    total_delta = stored_chi2 - stored_chi2[minimum_index]

    fig, ax = plt.subplots(figsize=(6.4, 4.4), dpi=180)
    colors = ("C0", "C4", "C1", "C2", "C3")
    labels = ("Radial position", "Transverse position", "Doppler", "8 ms phase", "16 ms phase")
    for name, color, label in zip(names, colors, labels):
        ax.plot(radius_um, delta[name], color=color, lw=1.1, label=label)
    ax.plot(radius_um, total_delta, color="black", lw=1.5, label="Total")
    ax.axhline(3.841458820694124, color="0.4", ls="--", lw=0.8, label="95% threshold")
    ax.set_xscale("log")
    ax.set_xlabel(r"Fixed initial radius $r_0$ ($\mu$m)")
    ax.set_ylabel(r"Contribution relative to profile minimum")
    ax.set_ylim(-10.0, min(100.0, max(70.0, float(np.nanpercentile(total_delta, 60)))))
    ax.grid(alpha=0.2, which="both", linewidth=0.5)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    args.output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_plot, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)

    args.output_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output_h5, "w") as output:
        output.attrs["schema"] = "pansy.mass_profile_contribution_audit.v1"
        output.attrs["sample_idx"] = int(observations["sample_idx"])
        output.attrs["maximum_profile_reconstruction_error"] = maximum_mismatch
        output.attrs["sigma_position_km"] = sigma_position
        output.attrs["sigma_doppler_km_s"] = sigma_doppler
        output["radius_um"] = radius_um
        output["stored_chi2"] = stored_chi2
        output["total_delta_chi2"] = total_delta
        for name in names:
            output[f"chi2/{name}"] = contribution[name]
            output[f"delta_chi2/{name}"] = delta[name]
        output["chi2/position_enu"] = coordinate_position_chi2
        output["delta_chi2/position_enu"] = (
            coordinate_position_chi2 - coordinate_position_chi2[:, minimum_index, None]
        )

    print(f"sigma position: {sigma_position:.6g} km")
    print(f"sigma Doppler: {sigma_doppler:.6g} km/s")
    print(f"maximum reconstruction mismatch: {maximum_mismatch:.3g}")
    for target in (10.0, 100.0, 1000.0, 10000.0):
        index = int(np.argmin(np.abs(np.log(radius_um) - np.log(target))))
        pieces = " ".join(f"{name}={delta[name][index]:.3f}" for name in names)
        enu_delta = coordinate_position_chi2[:, index] - coordinate_position_chi2[:, minimum_index]
        print(
            f"r0={radius_um[index]:g} total={total_delta[index]:.3f} {pieces} "
            f"east={enu_delta[0]:.3f} north={enu_delta[1]:.3f} up={enu_delta[2]:.3f}"
        )

    target_indices = {
        target: int(np.argmin(np.abs(np.log(radius_um) - np.log(target))))
        for target in (100.0, 1000.0)
    }
    target_models = {}
    for target, index in target_indices.items():
        target_models[target], *_ = physics.propagate_shrinking_radius_model(
            np.r_[parameters6[index], np.log10(radius_um[index] * 1e-6)], t_s, density
        )
    best_residual = points[echo_keep] - target_models[100.0][echo_keep]
    print("position residual component correlation:")
    print(np.array2string(np.corrcoef(best_residual.T), precision=3, suppress_small=True))
    for component, label in enumerate(("east", "north", "up")):
        correlations, correlation_time, effective_samples = residual_correlation_summary(
            best_residual[:, component]
        )
        print(
            f"{label} residual_rms={np.sqrt(np.mean(best_residual[:, component]**2)):.4f} km "
            f"lag1={correlations[0]:.3f} correlation_time={correlation_time:.2f} "
            f"effective_samples={effective_samples:.1f}/{len(best_residual)}"
        )
    model_difference = target_models[1000.0][echo_keep] - target_models[100.0][echo_keep]
    print(
        f"100-to-1000 model separation rms={np.sqrt(np.mean(model_difference**2)):.4f} km "
        f"maximum_3d={np.max(np.linalg.norm(model_difference, axis=1)):.4f} km"
    )
    print(f"wrote {args.output_plot}")
    print(f"wrote {args.output_h5}")


if __name__ == "__main__":
    main()
