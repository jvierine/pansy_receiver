#!/usr/bin/env python3
"""Audit inter-pulse acceleration aliases against range and Doppler."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

import pansy_config as pc


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--phase-profile", type=Path, required=True)
    parser.add_argument("--output-plot", type=Path, required=True)
    parser.add_argument("--output-h5", type=Path, required=True)
    return parser.parse_args()


def rms(values):
    return float(np.sqrt(np.mean(np.asarray(values, dtype=float) ** 2)))


def fit_free_acceleration(t_s, range_km, doppler_km_s):
    range_coefficients = np.polyfit(t_s, range_km, 2)
    doppler_coefficients = np.polyfit(t_s, doppler_km_s, 1)
    sigma_range = rms(range_km - np.polyval(range_coefficients, t_s))
    sigma_doppler = rms(doppler_km_s - np.polyval(doppler_coefficients, t_s))

    range_design = np.column_stack((np.ones_like(t_s), t_s, np.zeros_like(t_s), 0.5 * t_s**2))
    doppler_design = np.column_stack((np.zeros_like(t_s), np.zeros_like(t_s), np.ones_like(t_s), t_s))
    design = np.vstack((range_design / sigma_range, doppler_design / sigma_doppler))
    values = np.r_[range_km / sigma_range, doppler_km_s / sigma_doppler]
    parameters, *_ = np.linalg.lstsq(design, values, rcond=None)
    covariance = np.linalg.inv(design.T @ design)
    return parameters, covariance, sigma_range, sigma_doppler


def fit_fixed_acceleration(t_s, range_km, doppler_km_s, acceleration_km_s2):
    range_design = np.column_stack((np.ones_like(t_s), t_s))
    range_values = range_km - 0.5 * acceleration_km_s2 * t_s**2
    range_parameters, *_ = np.linalg.lstsq(range_design, range_values, rcond=None)
    range_model = range_design @ range_parameters + 0.5 * acceleration_km_s2 * t_s**2

    doppler_intercept = float(np.mean(doppler_km_s - acceleration_km_s2 * t_s))
    doppler_model = doppler_intercept + acceleration_km_s2 * t_s
    return range_model, doppler_model


def circular_residual(observed, predicted):
    return np.angle(np.exp(1j * (observed - predicted)))


def phase_chi2(group, acceleration_km_s2):
    selected = np.asarray(group["shared_inlier_mask"], dtype=bool)
    samples = np.asarray(group["samples"])
    observed = samples["observed_delta_phase_rad"]
    mean_dt = 0.5 * (samples["first_dt_s"] + samples["second_dt_s"])
    lag_s = float(group.attrs["acceleration_lag_s"])
    best_acceleration = np.asarray(group["best_model_radial_acceleration_mps2"], dtype=float)
    best_prediction = np.asarray(group["best_model_prediction_rad"], dtype=float)
    phase_shift = (
        4.0
        * np.pi
        * (acceleration_km_s2 * 1e3 - best_acceleration)
        * mean_dt
        * lag_s
        / pc.wavelength
    )
    prediction = best_prediction + phase_shift
    residual = circular_residual(observed[selected], prediction[selected])
    covariance = np.asarray(group["phase_difference_covariance_rad2"], dtype=float)
    covariance = covariance[np.ix_(selected, selected)]
    scale = float(group.attrs["covariance_scale"])
    return float(residual @ np.linalg.solve(covariance, residual) / scale**2)


def main():
    args = parse_args()
    with h5py.File(args.diagnostics, "r") as diagnostics, h5py.File(args.phase_profile, "r") as phase:
        hypothesis = str(diagnostics.attrs["selected_hypothesis"])
        observations = diagnostics[f"hypotheses/{hypothesis}"]
        t_all = np.asarray(observations["fit_t_s"], dtype=float)
        range_all = np.asarray(observations["range_km"], dtype=float)
        doppler_all = np.asarray(observations["doppler_mps"], dtype=float) / 1e3
        keep = np.asarray(phase["result/echo_shared_inlier_mask"], dtype=bool)
        finite = np.isfinite(t_all) & np.isfinite(range_all) & np.isfinite(doppler_all)
        keep &= finite
        t_s = t_all[keep]
        range_km = range_all[keep]
        doppler_km_s = doppler_all[keep]

        free, free_covariance, sigma_range, sigma_doppler = fit_free_acceleration(
            t_s, range_km, doppler_km_s
        )
        independent_acceleration = float(free[3])
        independent_acceleration_std = float(np.sqrt(free_covariance[3, 3]))

        phase8 = phase["phase_acceleration"]
        phase8_keep = np.asarray(phase8["shared_inlier_mask"], dtype=bool)
        principal = np.asarray(
            phase8["measured_radial_acceleration_principal_mps2"], dtype=float
        )[phase8_keep] / 1e3
        ambiguity = float(
            np.median(np.asarray(phase8["acceleration_ambiguity_period_mps2"])[phase8_keep])
            / 1e3
        )
        phase_branch_zero = float(np.median(principal))
        phase_branch_zero += round((independent_acceleration - phase_branch_zero) / ambiguity) * ambiguity
        branch_index = np.arange(-2, 3, dtype=int)
        branch_acceleration = phase_branch_zero + branch_index * ambiguity

        range_chi2 = np.empty(len(branch_index))
        doppler_chi2 = np.empty(len(branch_index))
        phase8_chi2 = np.empty(len(branch_index))
        phase16_chi2 = np.empty(len(branch_index))
        range_models = np.empty((len(branch_index), len(t_s)))
        doppler_models = np.empty((len(branch_index), len(t_s)))
        for index, acceleration in enumerate(branch_acceleration):
            range_model, doppler_model = fit_fixed_acceleration(
                t_s, range_km, doppler_km_s, acceleration
            )
            range_models[index] = range_model
            doppler_models[index] = doppler_model
            range_chi2[index] = np.sum(((range_km - range_model) / sigma_range) ** 2)
            doppler_chi2[index] = np.sum(((doppler_km_s - doppler_model) / sigma_doppler) ** 2)
            phase8_chi2[index] = phase_chi2(phase8, acceleration)
            phase16_chi2[index] = phase_chi2(phase["phase_acceleration_16ms"], acceleration)

        components = np.vstack((range_chi2, doppler_chi2, phase8_chi2, phase16_chi2))
        delta_components = components - components[:, branch_index == 0]
        total = np.sum(components, axis=0)
        delta_total = total - np.min(total)
        sample_idx = int(phase.attrs["sample_idx"])

    colors = {-1: "C1", 0: "C0", 1: "C3"}
    labels = {
        k: rf"$k={k:+d}$, $\ddot{{R}}={branch_acceleration[branch_index == k][0]:.1f}$ km s$^{{-2}}$"
        for k in (-1, 0, 1)
    }
    fig = plt.figure(figsize=(8.2, 6.2), dpi=180)
    grid = fig.add_gridspec(2, 2, height_ratios=(1.0, 0.9), hspace=0.28, wspace=0.27)
    ax_range = fig.add_subplot(grid[0, 0])
    ax_doppler = fig.add_subplot(grid[0, 1])
    ax_chi2 = fig.add_subplot(grid[1, :])

    ax_range.scatter(t_s, range_km, s=5, color="black", alpha=0.65, zorder=3)
    ax_doppler.scatter(t_s, doppler_km_s, s=5, color="black", alpha=0.65, zorder=3)
    for k in (-1, 0, 1):
        index = int(np.flatnonzero(branch_index == k)[0])
        ax_range.plot(t_s, range_models[index], color=colors[k], lw=1.2, label=labels[k])
        ax_doppler.plot(t_s, doppler_models[index], color=colors[k], lw=1.2)
    ax_range.set_xlabel("Time (s)")
    ax_range.set_ylabel("Range (km)")
    ax_doppler.set_xlabel("Time (s)")
    ax_doppler.set_ylabel(r"Doppler velocity (km s$^{-1}$)")
    ax_range.legend(frameon=False, fontsize=7, loc="best")

    bottoms = np.zeros(len(branch_index))
    component_labels = ("Range", "Doppler", "8 ms phase", "16 ms phase")
    component_colors = ("C0", "C1", "C2", "C4")
    for values, label, color in zip(delta_components, component_labels, component_colors):
        positive = np.maximum(values, 0.0)
        ax_chi2.bar(branch_index, positive, bottom=bottoms, width=0.65, color=color, label=label)
        bottoms += positive
    ax_chi2.plot(branch_index, np.maximum(delta_total, 1e-3), "ko", ms=3, label=r"Total $\Delta\chi^2$")
    ax_chi2.set_yscale("symlog", linthresh=1.0)
    ax_chi2.set_xticks(branch_index)
    ax_chi2.set_xticklabels(
        [f"{k:+d}\n{a:.1f}" for k, a in zip(branch_index, branch_acceleration)]
    )
    ax_chi2.set_xlabel(r"Common alias index $k$ and range acceleration $\ddot{R}$ (km s$^{-2}$)")
    ax_chi2.set_ylabel(r"Penalty relative to $k=0$")
    ax_chi2.legend(frameon=False, fontsize=7, ncol=5, loc="upper center")
    phase_inset = ax_chi2.inset_axes([0.69, 0.12, 0.28, 0.42])
    phase_delta = (phase8_chi2 + phase16_chi2) - (phase8_chi2 + phase16_chi2)[branch_index == 0]
    phase_inset.plot(branch_index, phase_delta, "o-", color="C2", ms=3, lw=0.9)
    phase_inset.axhline(0.0, color="0.3", lw=0.6)
    phase_inset.set_xticks(branch_index)
    phase_inset.set_title(r"Phase-only $\Delta\chi^2$", fontsize=7)
    phase_inset.tick_params(labelsize=6)
    phase_inset.grid(alpha=0.2, linewidth=0.4)
    for axis in (ax_range, ax_doppler, ax_chi2):
        axis.grid(alpha=0.2, linewidth=0.5)
    fig.text(
        0.5,
        0.995,
        rf"Event {sample_idx}: independent range acceleration "
        rf"$\ddot{{R}}={independent_acceleration:.2f}\pm{independent_acceleration_std:.2f}$ km s$^{{-2}}$",
        ha="center",
        va="top",
        fontsize=9,
    )
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.10, top=0.94)
    args.output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_plot, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)

    args.output_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output_h5, "w") as output:
        output.attrs["schema"] = "pansy.inter_pulse_acceleration_alias_audit.v1"
        output.attrs["sample_idx"] = sample_idx
        output.attrs["diagnostics_h5"] = str(args.diagnostics)
        output.attrs["phase_profile_h5"] = str(args.phase_profile)
        output.attrs["mask_policy"] = "one shared echo mask from the phase-aware best fit, fixed for every alias"
        output["branch_index"] = branch_index
        output["branch_acceleration_km_s2"] = branch_acceleration
        output["independent_range_acceleration_km_s2"] = independent_acceleration
        output["independent_range_acceleration_std_km_s2"] = independent_acceleration_std
        output["range_chi2"] = range_chi2
        output["doppler_chi2"] = doppler_chi2
        output["phase8_chi2"] = phase8_chi2
        output["phase16_chi2"] = phase16_chi2
        output["total_delta_chi2"] = delta_total
        output["echo_shared_inlier_mask"] = keep
    print(f"independent range acceleration: {independent_acceleration:.4f} +/- {independent_acceleration_std:.4f} km/s^2")
    for k, acceleration, delta_value in zip(branch_index, branch_acceleration, delta_total):
        print(f"k={k:+d} acceleration={acceleration:.3f} km/s^2 delta_chi2={delta_value:.3g}")
    print(f"wrote {args.output_plot}")
    print(f"wrote {args.output_h5}")


if __name__ == "__main__":
    main()
