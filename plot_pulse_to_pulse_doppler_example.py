#!/usr/bin/env python3
"""Plot a compact pulse-to-pulse Doppler example for the PANSY paper."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


METEOROID_DENSITY_KG_M3 = 3000.0


def robust_limits(values: np.ndarray, padding_fraction: float = 0.12) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 0.0, 1.0
    lo, hi = np.nanpercentile(finite, [0.5, 99.5])
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo, hi = np.nanmin(finite), np.nanmax(finite)
    pad = max((hi - lo) * padding_fraction, 1e-6)
    return float(lo - pad), float(hi + pad)


def polynomial_refit(
    time_ms: np.ndarray,
    velocity_mps: np.ndarray,
    reference_mps: np.ndarray,
    keep: np.ndarray,
    degree: int = 3,
):
    """Fit a low-order residual correction to the main velocity branch."""
    finite = (
        keep
        & np.isfinite(time_ms)
        & np.isfinite(velocity_mps)
        & np.isfinite(reference_mps)
    )
    if np.count_nonzero(finite) <= degree + 1:
        correction = np.zeros(1)
        center_ms = 0.0
        scale_ms = 1.0

        def evaluate(target_time_ms, target_reference_mps):
            return np.asarray(target_reference_mps, dtype=float)

        return evaluate, correction, center_ms, scale_ms

    center_ms = float(np.nanmedian(time_ms[finite]))
    scale_ms = float(np.nanpercentile(time_ms[finite], 84) - np.nanpercentile(time_ms[finite], 16))
    if not np.isfinite(scale_ms) or scale_ms <= 0.0:
        scale_ms = 1.0
    x = (time_ms[finite] - center_ms) / scale_ms
    residual = velocity_mps[finite] - reference_mps[finite]
    coefficients = np.polyfit(x, residual, min(degree, np.count_nonzero(finite) - 1))

    def evaluate(target_time_ms, target_reference_mps):
        target_x = (np.asarray(target_time_ms, dtype=float) - center_ms) / scale_ms
        return np.asarray(target_reference_mps, dtype=float) + np.polyval(coefficients, target_x)

    return evaluate, coefficients, center_ms, scale_ms


def radius_um_to_mass_kg(radius_um: np.ndarray | float) -> np.ndarray:
    radius_m = np.maximum(np.asarray(radius_um, dtype=float), 1e-300) * 1e-6
    return (4.0 * np.pi / 3.0) * METEOROID_DENSITY_KG_M3 * radius_m**3


def mass_kg_to_radius_um(mass_kg: np.ndarray | float) -> np.ndarray:
    mass_kg = np.maximum(np.asarray(mass_kg, dtype=float), 1e-300)
    return (3.0 * mass_kg / (4.0 * np.pi * METEOROID_DENSITY_KG_M3)) ** (1.0 / 3.0) * 1e6


def profile_interval(radius_um: np.ndarray, delta_chi2: np.ndarray, threshold: float = 3.841458820694124):
    radius_um = np.asarray(radius_um, dtype=float)
    delta_chi2 = np.asarray(delta_chi2, dtype=float)
    good = np.isfinite(radius_um) & np.isfinite(delta_chi2)
    radius_um = radius_um[good]
    delta_chi2 = delta_chi2[good]
    best_index = int(np.nanargmin(delta_chi2))
    lower = float(radius_um[0])
    upper = float(radius_um[-1])
    lower_bounded = False
    upper_bounded = False

    def crossing(index_a: int, index_b: int) -> float:
        log_a = np.log(radius_um[index_a])
        log_b = np.log(radius_um[index_b])
        delta_a = delta_chi2[index_a]
        delta_b = delta_chi2[index_b]
        if delta_a == delta_b:
            return float(radius_um[index_a])
        return float(np.exp(log_a + (threshold - delta_a) * (log_b - log_a) / (delta_b - delta_a)))

    for index in range(best_index, 0, -1):
        if delta_chi2[index] <= threshold and delta_chi2[index - 1] > threshold:
            lower = crossing(index, index - 1)
            lower_bounded = True
            break
    for index in range(best_index, len(radius_um) - 1):
        if delta_chi2[index] <= threshold and delta_chi2[index + 1] > threshold:
            upper = crossing(index, index + 1)
            upper_bounded = True
            break
    return {
        "best_radius_um": float(radius_um[best_index]),
        "lower_radius_um": lower,
        "upper_radius_um": upper,
        "lower_bounded": lower_bounded,
        "upper_bounded": upper_bounded,
    }


def normalized_profile_density(radius_um: np.ndarray, density_log_radius: np.ndarray) -> np.ndarray:
    density = np.asarray(density_log_radius, dtype=float).copy()
    density[~np.isfinite(density)] = 0.0
    peak = np.nanmax(density) if np.any(np.isfinite(density)) else 0.0
    if peak > 0.0:
        density /= peak
    return density


def sci_tex(value: float, precision: int = 1) -> str:
    if not np.isfinite(value) or value <= 0.0:
        return str(value)
    exponent = int(np.floor(np.log10(value)))
    mantissa = value / 10.0**exponent
    return rf"{mantissa:.{precision}f}\times10^{{{exponent}}}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot three-pulse velocity/acceleration against single-pulse velocity."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(
            "figs/fig3_same_event_1757458600402300/three_pulse_full_event_regen/"
            "three_pulse_full_event_1757458600402300.h5"
        ),
    )
    parser.add_argument(
        "--single-pulse",
        type=Path,
        default=Path(
            "figs/paper_event_1757458600402300_single_pulse_raw_voltage_20260721/"
            "single_pulse_raw_voltage_1757458600402300.h5"
        ),
    )
    parser.add_argument(
        "--beat-observables",
        type=Path,
        default=Path(
            "figs/fig3_same_event_1757458600402300/"
            "fit_observables_catalogue_1757458600402300.h5"
        ),
    )
    parser.add_argument(
        "--phase-profile",
        type=Path,
        default=Path(
            "figs/fig3_same_event_1757458600402300/"
            "mass_profile_phase_aware_1757458600402300.h5"
        ),
    )
    parser.add_argument(
        "--three-pulse-mass-profile",
        type=Path,
        default=Path(
            "figs/fig3_same_event_1757458600402300/three_pulse_full_event_regen/"
            "sideband_rejected_mass_profile_sigma14_accel_clean_1757458600402300.h5"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("/Users/j/src/pansy_paper/pulse_to_pulse_doppler_example.png"),
    )
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    with h5py.File(args.single_pulse, "r") as handle:
        single = handle["single_pulse_fit"][:]
        single_radius_um = np.asarray(handle["dynamics_refit/profile_radius_um"], dtype=float)
        single_delta_chi2 = np.asarray(handle["dynamics_refit/profile_delta_chi2"], dtype=float)
        single_density = normalized_profile_density(
            single_radius_um,
            np.asarray(handle["dynamics_refit/profile_probability_density_log_radius"], dtype=float),
        )
    with h5py.File(args.three_pulse_mass_profile, "r") as handle:
        three_profile_radius_um = np.asarray(handle["profile/radius_um"], dtype=float)
        three_profile_delta_chi2 = np.asarray(handle["profile/delta_chi2"], dtype=float)
        three_profile_density = normalized_profile_density(
            three_profile_radius_um,
            np.asarray(handle["profile/probability_density_log_radius"], dtype=float),
        )
        three_ci = {
            "best_radius_um": float(handle["result/best_radius_um"][()]),
            "lower_radius_um": float(handle["result/ci95_lower_radius_um"][()]),
            "upper_radius_um": float(handle["result/ci95_upper_radius_um"][()]),
            "lower_bounded": handle["result"].attrs.get("ci95_lower_status", "") == "bounded",
            "upper_bounded": handle["result"].attrs.get("ci95_upper_status", "") == "bounded",
        }
    single_ci = profile_interval(single_radius_um, single_delta_chi2)
    finite_single = (
        single["shared_inlier"].astype(bool)
        & np.isfinite(single["velocity_mps"])
        & np.isfinite(single["velocity_std_mps"])
        & np.isfinite(single["model_velocity_mps"])
    )
    sample_idx = -1
    have_triplet_table = args.input.exists()
    if have_triplet_table:
        with h5py.File(args.input, "r") as handle:
            result = handle["triplet_fit"][:]
            sample_idx = int(handle.attrs.get("sample_idx", -1))

        required_fields = {
            "time_s",
            "velocity_mps",
            "velocity_std_mps",
            "model_velocity_mps",
            "acceleration_mps2",
            "acceleration_std_mps2",
            "model_acceleration_mps2",
            "shared_inlier",
        }
        missing = required_fields.difference(result.dtype.names or ())
        if missing:
            raise KeyError(f"{args.input} is missing fields: {sorted(missing)}")

        time_zero = np.nanmin(
            np.concatenate(
                (
                    result["time_s"][np.isfinite(result["time_s"])],
                    single["time_s"][np.isfinite(single["time_s"])],
                )
            )
        )
        time_ms = 1e3 * (result["time_s"] - time_zero)
        single_time_ms = 1e3 * (single["time_s"] - time_zero)
        keep = result["shared_inlier"].astype(bool)
        if "acceleration_at_bound" in (result.dtype.names or ()):
            keep &= ~result["acceleration_at_bound"].astype(bool)
        finite_velocity = (
            keep
            & np.isfinite(result["velocity_mps"])
            & np.isfinite(result["velocity_std_mps"])
            & np.isfinite(result["model_velocity_mps"])
        )
        finite_acceleration = (
            keep
            & np.isfinite(result["acceleration_mps2"])
            & np.isfinite(result["acceleration_std_mps2"])
            & np.isfinite(result["model_acceleration_mps2"])
        )
    else:
        with h5py.File(args.beat_observables, "r") as handle:
            beat = {name: np.asarray(values) for name, values in handle["baud_averaged_beat_pairs"].items()}
            sample_idx = int(handle.attrs.get("sample_idx", -1))
        with h5py.File(args.phase_profile, "r") as handle:
            group = handle["phase_acceleration"]
            phase = {name: np.asarray(values) for name, values in group.items()}
            phase_keep = np.asarray(group["shared_inlier_mask"], dtype=bool)

        beat_time = beat["time_s"]
        phase_time = phase["time_s"]
        absolute_zero = np.nanmin(beat_time[np.isfinite(beat_time)])
        time_ms = 1e3 * (beat_time - absolute_zero)
        phase_time_ms = 1e3 * (phase_time - absolute_zero)
        single_time_ms = 1e3 * single["time_s"]
        finite_velocity = (
            np.isfinite(beat["phase_doppler_mps"])
            & np.isfinite(beat["phase_doppler_std_mps"])
            & np.isfinite(beat["time_s"])
        )
        finite_acceleration = (
            phase_keep
            & np.isfinite(phase["measured_radial_acceleration_display_mps2"])
            & np.isfinite(phase["best_model_radial_acceleration_mps2"])
            & np.isfinite(phase_time_ms)
        )
        phase_variance = np.diag(phase["phase_difference_covariance_rad2"])
        phase_acceleration_std = (
            np.sqrt(np.maximum(phase_variance, 0.0))
            * phase["acceleration_ambiguity_period_mps2"]
            / (2.0 * np.pi)
        )

    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "legend.fontsize": 6.8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(3.45, 5.05),
        sharex=False,
        constrained_layout=True,
        height_ratios=(1.0, 0.78, 0.72),
    )

    ax = axes[0]
    if have_triplet_table:
        three_velocity = result["velocity_mps"]
        three_velocity_std = result["velocity_std_mps"]
        three_model = result["model_velocity_mps"]
        initial_three_residual = three_velocity - three_model
        three_main = finite_velocity & (np.abs(initial_three_residual) < 200.0)
        three_sideband = finite_velocity & ~three_main
        evaluate_velocity_refit, _coefficients, _center_ms, _scale_ms = polynomial_refit(
            time_ms,
            three_velocity,
            three_model,
            three_main,
            degree=3,
        )
        three_refit = evaluate_velocity_refit(time_ms, three_model)
        single_refit = evaluate_velocity_refit(single_time_ms, single["model_velocity_mps"])
        three_refit_residual = three_velocity - three_refit
        single_refit_residual = single["velocity_mps"] - single_refit
        three_sigma = float(np.nanstd(three_refit_residual[three_main], ddof=1))
        single_sigma = float(np.nanstd(single_refit_residual[finite_single], ddof=1))
    else:
        three_velocity = beat["phase_doppler_mps"]
        three_velocity_std = beat["phase_doppler_std_mps"]
        three_model = single["model_velocity_mps"][finite_single]
        three_main = finite_velocity
        three_sideband = np.zeros_like(finite_velocity, dtype=bool)
        three_refit = None
        three_sigma = np.nan
        single_sigma = float(
            np.nanstd(
                single["velocity_mps"][finite_single]
                - single["model_velocity_mps"][finite_single],
                ddof=1,
            )
        )

    ax.scatter(
        single_time_ms[finite_single],
        single["velocity_mps"][finite_single] / 1e3,
        s=8,
        facecolors="none",
        edgecolors="0.55",
        linewidths=0.55,
        label="single-pulse Doppler",
        zorder=1,
    )
    if have_triplet_table and np.any(three_sideband):
        ax.errorbar(
            time_ms[three_sideband],
            three_velocity[three_sideband] / 1e3,
            yerr=three_velocity_std[three_sideband] / 1e3,
            fmt="o",
            markersize=2.1,
            color="#D55E00",
            ecolor="#D55E00",
            elinewidth=0.35,
            capsize=0.8,
            label=r"rejected side band",
            zorder=3,
        )
    ax.errorbar(
        time_ms[three_main],
        three_velocity[three_main] / 1e3,
        yerr=three_velocity_std[three_main] / 1e3,
        fmt="o",
        markersize=2.2,
        color="black",
        ecolor="black",
        elinewidth=0.35,
        capsize=1.0,
        label="three-pulse fit Doppler",
        zorder=4,
    )
    model_time_ms = time_ms[finite_velocity] if have_triplet_table else single_time_ms[finite_single]
    model_velocity = (
        three_refit[finite_velocity]
        if have_triplet_table
        else single["model_velocity_mps"][finite_single]
    )
    order = np.argsort(model_time_ms)
    ax.plot(
        model_time_ms[order],
        model_velocity[order] / 1e3,
        color="#0072B2",
        linewidth=1.25,
        label="sideband-rejected refit" if have_triplet_table else "trajectory fit",
        zorder=5,
    )
    velocity_limit_values = np.concatenate(
        (
            result["velocity_mps"][finite_velocity],
            model_velocity,
            single["velocity_mps"][finite_single],
        )
        if have_triplet_table
        else (
            beat["phase_doppler_mps"][finite_velocity],
            single["model_velocity_mps"][finite_single],
            single["velocity_mps"][finite_single],
        )
    )
    vmin, vmax = robust_limits(velocity_limit_values / 1e3, padding_fraction=0.08)
    ax.set_ylim(vmin, vmax)
    ax.set_ylabel(r"Range rate (km s$^{-1}$)")
    if have_triplet_table:
        ax.text(
            0.03,
            0.96,
            rf"$\sigma_{{3p}}={three_sigma:.0f}\,$m s$^{{-1}}$"
            "\n"
            rf"$\sigma_{{1p}}={single_sigma:.0f}\,$m s$^{{-1}}$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7,
            color="0.2",
        )
    ax.legend(
        frameon=False,
        loc="lower right",
        handlelength=1.5,
        borderaxespad=0.2,
        labelspacing=0.25,
    )

    ax = axes[1]
    ax.errorbar(
        (time_ms if have_triplet_table else phase_time_ms)[finite_acceleration],
        (
            result["acceleration_mps2"][finite_acceleration]
            if have_triplet_table
            else phase["measured_radial_acceleration_display_mps2"][finite_acceleration]
        )
        / 1e3,
        yerr=(
            result["acceleration_std_mps2"][finite_acceleration]
            if have_triplet_table
            else phase_acceleration_std[finite_acceleration]
        )
        / 1e3,
        fmt="o",
        markersize=1.8,
        color="black",
        ecolor="black",
        elinewidth=0.35,
        capsize=0.8,
        label="three-pulse acceleration",
        alpha=0.85,
        zorder=3,
    )
    acceleration_time_ms = (time_ms if have_triplet_table else phase_time_ms)[
        finite_acceleration
    ]
    model_acceleration = (
        result["model_acceleration_mps2"][finite_acceleration]
        if have_triplet_table
        else phase["best_model_radial_acceleration_mps2"][finite_acceleration]
    )
    order = np.argsort(acceleration_time_ms)
    ax.plot(
        acceleration_time_ms[order],
        model_acceleration[order] / 1e3,
        color="#0072B2",
        linewidth=1.25,
        label="trajectory fit",
        zorder=4,
    )
    acceleration_limit_values = np.concatenate(
        (
            result["acceleration_mps2"][finite_acceleration],
            result["model_acceleration_mps2"][finite_acceleration],
        )
        if have_triplet_table
        else (
            phase["measured_radial_acceleration_display_mps2"][finite_acceleration],
            phase["best_model_radial_acceleration_mps2"][finite_acceleration],
        )
    )
    amin, amax = robust_limits(acceleration_limit_values / 1e3, padding_fraction=0.14)
    ax.set_ylim(amin, amax)
    ax.set_xlabel("Time since first triplet (ms)")
    ax.set_ylabel(r"Acceleration (km s$^{-2}$)")
    ax.legend(
        frameon=False,
        loc="upper left",
        handlelength=1.5,
        borderaxespad=0.2,
        labelspacing=0.25,
    )

    ax = axes[2]
    single_mass_kg = radius_um_to_mass_kg(single_radius_um)
    three_mass_kg = radius_um_to_mass_kg(three_profile_radius_um)
    ax.plot(
        single_mass_kg,
        single_density,
        color="0.48",
        linewidth=1.25,
        label="single-pulse profile",
    )
    ax.plot(
        three_mass_kg,
        three_profile_density,
        color="#0072B2",
        linewidth=1.45,
        label="three-pulse profile",
    )

    single_lower_mass = float(radius_um_to_mass_kg(single_ci["lower_radius_um"]))
    single_upper_mass = float(radius_um_to_mass_kg(single_ci["upper_radius_um"]))
    three_lower_mass = float(radius_um_to_mass_kg(three_ci["lower_radius_um"]))
    three_upper_mass = float(radius_um_to_mass_kg(three_ci["upper_radius_um"]))
    three_best_mass = float(radius_um_to_mass_kg(three_ci["best_radius_um"]))
    single_best_mass = float(radius_um_to_mass_kg(single_ci["best_radius_um"]))
    ax.axvspan(single_lower_mass, single_upper_mass, color="0.6", alpha=0.16, linewidth=0)
    ax.axvspan(three_lower_mass, three_upper_mass, color="#0072B2", alpha=0.18, linewidth=0)
    ax.axvline(single_best_mass, color="0.42", linewidth=0.9, linestyle=":")
    ax.axvline(three_best_mass, color="#0072B2", linewidth=0.95, linestyle=":")
    ax.set_xscale("log")
    mass_min = float(radius_um_to_mass_kg(80.0))
    mass_max = float(radius_um_to_mass_kg(2500.0))
    ax.set_xlim(mass_min, mass_max)
    ax.set_ylim(0.0, 1.08)
    ax.set_xlabel(r"Initial mass $m_0$ (kg)")
    ax.set_ylabel("Relative density")
    ax.legend(
        frameon=False,
        loc="upper right",
        handlelength=1.5,
        borderaxespad=0.2,
        labelspacing=0.25,
    )
    radius_axis = ax.secondary_xaxis("top", functions=(mass_kg_to_radius_um, radius_um_to_mass_kg))
    radius_axis.set_xlabel(r"Initial radius $r_0$ ($\mu$m)")

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, color="0.88", linewidth=0.45)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi)
    fig.savefig(args.output.with_suffix(".pdf"))
    print(args.output)
    print(args.output.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
