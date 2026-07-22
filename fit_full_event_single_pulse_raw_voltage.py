#!/usr/bin/env python3
"""Fit independent raw-voltage Doppler measurements for one PANSY event."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import pansy_config as pc
import pansy_interferometry as pint
import pansy_modes as pmm
from aod_complex_doppler_fit import event_receiver_noise_variance
from fit_full_event_three_pulse_complex_envelope import (
    FS_HZ,
    beam_pixmap,
    refit_dynamics,
    uniform_ipp_rti,
)
from fit_three_pulse_complex_envelope import fit_single_pulse_raw_voltage
from inter_pulse_phase_deceleration import (
    decoded_pulse_responses,
    diagnostic_measurement_clock,
    diagnostic_to_raw_pulses,
    fractional_segment,
    load_selected,
    precise_matched_filter_estimates,
)
from interferometer_alias_diagnostics import load_cut, recompute_cut_observables
from pansy_coherent import (
    coherent_add_modules,
    enu_line_of_sight_to_arrival_uvw,
    estimate_event_module_voltage_gain,
)
from plot_interferometric_disambiguation import refine_coherence_peak
from run_catalogue_mass_profiles import load_selected_fit, profile_interval, radius_um_to_mass_kg


def write_mass_profile_compatible_h5(
    path: Path,
    args: argparse.Namespace,
    observations: dict,
    dynamics: dict,
    echo_keep: np.ndarray,
    fitted_points_km: np.ndarray,
    doppler_residual: np.ndarray,
    velocity_keep: np.ndarray,
    position_rms: np.ndarray,
    path_length_km: float,
    selected_hypothesis: str,
) -> None:
    """Write the single-pulse refit in the catalogue mass-profile summary schema."""
    radius_um = np.asarray(dynamics["profile_radius_um"], dtype=np.float64)
    delta_chi2 = np.asarray(dynamics["profile_delta_chi2"], dtype=np.float64)
    threshold = 3.841458820694124
    lower_um, upper_um, lower_bounded, upper_bounded, lower_status, upper_status = profile_interval(
        radius_um, delta_chi2, threshold
    )
    best_radius_um = float(1e6 * 10.0 ** np.asarray(dynamics["parameters"], dtype=float)[6])
    velocity = np.asarray(dynamics["velocity_km_s"], dtype=float)
    temporary = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(temporary, "w") as handle:
        handle.attrs["schema"] = "pansy.catalogue_mass_profile.single_pulse_raw_voltage.v1"
        handle.attrs["sample_idx"] = int(args.sample_idx)
        handle.attrs["sample_epoch_unix"] = float(observations["sample_epoch_unix"])
        handle.attrs["diagnostics_h5"] = str(args.diagnostics_h5)
        handle.attrs["selected_hypothesis"] = str(selected_hypothesis)
        handle.attrs["meteoroid_density_kg_m3"] = 3000.0
        handle.attrs["measurement"] = "single-pulse raw-voltage Doppler refit used for Figure 3"
        profile = handle.create_group("profile")
        profile.create_dataset("radius_um", data=radius_um)
        profile.create_dataset("mass_kg", data=radius_um_to_mass_kg(radius_um))
        profile.create_dataset("chi2", data=np.asarray(dynamics["profile_chi2"], dtype=np.float64))
        profile.create_dataset("delta_chi2", data=delta_chi2)
        profile.create_dataset(
            "relative_probability_density_log_radius",
            data=np.asarray(dynamics["profile_probability_density_log_radius"], dtype=np.float64),
        )
        profile.create_dataset(
            "probability_weight",
            data=np.asarray(dynamics["profile_probability_weights"], dtype=np.float64),
        )
        profile.create_dataset("success", data=np.asarray(dynamics["profile_success"], dtype=bool))
        profile.create_dataset("parameters6", data=np.asarray(dynamics["profile_parameters6"], dtype=np.float64))

        result = handle.create_group("result")
        result.create_dataset("free_best_parameters7", data=np.asarray(dynamics["parameters"], dtype=np.float64))
        result.create_dataset("free_best_radius_um", data=best_radius_um)
        result.create_dataset("free_best_mass_kg", data=radius_um_to_mass_kg(best_radius_um))
        result.create_dataset("profile_grid_best_radius_um", data=float(radius_um[int(np.nanargmin(delta_chi2))]))
        result.create_dataset("profile_grid_best_mass_kg", data=radius_um_to_mass_kg(radius_um[int(np.nanargmin(delta_chi2))]))
        result.create_dataset("profile_ci95_lower_radius_um", data=lower_um)
        result.create_dataset("profile_ci95_upper_radius_um", data=upper_um)
        result.create_dataset("profile_ci95_lower_mass_kg", data=radius_um_to_mass_kg(lower_um))
        result.create_dataset("profile_ci95_upper_mass_kg", data=radius_um_to_mass_kg(upper_um))
        result.create_dataset("profile_ci95_lower_bounded", data=lower_bounded)
        result.create_dataset("profile_ci95_upper_bounded", data=upper_bounded)
        result.create_dataset("profile_ci95_delta_chi2_threshold", data=threshold)
        result.create_dataset(
            "marginal_radius_quantiles_um",
            data=np.asarray(dynamics["marginal_radius_quantiles_um"], dtype=np.float64),
        )
        result.create_dataset(
            "marginal_mass_quantiles_kg",
            data=np.asarray(dynamics["marginal_mass_quantiles_kg"], dtype=np.float64),
        )
        result.create_dataset("minimum_chi2", data=float(np.nanmin(np.asarray(dynamics["profile_chi2"], dtype=float))))
        result.create_dataset("free_chi2", data=float(np.nanmin(np.asarray(dynamics["profile_chi2"], dtype=float))))
        result.attrs["profile_ci95_lower_status"] = lower_status
        result.attrs["profile_ci95_upper_status"] = upper_status

        quality = handle.create_group("quality")
        quality.create_dataset("n_measurements", data=int(np.count_nonzero(echo_keep)))
        quality.create_dataset(
            "sigma_position_component_km",
            data=float(np.sqrt(np.nanmean(np.asarray(position_rms, dtype=float) ** 2))),
        )
        quality.create_dataset(
            "sigma_doppler_km_s",
            data=float(np.sqrt(np.nanmean(np.asarray(doppler_residual[velocity_keep], dtype=float) ** 2))) / 1e3,
        )
        quality.create_dataset(
            "free_position_3d_rms_km",
            data=float(
                np.sqrt(
                    np.nanmean(
                        np.sum((fitted_points_km[echo_keep] - np.asarray(dynamics["position_km"])[echo_keep]) ** 2, axis=1)
                    )
                )
            ),
        )
        quality.create_dataset(
            "free_doppler_rms_km_s",
            data=float(np.sqrt(np.nanmean(np.asarray(doppler_residual[velocity_keep], dtype=float) ** 2))) / 1e3,
        )
        quality.create_dataset("path_length_km", data=float(path_length_km))
        quality.create_dataset("initial_speed_km_s", data=float(np.linalg.norm(np.asarray(dynamics["parameters"])[3:6]) / 1e3))
        quality.create_dataset("fitted_speed_mean_km_s", data=float(np.nanmean(np.linalg.norm(velocity[echo_keep], axis=1))))
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-idx", type=int, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--diagnostics-h5", type=Path, required=True)
    parser.add_argument("--initial-fit-h5", type=Path)
    parser.add_argument("--prior-profile-h5", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    parser.add_argument("--frequency-half-width-hz", type=float, default=1.5e4)
    parser.add_argument("--frequency-grid-size", type=int, default=2001)
    parser.add_argument("--refine-interferometry-subgrid", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--mass-profile-output-dir", type=Path)
    args = parser.parse_args()

    cut = load_cut(args.base / "metadata/cut", args.sample_idx)
    hypothesis = load_selected(args.diagnostics_h5)
    if args.initial_fit_h5 is None:
        precise = precise_matched_filter_estimates(cut, hypothesis, args.snr_threshold)
        precise_range_km = np.asarray(precise["range_km"], dtype=float)
        precise_doppler_mps = np.asarray(precise["doppler_mps"], dtype=float)
        precise_raw_idx = np.asarray(precise["raw_idx"], dtype=int)
    else:
        with h5py.File(args.initial_fit_h5, "r") as handle:
            precise_range_km = np.asarray(handle["range_km"], dtype=float)
            precise_doppler_mps = np.asarray(handle["doppler_mps"], dtype=float)
            precise_raw_idx = (
                np.asarray(handle["raw_idx"], dtype=int) if "raw_idx" in handle else None
            )

    precise_selection = None
    if precise_raw_idx is not None:
        clock = diagnostic_measurement_clock(cut, hypothesis, args.snr_threshold)
        target_raw_idx = diagnostic_to_raw_pulses(cut, clock, hypothesis["t_rel_s"])
        raw_to_precise = {int(raw): index for index, raw in enumerate(precise_raw_idx)}
        missing = [int(raw) for raw in target_raw_idx if int(raw) not in raw_to_precise]
        if missing:
            raise RuntimeError(f"precise sidecar is missing raw pulse indices: {missing}")
        precise_selection = np.asarray(
            [raw_to_precise[int(raw)] for raw in target_raw_idx], dtype=int
        )
        precise_range_km = precise_range_km[precise_selection]
        precise_doppler_mps = precise_doppler_mps[precise_selection]
    decoded = decoded_pulse_responses(
        cut,
        hypothesis,
        args.snr_threshold,
        precise_range_km=precise_range_km,
        precise_doppler_mps=precise_doppler_mps,
    )
    observations, stored_fit = load_selected_fit(args.diagnostics_h5)
    with h5py.File(args.prior_profile_h5, "r") as handle:
        if "result/echo_shared_inlier_mask" in handle:
            prior_echo_keep = np.asarray(handle["result/echo_shared_inlier_mask"], dtype=bool)
        else:
            prior_echo_keep = np.asarray(stored_fit["keep"], dtype=bool)
        if "position_residual_covariance_km2" in handle:
            position_covariance_km2 = np.asarray(
                handle["position_residual_covariance_km2"], dtype=float
            )
        else:
            residual = np.asarray(observations["points_km"], dtype=float) - np.asarray(
                stored_fit["model_km"], dtype=float
            )
            covariance_keep = prior_echo_keep & np.all(np.isfinite(residual), axis=1)
            if np.count_nonzero(covariance_keep) >= 3:
                position_covariance_km2 = np.cov(residual[covariance_keep].T)
            else:
                position_covariance_km2 = np.eye(3) * 0.1**2
            position_covariance_km2 = np.asarray(position_covariance_km2, dtype=float)
            position_covariance_km2 += np.eye(3) * 1e-8
        delta = np.asarray(handle["profile/delta_chi2"], dtype=float)
        best_index = int(np.nanargmin(delta))
        parameters6 = np.asarray(handle["profile/parameters6"][best_index], dtype=float)
        radius_um = float(handle["profile/radius_um"][best_index])
        profile_grid_radius_um = np.asarray(handle["profile/radius_um"], dtype=float)
    if precise_selection is not None:
        echo_keep = np.zeros(len(precise_selection), dtype=bool)
        covered = precise_selection < len(prior_echo_keep)
        echo_keep[covered] = prior_echo_keep[precise_selection[covered]]
    else:
        echo_keep = prior_echo_keep.copy()

    trajectory_time = np.asarray(observations["t_s"], dtype=float)
    import fit_best_alias_physics_models as physics

    steering_points_km = np.asarray(stored_fit["model_km"], dtype=float)
    steering_uvw = enu_line_of_sight_to_arrival_uvw(steering_points_km)
    module_noise_variance = event_receiver_noise_variance(
        decoded["z_rx"],
        decoded["raw_idx"],
        decoded["range_gate"],
        decoded["z_tx"].shape[1],
    )
    module_voltage_gain = estimate_event_module_voltage_gain(decoded["response"])
    beamformed_echo = np.stack(
        [
            coherent_add_modules(
                decoded["z_rx"][raw_index],
                steering_uvw[observation],
                int(decoded["beam_id"][observation]),
                module_noise_variance=module_noise_variance[raw_index],
                module_voltage_gain=module_voltage_gain,
                normalization="weighted_mean",
            )
            for observation, raw_index in enumerate(decoded["raw_idx"])
        ]
    )
    absolute_zero = decoded["tx_idx"][0] / FS_HZ - trajectory_time[0]
    n_fast = decoded["z_tx"].shape[1]
    dtype = [
        ("observation", "i4"),
        ("time_s", "f8"),
        ("velocity_mps", "f8"),
        ("velocity_std_mps", "f8"),
        ("acceleration_mps2", "f8"),
        ("acceleration_std_mps2", "f8"),
        ("model_acceleration_mps2", "f8"),
        ("initial_fft_velocity_mps", "f8"),
        ("model_velocity_mps", "f8"),
        ("weighted_sse", "f8"),
        ("degrees_of_freedom", "i4"),
        ("shared_inlier", "?"),
    ]
    initial_model_doppler_mps = (
        physics.predicted_doppler(
            np.asarray(stored_fit["model_km"], dtype=float),
            np.asarray(stored_fit["velocity_km_s"], dtype=float),
        )
        * 1e3
    )
    rows = []
    fits = []
    for observation, raw_index in enumerate(decoded["raw_idx"]):
        raw_pulse = fractional_segment(
            beamformed_echo[observation], decoded["range_gate"][observation], n_fast
        )
        template = decoded["z_tx"][raw_index].astype(np.complex128)
        fit = fit_single_pulse_raw_voltage(
            raw_pulse,
            template,
            2.0 * float(decoded["coarse_doppler_mps"][observation]) / pc.wavelength,
            pc.wavelength,
            frequency_half_width_hz=args.frequency_half_width_hz,
            frequency_grid_size=args.frequency_grid_size,
        )
        fit_time = (
            decoded["tx_idx"][observation] / FS_HZ
            - absolute_zero
            + fit["time_center_s"]
        )
        rows.append(
            (
                observation,
                fit_time,
                fit["velocity_mps"],
                fit["velocity_std_mps"],
                np.nan,
                np.nan,
                np.nan,
                decoded["coarse_doppler_mps"][observation],
                float(np.interp(fit_time, trajectory_time, initial_model_doppler_mps)),
                fit["weighted_sse"],
                fit["degrees_of_freedom"],
                bool(echo_keep[observation]),
            )
        )
        fits.append(fit)
    result = np.asarray(rows, dtype=dtype)

    fitted_points_km = np.asarray(observations["points_km"], dtype=float).copy()
    if args.refine_interferometry_subgrid:
        ch_pairs = np.asarray(pint.ch_pairs, dtype=int)
        dmat = pint.pair_mat(ch_pairs, pint.get_antpos())
        phasecal = pint.get_phasecal()
        refined_direction = np.asarray(hypothesis["direction_uvw"], dtype=float).copy()
        for observation in range(len(refined_direction)):
            refined_direction[observation] = refine_coherence_peak(
                decoded["xc"][observation],
                int(decoded["beam_id"][observation]),
                phasecal,
                ch_pairs,
                dmat,
                float(refined_direction[observation, 0]),
                float(refined_direction[observation, 1]),
                2.0 / 500.0,
            )[:3]
        fitted_points_km = precise_range_km[:, None] * np.column_stack(
            (refined_direction[:, 0], refined_direction[:, 1], -refined_direction[:, 2])
        )
    catalogue_range_km = np.linalg.norm(fitted_points_km, axis=1)
    valid_range = np.isfinite(precise_range_km) & (catalogue_range_km > 0.0)
    fitted_points_km[valid_range] *= (
        precise_range_km[valid_range] / catalogue_range_km[valid_range]
    )[:, None]

    density, _metadata = physics.pbal.density_interpolator(
        observations["sample_epoch_unix"]
    )
    velocity_keep = np.asarray(result["shared_inlier"], dtype=bool)
    acceleration_keep = np.zeros(len(result), dtype=bool)
    correlation = np.eye(2 * len(result))
    dynamics = refit_dynamics(
        physics,
        trajectory_time,
        fitted_points_km,
        position_covariance_km2,
        echo_keep,
        result,
        np.r_[parameters6, np.log10(radius_um * 1e-6)],
        density,
        profile_grid_radius_um,
        triplet_correlation=correlation,
        velocity_keep=velocity_keep,
        acceleration_keep=acceleration_keep,
        compute_profile=True,
        velocity_sigma_override_mps=None,
        acceleration_sigma_override_mps2=None,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_h5 = args.output_dir / f"single_pulse_raw_voltage_{args.sample_idx}.h5"
    event_png = args.output_dir / f"single_pulse_raw_voltage_plot_{args.sample_idx}.png"
    with h5py.File(output_h5, "w") as handle:
        handle.attrs["schema"] = "pansy.single_pulse_raw_voltage_event.v1"
        handle.attrs["sample_idx"] = args.sample_idx
        handle.attrs["measurement"] = (
            "independent amplitude, phase, and Doppler fit to each fractional-range-aligned raw-voltage pulse"
        )
        handle.attrs["voltage_source"] = (
            "AoA-steered inverse-noise coherent sum of seven receiver modules"
        )
        handle.attrs["acceleration_measurement_used"] = False
        handle.create_dataset("single_pulse_fit", data=result)
        handle.create_dataset("event_module_voltage_gain", data=module_voltage_gain)
        fit_group = handle.create_group("dynamics_refit")
        for name, value in dynamics.items():
            if isinstance(value, np.ndarray):
                fit_group.create_dataset(name, data=value)
            else:
                fit_group.attrs[name] = value

    refit_position = dynamics["position_km"]
    refit_doppler = dynamics["doppler_mps"]
    refit_speed = np.linalg.norm(dynamics["velocity_km_s"], axis=1)
    time_origin = float(trajectory_time[echo_keep][0])
    observation_time = trajectory_time - time_origin
    measurement_time = result["time_s"] - time_origin
    refit_range = np.linalg.norm(refit_position, axis=1)
    fitted_residual = fitted_points_km - refit_position
    position_rms = np.sqrt(np.mean(fitted_residual[echo_keep] ** 2, axis=0))
    range_residual = precise_range_km[echo_keep] - refit_range[echo_keep]
    measured_model = np.interp(result["time_s"], trajectory_time, refit_doppler)
    doppler_residual = result["velocity_mps"] - measured_model
    rms = lambda values: float(np.sqrt(np.nanmean(np.asarray(values) ** 2)))
    path_length_km = float(
        np.sum(np.linalg.norm(np.diff(refit_position[echo_keep], axis=0), axis=1))
    )
    snr_db = 10.0 * np.log10(
        np.maximum(np.asarray(observations["snr"], dtype=float), 1e-12)
    )
    if args.mass_profile_output_dir is not None:
        mass_profile_path = args.mass_profile_output_dir / f"mass_profile_{args.sample_idx}.h5"
        write_mass_profile_compatible_h5(
            mass_profile_path,
            args,
            observations,
            dynamics,
            echo_keep,
            fitted_points_km,
            doppler_residual,
            velocity_keep,
            position_rms,
            path_length_km,
            stored_fit["label"],
        )
    if args.no_plot:
        print(
            f"single_pulses={len(result)} retained={np.count_nonzero(velocity_keep)} "
            f"doppler_rms_mps={rms(doppler_residual[velocity_keep]):.6f} "
            f"radius_um={1e6 * 10.0 ** dynamics['parameters'][6]:.6f}",
            flush=True,
        )
        return 0

    observation_rti = recompute_cut_observables(cut, interp=1)
    rti_tx_idx = np.asarray(observation_rti["tx_idx"], dtype=np.int64)
    measurement_tx_idx = np.asarray(decoded["tx_idx"], dtype=np.int64)
    rti_lookup = {int(value): row for row, value in enumerate(rti_tx_idx)}
    rti_rows = np.asarray([rti_lookup[int(value)] for value in measurement_tx_idx])
    native_ipp_samples = int(round(float(pmm.get_m_mode()["ipp_us"]) * 1e-6 * FS_HZ))
    uniform_tx_idx, rti_linear = uniform_ipp_rti(
        measurement_tx_idx,
        np.asarray(observation_rti["rti_snr"], dtype=float)[rti_rows],
        native_ipp_samples,
    )
    rti_db = 10.0 * np.log10(np.maximum(rti_linear, 1e-12))
    rti_db = np.nan_to_num(rti_db, nan=0.0, neginf=0.0, posinf=0.0)
    rti_plot_time = (
        trajectory_time[0] - time_origin + (uniform_tx_idx - measurement_tx_idx[0]) / FS_HZ
    )
    rti_range_km = np.asarray(observation_rti["range_grid_km"], dtype=float)
    snr_vmin = 0.0
    snr_vmax = max(
        18.0,
        float(
            np.nanpercentile(
                np.r_[snr_db[np.isfinite(snr_db)], rti_db[np.isfinite(rti_db)]], 99.7
            )
        ),
    )

    fig, axes = plt.subplots(2, 3, figsize=(13.0, 8.0), constrained_layout=True)
    fixed_colors = ["tab:purple", "tab:orange", "tab:red"]
    axis = axes[0, 0]
    beam_east_km, beam_north_km, beam_map = beam_pixmap()
    axis.pcolormesh(beam_east_km, beam_north_km, beam_map, cmap="gist_yarg", vmax=5.0, shading="auto")
    axis.scatter(fitted_points_km[~echo_keep, 0], fitted_points_km[~echo_keep, 1], color="0.75", s=5, edgecolors="none")
    scatter = axis.scatter(fitted_points_km[echo_keep, 0], fitted_points_km[echo_keep, 1], c=snr_db[echo_keep], cmap="plasma", vmin=snr_vmin, vmax=snr_vmax, s=5, edgecolors="none")
    axis.plot(refit_position[echo_keep, 0], refit_position[echo_keep, 1], color="C0", lw=1.0)
    axis.text(fitted_points_km[0, 0], fitted_points_km[0, 1], r"$t_0$", fontsize=7)
    axis.set(aspect="equal", xlabel="EW (km)", ylabel="NS (km)", xlim=(-35, 35), ylim=(-35, 35))
    axis.text(0.03, 0.97, f"Path {path_length_km:.1f} km\nEW RMS {position_rms[0] * 1e3:.0f} m\nNS RMS {position_rms[1] * 1e3:.0f} m", transform=axis.transAxes, va="top")
    color_axis = inset_axes(axis, width="4%", height="35%", loc="lower right", borderpad=2.0)
    colorbar = fig.colorbar(scatter, cax=color_axis)
    colorbar.set_label("SNR (dB)", fontsize=7)
    colorbar.ax.yaxis.set_label_position("left")

    axis = axes[0, 1]
    axis.fill_between(observation_time[echo_keep], dynamics["position_interval_km"][0, echo_keep, 2], dynamics["position_interval_km"][1, echo_keep, 2], color="C0", alpha=0.18, linewidth=0)
    axis.plot(observation_time[echo_keep], fitted_points_km[echo_keep, 2], ".", color="black")
    axis.plot(observation_time[echo_keep], refit_position[echo_keep, 2], color="C0")
    axis.set(xlabel="Time (s)", ylabel="Up (km)")
    axis.text(0.97, 0.97, f"Up RMS {position_rms[2] * 1e3:.0f} m", transform=axis.transAxes, ha="right", va="top")

    axis = axes[0, 2]
    axis.fill_between(observation_time[echo_keep], dynamics["doppler_interval_mps"][0, echo_keep] / 1e3, dynamics["doppler_interval_mps"][1, echo_keep] / 1e3, color="C0", alpha=0.18, linewidth=0)
    axis.errorbar(measurement_time[velocity_keep], result["velocity_mps"][velocity_keep] / 1e3, yerr=result["velocity_std_mps"][velocity_keep] / 1e3, fmt=".", color="black", zorder=2)
    for fixed_doppler, color in zip(dynamics["fixed_doppler_mps"], fixed_colors):
        axis.plot(observation_time[echo_keep], fixed_doppler[echo_keep] / 1e3, color=color, ls="--", lw=0.9, zorder=4)
    axis.plot(observation_time, refit_doppler / 1e3, color="C0", zorder=5)
    axis.set(xlabel="Time (s)", ylabel="Doppler (km/s)")
    axis.text(0.03, 0.97, f"RMS {rms(doppler_residual[velocity_keep]):.0f} m/s", transform=axis.transAxes, va="top")

    axis = axes[1, 0]
    rti_cmap = plt.get_cmap("plasma").copy()
    axis.pcolormesh(rti_plot_time, rti_range_km, rti_db.T, cmap=rti_cmap, shading="auto", vmin=snr_vmin, vmax=snr_vmax)
    axis.plot(observation_time[echo_keep], precise_range_km[echo_keep], ".", color="black")
    axis.fill_between(observation_time[echo_keep], dynamics["range_interval_km"][0, echo_keep], dynamics["range_interval_km"][1, echo_keep], color="C0", alpha=0.18, linewidth=0)
    axis.plot(observation_time, refit_range, color="C0")
    range_values = np.r_[refit_range, precise_range_km[echo_keep]]
    axis.set_ylim(max(0.0, float(np.nanmin(range_values)) - 4.0), float(np.nanmax(range_values)) + 4.0)
    axis.set(xlabel="Time (s)", ylabel="Range (km)")
    axis.text(0.03, 0.97, f"RMS {rms(range_residual) * 1e3:.0f} m", transform=axis.transAxes, va="top")

    axis = axes[1, 1]
    probability = dynamics["profile_probability_density_log_radius"].copy()
    probability /= max(float(np.nanmax(probability)), 1e-30)
    profile_radius_um = dynamics["profile_radius_um"]
    radius_quantiles_um = dynamics["marginal_radius_quantiles_um"]
    mass_quantiles_kg = dynamics["marginal_mass_quantiles_kg"]
    upper_constrained = bool(dynamics["radius_upper_limit_data_constrained"])
    axis.plot(profile_radius_um, probability, color="black")
    axis.fill_between(profile_radius_um, 0, probability, color="0.75", alpha=0.7)
    axis.axvline(1e6 * 10.0 ** dynamics["parameters"][6], color="C0")
    axis.axvline(radius_quantiles_um[0], color="C0", ls="--", lw=1)
    if upper_constrained:
        axis.axvspan(radius_quantiles_um[0], radius_quantiles_um[-1], color="C0", alpha=0.12)
        interval = rf"95% $r_0$ {radius_quantiles_um[0]:.0f}--{radius_quantiles_um[-1]:.0f} $\mu$m" + "\n" + rf"95% $m_0$ {mass_quantiles_kg[0]:.1e}--{mass_quantiles_kg[-1]:.1e} kg"
    else:
        interval = rf"95% lower $r_0>{radius_quantiles_um[0]:.0f}$ $\mu$m; no upper bound" + "\n" + rf"95% lower $m_0>{mass_quantiles_kg[0]:.1e}$ kg; no upper bound"
    axis.text(0.04, 0.94, interval, transform=axis.transAxes, va="top", fontsize=8)
    axis.set(xscale="log", xlim=(1, 1e4), ylim=(0, 1.05), xlabel=r"Initial radius $r_0$ ($\mu$m)", ylabel="Relative probability")
    secondary = axis.secondary_xaxis("top", functions=(radius_um_to_mass_kg, lambda mass: 1e6 * (3.0 * mass / (4.0 * np.pi * 3000.0)) ** (1.0 / 3.0)))
    secondary.set_xscale("log")
    secondary.set_xlabel(r"Initial mass $m_0$ (kg)")

    axis = axes[1, 2]
    axis.fill_between(observation_time[echo_keep], dynamics["speed_interval_km_s"][0, echo_keep], dynamics["speed_interval_km_s"][1, echo_keep], color="C0", alpha=0.18, linewidth=0)
    axis.plot(observation_time[echo_keep], refit_speed[echo_keep], color="C0", label="fit")
    speed_values = [refit_speed[echo_keep], dynamics["speed_interval_km_s"][:, echo_keep].ravel()]
    for radius, fixed_speed, color in zip(dynamics["fixed_radius_um"], dynamics["fixed_speed_km_s"], fixed_colors):
        axis.plot(observation_time[echo_keep], fixed_speed[echo_keep], color=color, ls="--", lw=0.9, label=rf"$r_0={radius:g}\,\mu$m")
    finite_speed = np.concatenate(speed_values)
    finite_speed = finite_speed[np.isfinite(finite_speed)]
    margin = max(0.2, 0.05 * float(np.ptp(finite_speed)))
    axis.set_ylim(float(np.min(finite_speed)) - margin, float(np.max(finite_speed)) + margin)
    axis.set(xlabel="Time (s)", ylabel="Speed (km/s)")
    axis.legend(frameon=False, loc="lower left", fontsize=7)

    for axis in axes.ravel():
        axis.grid(alpha=0.2, lw=0.5)
    fig.savefig(event_png, dpi=190)
    plt.close(fig)
    print(
        f"single_pulses={len(result)} retained={np.count_nonzero(velocity_keep)} "
        f"doppler_rms_mps={rms(doppler_residual[velocity_keep]):.6f} "
        f"radius_um={1e6 * 10.0 ** dynamics['parameters'][6]:.6f}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
