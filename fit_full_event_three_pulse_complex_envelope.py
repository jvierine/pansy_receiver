#!/usr/bin/env python3
"""Fit all valid three-pulse complex envelopes in one PANSY event."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

import pansy_config as pc
from aod_complex_doppler_fit import (
    ESTABLISHED_AOD_GEOMETRY_SIGN,
    ESTABLISHED_AOD_PHASECAL_SIGN,
    beamform_echo,
    event_module_voltage_gain,
    event_receiver_noise_variance,
)
import pansy_interferometry as pint
import pansy_modes as pmm
from fit_three_pulse_complex_envelope import (
    FS_HZ,
    PAIR_SPACING_S,
    PAIR_SPACING_TOLERANCE_S,
    fit_three_pulse_envelope,
    rms_baud_amplitudes,
)
from interferometer_alias_diagnostics import load_cut, recompute_cut_observables
from inter_pulse_phase_deceleration import (
    baud_averaged_beat_pairs,
    decoded_pulse_responses,
    fractional_segment,
    load_selected,
)
from run_catalogue_mass_profiles import load_selected_fit
from run_catalogue_mass_profiles import log_grid_probability
from run_catalogue_mass_profiles import radius_um_to_mass_kg
from run_catalogue_mass_profiles import weighted_quantile


def phase_acceleration_lookup(profile_h5: Path) -> dict[tuple[int, int, int, int], dict]:
    lookup = {}
    with h5py.File(profile_h5, "r") as handle:
        group = handle["phase_acceleration"]
        support = np.asarray(group["support_observation_indices"], dtype=int)
        model = np.asarray(group["best_model_radial_acceleration_mps2"], dtype=float)
        measured = np.asarray(group["measured_radial_acceleration_display_mps2"], dtype=float)
        keep = np.asarray(group["shared_inlier_mask"], dtype=bool)
        for index, row in enumerate(support):
            lookup[tuple(int(value) for value in row)] = {
                "model_mps2": float(model[index]),
                "measured_mps2": float(measured[index]),
                "keep": bool(keep[index]),
            }
    return lookup


def beam_pixmap() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mode = pmm.get_m_mode()
    u, v, w = pint.uv_coverage(N=450, max_zenith_angle=20.0)
    bitmap = np.zeros(u.shape)
    for beam in range(5):
        azimuth = mode["beam_pos_az_za"][beam][0]
        elevation = 90.0 - mode["beam_pos_az_za"][beam][1]
        horizontal = np.cos(np.deg2rad(elevation))
        vertical = -np.sin(np.deg2rad(elevation))
        east = -horizontal * np.sin(np.deg2rad(-azimuth))
        north = horizontal * np.cos(np.deg2rad(-azimuth))
        separation = np.rad2deg(
            np.arccos(np.clip(u * east + v * north + vertical * w, -1.0, 1.0))
        )
        bitmap[separation < 10.0] += 1.0
    return u * 100.0, v * 100.0, bitmap


def valid_triplets(pairs: dict) -> list[tuple[int, int, int]]:
    previous = np.asarray(pairs["previous_index"], dtype=int)
    current = np.asarray(pairs["current_index"], dtype=int)
    delta_t = np.asarray(pairs["delta_t_s"], dtype=float)
    valid = np.isfinite(delta_t) & (np.abs(delta_t - PAIR_SPACING_S) <= PAIR_SPACING_TOLERANCE_S)
    rows = []
    for first in np.flatnonzero(valid):
        second_rows = np.flatnonzero(valid & (previous == current[first]))
        for second in second_rows:
            rows.append((int(previous[first]), int(current[first]), int(current[second])))
    return rows


def pair_velocity_at_middle(pairs: dict, previous: int, middle: int, current: int) -> float:
    pair_previous = np.asarray(pairs["previous_index"], dtype=int)
    pair_current = np.asarray(pairs["current_index"], dtype=int)
    first = np.flatnonzero((pair_previous == previous) & (pair_current == middle))
    second = np.flatnonzero((pair_previous == middle) & (pair_current == current))
    if len(first) != 1 or len(second) != 1:
        return np.nan
    velocity = np.asarray(pairs["phase_doppler_mps"], dtype=float)
    return float(0.5 * (velocity[first[0]] + velocity[second[0]]))


def refit_dynamics(
    physics,
    trajectory_time,
    points_km,
    position_covariance_km2,
    echo_keep,
    result,
    initial_parameters,
    density,
    profile_radius_um,
    joint_correlation=None,
) -> dict:
    velocity_keep = np.asarray(result["shared_inlier"], dtype=bool)
    acceleration_keep = velocity_keep & ~np.asarray(result["acceleration_at_bound"], dtype=bool)
    velocity_sigma_mps = max(
        20.0,
        float(
            np.sqrt(
                np.mean(
                    (
                        result["velocity_mps"][velocity_keep]
                        - result["model_velocity_mps"][velocity_keep]
                    )
                    ** 2
                )
            )
        ),
    )
    acceleration_sigma_mps2 = max(
        1e3,
        float(
            np.sqrt(
                np.mean(
                    (
                        result["acceleration_mps2"][acceleration_keep]
                        - result["model_acceleration_mps2"][acceleration_keep]
                    )
                    ** 2
                )
            )
        ),
    )
    position_sigma_km = np.sqrt(np.maximum(np.diag(position_covariance_km2), 1e-8))
    measurement_count = (
        3 * np.count_nonzero(echo_keep)
        + np.count_nonzero(velocity_keep)
        + np.count_nonzero(acceleration_keep)
    )
    if joint_correlation is None:
        joint_correlation = np.eye(measurement_count)
    joint_correlation = np.asarray(joint_correlation, dtype=float)
    if joint_correlation.shape != (measurement_count, measurement_count):
        raise ValueError(
            f"joint correlation shape {joint_correlation.shape} does not match "
            f"the {measurement_count}-element fit vector"
        )

    def covariance_cholesky():
        measurement_scale = np.r_[
            np.tile(position_sigma_km, np.count_nonzero(echo_keep)),
            np.full(np.count_nonzero(velocity_keep), velocity_sigma_mps),
            np.full(np.count_nonzero(acceleration_keep), acceleration_sigma_mps2),
        ]
        covariance = joint_correlation * np.outer(measurement_scale, measurement_scale)
        covariance += np.eye(len(covariance)) * max(
            1e-12, 1e-10 * float(np.median(np.diag(covariance)))
        )
        return np.linalg.cholesky(covariance)

    joint_cholesky = covariance_cholesky()
    fit_time = np.asarray(result["time_s"], dtype=float)

    def predict(parameters):
        position, velocity, radius, mass, success, message = physics.propagate_shrinking_radius_model(
            parameters, trajectory_time, density
        )
        if not success or not np.all(np.isfinite(position)) or not np.all(np.isfinite(velocity)):
            raise RuntimeError(message)
        doppler_mps = physics.predicted_doppler(position, velocity) * 1e3
        acceleration_mps2 = np.gradient(doppler_mps, trajectory_time)
        return position, velocity, doppler_mps, acceleration_mps2, radius, mass

    def residual(parameters):
        try:
            position, _velocity, doppler_mps, acceleration_mps2, _radius, _mass = predict(parameters)
        except Exception:
            count = 3 * np.count_nonzero(echo_keep) + np.count_nonzero(velocity_keep) + np.count_nonzero(acceleration_keep)
            return np.full(count, 1e6)
        position_residual = (points_km[echo_keep] - position[echo_keep]).ravel()
        velocity_prediction = np.interp(fit_time, trajectory_time, doppler_mps)
        acceleration_prediction = np.interp(fit_time, trajectory_time, acceleration_mps2)
        raw_residual = np.r_[
            position_residual,
            result["velocity_mps"][velocity_keep] - velocity_prediction[velocity_keep],
            result["acceleration_mps2"][acceleration_keep]
            - acceleration_prediction[acceleration_keep],
        ]
        return solve_triangular(joint_cholesky, raw_residual, lower=True, check_finite=False)

    lower = np.asarray([-np.inf, -np.inf, 20e3, -90e3, -90e3, -90e3, -6.0])
    upper = np.asarray([np.inf, np.inf, 220e3, 90e3, 90e3, 90e3, -2.0])
    scale = np.asarray([1e5, 1e5, 1e5, 7e4, 7e4, 7e4, 1.0])
    robust = least_squares(
        residual,
        np.clip(initial_parameters, lower + 1e-10, upper - 1e-10),
        bounds=(lower, upper),
        x_scale=scale,
        loss="soft_l1",
        max_nfev=160,
    )
    final = least_squares(
        residual,
        robust.x,
        bounds=(lower, upper),
        x_scale=scale,
        loss="linear",
        max_nfev=200,
    )
    preliminary = predict(final.x)
    preliminary_velocity = np.interp(fit_time, trajectory_time, preliminary[2])
    preliminary_acceleration = np.interp(fit_time, trajectory_time, preliminary[3])
    position_sigma_km = np.maximum(
        np.sqrt(np.mean((points_km[echo_keep] - preliminary[0][echo_keep]) ** 2, axis=0)),
        0.01,
    )
    velocity_sigma_mps = max(
        20.0,
        float(
            np.sqrt(
                np.mean(
                    (result["velocity_mps"][velocity_keep] - preliminary_velocity[velocity_keep])
                    ** 2
                )
            )
        ),
    )
    acceleration_sigma_mps2 = max(
        1e3,
        float(
            np.sqrt(
                np.mean(
                    (
                        result["acceleration_mps2"][acceleration_keep]
                        - preliminary_acceleration[acceleration_keep]
                    )
                    ** 2
                )
            )
        ),
    )
    joint_cholesky = covariance_cholesky()
    final = least_squares(
        residual,
        final.x,
        bounds=(lower, upper),
        x_scale=scale,
        loss="linear",
        max_nfev=240,
    )
    covariance_degrees_of_freedom = max(measurement_count - len(final.x), 1)
    covariance_variance_inflation = max(
        1.0, float(np.sum(residual(final.x) ** 2) / covariance_degrees_of_freedom)
    )
    joint_cholesky *= np.sqrt(covariance_variance_inflation)
    final = least_squares(
        residual,
        final.x,
        bounds=(lower, upper),
        x_scale=scale,
        loss="linear",
        max_nfev=160,
    )
    position, velocity, doppler_mps, acceleration_mps2, radius, mass = predict(final.x)

    profile_radius_um = np.unique(
        np.r_[np.asarray(profile_radius_um, dtype=float), 10.0, 100.0, 1000.0]
    )
    profile_log_radius_m = np.log10(profile_radius_um * 1e-6)
    profile_parameters6 = np.full((len(profile_radius_um), 6), np.nan)
    profile_chi2 = np.full(len(profile_radius_um), np.nan)
    profile_success = np.zeros(len(profile_radius_um), dtype=bool)
    lower6 = lower[:6]
    upper6 = upper[:6]
    scale6 = scale[:6]
    center = int(np.argmin(np.abs(profile_log_radius_m - final.x[6])))

    def fit_profile_index(index, seed):
        fixed_log_radius = profile_log_radius_m[index]

        def fixed_residual(parameters6):
            return residual(np.r_[parameters6, fixed_log_radius])

        fitted = least_squares(
            fixed_residual,
            np.clip(seed, lower6 + 1e-10, upper6 - 1e-10),
            bounds=(lower6, upper6),
            x_scale=scale6,
            loss="linear",
            max_nfev=180,
        )
        profile_parameters6[index] = fitted.x
        profile_chi2[index] = float(np.sum(fixed_residual(fitted.x) ** 2))
        profile_success[index] = bool(fitted.success)
        return fitted.x

    center_seed = fit_profile_index(center, final.x[:6])
    seed = center_seed
    for index in range(center - 1, -1, -1):
        seed = fit_profile_index(index, seed)
    seed = center_seed
    for index in range(center + 1, len(profile_radius_um)):
        seed = fit_profile_index(index, seed)

    minimum_chi2 = min(float(np.sum(residual(final.x) ** 2)), float(np.nanmin(profile_chi2)))
    profile_delta_chi2 = profile_chi2 - minimum_chi2
    profile_probability, profile_weights = log_grid_probability(
        profile_radius_um, profile_delta_chi2
    )
    radius_quantiles_um = weighted_quantile(
        profile_radius_um, profile_weights, [0.025, 0.5, 0.975]
    )

    profile_position = []
    profile_doppler_mps = []
    profile_acceleration_mps2 = []
    profile_speed_km_s = []
    model_weights = []
    for index in np.flatnonzero(profile_success & np.isfinite(profile_weights)):
        try:
            predicted = predict(
                np.r_[profile_parameters6[index], profile_log_radius_m[index]]
            )
        except Exception:
            continue
        profile_position.append(predicted[0])
        profile_doppler_mps.append(predicted[2])
        profile_acceleration_mps2.append(predicted[3])
        profile_speed_km_s.append(np.linalg.norm(predicted[1], axis=1))
        model_weights.append(profile_weights[index])

    def pointwise_interval(values, weights):
        values = np.asarray(values, dtype=float)
        weights = np.asarray(weights, dtype=float)
        flattened = values.reshape(values.shape[0], -1)
        interval = np.asarray(
            [weighted_quantile(flattened[:, index], weights, [0.025, 0.975]) for index in range(flattened.shape[1])]
        ).T
        return interval.reshape((2,) + values.shape[1:])

    model_weights = np.asarray(model_weights, dtype=float)
    model_weights /= np.sum(model_weights)
    position_interval_km = pointwise_interval(profile_position, model_weights)
    doppler_interval_mps = pointwise_interval(profile_doppler_mps, model_weights)
    acceleration_interval_mps2 = pointwise_interval(profile_acceleration_mps2, model_weights)
    speed_interval_km_s = pointwise_interval(profile_speed_km_s, model_weights)
    range_interval_km = pointwise_interval(
        np.linalg.norm(np.asarray(profile_position), axis=2), model_weights
    )

    fixed_radius_um = np.asarray([10.0, 100.0, 1000.0])
    fixed_doppler_mps = []
    fixed_acceleration_mps2 = []
    fixed_speed_km_s = []
    for target_um in fixed_radius_um:
        index = int(np.flatnonzero(profile_radius_um == target_um)[0])
        predicted = predict(np.r_[profile_parameters6[index], profile_log_radius_m[index]])
        fixed_doppler_mps.append(predicted[2])
        fixed_acceleration_mps2.append(predicted[3])
        fixed_speed_km_s.append(np.linalg.norm(predicted[1], axis=1))

    return {
        "parameters": final.x,
        "position_km": position,
        "velocity_km_s": velocity,
        "doppler_mps": doppler_mps,
        "acceleration_mps2": acceleration_mps2,
        "radius_m": radius,
        "mass_kg": mass,
        "velocity_sigma_mps": velocity_sigma_mps,
        "acceleration_sigma_mps2": acceleration_sigma_mps2,
        "position_sigma_km": position_sigma_km,
        "joint_correlation": joint_correlation,
        "covariance_degrees_of_freedom": covariance_degrees_of_freedom,
        "covariance_variance_inflation": covariance_variance_inflation,
        "cost": float(final.cost),
        "success": bool(final.success),
        "profile_radius_um": profile_radius_um,
        "profile_parameters6": profile_parameters6,
        "profile_chi2": profile_chi2,
        "profile_delta_chi2": profile_delta_chi2,
        "profile_probability_density_log_radius": profile_probability,
        "profile_probability_weights": profile_weights,
        "profile_success": profile_success,
        "marginal_radius_quantiles_um": radius_quantiles_um,
        "marginal_mass_quantiles_kg": radius_um_to_mass_kg(radius_quantiles_um),
        "position_interval_km": position_interval_km,
        "doppler_interval_mps": doppler_interval_mps,
        "acceleration_interval_mps2": acceleration_interval_mps2,
        "speed_interval_km_s": speed_interval_km_s,
        "range_interval_km": range_interval_km,
        "fixed_radius_um": fixed_radius_um,
        "fixed_doppler_mps": np.asarray(fixed_doppler_mps),
        "fixed_acceleration_mps2": np.asarray(fixed_acceleration_mps2),
        "fixed_speed_km_s": np.asarray(fixed_speed_km_s),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-idx", type=int, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--diagnostics-h5", type=Path, required=True)
    parser.add_argument("--initial-fit-h5", type=Path, required=True)
    parser.add_argument("--prior-profile-h5", type=Path, required=True)
    parser.add_argument("--joint-covariance-h5", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    args = parser.parse_args()

    cut = load_cut(args.base / "metadata/cut", args.sample_idx)
    hypothesis = load_selected(args.diagnostics_h5)
    with h5py.File(args.initial_fit_h5, "r") as handle:
        precise_range_km = np.asarray(handle["range_km"], dtype=float)
        precise_doppler_mps = np.asarray(handle["doppler_mps"], dtype=float)
    decoded = decoded_pulse_responses(
        cut,
        hypothesis,
        args.snr_threshold,
        precise_range_km=precise_range_km,
        precise_doppler_mps=precise_doppler_mps,
    )
    pairs = baud_averaged_beat_pairs(decoded, same_beam=True)
    triplets = valid_triplets(pairs)
    prior_acceleration = phase_acceleration_lookup(args.prior_profile_h5)
    with h5py.File(args.prior_profile_h5, "r") as handle:
        echo_keep = np.asarray(handle["result/echo_shared_inlier_mask"], dtype=bool)

    observations, stored_fit = load_selected_fit(args.diagnostics_h5)
    trajectory_time = np.asarray(observations["t_s"], dtype=float)
    steering_points_km = np.asarray(stored_fit["model_km"], dtype=float)
    steering_uvw = steering_points_km / np.linalg.norm(
        steering_points_km, axis=1, keepdims=True
    )
    module_noise_variance = event_receiver_noise_variance(
        decoded["z_rx"],
        decoded["raw_idx"],
        decoded["range_gate"],
        decoded["z_tx"].shape[1],
    )
    module_voltage_gain = event_module_voltage_gain(decoded["response"])
    steering_signs = (
        ESTABLISHED_AOD_PHASECAL_SIGN,
        ESTABLISHED_AOD_GEOMETRY_SIGN,
    )
    beamformed_echo = np.stack(
        [
            beamform_echo(
                decoded["z_rx"][raw_index],
                steering_uvw[observation],
                int(decoded["beam_id"][observation]),
                steering_signs,
                module_noise_variance[raw_index],
                module_voltage_gain,
            )
            for observation, raw_index in enumerate(decoded["raw_idx"])
        ]
    )
    absolute_zero = decoded["tx_idx"][0] / FS_HZ - trajectory_time[0]
    with h5py.File(args.prior_profile_h5, "r") as handle:
        delta = np.asarray(handle["profile/delta_chi2"], dtype=float)
        best_index = int(np.nanargmin(delta))
        parameters6 = np.asarray(handle["profile/parameters6"][best_index], dtype=float)
        radius_um = float(handle["profile/radius_um"][best_index])
        profile_grid_radius_um = np.asarray(handle["profile/radius_um"], dtype=float)
    import fit_best_alias_physics_models as physics

    density, _metadata = physics.pbal.density_interpolator(observations["sample_epoch_unix"])
    position, velocity, _radius, _mass, success, message = physics.propagate_shrinking_radius_model(
        np.r_[parameters6, np.log10(radius_um * 1e-6)], trajectory_time, density
    )
    if not success:
        raise RuntimeError(f"prior physical trajectory failed: {message}")
    model_doppler_mps = physics.predicted_doppler(position, velocity) * 1e3
    model_acceleration_mps2 = np.gradient(model_doppler_mps, trajectory_time)

    dtype = [
        ("previous", "i4"),
        ("middle", "i4"),
        ("current", "i4"),
        ("channel", "i4"),
        ("time_s", "f8"),
        ("velocity_mps", "f8"),
        ("velocity_std_mps", "f8"),
        ("fft_velocity_mps", "f8"),
        ("beat_velocity_mps", "f8"),
        ("model_velocity_mps", "f8"),
        ("acceleration_mps2", "f8"),
        ("acceleration_std_mps2", "f8"),
        ("prior_measured_acceleration_mps2", "f8"),
        ("model_acceleration_mps2", "f8"),
        ("frequency_acceleration_correlation", "f8"),
        ("shared_inlier", "?"),
        ("acceleration_at_bound", "?"),
    ]
    rows = []
    n_fast = decoded["z_tx"].shape[1]
    for previous, middle, current in triplets:
        observation_indices = np.asarray([previous, middle, current], dtype=int)
        raw_indices = np.asarray(decoded["raw_idx"][observation_indices], dtype=int)
        channel = -1
        raw_pulses = np.stack(
            [
                fractional_segment(
                    beamformed_echo[observation],
                    decoded["range_gate"][observation],
                    n_fast,
                )
                for raw_index, observation in zip(raw_indices, observation_indices)
            ]
        )
        templates = decoded["z_tx"][raw_indices].astype(np.complex128)
        amplitude_prior = rms_baud_amplitudes(raw_pulses, templates)
        support = (previous, middle, middle, current)
        local_prior = prior_acceleration.get(support)
        middle_time = decoded["tx_idx"][middle] / FS_HZ - absolute_zero
        acceleration_guess = (
            local_prior["model_mps2"]
            if local_prior is not None
            else float(np.interp(middle_time, trajectory_time, model_acceleration_mps2))
        )
        velocity_guess = float(np.interp(middle_time, trajectory_time, model_doppler_mps))
        fit = fit_three_pulse_envelope(
            raw_pulses,
            templates,
            PAIR_SPACING_S,
            2.0 * velocity_guess / pc.wavelength,
            acceleration_guess,
            pc.wavelength,
            acceleration_half_width_mps2=2.0e4,
            frequency_half_width_hz=1.0e3,
            pulse_snr=np.asarray(decoded["snr"][observation_indices], dtype=float),
            matched_filter_amplitudes=amplitude_prior,
        )
        covariance = np.asarray(fit["parameter_covariance"], dtype=float)
        fit_time = decoded["tx_idx"][previous] / FS_HZ - absolute_zero + fit["time_center_s"]
        model_velocity = float(np.interp(fit_time, trajectory_time, model_doppler_mps))
        model_acceleration = (
            local_prior["model_mps2"]
            if local_prior is not None
            else float(np.interp(fit_time, trajectory_time, model_acceleration_mps2))
        )
        rows.append(
            (
                previous,
                middle,
                current,
                channel,
                fit_time,
                fit["parameters"][4] * pc.wavelength / 2.0,
                np.sqrt(covariance[4, 4]) * pc.wavelength / 2.0,
                decoded["coarse_doppler_mps"][middle],
                pair_velocity_at_middle(pairs, previous, middle, current),
                model_velocity,
                fit["parameters"][5],
                np.sqrt(covariance[5, 5]),
                local_prior["measured_mps2"] if local_prior is not None else np.nan,
                model_acceleration,
                covariance[4, 5] / np.sqrt(covariance[4, 4] * covariance[5, 5]),
                bool(np.all(echo_keep[observation_indices]))
                and (local_prior is None or local_prior["keep"]),
                bool(abs(fit["parameters"][5] - acceleration_guess) > 1.99e4),
            )
        )
    result = np.asarray(rows, dtype=dtype)
    if len(result) == 0:
        raise RuntimeError("no valid three-pulse fits")

    fitted_points_km = np.asarray(observations["points_km"], dtype=float).copy()
    catalogue_range_km = np.linalg.norm(fitted_points_km, axis=1)
    valid_range = (
        np.isfinite(precise_range_km)
        & np.isfinite(catalogue_range_km)
        & (catalogue_range_km > 0.0)
    )
    fitted_points_km[valid_range] *= (
        precise_range_km[valid_range] / catalogue_range_km[valid_range]
    )[:, None]
    with h5py.File(args.prior_profile_h5, "r") as handle:
        position_covariance_km2 = np.asarray(handle["position_residual_covariance_km2"], dtype=float)
    joint_correlation = None
    if args.joint_covariance_h5 is not None:
        with h5py.File(args.joint_covariance_h5, "r") as handle:
            covariance_echo_keep = np.asarray(handle["echo_keep"], dtype=bool)
            covariance_velocity_keep = np.asarray(handle["velocity_keep"], dtype=bool)
            covariance_acceleration_keep = np.asarray(handle["acceleration_keep"], dtype=bool)
            joint_correlation = np.asarray(handle["correlation"], dtype=float)
        expected_velocity_keep = np.asarray(result["shared_inlier"], dtype=bool)
        expected_acceleration_keep = expected_velocity_keep & ~np.asarray(
            result["acceleration_at_bound"], dtype=bool
        )
        if not np.array_equal(covariance_echo_keep, echo_keep):
            raise RuntimeError("joint covariance echo mask does not match this fit")
        if not np.array_equal(covariance_velocity_keep, expected_velocity_keep):
            raise RuntimeError("joint covariance velocity mask does not match this fit")
        if not np.array_equal(covariance_acceleration_keep, expected_acceleration_keep):
            raise RuntimeError("joint covariance acceleration mask does not match this fit")
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
        joint_correlation=joint_correlation,
    )

    output_h5 = args.output_dir / f"three_pulse_full_event_{args.sample_idx}.h5"
    output_png = args.output_dir / f"three_pulse_full_event_{args.sample_idx}.png"
    event_png = args.output_dir / f"three_pulse_full_event_plot_{args.sample_idx}.png"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_h5, "w") as handle:
        handle.attrs["schema"] = "pansy.three_pulse_full_event.v1"
        handle.attrs["sample_idx"] = args.sample_idx
        handle.attrs["voltage_source"] = (
            "AoA-steered inverse-noise coherent sum of seven receiver modules"
        )
        handle.attrs["steering_direction_source"] = (
            "selected first-stage trajectory model"
        )
        handle.attrs["receiver_phase_calibration"] = "stored calibration applied"
        handle.attrs["static_amplitude_calibration"] = "mesocal amp_scale applied"
        handle.attrs["event_amplitude_calibration"] = (
            "robust per-module voltage gain from matched event responses"
        )
        handle.attrs["phasecal_sign"] = ESTABLISHED_AOD_PHASECAL_SIGN
        handle.attrs["geometry_sign"] = ESTABLISHED_AOD_GEOMETRY_SIGN
        handle.create_dataset("triplet_fit", data=result)
        handle.create_dataset("event_module_voltage_gain", data=module_voltage_gain)
        fit_group = handle.create_group("dynamics_refit")
        for name, value in dynamics.items():
            if isinstance(value, np.ndarray):
                fit_group.create_dataset(name, data=value)
            else:
                fit_group.attrs[name] = value

    time = result["time_s"] - trajectory_time[0]
    fft_residual = result["fft_velocity_mps"] - result["model_velocity_mps"]
    beat_residual = result["beat_velocity_mps"] - result["model_velocity_mps"]
    fit_residual = result["velocity_mps"] - result["model_velocity_mps"]
    old_acceleration_residual = (
        result["prior_measured_acceleration_mps2"] - result["model_acceleration_mps2"]
    )
    fit_acceleration_residual = result["acceleration_mps2"] - result["model_acceleration_mps2"]
    rms = lambda values: float(np.sqrt(np.nanmean(np.asarray(values) ** 2)))
    compare = np.asarray(result["shared_inlier"], dtype=bool)
    acceleration_compare = compare & ~np.asarray(result["acceleration_at_bound"], dtype=bool)

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.0), sharex=True, constrained_layout=True)
    axes[0].plot(time[compare], fft_residual[compare], ".", color="0.65", label=f"FFT Doppler, RMS {rms(fft_residual[compare]):.0f} m/s")
    axes[0].plot(time[compare], beat_residual[compare], ".", color="C1", label=f"beat phase, RMS {rms(beat_residual[compare]):.0f} m/s")
    axes[0].errorbar(
        time[compare],
        fit_residual[compare],
        yerr=result["velocity_std_mps"][compare],
        fmt=".",
        color="C0",
        ecolor="C0",
        alpha=0.8,
        label=f"three-pulse complex fit, RMS {rms(fit_residual[compare]):.0f} m/s",
    )
    axes[0].axhline(0.0, color="black", lw=0.8)
    axes[0].set_ylabel("Velocity residual (m/s)")
    axes[0].legend(frameon=False)

    axes[1].plot(
        time[acceleration_compare],
        old_acceleration_residual[acceleration_compare] / 1e3,
        ".",
        color="C1",
        label=f"existing phase acceleration, RMS {rms(old_acceleration_residual[acceleration_compare]) / 1e3:.2f} km/s2",
    )
    axes[1].errorbar(
        time[acceleration_compare],
        fit_acceleration_residual[acceleration_compare] / 1e3,
        yerr=result["acceleration_std_mps2"][acceleration_compare] / 1e3,
        fmt=".",
        color="C0",
        ecolor="C0",
        alpha=0.8,
        label=f"three-pulse complex fit, RMS {rms(fit_acceleration_residual[acceleration_compare]) / 1e3:.2f} km/s2",
    )
    at_bound = compare & np.asarray(result["acceleration_at_bound"], dtype=bool)
    axes[1].plot(
        time[at_bound],
        fit_acceleration_residual[at_bound] / 1e3,
        "x",
        color="C3",
        label=f"complex fit at acceleration bound ({np.count_nonzero(at_bound)})",
    )
    axes[1].axhline(0.0, color="black", lw=0.8)
    axes[1].set_ylabel("Acceleration residual (km/s2)")
    axes[1].set_xlabel("Time (s)")
    axes[1].legend(frameon=False)
    for axis in axes:
        axis.grid(alpha=0.2, lw=0.5)
    fig.savefig(output_png, dpi=190)
    plt.close(fig)

    refit_position = dynamics["position_km"]
    refit_velocity = dynamics["velocity_km_s"]
    refit_doppler = dynamics["doppler_mps"]
    refit_acceleration = dynamics["acceleration_mps2"]
    time_origin = float(trajectory_time[echo_keep][0])
    observation_time = trajectory_time - time_origin
    triplet_time = result["time_s"] - time_origin
    refit_range = np.linalg.norm(refit_position, axis=1)
    refit_speed = np.linalg.norm(refit_velocity, axis=1)
    snr_db = 10.0 * np.log10(np.maximum(np.asarray(observations["snr"], dtype=float), 1e-12))
    position_residual = fitted_points_km - refit_position
    position_rms = np.sqrt(np.mean(position_residual[echo_keep] ** 2, axis=0))
    range_residual = precise_range_km[echo_keep] - refit_range[echo_keep]
    velocity_refit_at_triplet = np.interp(result["time_s"], trajectory_time, refit_doppler)
    acceleration_refit_at_triplet = np.interp(
        result["time_s"], trajectory_time, refit_acceleration
    )
    velocity_residual_refit = result["velocity_mps"] - velocity_refit_at_triplet
    acceleration_residual_refit = result["acceleration_mps2"] - acceleration_refit_at_triplet
    fit_velocity_keep = np.asarray(result["shared_inlier"], dtype=bool)
    fit_acceleration_keep = fit_velocity_keep & ~np.asarray(result["acceleration_at_bound"], dtype=bool)
    path_length_km = float(np.sum(np.linalg.norm(np.diff(refit_position[echo_keep], axis=0), axis=1)))

    profile_radius_um = dynamics["profile_radius_um"]
    profile_probability = dynamics["profile_probability_density_log_radius"]
    radius_quantiles_um = dynamics["marginal_radius_quantiles_um"]
    mass_quantiles_kg = dynamics["marginal_mass_quantiles_kg"]
    profile_probability /= max(float(np.nanmax(profile_probability)), 1e-30)
    observation_rti = recompute_cut_observables(cut, interp=1)
    rti_time_absolute = np.asarray(observation_rti["tx_idx"], dtype=float) / FS_HZ
    measurement_time_absolute = np.asarray(decoded["tx_idx"], dtype=float) / FS_HZ
    rti_rows = np.asarray(
        [np.argmin(np.abs(rti_time_absolute - value)) for value in measurement_time_absolute],
        dtype=int,
    )
    rti_db = 10.0 * np.log10(
        np.maximum(np.asarray(observation_rti["rti_snr"], dtype=float)[rti_rows], 1e-12)
    )
    rti_range_km = np.asarray(observation_rti["range_grid_km"], dtype=float)
    snr_vmin = 0.0
    snr_vmax = max(
        18.0,
        float(np.nanpercentile(np.r_[snr_db[np.isfinite(snr_db)], rti_db[np.isfinite(rti_db)]], 99.7)),
    )

    event_fig, event_axes = plt.subplots(2, 4, figsize=(16.0, 8.0), constrained_layout=True)
    fixed_colors = ["tab:purple", "tab:orange", "tab:red"]
    axis = event_axes[0, 0]
    beam_east_km, beam_north_km, beam_map = beam_pixmap()
    axis.pcolormesh(
        beam_east_km,
        beam_north_km,
        beam_map,
        cmap="gist_yarg",
        vmax=5.0,
        shading="auto",
    )
    if np.any(~echo_keep):
        axis.scatter(
            fitted_points_km[~echo_keep, 0],
            fitted_points_km[~echo_keep, 1],
            color="0.75",
            s=5,
            edgecolors="none",
            zorder=2,
        )
    scatter = axis.scatter(
        fitted_points_km[echo_keep, 0],
        fitted_points_km[echo_keep, 1],
        c=snr_db[echo_keep],
        cmap="plasma",
        vmin=snr_vmin,
        vmax=snr_vmax,
        s=5,
        edgecolors="none",
        zorder=3,
    )
    axis.plot(
        refit_position[echo_keep, 0],
        refit_position[echo_keep, 1],
        color="C0",
        lw=1.0,
        zorder=4,
    )
    axis.text(fitted_points_km[0, 0], fitted_points_km[0, 1], r"$t_0$", fontsize=7)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("EW (km)")
    axis.set_ylabel("NS (km)")
    axis.set_xlim(-35.0, 35.0)
    axis.set_ylim(-35.0, 35.0)
    axis.text(
        0.03,
        0.97,
        f"Path {path_length_km:.1f} km\nEW RMS {position_rms[0]:.2f} km\nNS RMS {position_rms[1]:.2f} km",
        transform=axis.transAxes,
        va="top",
    )
    color_axis = inset_axes(axis, width="4%", height="35%", loc="lower right", borderpad=2.0)
    colorbar = event_fig.colorbar(scatter, cax=color_axis)
    colorbar.set_label("SNR (dB)", fontsize=7)
    colorbar.ax.yaxis.set_label_position("left")
    colorbar.ax.yaxis.set_ticks_position("right")
    colorbar.ax.tick_params(labelsize=7, length=2)

    axis = event_axes[0, 1]
    axis.fill_between(
        observation_time[echo_keep],
        dynamics["position_interval_km"][0, echo_keep, 2],
        dynamics["position_interval_km"][1, echo_keep, 2],
        color="C0",
        alpha=0.18,
        linewidth=0.0,
    )
    axis.plot(observation_time[echo_keep], fitted_points_km[echo_keep, 2], ".", color="black")
    axis.plot(observation_time[echo_keep], refit_position[echo_keep, 2], color="C0")
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Up (km)")
    axis.text(0.97, 0.97, f"Up RMS {position_rms[2]:.2f} km", transform=axis.transAxes, ha="right", va="top")

    axis = event_axes[0, 2]
    axis.fill_between(
        observation_time[echo_keep],
        dynamics["doppler_interval_mps"][0, echo_keep] / 1e3,
        dynamics["doppler_interval_mps"][1, echo_keep] / 1e3,
        color="C0",
        alpha=0.18,
        linewidth=0.0,
    )
    axis.plot(observation_time[echo_keep], precise_doppler_mps[echo_keep] / 1e3, ".", color="0.7", label="FFT")
    axis.errorbar(
        triplet_time[fit_velocity_keep],
        result["velocity_mps"][fit_velocity_keep] / 1e3,
        yerr=result["velocity_std_mps"][fit_velocity_keep] / 1e3,
        fmt=".",
        color="black",
        label="three-pulse complex",
    )
    axis.plot(observation_time, refit_doppler / 1e3, color="C0", label="refit")
    for target_um, fixed_doppler, color in zip(
        dynamics["fixed_radius_um"], dynamics["fixed_doppler_mps"], fixed_colors
    ):
        axis.plot(
            observation_time[echo_keep],
            fixed_doppler[echo_keep] / 1e3,
            color=color,
            linestyle="--",
            linewidth=0.9,
            label=rf"$r_0={target_um:g}\,\mu$m",
        )
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Doppler (km/s)")
    measured_doppler_km_s = result["velocity_mps"][fit_velocity_keep] / 1e3
    measured_doppler_km_s = measured_doppler_km_s[np.isfinite(measured_doppler_km_s)]
    if len(measured_doppler_km_s):
        doppler_span = float(np.ptp(measured_doppler_km_s))
        doppler_padding = max(0.5, 0.08 * doppler_span)
        axis.set_ylim(
            float(np.min(measured_doppler_km_s) - doppler_padding),
            float(np.max(measured_doppler_km_s) + doppler_padding),
        )
    axis.text(
        0.03,
        0.97,
        f"RMS {rms(velocity_residual_refit[fit_velocity_keep]):.0f} m/s",
        transform=axis.transAxes,
        va="top",
    )

    axis = event_axes[0, 3]
    axis.fill_between(
        observation_time[echo_keep],
        dynamics["acceleration_interval_mps2"][0, echo_keep] / 1e3,
        dynamics["acceleration_interval_mps2"][1, echo_keep] / 1e3,
        color="C0",
        alpha=0.18,
        linewidth=0.0,
    )
    axis.errorbar(
        triplet_time[fit_acceleration_keep],
        result["acceleration_mps2"][fit_acceleration_keep] / 1e3,
        yerr=result["acceleration_std_mps2"][fit_acceleration_keep] / 1e3,
        fmt=".",
        color="black",
        label="three-pulse complex",
    )
    axis.plot(observation_time, refit_acceleration / 1e3, color="C0", label="refit")
    for target_um, fixed_acceleration, color in zip(
        dynamics["fixed_radius_um"], dynamics["fixed_acceleration_mps2"], fixed_colors
    ):
        axis.plot(
            observation_time[echo_keep],
            fixed_acceleration[echo_keep] / 1e3,
            color=color,
            linestyle="--",
            linewidth=0.9,
            label=rf"$r_0={target_um:g}\,\mu$m",
        )
    axis.set_xlabel("Time (s)")
    axis.set_ylabel(r"Radial acceleration (km s$^{-2}$)")
    axis.text(
        0.03,
        0.97,
        f"RMS {rms(acceleration_residual_refit[fit_acceleration_keep]) / 1e3:.2f} km/s2",
        transform=axis.transAxes,
        va="top",
    )

    axis = event_axes[1, 0]
    axis.pcolormesh(
        observation_time,
        rti_range_km,
        rti_db.T,
        cmap="plasma",
        shading="auto",
        vmin=snr_vmin,
        vmax=snr_vmax,
    )
    axis.plot(observation_time[echo_keep], precise_range_km[echo_keep], ".", color="black", label="measurement")
    axis.fill_between(
        observation_time[echo_keep],
        dynamics["range_interval_km"][0, echo_keep],
        dynamics["range_interval_km"][1, echo_keep],
        color="C0",
        alpha=0.18,
        linewidth=0.0,
    )
    axis.plot(observation_time, refit_range, color="C0", label="refit")
    range_panel_values = np.concatenate(
        (refit_range[np.isfinite(refit_range)], precise_range_km[echo_keep])
    )
    axis.set_ylim(
        max(0.0, float(np.nanmin(range_panel_values)) - 4.0),
        float(np.nanmax(range_panel_values)) + 4.0,
    )
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Range (km)")
    axis.text(0.03, 0.97, f"RMS {rms(range_residual):.2f} km", transform=axis.transAxes, va="top")

    axis = event_axes[1, 1]
    axis.plot(profile_radius_um, profile_probability, color="black")
    axis.fill_between(profile_radius_um, 0.0, profile_probability, color="0.75", alpha=0.7)
    axis.axvline(1e6 * 10.0 ** dynamics["parameters"][6], color="C0", label="new best fit")
    axis.axvspan(radius_quantiles_um[0], radius_quantiles_um[-1], color="C0", alpha=0.12)
    axis.text(
        0.04,
        0.94,
        rf"95% $r_0$ {radius_quantiles_um[0]:.0f}--{radius_quantiles_um[-1]:.0f} $\mu$m"
        + "\n"
        + rf"95% $m_0$ {mass_quantiles_kg[0]:.1e}--{mass_quantiles_kg[-1]:.1e} kg",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    axis.set_xscale("log")
    axis.set_xlim(1.0, 1e4)
    axis.set_ylim(0.0, 1.05)
    axis.set_xlabel(r"Initial radius $r_0$ ($\mu$m)")
    axis.set_ylabel("Relative probability")
    secondary = axis.secondary_xaxis(
        "top",
        functions=(radius_um_to_mass_kg, lambda mass: 1e6 * (3.0 * mass / (4.0 * np.pi * 3000.0)) ** (1.0 / 3.0)),
    )
    secondary.set_xscale("log")
    secondary.set_xlabel(r"Initial mass $m_0$ (kg)")

    axis = event_axes[1, 2]
    axis.fill_between(
        observation_time[echo_keep],
        dynamics["speed_interval_km_s"][0, echo_keep],
        dynamics["speed_interval_km_s"][1, echo_keep],
        color="C0",
        alpha=0.18,
        linewidth=0.0,
    )
    axis.plot(observation_time[echo_keep], refit_speed[echo_keep], color="C0", label="fit")
    speed_limits = [
        refit_speed[echo_keep],
        dynamics["speed_interval_km_s"][:, echo_keep].ravel(),
    ]
    for target_um, fixed_speed, color in zip(
        dynamics["fixed_radius_um"], dynamics["fixed_speed_km_s"], fixed_colors
    ):
        axis.plot(
            observation_time[echo_keep],
            fixed_speed[echo_keep],
            color=color,
            linestyle="--",
            linewidth=0.9,
            label=rf"$r_0={target_um:g}\,\mu$m",
        )
    finite_speed = np.concatenate(speed_limits)
    finite_speed = finite_speed[np.isfinite(finite_speed)]
    speed_margin = max(0.2, 0.05 * float(np.ptp(finite_speed)))
    axis.set_ylim(float(np.min(finite_speed)) - speed_margin, float(np.max(finite_speed)) + speed_margin)
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Speed (km/s)")
    axis.legend(frameon=False, loc="lower left", fontsize=7)

    axis = event_axes[1, 3]
    doppler_artist = axis.errorbar(
        triplet_time[fit_velocity_keep],
        velocity_residual_refit[fit_velocity_keep],
        yerr=result["velocity_std_mps"][fit_velocity_keep],
        fmt=".",
        color="C0",
        label="Doppler",
    )
    acceleration_axis = axis.twinx()
    acceleration_artist = acceleration_axis.plot(
        triplet_time[fit_acceleration_keep],
        acceleration_residual_refit[fit_acceleration_keep] / 1e3,
        ".",
        color="C1",
        label="Acceleration",
    )[0]
    axis.axhline(0.0, color="black", lw=0.8)
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Doppler residual (m/s)", color="C0")
    acceleration_axis.set_ylabel(r"Acceleration residual (km s$^{-2}$)", color="C1")
    axis.tick_params(axis="y", colors="C0")
    acceleration_axis.tick_params(axis="y", colors="C1")
    axis.legend(
        [doppler_artist, acceleration_artist],
        ["Doppler", "Acceleration"],
        frameon=False,
        loc="upper left",
    )
    for axis in event_axes.ravel():
        axis.grid(alpha=0.2, lw=0.5)
    event_fig.savefig(event_png, dpi=190)
    plt.close(event_fig)
    print(
        f"triplets={len(result)} shared_inliers={np.count_nonzero(compare)} "
        f"fft_velocity_rms_mps={rms(fft_residual[compare]):.6f} "
        f"beat_velocity_rms_mps={rms(beat_residual[compare]):.6f} "
        f"complex_velocity_rms_mps={rms(fit_residual[compare]):.6f} "
        f"acceleration_comparison={np.count_nonzero(acceleration_compare)} "
        f"phase_acceleration_rms_mps2={rms(old_acceleration_residual[acceleration_compare]):.6f} "
        f"complex_acceleration_rms_mps2={rms(fit_acceleration_residual[acceleration_compare]):.6f}"
        f" refit_radius_um={1e6 * 10.0 ** dynamics['parameters'][6]:.6f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
