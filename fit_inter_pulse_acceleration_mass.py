#!/usr/bin/env python3
"""Compare shrinking-radius profiles with and without decoded beat-phase acceleration."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la
import scipy.optimize as opt

import fit_best_alias_physics_models as physics
from mass_profile_grid import adaptive_profile_radii
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


def lower_banded(matrix: np.ndarray, tolerance: float = 0.0) -> tuple[np.ndarray, int]:
    rows, columns = np.nonzero(np.abs(matrix) > tolerance)
    bandwidth = int(np.max(rows - columns, initial=0))
    banded = np.zeros((bandwidth + 1, len(matrix)), dtype=float)
    for offset in range(bandwidth + 1):
        diagonal = np.diag(matrix, k=-offset)
        banded[offset, : len(diagonal)] = diagonal
    return banded, bandwidth


def covariance_weighted_residual(covariance: np.ndarray, residual: np.ndarray) -> np.ndarray:
    banded, bandwidth = lower_banded(covariance)
    cholesky = la.cholesky_banded(banded, lower=True, check_finite=False)
    return la.solve_banded(
        (bandwidth, 0), cholesky, residual, check_finite=False
    )


def map_measurements_to_observations(measurement_tx_s: np.ndarray, t_s: np.ndarray) -> np.ndarray:
    absolute_zero = measurement_tx_s[0] - t_s[0]
    measurement_t_s = measurement_tx_s - absolute_zero
    nearest = np.abs(measurement_t_s[:, None] - t_s[None, :]).argmin(axis=1)
    error_s = np.abs(measurement_t_s - t_s[nearest])
    if np.any(error_s > 5e-6):
        raise RuntimeError(f"measurement/trajectory time alignment error is {np.max(error_s):.3g} s")
    return nearest


def load_nonoverlapping_phase_acceleration(path: Path, acceleration_lag_s: float = 0.008) -> dict:
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
        for first_index, first in enumerate(indices[:-1]):
            later = indices[first_index + 1 :]
            separation = pair["time_s"][later] - pair["time_s"][first]
            matches = later[np.abs(separation - acceleration_lag_s) <= 5e-6]
            if len(matches) == 0:
                continue
            second = int(matches[0])
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
    samples = np.asarray(rows, dtype=dtype)
    incidence = np.zeros((len(samples), len(beat)), dtype=float)
    row_index = np.arange(len(samples))
    incidence[row_index, samples["first_pair"]] = -1.0
    incidence[row_index, samples["second_pair"]] = 1.0
    beat_variance = np.asarray(pair["beat_phase_std_rad"], dtype=float) ** 2
    finite_variance = beat_variance[np.isfinite(beat_variance) & (beat_variance > 0.0)]
    variance_floor = float(np.nanmin(finite_variance)) if len(finite_variance) else 1.0
    beat_variance = np.where(
        np.isfinite(beat_variance) & (beat_variance > 0.0),
        beat_variance,
        variance_floor,
    )
    covariance = (incidence * beat_variance[None, :]) @ incidence.T
    return {
        "samples": samples,
        "measurement_tx_s": measurement_tx_s,
        "pair_previous_index": previous,
        "pair_current_index": current,
        "acceleration_lag_s": float(acceleration_lag_s),
        "phase_difference_covariance_rad2": covariance,
        "phase_incidence": incidence,
        "beat_phase_variance_rad2": beat_variance,
    }


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


def fit_profile(
    diagnostics_h5: Path,
    baseline_h5: Path,
    beat_h5: Path,
    output_h5: Path,
    output_png: Path | None,
    max_nfev: int = 60,
    global_starts: bool = False,
    seed_profile_h5: Path | None = None,
    adaptive_profile: bool = False,
    adaptive_spacing_dex: float = 0.02,
    adaptive_max_nfev: int = 120,
) -> dict:
    observations, stored = load_selected_fit(diagnostics_h5)
    t_s = observations["t_s"]
    points = np.asarray(observations["points_km"], dtype=float).copy()
    doppler = observations["doppler_km_s"]
    with h5py.File(beat_h5, "r") as beat_handle:
        fft_range_km = np.asarray(beat_handle["measurement/range_km"], dtype=float)
        range_estimator = str(beat_handle.attrs.get("range_estimator", "unknown"))
    catalogue_range_km = np.linalg.norm(points, axis=1)
    if len(fft_range_km) != len(points):
        raise RuntimeError(
            f"FFT range count {len(fft_range_km)} does not match trajectory count {len(points)}"
        )
    valid_range = (
        np.isfinite(fft_range_km)
        & np.isfinite(catalogue_range_km)
        & (catalogue_range_km > 0.0)
    )
    points[valid_range] *= (
        fft_range_km[valid_range] / catalogue_range_km[valid_range]
    )[:, None]
    finite = np.isfinite(t_s) & np.all(np.isfinite(points), axis=1) & np.isfinite(doppler)
    keep = stored["keep"] & finite
    if np.count_nonzero(keep) < 10:
        keep = finite
    stored_prediction = physics.predicted_doppler(stored["model_km"], stored["velocity_km_s"])
    stored_position_residual = points[keep] - stored["model_km"][keep]
    sigma_position = float(np.sqrt(np.mean(stored_position_residual**2)))
    position_covariance = np.cov(stored_position_residual, rowvar=False)
    position_variance_floor = max(
        1e-10, 1e-6 * float(np.trace(position_covariance)) / 3.0
    )
    position_covariance += np.eye(3) * position_variance_floor
    position_cholesky = np.linalg.cholesky(position_covariance)
    sigma_doppler = float(np.sqrt(np.mean((doppler[keep] - stored_prediction[keep]) ** 2)))
    def prepare_phase_modality(acceleration_lag_s):
        data = load_nonoverlapping_phase_acceleration(
            beat_h5, acceleration_lag_s=acceleration_lag_s
        )
        samples = data["samples"]
        if len(samples) < 3:
            raise RuntimeError(
                f"fewer than three non-overlapping {1e3 * acceleration_lag_s:.0f}-ms acceleration samples"
            )
        stored_phase_prediction = predicted_delta_phase(
            stored_prediction * 1e3, t_s, data
        )
        stored_phase_residual = circular_residual(
            samples["observed_delta_phase_rad"], stored_phase_prediction
        )
        covariance = np.asarray(
            data["phase_difference_covariance_rad2"], dtype=float
        )
        covariance += np.eye(len(covariance)) * max(
            1e-12, 1e-9 * float(np.nanmax(np.diag(covariance)))
        )
        covariance_banded, covariance_bandwidth = lower_banded(covariance)
        stored_covariance_weighted_residual = covariance_weighted_residual(
            covariance, stored_phase_residual
        )
        covariance_scale = float(
            np.sqrt(np.mean(stored_covariance_weighted_residual**2))
        )
        if not np.isfinite(covariance_scale) or covariance_scale <= 0.0:
            raise RuntimeError(
                f"invalid {1e3 * acceleration_lag_s:.0f}-ms covariance scale"
            )
        absolute_zero = data["measurement_tx_s"][0] - t_s[0]
        midpoint_t = (
            0.5 * (samples["first_time_s"] + samples["second_time_s"])
            - absolute_zero
        )
        snr = np.interp(midpoint_t, t_s, observations["snr"])
        weight = np.clip(snr, 0.0, 100.0)
        finite_positive_weight = np.isfinite(weight) & (weight > 0.0)
        if not np.any(finite_positive_weight):
            weight = np.ones_like(stored_phase_residual)
        else:
            floor = float(np.nanmin(weight[finite_positive_weight]))
            weight = np.where(finite_positive_weight, weight, floor)
        weight /= np.mean(weight)
        sigma_phase = float(
            np.sqrt(np.sum(weight * stored_phase_residual**2) / np.sum(weight))
        )
        if not np.isfinite(sigma_phase) or sigma_phase <= 0.0:
            raise RuntimeError(
                f"invalid {1e3 * acceleration_lag_s:.0f}-ms beat-phase residual RMS"
            )
        mean_dt = 0.5 * (samples["first_dt_s"] + samples["second_dt_s"])
        measured_principal = (
            pc.wavelength
            * samples["observed_delta_phase_rad"]
            / (4.0 * np.pi * mean_dt * acceleration_lag_s)
        )
        stored_acceleration = predicted_radial_acceleration(
            stored_prediction * 1e3, t_s, data
        )
        ambiguity = pc.wavelength / (2.0 * mean_dt * acceleration_lag_s)
        acceleration_residual = measured_principal - stored_acceleration
        acceleration_residual -= (
            np.rint(acceleration_residual / ambiguity) * ambiguity
        )
        sigma_acceleration = float(np.sqrt(np.mean(acceleration_residual**2)))
        measurement_observation = map_measurements_to_observations(
            data["measurement_tx_s"], t_s
        )
        first_pair = samples["first_pair"]
        second_pair = samples["second_pair"]
        support = measurement_observation[
            np.column_stack(
                (
                    data["pair_previous_index"][first_pair],
                    data["pair_current_index"][first_pair],
                    data["pair_previous_index"][second_pair],
                    data["pair_current_index"][second_pair],
                )
            )
        ]
        return {
            "data": data,
            "samples": samples,
            "stored_prediction": stored_phase_prediction,
            "stored_residual": stored_phase_residual,
            "snr": snr,
            "weight": weight,
            "sigma_phase": sigma_phase,
            "covariance": covariance,
            "covariance_banded": covariance_banded,
            "covariance_bandwidth": covariance_bandwidth,
            "covariance_scale": covariance_scale,
            "marginal_outlier_sigma_rad": np.sqrt(np.diag(covariance))
            * covariance_scale,
            "mean_dt": mean_dt,
            "measured_principal": measured_principal,
            "stored_acceleration": stored_acceleration,
            "acceleration_ambiguity": ambiguity,
            "acceleration_residual": acceleration_residual,
            "sigma_acceleration": sigma_acceleration,
            "support": support,
        }

    phase8 = prepare_phase_modality(0.008)
    phase16 = prepare_phase_modality(0.016)
    cross_phase_covariance = (
        phase8["data"]["phase_incidence"]
        * phase8["data"]["beat_phase_variance_rad2"][None, :]
    ) @ phase16["data"]["phase_incidence"].T
    combined_phase_covariance = np.block(
        [
            [phase8["covariance"], cross_phase_covariance],
            [cross_phase_covariance.T, phase16["covariance"]],
        ]
    )
    combined_phase_covariance += np.eye(len(combined_phase_covariance)) * max(
        1e-12, 1e-9 * float(np.nanmax(np.diag(combined_phase_covariance)))
    )
    combined_stored_phase_residual = np.r_[
        phase8["stored_residual"], phase16["stored_residual"]
    ]
    combined_phase_covariance_scale = float(
        np.sqrt(
            np.mean(
                covariance_weighted_residual(
                    combined_phase_covariance, combined_stored_phase_residual
                )
                ** 2
            )
        )
    )
    if not np.isfinite(combined_phase_covariance_scale) or combined_phase_covariance_scale <= 0.0:
        raise RuntimeError("invalid combined beat-phase covariance scale")
    phase_offset = 0
    for modality in (phase8, phase16):
        count = len(modality["samples"])
        modality["combined_covariance_slice"] = slice(phase_offset, phase_offset + count)
        modality["covariance_scale"] = combined_phase_covariance_scale
        modality["marginal_outlier_sigma_rad"] = np.sqrt(
            np.diag(combined_phase_covariance)[phase_offset : phase_offset + count]
        ) * combined_phase_covariance_scale
        phase_offset += count
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
    seed_parameters6 = None
    if seed_profile_h5 is not None:
        with h5py.File(seed_profile_h5, "r") as seed_profile:
            seed_radius_um = np.asarray(seed_profile["profile/radius_um"], dtype=float)
            source_parameters6 = np.asarray(seed_profile["profile/parameters6"], dtype=float)
        seed_parameters6 = np.column_stack(
            [
                np.interp(np.log(radius_um), np.log(seed_radius_um), source_parameters6[:, column])
                for column in range(source_parameters6.shape[1])
            ]
        )
    lower6 = np.asarray([-np.inf, -np.inf, 20e3, -90e3, -90e3, -90e3])
    upper6 = np.asarray([np.inf, np.inf, 220e3, 90e3, 90e3, 90e3])
    scale6 = np.asarray([1e5, 1e5, 1e5, 7e4, 7e4, 7e4])

    def normalized_residuals(parameters6, fixed_log_radius):
        position, velocity, _radius, _mass, _success, _message = physics.propagate_shrinking_radius_model(
            np.r_[parameters6, fixed_log_radius], t_s, density
        )
        prediction = physics.predicted_doppler(position, velocity)
        phase_residuals = []
        phase_outlier_residuals = []
        for modality in (phase8, phase16):
            phase_prediction = predicted_delta_phase(
                prediction * 1e3, t_s, modality["data"]
            )
            raw_residual = circular_residual(
                modality["samples"]["observed_delta_phase_rad"],
                phase_prediction,
            )
            phase_residuals.append(raw_residual)
            phase_outlier_residuals.append(
                raw_residual / modality["marginal_outlier_sigma_rad"]
            )
        return (
            la.solve_triangular(
                position_cholesky,
                (points - position).T,
                lower=True,
                check_finite=False,
            ).T,
            (doppler - prediction) / sigma_doppler,
            phase_residuals[0],
            phase_residuals[1],
            phase_outlier_residuals[0],
            phase_outlier_residuals[1],
        )

    def residual6(parameters6, fixed_log_radius, echo_keep=None):
        (
            position_residual,
            doppler_residual,
            phase8_residual,
            phase16_residual,
            _phase8_outlier_residual,
            _phase16_outlier_residual,
        ) = normalized_residuals(parameters6, fixed_log_radius)
        if echo_keep is None:
            point_keep = finite
            phase8_keep = np.ones(len(phase8["support"]), dtype=bool)
            phase16_keep = np.ones(len(phase16["support"]), dtype=bool)
        else:
            point_keep = finite & echo_keep
            phase8_keep = np.all(echo_keep[phase8["support"]], axis=1)
            phase16_keep = np.all(echo_keep[phase16["support"]], axis=1)

        phase8_indices = np.flatnonzero(phase8_keep)
        phase16_indices = np.flatnonzero(phase16_keep) + len(phase8_residual)
        combined_phase_indices = np.r_[phase8_indices, phase16_indices]
        combined_phase_residual = np.r_[
            phase8_residual[phase8_keep], phase16_residual[phase16_keep]
        ]
        selected_phase_covariance = combined_phase_covariance[
            np.ix_(combined_phase_indices, combined_phase_indices)
        ]
        weighted_phase_residual = (
            covariance_weighted_residual(
                selected_phase_covariance, combined_phase_residual
            )
            / combined_phase_covariance_scale
        )

        return np.r_[
            position_residual[point_keep].ravel(),
            doppler_residual[point_keep],
            weighted_phase_residual,
        ]

    robust_objective = np.full(len(radius_um), np.nan)
    robust_parameters6 = np.full((len(radius_um), 6), np.nan)
    for index in np.argsort(np.abs(np.log(radius_um) - np.log(baseline_best))):
        starts = [baseline_params6[index]]
        if seed_parameters6 is not None:
            starts.append(seed_parameters6[index])
        if index > 0 and np.all(np.isfinite(robust_parameters6[index - 1])):
            starts.append(robust_parameters6[index - 1])
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
                    loss="soft_l1",
                    f_scale=1.0,
                    max_nfev=max_nfev,
                )
                candidates.append(result)
            except Exception:
                pass
        if candidates:
            best = min(candidates, key=lambda candidate: candidate.cost)
            robust_parameters6[index] = best.x
            robust_objective[index] = float(2.0 * best.cost)

    robust_best_index = int(np.nanargmin(robust_objective))
    robust_residuals = normalized_residuals(
        robust_parameters6[robust_best_index], log_radius_m[robust_best_index]
    )
    echo_residual_score = np.full(len(t_s), np.nan)
    (
        position_residual,
        doppler_residual,
        _phase8_residual,
        _phase16_residual,
        phase8_outlier_residual,
        phase16_outlier_residual,
    ) = robust_residuals
    echo_residual_score[finite] = np.sqrt(
        np.sum(position_residual[finite] ** 2, axis=1) + doppler_residual[finite] ** 2
    )
    for modality, residuals in (
        (phase8, phase8_outlier_residual),
        (phase16, phase16_outlier_residual),
    ):
        for support, residual in zip(modality["support"], residuals):
            unique_support = np.unique(support)
            echo_residual_score[unique_support] = np.fmax(
                echo_residual_score[unique_support], abs(residual)
            )
    scored = finite & np.isfinite(echo_residual_score)
    echo_keep = scored & (echo_residual_score < 3.5)
    if np.count_nonzero(echo_keep) < 10:
        raise RuntimeError("joint outlier rejection retained fewer than ten echoes")
    phase8_keep = np.all(echo_keep[phase8["support"]], axis=1)
    phase16_keep = np.all(echo_keep[phase16["support"]], axis=1)

    chi2 = np.full(len(radius_um), np.nan)
    parameters6 = np.full((len(radius_um), 6), np.nan)
    success = np.zeros(len(radius_um), dtype=bool)
    for index in np.argsort(np.abs(np.log(radius_um) - np.log(radius_um[robust_best_index]))):
        starts = [robust_parameters6[index]]
        if seed_parameters6 is not None:
            starts.append(seed_parameters6[index])
        if index > 0 and np.all(np.isfinite(parameters6[index - 1])):
            starts.append(parameters6[index - 1])
        if global_starts and np.all(np.isfinite(robust_parameters6[robust_best_index])):
            starts.append(robust_parameters6[robust_best_index])
        candidates = []
        for start in starts:
            if not np.all(np.isfinite(start)):
                continue
            try:
                result = opt.least_squares(
                    lambda values: residual6(values, log_radius_m[index], echo_keep),
                    start,
                    bounds=(lower6, upper6),
                    x_scale=scale6,
                    loss="linear",
                    max_nfev=max_nfev,
                )
                candidates.append(result)
            except Exception:
                pass
        if candidates:
            best = min(candidates, key=lambda candidate: candidate.cost)
            parameters6[index] = best.x
            chi2[index] = float(2.0 * best.cost)
            success[index] = bool(best.success)
    minimum = float(np.nanmin(chi2))
    delta = chi2 - minimum
    coarse_radius_um = radius_um.copy()
    coarse_delta = delta.copy()
    adaptive_points_added = 0
    continuation_points_improved = 0
    if adaptive_profile:
        extra_radius_um = adaptive_profile_radii(
            radius_um,
            delta,
            spacing_dex=adaptive_spacing_dex,
        )
        if len(extra_radius_um):
            existing_radius_um = radius_um.copy()
            existing_parameters6 = parameters6.copy()
            existing_chi2 = chi2.copy()
            existing_success = success.copy()
            best_existing = int(np.nanargmin(existing_chi2))
            extra_parameters6 = np.full((len(extra_radius_um), 6), np.nan)
            extra_chi2 = np.full(len(extra_radius_um), np.nan)
            extra_success = np.zeros(len(extra_radius_um), dtype=bool)
            finite_parameter_rows = np.all(np.isfinite(existing_parameters6), axis=1)
            for extra_index in np.argsort(
                np.abs(np.log(extra_radius_um) - np.log(existing_radius_um[best_existing]))
            ):
                target_radius = extra_radius_um[extra_index]
                target_log_radius_m = np.log10(target_radius * 1e-6)
                starts = []
                if np.count_nonzero(finite_parameter_rows) >= 2:
                    starts.append(
                        np.asarray(
                            [
                                np.interp(
                                    np.log(target_radius),
                                    np.log(existing_radius_um[finite_parameter_rows]),
                                    existing_parameters6[finite_parameter_rows, column],
                                )
                                for column in range(6)
                            ]
                        )
                    )
                nearest = int(np.nanargmin(np.abs(np.log(existing_radius_um) - np.log(target_radius))))
                starts.append(existing_parameters6[nearest])
                starts.append(existing_parameters6[best_existing])
                if np.any(np.all(np.isfinite(extra_parameters6), axis=1)):
                    fitted_extra = np.flatnonzero(np.all(np.isfinite(extra_parameters6), axis=1))
                    nearest_extra = fitted_extra[
                        np.argmin(np.abs(np.log(extra_radius_um[fitted_extra]) - np.log(target_radius)))
                    ]
                    starts.append(extra_parameters6[nearest_extra])
                candidates = []
                for start in starts:
                    if not np.all(np.isfinite(start)):
                        continue
                    try:
                        result = opt.least_squares(
                            lambda values: residual6(values, target_log_radius_m, echo_keep),
                            start,
                            bounds=(lower6, upper6),
                            x_scale=scale6,
                            loss="linear",
                            max_nfev=adaptive_max_nfev,
                        )
                        candidates.append(result)
                    except Exception:
                        pass
                if candidates:
                    best = min(candidates, key=lambda candidate: candidate.cost)
                    extra_parameters6[extra_index] = best.x
                    extra_chi2[extra_index] = float(2.0 * best.cost)
                    extra_success[extra_index] = bool(best.success)
            radius_um = np.r_[existing_radius_um, extra_radius_um]
            parameters6 = np.vstack((existing_parameters6, extra_parameters6))
            chi2 = np.r_[existing_chi2, extra_chi2]
            success = np.r_[existing_success, extra_success]
            sort_order = np.argsort(radius_um)
            radius_um = radius_um[sort_order]
            parameters6 = parameters6[sort_order]
            chi2 = chi2[sort_order]
            success = success[sort_order]
            log_radius_m = np.log10(radius_um * 1e-6)
            minimum = float(np.nanmin(chi2))
            delta = chi2 - minimum
            adaptive_points_added = int(len(extra_radius_um))

            # A newly discovered optimizer branch must continue smoothly in
            # radius. Walk its solution toward both sides of the low-chi-square
            # basin so an isolated numerical dip cannot define the interval.
            continuation_limit = 50.0
            best_index = int(np.nanargmin(chi2))
            selected = np.isfinite(delta) & (delta <= continuation_limit)
            selected_indices = np.flatnonzero(selected)
            if len(selected_indices):
                continuation_left = max(0, int(selected_indices[0]) - 1)
                continuation_right = min(len(radius_um) - 1, int(selected_indices[-1]) + 1)
                for indices in (
                    range(best_index - 1, continuation_left - 1, -1),
                    range(best_index + 1, continuation_right + 1),
                ):
                    previous = best_index
                    for index in indices:
                        start = parameters6[previous]
                        if not np.all(np.isfinite(start)):
                            break
                        try:
                            result = opt.least_squares(
                                lambda values: residual6(values, log_radius_m[index], echo_keep),
                                start,
                                bounds=(lower6, upper6),
                                x_scale=scale6,
                                loss="linear",
                                max_nfev=adaptive_max_nfev,
                            )
                        except Exception:
                            break
                        candidate_chi2 = float(2.0 * result.cost)
                        if not np.isfinite(candidate_chi2):
                            break
                        if not np.isfinite(chi2[index]) or candidate_chi2 < chi2[index]:
                            parameters6[index] = result.x
                            chi2[index] = candidate_chi2
                            success[index] = bool(result.success)
                            continuation_points_improved += 1
                        previous = index
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
    for modality in (phase8, phase16):
        modality["best_prediction"] = predicted_delta_phase(
            best_prediction * 1e3, t_s, modality["data"]
        )
        modality["best_acceleration"] = predicted_radial_acceleration(
            best_prediction * 1e3, t_s, modality["data"]
        )
        modality["display_acceleration"] = modality["measured_principal"] + np.rint(
            (modality["best_acceleration"] - modality["measured_principal"])
            / modality["acceleration_ambiguity"]
        ) * modality["acceleration_ambiguity"]

    def write_phase_group(handle, name, modality, shared_keep):
        group = handle.create_group(name)
        samples = modality["samples"]
        group["samples"] = samples
        group["time_s"] = 0.5 * (
            samples["first_time_s"] + samples["second_time_s"]
        )
        group["snr_linear"] = modality["snr"]
        group["normalized_snr_weight"] = modality["weight"]
        group["shared_inlier_mask"] = shared_keep
        group["support_observation_indices"] = modality["support"]
        group["phase_difference_covariance_rad2"] = modality["covariance"]
        group["phase_difference_covariance_lower_banded_rad2"] = modality[
            "covariance_banded"
        ]
        group["stored_model_prediction_rad"] = modality["stored_prediction"]
        group["stored_model_residual_rad"] = modality["stored_residual"]
        group["best_model_prediction_rad"] = modality["best_prediction"]
        group["measured_radial_acceleration_principal_mps2"] = modality[
            "measured_principal"
        ]
        group["measured_radial_acceleration_display_mps2"] = modality[
            "display_acceleration"
        ]
        group["acceleration_ambiguity_period_mps2"] = modality[
            "acceleration_ambiguity"
        ]
        group["best_model_radial_acceleration_mps2"] = modality[
            "best_acceleration"
        ]
        group["stored_model_radial_acceleration_mps2"] = modality[
            "stored_acceleration"
        ]
        group["stored_model_radial_acceleration_residual_mps2"] = modality[
            "acceleration_residual"
        ]
        group.attrs["acceleration_lag_s"] = modality["data"]["acceleration_lag_s"]
        group.attrs["time_reference"] = "midpoint of the two beat-phase measurements"
        group.attrs["sigma_phase_rad"] = modality["sigma_phase"]
        group.attrs["covariance_scale"] = modality["covariance_scale"]
        group.attrs["covariance_lower_bandwidth"] = modality[
            "covariance_bandwidth"
        ]
        group.attrs["sigma_radial_acceleration_mps2"] = modality[
            "sigma_acceleration"
        ]

    output_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_h5, "w") as handle:
        handle.attrs["sample_idx"] = int(observations["sample_idx"])
        handle.attrs["phase_weighting"] = "linear SNR interpolated to each phase sample, capped at 100 and normalized to unit mean"
        handle.attrs["joint_likelihood"] = "initial soft-L1 fit to all finite position, Doppler, 8-ms radial beat-phase, and 16-ms radial beat-phase measurements; shared echo rejection is the union of modality scores above 3.5; final linear profile fit uses that mask in every modality"
        handle.attrs["joint_echo_outlier_threshold"] = 3.5
        handle.attrs["phase_likelihood"] = "beat-phase residual wrapped to [-pi, pi) independently for each model"
        handle.attrs["profile_strategy"] = "adaptive_log_radius" if adaptive_profile else "baseline_grid"
        handle.attrs["outlier_mask_strategy"] = "locked once from robust coarse-grid best fit before all fixed-radius profile fits"
        handle.attrs["position_likelihood"] = "fixed full 3x3 ENU residual covariance estimated from the stored best fit"
        handle.attrs["phase_cross_likelihood"] = "single covariance over 8-ms and 16-ms phase differences, including covariance from reused decoded beat phases"
        handle.attrs["fitted_range_estimator"] = range_estimator
        handle.attrs["fitted_position_direction"] = "catalogue interferometric ENU direction rescaled to the FFT fractional-delay range"
        handle["position_residual_covariance_km2"] = position_covariance
        handle["phase_acceleration_cross_covariance_8ms_16ms_rad2"] = cross_phase_covariance
        profile = handle.create_group("profile")
        profile["radius_um"] = radius_um
        profile["mass_kg"] = radius_um_to_mass_kg(radius_um)
        profile["chi2"] = chi2
        profile["delta_chi2"] = delta
        profile["probability_density_log_radius"] = probability
        profile["parameters6"] = parameters6
        profile["success"] = success
        profile["coarse_radius_um"] = coarse_radius_um
        profile["coarse_delta_chi2"] = coarse_delta
        profile.attrs["adaptive_points_added"] = adaptive_points_added
        profile.attrs["continuation_points_improved"] = continuation_points_improved
        profile.attrs["adaptive_spacing_dex"] = adaptive_spacing_dex
        profile.attrs["adaptive_max_nfev"] = adaptive_max_nfev
        result = handle.create_group("result")
        result["best_radius_um"] = best_radius
        result["ci95_lower_radius_um"] = lower
        result["ci95_upper_radius_um"] = upper
        result["ci95_lower_bounded"] = lower_bounded
        result["ci95_upper_bounded"] = upper_bounded
        result.attrs["ci95_lower_status"] = lower_status
        result.attrs["ci95_upper_status"] = upper_status
        result["marginal_radius_quantiles_um"] = quantiles
        result["echo_normalized_joint_residual_score"] = echo_residual_score
        result["echo_shared_inlier_mask"] = echo_keep
        write_phase_group(handle, "phase_acceleration", phase8, phase8_keep)
        write_phase_group(handle, "phase_acceleration_16ms", phase16, phase16_keep)

    summary = {
        "sample_idx": int(observations["sample_idx"]),
        "baseline_best_radius_um": baseline_best,
        "baseline_interval_um": baseline_interval,
        "new_best_radius_um": best_radius,
        "new_interval_um": (lower, upper),
        "new_quantiles_um": quantiles,
        "sigma_phase_8ms_rad": phase8["sigma_phase"],
        "sigma_phase_16ms_rad": phase16["sigma_phase"],
        "n_phase_8ms": len(phase8["samples"]),
        "n_phase_16ms": len(phase16["samples"]),
        "adaptive_points_added": adaptive_points_added,
        "continuation_points_improved": continuation_points_improved,
    }
    if output_png is None:
        return summary

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3), constrained_layout=True)
    ax = axes[0]
    baseline_plot_delta = np.interp(
        np.log(radius_um), np.log(coarse_radius_um), baseline_delta
    )
    ax.plot(radius_um, baseline_plot_delta, color="0.5", lw=1.2, label="position + Doppler")
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
        + f"\nN accel={len(phase8['samples'])}+{len(phase16['samples'])}",
        transform=ax.transAxes,
        va="top",
    )
    fig.savefig(output_png, dpi=190)
    plt.close(fig)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics", type=Path)
    parser.add_argument("--baseline-profile", type=Path)
    parser.add_argument("--beat-h5", type=Path)
    parser.add_argument("--output-h5", type=Path)
    parser.add_argument("--output-png", type=Path)
    parser.add_argument("--max-nfev", type=int, default=60)
    parser.add_argument("--global-starts", action="store_true")
    parser.add_argument("--seed-profile", type=Path)
    parser.add_argument("--adaptive-profile", action="store_true")
    parser.add_argument("--adaptive-spacing-dex", type=float, default=0.02)
    parser.add_argument("--adaptive-max-nfev", type=int, default=120)
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
    summary = fit_profile(
        diagnostics,
        baseline_profile,
        beat_h5,
        output_h5,
        output_png,
        max_nfev=args.max_nfev,
        global_starts=args.global_starts,
        seed_profile_h5=args.seed_profile,
        adaptive_profile=args.adaptive_profile,
        adaptive_spacing_dex=args.adaptive_spacing_dex,
        adaptive_max_nfev=args.adaptive_max_nfev,
    )
    print(summary, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
