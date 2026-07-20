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
    event_receiver_noise_variance,
)
from pansy_coherent import (
    coherent_add_modules,
    enu_line_of_sight_to_arrival_uvw,
    estimate_event_module_voltage_gain,
)
import pansy_interferometry as pint
import pansy_modes as pmm
from fit_three_pulse_complex_envelope import (
    FS_HZ,
    PAIR_SPACING_S,
    PAIR_SPACING_TOLERANCE_S,
    fit_three_pulse_acceleration_aliases,
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
        principal = np.asarray(
            group["measured_radial_acceleration_principal_mps2"], dtype=float
        )
        covariance = np.asarray(group["phase_difference_covariance_rad2"], dtype=float)
        keep = np.asarray(group["shared_inlier_mask"], dtype=bool)
        for index, row in enumerate(support):
            lookup[tuple(int(value) for value in row)] = {
                "model_mps2": float(model[index]),
                "measured_mps2": float(measured[index]),
                "principal_mps2": float(principal[index]),
                "phase_std_rad": float(np.sqrt(covariance[index, index])),
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


def plot_triplet_components(
    fit: dict,
    observation_indices: np.ndarray,
    event_time_s: np.ndarray,
    output: Path,
) -> None:
    """Plot the exact baud-averaged beamformed voltages used by one fit."""
    pulse_index = np.asarray(fit["pulse_index"], dtype=int)
    use = np.asarray(fit["use"], dtype=bool)
    local_time_us = 1e6 * (
        np.asarray(fit["absolute_time_s"], dtype=float) - pulse_index * PAIR_SPACING_S
    )
    fig, axes = plt.subplots(
        1, 3, figsize=(11.5, 3.6), sharex=True, sharey=True, constrained_layout=True
    )
    for pulse, axis in enumerate(axes):
        select = use & (pulse_index == pulse)
        order = np.argsort(local_time_us[select])
        for component, color, label in (
            (np.real, "C0", "real"),
            (np.imag, "C1", "imaginary"),
        ):
            axis.plot(
                local_time_us[select],
                component(fit["data"][select]),
                ".",
                color=color,
                markersize=5,
                label=f"{label} measurement",
                zorder=3,
            )
            axis.plot(
                local_time_us[select][order],
                component(fit["model"][select])[order],
                color=color,
                linewidth=1.2,
                label=f"{label} model",
            )
        observation = int(observation_indices[pulse])
        axis.set_title(
            (r"Pulse $n$", r"Pulse $n+1$", r"Pulse $n+2$")[pulse]
            + rf", $t={event_time_s[observation] - event_time_s[0]:.3f}$ s"
        )
        axis.set_xlabel(r"Fast time within pulse ($\mu$s)")
        axis.grid(alpha=0.2, linewidth=0.5)
    axes[0].set_ylabel("Normalized decoded complex voltage")
    axes[0].legend(loc="best", fontsize=7, frameon=False)
    axes[-1].text(
        0.97,
        0.97,
        rf"$v_r={fit['parameters'][4] * pc.wavelength / 2e3:.3f}$ km s$^{{-1}}$"
        + "\n"
        + rf"$a_r={fit['parameters'][5] / 1e3:.2f}$ km s$^{{-2}}$",
        transform=axes[-1].transAxes,
        ha="right",
        va="top",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def overlapping_triplet_correlation(triplet_fits, observation_indices):
    """Propagate shared-pulse noise into a banded (velocity, acceleration) correlation."""
    observation_indices = np.asarray(observation_indices, dtype=int)
    if observation_indices.shape != (len(triplet_fits), 3):
        raise ValueError("triplet observation indices must have shape (n_triplets, 3)")
    column_lookup = {}
    fit_columns = []
    for triplet, (fit, observations) in enumerate(zip(triplet_fits, observation_indices)):
        pulse_index = np.asarray(fit["residual_pulse_index"], dtype=int)
        local_columns = []
        for column, local_pulse in enumerate(pulse_index):
            if local_pulse < 0:
                key = ("independent", triplet, column)
            else:
                pulse_columns = np.flatnonzero(pulse_index == local_pulse)
                ordinal = int(np.flatnonzero(pulse_columns == column)[0])
                key = ("pulse", int(observations[local_pulse]), ordinal)
            local_columns.append(column_lookup.setdefault(key, len(column_lookup)))
        fit_columns.append(np.asarray(local_columns, dtype=int))

    loading = np.zeros((2 * len(triplet_fits), len(column_lookup)), dtype=float)
    for triplet, (fit, columns) in enumerate(zip(triplet_fits, fit_columns)):
        influence = np.asarray(fit["velocity_acceleration_noise_influence"], dtype=float)
        if influence.shape != (2, len(columns)):
            raise ValueError("unexpected three-pulse influence matrix shape")
        loading[2 * triplet : 2 * triplet + 2, columns] += influence
    covariance = loading @ loading.T
    diagonal = np.sqrt(np.maximum(np.diag(covariance), np.finfo(float).tiny))
    correlation = covariance / np.outer(diagonal, diagonal)
    correlation = 0.5 * (correlation + correlation.T)
    np.fill_diagonal(correlation, 1.0)
    return correlation


def smooth_global_branch_indices(
    branch_candidates,
    predicted_velocity_mps,
    predicted_acceleration_mps2,
    velocity_sigma_mps,
    acceleration_sigma_mps2,
    smoothness_sigma_mps2=2.5e4,
    chain_id=None,
):
    """Select equally likely k=-1,0,+1 branches with weak acceleration smoothness."""
    count = len(branch_candidates)
    if count == 0:
        return np.empty(0, dtype=int)
    candidate_velocity = np.asarray(
        [row["velocity_mps"] for row in branch_candidates], dtype=float
    )
    candidate_acceleration = np.asarray(
        [row["acceleration_mps2"] for row in branch_candidates], dtype=float
    )
    if candidate_velocity.shape != (count, 3) or candidate_acceleration.shape != (
        count,
        3,
    ):
        raise ValueError("each triplet must retain exactly three branch candidates")
    emission = (
        (candidate_velocity - np.asarray(predicted_velocity_mps)[:, None])
        / float(velocity_sigma_mps)
    ) ** 2 + (
        (candidate_acceleration - np.asarray(predicted_acceleration_mps2)[:, None])
        / float(acceleration_sigma_mps2)
    ) ** 2

    if chain_id is None:
        chain_id = np.zeros(count, dtype=int)
    chain_id = np.asarray(chain_id, dtype=int)
    if chain_id.shape != (count,):
        raise ValueError("chain_id must contain one value per triplet")

    selected = np.empty(count, dtype=int)
    for chain in np.unique(chain_id):
        rows = np.flatnonzero(chain_id == chain)
        cost = emission[rows[0]].copy()
        predecessor = np.full((len(rows), 3), -1, dtype=int)
        for local_row in range(1, len(rows)):
            row = rows[local_row]
            previous = rows[local_row - 1]
            smoothness = (
                (
                    candidate_acceleration[row][None, :]
                    - candidate_acceleration[previous][:, None]
                )
                / float(smoothness_sigma_mps2)
            ) ** 2
            transition_cost = cost[:, None] + smoothness
            predecessor[local_row] = np.argmin(transition_cost, axis=0)
            cost = emission[row] + np.min(transition_cost, axis=0)

        selected[rows[-1]] = int(np.argmin(cost))
        for local_row in range(len(rows) - 1, 0, -1):
            selected[rows[local_row - 1]] = predecessor[
                local_row, selected[rows[local_row]]
            ]
    return selected


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
    triplet_correlation=None,
    velocity_keep=None,
    acceleration_keep=None,
    branch_candidates=None,
    branch_observation_indices=None,
) -> dict:
    if velocity_keep is None:
        velocity_keep = np.asarray(result["shared_inlier"], dtype=bool)
    else:
        velocity_keep = np.asarray(velocity_keep, dtype=bool)
    if acceleration_keep is None:
        acceleration_keep = velocity_keep & ~np.asarray(
            result["acceleration_at_bound"], dtype=bool
        )
    else:
        acceleration_keep = np.asarray(acceleration_keep, dtype=bool)
    if velocity_keep.shape != result.shape or acceleration_keep.shape != result.shape:
        raise ValueError("triplet measurement masks must match the fit result")
    if np.any(acceleration_keep & ~velocity_keep):
        raise ValueError("acceleration mask cannot retain a rejected velocity triplet")
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
    triplet_correlation = np.asarray(triplet_correlation, dtype=float)
    if triplet_correlation.shape != (2 * len(result), 2 * len(result)):
        raise ValueError(
            f"triplet correlation shape {triplet_correlation.shape} does not match "
            f"the {len(result)} three-pulse estimates"
        )
    fit_time = np.asarray(result["time_s"], dtype=float)
    position_count = np.count_nonzero(echo_keep)
    velocity_count = np.count_nonzero(velocity_keep)
    position_correlation = position_covariance_km2 / np.outer(
        position_sigma_km, position_sigma_km
    )
    measurement_correlation = np.eye(measurement_count)
    for echo in range(position_count):
        rows = slice(3 * echo, 3 * echo + 3)
        measurement_correlation[rows, rows] = position_correlation
    selected_triplet_rows = np.r_[
        2 * np.flatnonzero(velocity_keep),
        2 * np.flatnonzero(acceleration_keep) + 1,
    ]
    triplet_start = 3 * position_count
    measurement_correlation[triplet_start:, triplet_start:] = triplet_correlation[
        np.ix_(selected_triplet_rows, selected_triplet_rows)
    ]

    def covariance_cholesky():
        measurement_scale = np.r_[
            np.tile(position_sigma_km, np.count_nonzero(echo_keep)),
            np.full(np.count_nonzero(velocity_keep), velocity_sigma_mps),
            np.full(np.count_nonzero(acceleration_keep), acceleration_sigma_mps2),
        ]
        covariance = measurement_correlation * np.outer(
            measurement_scale, measurement_scale
        )
        covariance += np.eye(len(covariance)) * max(
            1e-12, 1e-10 * float(np.median(np.diag(covariance)))
        )
        return np.linalg.cholesky(covariance)

    joint_cholesky = covariance_cholesky()

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

    selected_branch_k = np.asarray(result["selected_alias_number"], dtype=int).copy()
    branch_iterations = 0
    if branch_candidates is not None:
        if len(branch_candidates) != len(result):
            raise ValueError("branch candidates must match the triplet result rows")
        branch_observation_indices = np.asarray(branch_observation_indices, dtype=int)
        if branch_observation_indices.shape != (len(result), 3):
            raise ValueError("branch observation indices must have shape (n, 3)")

        for branch_iterations in range(1, 7):
            prediction = predict(final.x)
            predicted_velocity = np.interp(fit_time, trajectory_time, prediction[2])
            predicted_acceleration = np.interp(fit_time, trajectory_time, prediction[3])
            selected_indices = smooth_global_branch_indices(
                branch_candidates,
                predicted_velocity,
                predicted_acceleration,
                velocity_sigma_mps,
                acceleration_sigma_mps2,
                chain_id=np.asarray(
                    [row["chain_id"] for row in branch_candidates], dtype=int
                ),
            )
            changed = False
            selected_fits = []
            for row_index, (candidates, selected) in enumerate(
                zip(branch_candidates, selected_indices)
            ):
                candidate_velocity = np.asarray(candidates["velocity_mps"], dtype=float)
                candidate_acceleration = np.asarray(
                    candidates["acceleration_mps2"], dtype=float
                )
                candidate_delta = np.asarray(candidates["delta_chi2"], dtype=float)
                selected = int(selected)
                fit = candidates["fits"][selected]
                covariance = np.asarray(fit["parameter_covariance"], dtype=float)
                selected_fits.append(fit)
                new_k = int(candidates["k"][selected])
                changed |= new_k != selected_branch_k[row_index]
                selected_branch_k[row_index] = new_k
                result["selected_alias_number"][row_index] = new_k
                result["selected_alias_delta_chi2"][row_index] = candidate_delta[selected]
                result["velocity_mps"][row_index] = candidate_velocity[selected]
                result["velocity_std_mps"][row_index] = (
                    np.sqrt(covariance[4, 4]) * pc.wavelength / 2.0
                )
                result["acceleration_mps2"][row_index] = candidate_acceleration[selected]
                result["acceleration_std_mps2"][row_index] = np.sqrt(covariance[5, 5])
                result["frequency_acceleration_correlation"][row_index] = (
                    covariance[4, 5]
                    / np.sqrt(covariance[4, 4] * covariance[5, 5])
                )

            triplet_correlation = overlapping_triplet_correlation(
                selected_fits, branch_observation_indices
            )
            measurement_correlation[triplet_start:, triplet_start:] = triplet_correlation[
                np.ix_(selected_triplet_rows, selected_triplet_rows)
            ]
            preliminary_velocity = np.interp(fit_time, trajectory_time, prediction[2])
            preliminary_acceleration = np.interp(fit_time, trajectory_time, prediction[3])
            velocity_sigma_mps = max(
                20.0,
                float(
                    np.sqrt(
                        np.mean(
                            (
                                result["velocity_mps"][velocity_keep]
                                - preliminary_velocity[velocity_keep]
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
                max_nfev=180,
            )
            if not changed:
                break

        preliminary = predict(final.x)
    preliminary_velocity = np.interp(fit_time, trajectory_time, preliminary[2])
    preliminary_acceleration = np.interp(fit_time, trajectory_time, preliminary[3])
    # Keep the position error model fixed from the established trajectory
    # residuals.  Re-estimating it from this candidate fit lets a constant
    # position offset inflate its own uncertainty and become self-consistent.
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

    retained_time = trajectory_time[echo_keep] - float(np.min(trajectory_time))
    kinematic_seed = np.empty(6, dtype=float)
    for component in range(3):
        slope, intercept = np.polyfit(
            retained_time, points_km[echo_keep, component] * 1e3, 1
        )
        kinematic_seed[component] = intercept
        kinematic_seed[3 + component] = slope

    def fit_profile_index(index, seeds):
        fixed_log_radius = profile_log_radius_m[index]

        def fixed_residual(parameters6):
            return residual(np.r_[parameters6, fixed_log_radius])

        fitted_candidates = []
        for seed in seeds:
            fitted_candidates.append(
                least_squares(
                    fixed_residual,
                    np.clip(seed, lower6 + 1e-10, upper6 - 1e-10),
                    bounds=(lower6, upper6),
                    x_scale=scale6,
                    loss="linear",
                    max_nfev=400,
                )
            )
        fitted = min(
            fitted_candidates,
            key=lambda candidate: float(np.sum(fixed_residual(candidate.x) ** 2)),
        )
        fitted_chi2 = float(np.sum(fixed_residual(fitted.x) ** 2))
        if np.isfinite(profile_chi2[index]) and profile_chi2[index] <= fitted_chi2:
            return profile_parameters6[index]
        profile_parameters6[index] = fitted.x
        profile_chi2[index] = fitted_chi2
        profile_success[index] = bool(fitted.success)
        return profile_parameters6[index]

    center_seed = fit_profile_index(center, [final.x[:6], kinematic_seed])
    seed = center_seed
    for index in range(center - 1, -1, -1):
        seed = fit_profile_index(index, [seed, final.x[:6]])
    seed = center_seed
    for index in range(center + 1, len(profile_radius_um)):
        seed = fit_profile_index(index, [seed, kinematic_seed, final.x[:6]])

    # Propagate any better basin found at the large-radius asymptote back toward
    # the center.  The radius profile must be the lower envelope of converged
    # fixed-radius fits, not the result of one continuation direction.
    seed = profile_parameters6[-1]
    for index in range(len(profile_radius_um) - 2, center - 1, -1):
        seed = fit_profile_index(index, [profile_parameters6[index], seed])

    minimum_chi2 = min(float(np.sum(residual(final.x) ** 2)), float(np.nanmin(profile_chi2)))
    profile_delta_chi2 = profile_chi2 - minimum_chi2
    profile_probability, profile_weights = log_grid_probability(
        profile_radius_um, profile_delta_chi2
    )
    radius_quantiles_um = weighted_quantile(
        profile_radius_um, profile_weights, [0.025, 0.5, 0.975]
    )
    marginal_best_index = int(np.nanargmin(profile_delta_chi2))
    upper_tail = profile_delta_chi2[marginal_best_index:]
    upper_limit_data_constrained = bool(
        len(upper_tail) >= 3
        and np.all(np.isfinite(upper_tail[-3:]))
        and np.min(upper_tail[-3:]) > 3.841458820694124
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
    fixed_position_km = []
    fixed_doppler_mps = []
    fixed_acceleration_mps2 = []
    fixed_speed_km_s = []
    for target_um in fixed_radius_um:
        index = int(np.flatnonzero(profile_radius_um == target_um)[0])
        predicted = predict(np.r_[profile_parameters6[index], profile_log_radius_m[index]])
        fixed_position_km.append(predicted[0])
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
        "measurement_correlation": measurement_correlation,
        "triplet_correlation": triplet_correlation,
        "covariance_degrees_of_freedom": covariance_degrees_of_freedom,
        "covariance_variance_inflation": covariance_variance_inflation,
        "selected_branch_k": selected_branch_k,
        "branch_assignment_iterations": branch_iterations,
        "cost": float(final.cost),
        "success": bool(final.success),
        "profile_radius_um": profile_radius_um,
        "profile_parameters6": profile_parameters6,
        "profile_chi2": profile_chi2,
        "profile_delta_chi2": profile_delta_chi2,
        "profile_probability_method": "standard fixed-covariance profile likelihood with shared-pulse banded covariance",
        "profile_probability_density_log_radius": profile_probability,
        "profile_probability_weights": profile_weights,
        "profile_success": profile_success,
        "marginal_radius_quantiles_um": radius_quantiles_um,
        "marginal_mass_quantiles_kg": radius_um_to_mass_kg(radius_quantiles_um),
        "radius_upper_limit_data_constrained": upper_limit_data_constrained,
        "position_interval_km": position_interval_km,
        "doppler_interval_mps": doppler_interval_mps,
        "acceleration_interval_mps2": acceleration_interval_mps2,
        "speed_interval_km_s": speed_interval_km_s,
        "range_interval_km": range_interval_km,
        "fixed_radius_um": fixed_radius_um,
        "fixed_position_km": np.asarray(fixed_position_km),
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
    parser.add_argument("--debug-triplet-dir", type=Path)
    parser.add_argument("--debug-triplet-count", type=int, default=0)
    parser.add_argument("--debug-only", action="store_true")
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
    triplets.sort(
        key=lambda row: (
            float(decoded["tx_idx"][row[1]]),
            float(decoded["tx_idx"][row[0]]),
            float(decoded["tx_idx"][row[2]]),
        )
    )
    prior_acceleration = phase_acceleration_lookup(args.prior_profile_h5)
    with h5py.File(args.prior_profile_h5, "r") as handle:
        echo_keep = np.asarray(handle["result/echo_shared_inlier_mask"], dtype=bool)

    observations, stored_fit = load_selected_fit(args.diagnostics_h5)
    trajectory_time = np.asarray(observations["t_s"], dtype=float)
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
        ("selected_alias_number", "i4"),
        ("complex_best_alias_number", "i4"),
        ("selected_alias_delta_chi2", "f8"),
        ("second_alias_delta_chi2", "f8"),
    ]
    rows = []
    alias_rows = []
    selected_triplet_fits = []
    branch_candidates = []
    triplet_debug = []
    n_fast = decoded["z_tx"].shape[1]
    chain_id = -1
    previous_triplet = None
    for previous, middle, current in triplets:
        observation_indices = np.asarray([previous, middle, current], dtype=int)
        beam = int(decoded["beam_id"][middle])
        if (
            previous_triplet is None
            or previous_triplet[1] != previous
            or previous_triplet[2] != middle
            or previous_triplet[3] != beam
        ):
            chain_id += 1
        previous_triplet = (previous, middle, current, beam)
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
            local_prior["measured_mps2"]
            if local_prior is not None
            and np.isfinite(local_prior["measured_mps2"])
            else float(np.interp(middle_time, trajectory_time, model_acceleration_mps2))
        )
        velocity_guess = float(np.interp(middle_time, trajectory_time, model_doppler_mps))
        pulse_spacing_s = float(
            np.mean(np.diff(decoded["tx_idx"][observation_indices])) / FS_HZ
        )
        if local_prior is not None and np.isfinite(local_prior["phase_std_rad"]):
            phase_acceleration_sigma = (
                pc.wavelength
                * local_prior["phase_std_rad"]
                / (4.0 * np.pi * pulse_spacing_s**2)
            )
            branch_half_width = float(
                np.clip(4.0 * phase_acceleration_sigma, 2.0e3, 1.0e4)
            )
        else:
            branch_half_width = 1.0e4
        alias_fits = fit_three_pulse_acceleration_aliases(
            raw_pulses,
            templates,
            pulse_spacing_s,
            2.0 * velocity_guess / pc.wavelength,
            acceleration_guess,
            pc.wavelength,
            pulse_snr=10.0
            ** (np.asarray(decoded["snr"][observation_indices], dtype=float) / 10.0),
            matched_filter_amplitudes=amplitude_prior,
            acceleration_branch_half_width_mps2=branch_half_width,
        )
        alias_number = np.asarray(alias_fits["alias_number"], dtype=int)
        branch_indices = np.asarray(
            [int(np.flatnonzero(alias_number == k)[0]) for k in (-1, 0, 1)],
            dtype=int,
        )
        branch_candidates.append(
            {
                "chain_id": chain_id,
                "k": alias_number[branch_indices].copy(),
                "velocity_mps": np.asarray(alias_fits["fit_velocity_mps"])[
                    branch_indices
                ].copy(),
                "acceleration_mps2": np.asarray(
                    alias_fits["fit_acceleration_mps2"]
                )[branch_indices].copy(),
                "delta_chi2": np.asarray(alias_fits["delta_chi2"])[
                    branch_indices
                ].copy(),
                "fits": [alias_fits["fits"][index] for index in branch_indices],
            }
        )
        # Start at k=0. The event-wide trajectory fit chooses among all three.
        selected_alias_index = int(branch_indices[1])
        fit = alias_fits["fits"][selected_alias_index]
        selected_triplet_fits.append(fit)
        alias_rows.append(
            {
                "observation_indices": observation_indices.copy(),
                **{
                    name: np.asarray(alias_fits[name]).copy()
                    for name in (
                        "alias_number",
                        "fit_velocity_mps",
                        "fit_velocity_std_mps",
                        "fit_acceleration_mps2",
                        "fit_acceleration_std_mps2",
                        "weighted_sse",
                        "delta_chi2",
                    )
                },
            }
        )
        if args.debug_triplet_count > 0:
            triplet_debug.append((fit, observation_indices.copy()))
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
                False,
                int(alias_fits["alias_number"][selected_alias_index]),
                int(alias_fits["alias_number"][alias_fits["selected_index"]]),
                float(alias_fits["delta_chi2"][selected_alias_index]),
                float(
                    np.partition(alias_fits["delta_chi2"], 1)[1]
                    if len(alias_fits["delta_chi2"]) > 1
                    else np.inf
                ),
            )
        )
    result = np.asarray(rows, dtype=dtype)
    if len(result) == 0:
        raise RuntimeError("no valid three-pulse fits")

    if args.debug_triplet_count > 0:
        if args.debug_triplet_dir is None:
            raise ValueError("--debug-triplet-dir is required with --debug-triplet-count")
        selected = np.unique(
            np.linspace(
                0,
                len(triplet_debug) - 1,
                min(args.debug_triplet_count, len(triplet_debug)),
                dtype=int,
            )
        )
        for debug_rank, triplet_index in enumerate(selected, start=1):
            fit, observation_indices = triplet_debug[int(triplet_index)]
            plot_triplet_components(
                fit,
                observation_indices,
                trajectory_time,
                args.debug_triplet_dir
                / (
                    f"triplet_{debug_rank:02d}_index_{triplet_index:03d}_"
                    f"obs_{observation_indices[0]:03d}_{observation_indices[1]:03d}_"
                    f"{observation_indices[2]:03d}.png"
                ),
            )
        if args.debug_only:
            print(
                f"debug_triplets={len(selected)} output={args.debug_triplet_dir}",
                flush=True,
            )
            return 0

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
    fit_velocity_keep = np.asarray(result["shared_inlier"], dtype=bool)
    fit_acceleration_keep = fit_velocity_keep & ~np.asarray(
        result["acceleration_at_bound"], dtype=bool
    )
    triplet_correlation = overlapping_triplet_correlation(
        selected_triplet_fits,
        np.column_stack((result["previous"], result["middle"], result["current"])),
    )
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
        triplet_correlation=triplet_correlation,
        velocity_keep=fit_velocity_keep,
        acceleration_keep=fit_acceleration_keep,
        branch_candidates=branch_candidates,
        branch_observation_indices=np.column_stack(
            (result["previous"], result["middle"], result["current"])
        ),
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
        handle.attrs["steering_convention"] = "catalogue_arrival_uvw"
        handle.create_dataset("triplet_fit", data=result)
        handle.create_dataset("dynamics_velocity_keep", data=fit_velocity_keep)
        handle.create_dataset("dynamics_acceleration_keep", data=fit_acceleration_keep)
        alias_group = handle.create_group("triplet_aliases")
        alias_group.create_dataset(
            "observation_indices",
            data=np.asarray([row["observation_indices"] for row in alias_rows]),
        )
        for name in (
            "alias_number",
            "fit_velocity_mps",
            "fit_velocity_std_mps",
            "fit_acceleration_mps2",
            "fit_acceleration_std_mps2",
            "weighted_sse",
            "delta_chi2",
        ):
            base_dtype = np.int32 if name == "alias_number" else np.float64
            dataset = alias_group.create_dataset(
                name,
                shape=(len(alias_rows),),
                dtype=h5py.vlen_dtype(np.dtype(base_dtype)),
            )
            for row_index, row in enumerate(alias_rows):
                dataset[row_index] = np.asarray(row[name], dtype=base_dtype)
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
    compare = fit_velocity_keep
    acceleration_compare = fit_acceleration_keep

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
    path_length_km = float(np.sum(np.linalg.norm(np.diff(refit_position[echo_keep], axis=0), axis=1)))

    profile_radius_um = dynamics["profile_radius_um"]
    profile_probability = dynamics["profile_probability_density_log_radius"]
    radius_quantiles_um = dynamics["marginal_radius_quantiles_um"]
    mass_quantiles_kg = dynamics["marginal_mass_quantiles_kg"]
    radius_upper_limit_data_constrained = bool(
        dynamics["radius_upper_limit_data_constrained"]
    )
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
    axis.axvline(1e6 * 10.0 ** dynamics["parameters"][6], color="C0", label="best fit")
    axis.axvline(radius_quantiles_um[0], color="C0", ls="--", lw=1.0)
    if radius_upper_limit_data_constrained:
        axis.axvspan(radius_quantiles_um[0], radius_quantiles_um[-1], color="C0", alpha=0.12)
        interval_text = (
            rf"95% $r_0$ {radius_quantiles_um[0]:.0f}--{radius_quantiles_um[-1]:.0f} $\mu$m"
            + "\n"
            + rf"95% $m_0$ {mass_quantiles_kg[0]:.1e}--{mass_quantiles_kg[-1]:.1e} kg"
        )
    else:
        interval_text = (
            rf"95% lower $r_0>{radius_quantiles_um[0]:.0f}$ $\mu$m; no upper bound"
            + "\n"
            + rf"95% lower $m_0>{mass_quantiles_kg[0]:.1e}$ kg; no upper bound"
        )
    axis.text(
        0.04,
        0.94,
        interval_text,
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
    along_track = refit_velocity / np.maximum(
        np.linalg.norm(refit_velocity, axis=1, keepdims=True), 1e-12
    )

    def along_track_residual(model_position_km):
        return np.einsum(
            "ij,ij->i", fitted_points_km - model_position_km, along_track
        )

    axis.plot(
        observation_time[echo_keep],
        along_track_residual(refit_position)[echo_keep],
        ".",
        color="C0",
        markersize=3.0,
        label="fit",
    )
    for target_um, fixed_position, color in zip(
        dynamics["fixed_radius_um"], dynamics["fixed_position_km"], fixed_colors
    ):
        axis.plot(
            observation_time[echo_keep],
            along_track_residual(fixed_position)[echo_keep],
            linestyle="none",
            marker=".",
            color=color,
            markersize=3.0,
            label=rf"$r_0={target_um:g}\,\mu$m",
        )
    axis.axhline(0.0, color="black", lw=0.8)
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Along-track position residual (km)")
    axis.legend(frameon=False, loc="best", fontsize=7)
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
