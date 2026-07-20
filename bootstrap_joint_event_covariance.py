#!/usr/bin/env python3
"""Estimate joint position/Doppler/acceleration covariance from pulse voltage."""

from __future__ import annotations

import argparse
import multiprocessing as mp
from pathlib import Path

import h5py
import numpy as np

import pansy_config as pc
import pansy_interferometry as pint
from aod_complex_doppler_fit import (
    event_receiver_noise_variance,
)
from pansy_coherent import (
    coherent_add_modules,
    enu_line_of_sight_to_arrival_uvw,
    estimate_event_module_voltage_gain,
)
from fit_full_event_three_pulse_complex_envelope import (
    phase_acceleration_lookup,
    valid_triplets,
)
from fit_three_pulse_complex_envelope import (
    FS_HZ,
    PAIR_SPACING_S,
    fit_three_pulse_envelope,
    rms_baud_amplitudes,
)
from interferometer_alias_diagnostics import load_cut
from inter_pulse_phase_deceleration import (
    baud_averaged_beat_pairs,
    decoded_pulse_responses,
    fractional_segment,
    load_selected,
)
from run_catalogue_mass_profiles import load_selected_fit


_STATE = None


def tangent_basis(unit_vector: np.ndarray) -> np.ndarray:
    reference = np.asarray([0.0, 0.0, 1.0])
    if abs(float(np.dot(unit_vector, reference))) > 0.9:
        reference = np.asarray([1.0, 0.0, 0.0])
    first = np.cross(unit_vector, reference)
    first /= np.linalg.norm(first)
    second = np.cross(unit_vector, first)
    return np.column_stack((first, second))


def perturbed_positions(response: np.ndarray, state: dict) -> np.ndarray:
    nominal_response = state["nominal_response"]
    nominal_position = state["points_km"]
    antpos = state["antpos_m"]
    position = nominal_position.copy()
    wave_number = 2.0 * np.pi / pc.wavelength
    for pulse in range(len(position)):
        radius = float(np.linalg.norm(nominal_position[pulse]))
        direction = nominal_position[pulse] / radius
        basis = tangent_basis(direction)
        phase_change = np.angle(response[pulse] * np.conj(nominal_response[pulse]))
        weight = np.abs(nominal_response[pulse]) ** 2
        good = np.isfinite(phase_change) & np.isfinite(weight) & (weight > 0.0)
        if np.count_nonzero(good) < 4:
            continue
        design = np.column_stack(
            (
                np.ones(np.count_nonzero(good)),
                -wave_number * antpos[good] @ basis,
            )
        )
        sqrt_weight = np.sqrt(weight[good] / np.nanmax(weight[good]))
        solution = np.linalg.lstsq(
            sqrt_weight[:, None] * design,
            sqrt_weight * phase_change[good],
            rcond=None,
        )[0]
        changed_direction = direction + basis @ solution[1:]
        changed_direction /= np.linalg.norm(changed_direction)
        position[pulse] = radius * changed_direction
    return position


def decoded_response(z_rx: np.ndarray, state: dict) -> np.ndarray:
    decoded = state["decoded"]
    response = np.full_like(decoded["response"], np.nan + 1j * np.nan)
    n_fast = decoded["z_tx"].shape[1]
    fast_time = decoded["fast_time_s"]
    for observation, raw_index in enumerate(decoded["raw_idx"]):
        template = decoded["z_tx"][raw_index].astype(np.complex128)
        envelope = np.abs(template) ** 2
        use = envelope > 0.05 * np.nanmax(envelope)
        frequency_hz = 2.0 * decoded["coarse_doppler_mps"][observation] / pc.wavelength
        derotation = np.exp(-1j * 2.0 * np.pi * frequency_hz * fast_time)
        for channel in range(z_rx.shape[1]):
            echo = fractional_segment(
                z_rx[raw_index, channel], decoded["range_gate"][observation], n_fast
            )
            response[observation, channel] = np.mean(
                echo[use] * np.conj(template[use]) * derotation[use]
            )
    return response


def beamformed_echoes(
    z_rx: np.ndarray, state: dict, response: np.ndarray | None = None
) -> np.ndarray:
    """Apply stored phase calibration and fixed first-stage AoA steering."""
    decoded = state["decoded"]
    if response is None:
        response = decoded["response"]
    module_gain = estimate_event_module_voltage_gain(response)
    return np.stack(
        [
            coherent_add_modules(
                z_rx[raw_index],
                state["steering_uvw"][observation],
                int(decoded["beam_id"][observation]),
                module_noise_variance=state["noise_variance"][raw_index],
                module_voltage_gain=module_gain,
                normalization="weighted_mean",
            )
            for observation, raw_index in enumerate(decoded["raw_idx"])
        ]
    )


def fit_triplets(
    z_rx: np.ndarray, state: dict, response: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    decoded = state["decoded"]
    beamformed = beamformed_echoes(z_rx, state, response=response)
    velocities = np.full(len(state["triplet_specs"]), np.nan)
    accelerations = np.full(len(state["triplet_specs"]), np.nan)
    n_fast = decoded["z_tx"].shape[1]
    for row, spec in enumerate(state["triplet_specs"]):
        observations = np.asarray(spec["observations"], dtype=int)
        raw_indices = np.asarray(decoded["raw_idx"][observations], dtype=int)
        raw_pulses = np.stack(
            [
                fractional_segment(
                    beamformed[observation], decoded["range_gate"][observation], n_fast
                )
                for raw_index, observation in zip(raw_indices, observations)
            ]
        )
        templates = decoded["z_tx"][raw_indices].astype(np.complex128)
        amplitudes = rms_baud_amplitudes(raw_pulses, templates)
        try:
            fit = fit_three_pulse_envelope(
                raw_pulses,
                templates,
                PAIR_SPACING_S,
                spec["frequency_guess_hz"],
                spec["acceleration_guess_mps2"],
                pc.wavelength,
                acceleration_half_width_mps2=2.0e4,
                frequency_half_width_hz=1.0e3,
                pulse_snr=decoded["snr"][observations],
                matched_filter_amplitudes=amplitudes,
                initial_parameters=spec["parameters"],
                max_nfev=160,
            )
        except Exception:
            continue
        velocities[row] = fit["parameters"][4] * pc.wavelength / 2.0
        accelerations[row] = fit["parameters"][5]
    return velocities, accelerations


def bootstrap_one(seed: int) -> np.ndarray:
    state = _STATE
    rng = np.random.default_rng(seed)
    noise = (
        rng.standard_normal(state["z_rx"].shape)
        + 1j * rng.standard_normal(state["z_rx"].shape)
    ) * np.sqrt(state["noise_variance"][:, :, None] / 2.0)
    z_rx = state["z_rx"] + noise
    response = decoded_response(z_rx, state)
    position = perturbed_positions(response, state)
    velocity, acceleration = fit_triplets(z_rx, state, response=response)
    return np.r_[
        position[state["echo_keep"]].ravel(),
        velocity[state["velocity_keep"]],
        acceleration[state["acceleration_keep"]],
    ]


def initialize_worker(state):
    global _STATE
    _STATE = state


def oas_correlation(samples: np.ndarray) -> tuple[np.ndarray, float]:
    centered = samples - np.mean(samples, axis=0)
    scale = np.std(centered, axis=0, ddof=1)
    scale = np.maximum(scale, np.nanmedian(scale[scale > 0.0]) * 1e-6)
    standardized = centered / scale
    empirical = standardized.T @ standardized / len(standardized)
    dimension = empirical.shape[0]
    mu = float(np.trace(empirical) / dimension)
    alpha = float(np.mean(empirical**2))
    denominator = (len(samples) + 1.0) * (alpha - mu**2 / dimension)
    shrinkage = 1.0 if denominator <= 0.0 else min((alpha + mu**2) / denominator, 1.0)
    correlation = (1.0 - shrinkage) * empirical + shrinkage * mu * np.eye(dimension)
    diagonal = np.sqrt(np.diag(correlation))
    correlation /= diagonal[:, None] * diagonal[None, :]
    return correlation, float(shrinkage)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-idx", type=int, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--diagnostics-h5", type=Path, required=True)
    parser.add_argument("--initial-fit-h5", type=Path, required=True)
    parser.add_argument("--prior-profile-h5", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=256)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260719)
    args = parser.parse_args()

    cut = load_cut(args.base / "metadata/cut", args.sample_idx)
    hypothesis = load_selected(args.diagnostics_h5)
    with h5py.File(args.initial_fit_h5, "r") as handle:
        precise_range_km = np.asarray(handle["range_km"], dtype=float)
        precise_doppler_mps = np.asarray(handle["doppler_mps"], dtype=float)
    decoded = decoded_pulse_responses(
        cut,
        hypothesis,
        7.0,
        precise_range_km=precise_range_km,
        precise_doppler_mps=precise_doppler_mps,
    )
    pairs = baud_averaged_beat_pairs(decoded, same_beam=True)
    triplets = valid_triplets(pairs)
    observations, stored_fit = load_selected_fit(args.diagnostics_h5)
    points_km = np.asarray(observations["points_km"], dtype=float)
    steering_points_km = np.asarray(stored_fit["model_km"], dtype=float)
    steering_uvw = enu_line_of_sight_to_arrival_uvw(steering_points_km)
    trajectory_time = np.asarray(observations["t_s"], dtype=float)
    absolute_zero = decoded["tx_idx"][0] / FS_HZ - trajectory_time[0]
    prior_acceleration = phase_acceleration_lookup(args.prior_profile_h5)
    with h5py.File(args.prior_profile_h5, "r") as handle:
        echo_keep = np.asarray(handle["result/echo_shared_inlier_mask"], dtype=bool)
        delta = np.asarray(handle["profile/delta_chi2"], dtype=float)
        best_index = int(np.nanargmin(delta))
        parameters6 = np.asarray(handle["profile/parameters6"][best_index], dtype=float)
        radius_um = float(handle["profile/radius_um"][best_index])
    import fit_best_alias_physics_models as physics

    density, _metadata = physics.pbal.density_interpolator(observations["sample_epoch_unix"])
    position, velocity, _radius, _mass, success, message = physics.propagate_shrinking_radius_model(
        np.r_[parameters6, np.log10(radius_um * 1e-6)], trajectory_time, density
    )
    if not success:
        raise RuntimeError(message)
    model_doppler_mps = physics.predicted_doppler(position, velocity) * 1e3
    model_acceleration_mps2 = np.gradient(model_doppler_mps, trajectory_time)

    specs = []
    nominal_velocity = []
    nominal_acceleration = []
    shared = []
    at_bound = []
    n_fast = decoded["z_tx"].shape[1]
    z_rx = np.asarray(decoded["z_rx"], dtype=np.complex128)
    module_noise_variance = event_receiver_noise_variance(
        z_rx,
        decoded["raw_idx"],
        decoded["range_gate"],
        decoded["z_tx"].shape[1],
    )
    nominal_beamformed = beamformed_echoes(
        z_rx,
        {
            "decoded": decoded,
            "steering_uvw": steering_uvw,
            "noise_variance": module_noise_variance,
        },
    )
    for previous, middle, current in triplets:
        indices = np.asarray([previous, middle, current], dtype=int)
        raw_indices = np.asarray(decoded["raw_idx"][indices], dtype=int)
        raw_pulses = np.stack(
            [
                fractional_segment(
                    nominal_beamformed[index], decoded["range_gate"][index], n_fast
                )
                for raw_index, index in zip(raw_indices, indices)
            ]
        )
        templates = decoded["z_tx"][raw_indices].astype(np.complex128)
        amplitudes = rms_baud_amplitudes(raw_pulses, templates)
        support = (previous, middle, middle, current)
        local = prior_acceleration.get(support)
        fit_time = decoded["tx_idx"][middle] / FS_HZ - absolute_zero
        acceleration_guess = (
            local["model_mps2"]
            if local is not None
            else float(np.interp(fit_time, trajectory_time, model_acceleration_mps2))
        )
        velocity_guess = float(np.interp(fit_time, trajectory_time, model_doppler_mps))
        fit = fit_three_pulse_envelope(
            raw_pulses,
            templates,
            PAIR_SPACING_S,
            2.0 * velocity_guess / pc.wavelength,
            acceleration_guess,
            pc.wavelength,
            acceleration_half_width_mps2=2.0e4,
            frequency_half_width_hz=1.0e3,
            pulse_snr=decoded["snr"][indices],
            matched_filter_amplitudes=amplitudes,
        )
        specs.append(
            {
                "observations": (previous, middle, current),
                "channel": -1,
                "frequency_guess_hz": 2.0 * velocity_guess / pc.wavelength,
                "acceleration_guess_mps2": acceleration_guess,
                "parameters": np.asarray(fit["parameters"], dtype=float),
            }
        )
        nominal_velocity.append(fit["parameters"][4] * pc.wavelength / 2.0)
        nominal_acceleration.append(fit["parameters"][5])
        shared.append(bool(np.all(echo_keep[indices])) and (local is None or local["keep"]))
        at_bound.append(abs(fit["parameters"][5] - acceleration_guess) > 1.99e4)

    velocity_keep = np.asarray(shared, dtype=bool)
    acceleration_keep = velocity_keep & ~np.asarray(at_bound, dtype=bool)
    state = {
        "z_rx": z_rx,
        "noise_variance": module_noise_variance,
        "decoded": decoded,
        "steering_uvw": steering_uvw,
        "triplet_specs": specs,
        "nominal_response": np.asarray(decoded["response"], dtype=np.complex128),
        "points_km": points_km,
        "antpos_m": np.asarray(pint.get_antpos(), dtype=float),
        "echo_keep": echo_keep,
        "velocity_keep": velocity_keep,
        "acceleration_keep": acceleration_keep,
    }
    seeds = np.random.SeedSequence(args.seed).generate_state(args.replicates)
    context = mp.get_context("fork")
    with context.Pool(args.workers, initializer=initialize_worker, initargs=(state,)) as pool:
        rows = list(pool.imap_unordered(bootstrap_one, seeds, chunksize=1))
    samples = np.asarray(rows, dtype=float)
    finite = np.all(np.isfinite(samples), axis=1)
    samples = samples[finite]
    if len(samples) < 32:
        raise RuntimeError(f"only {len(samples)} complete bootstrap replicates")
    covariance = np.cov(samples, rowvar=False, ddof=1)
    correlation, shrinkage = oas_correlation(samples)
    nominal_vector = np.r_[
        points_km[echo_keep].ravel(),
        np.asarray(nominal_velocity)[velocity_keep],
        np.asarray(nominal_acceleration)[acceleration_keep],
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output, "w") as handle:
        handle.attrs["schema"] = "pansy.joint_event_measurement_covariance.v1"
        handle.attrs["sample_idx"] = args.sample_idx
        handle.attrs["requested_replicates"] = args.replicates
        handle.attrs["complete_replicates"] = len(samples)
        handle.attrs["oas_correlation_shrinkage"] = shrinkage
        handle.attrs["voltage_source"] = (
            "AoA-steered inverse-noise coherent sum of seven receiver modules"
        )
        handle.attrs["steering_direction_source"] = (
            "selected first-stage trajectory model"
        )
        handle.attrs["receiver_phase_calibration"] = "stored calibration applied"
        handle.attrs["static_amplitude_calibration"] = "mesocal amp_scale applied"
        handle.attrs["event_amplitude_calibration"] = (
            "robust per-module voltage gain re-estimated in each bootstrap replicate"
        )
        handle.attrs["steering_convention"] = "catalogue_arrival_uvw"
        handle.attrs["measurement_order"] = "position_enu_km[echo_keep].ravel, triplet_velocity_mps[velocity_keep], triplet_acceleration_mps2[acceleration_keep]"
        handle.create_dataset("bootstrap_samples", data=samples, compression="gzip")
        handle.create_dataset("nominal_measurement", data=nominal_vector)
        handle.create_dataset("covariance", data=covariance)
        handle.create_dataset("correlation", data=correlation)
        handle.create_dataset("echo_keep", data=echo_keep)
        handle.create_dataset("velocity_keep", data=velocity_keep)
        handle.create_dataset("acceleration_keep", data=acceleration_keep)
        handle.create_dataset("triplet_observation_indices", data=np.asarray([s["observations"] for s in specs]))
        handle.create_dataset("triplet_channels", data=np.asarray([s["channel"] for s in specs]))
    print(
        f"replicates={len(samples)} dimension={samples.shape[1]} "
        f"shrinkage={shrinkage:.4f} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
