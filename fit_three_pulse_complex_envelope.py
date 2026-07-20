#!/usr/bin/env python3
"""Fit a phase-continuous accelerating sinusoid to three same-beam pulses."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

import pansy_config as pc


FS_HZ = 1e6
PAIR_SPACING_S = 0.008
PAIR_SPACING_TOLERANCE_S = 5e-6


def complex_envelope_model(
    parameters: np.ndarray,
    time_s: np.ndarray,
    envelope: np.ndarray,
    pulse_index: np.ndarray,
    wavelength_m: float,
) -> np.ndarray:
    """Evaluate baud amplitudes with one phase-continuous quadratic sinusoid."""
    n_pulses = int(np.max(pulse_index)) + 1
    log_amplitude = np.asarray(parameters[:n_pulses], dtype=float)
    phase0, frequency_hz, acceleration_mps2 = np.asarray(
        parameters[n_pulses : n_pulses + 3], dtype=float
    )
    amplitude = np.exp(log_amplitude[np.asarray(pulse_index, dtype=int)])
    frequency_rate_hz_s = 2.0 * acceleration_mps2 / float(wavelength_m)
    phase = phase0 + 2.0 * np.pi * (
        frequency_hz * time_s + 0.5 * frequency_rate_hz_s * time_s**2
    )
    return amplitude * envelope * np.exp(1j * phase)


def fractional_delay_fft(values: np.ndarray, delay_samples: float) -> np.ndarray:
    """Shift a complex pulse template by a fractional sample using a phase ramp."""
    values = np.asarray(values, dtype=np.complex128)
    frequency = np.fft.fftfreq(len(values))
    phase = np.exp(-1j * 2.0 * np.pi * frequency * float(delay_samples))
    return np.fft.ifft(np.fft.fft(values) * phase)


def envelope_segments(template: np.ndarray, threshold: float = 0.05) -> list[np.ndarray]:
    power = np.abs(np.asarray(template)) ** 2
    power /= max(float(np.nanmax(power)), 1e-30)
    indices = np.flatnonzero(power > threshold)
    if len(indices) == 0:
        return []
    breaks = np.flatnonzero(np.diff(indices) > 1) + 1
    return [segment for segment in np.split(indices, breaks) if len(segment) >= 2]


def baud_measurements(
    raw_pulses: np.ndarray,
    pulse_templates: np.ndarray,
    range_offsets_samples: np.ndarray,
    pulse_spacing_s: float,
    pulse_snr: np.ndarray,
) -> dict:
    """Range-align, decode, and form one SNR-weighted complex value per baud."""
    times = []
    values = []
    envelopes = []
    weights = []
    pulse_indices = []
    for pulse in range(len(raw_pulses)):
        aligned = fractional_delay_fft(raw_pulses[pulse], -float(range_offsets_samples[pulse]))
        template = pulse_templates[pulse]
        power = np.abs(template) ** 2
        decoded = aligned * np.conj(template)
        for segment in envelope_segments(template):
            sample_weight = power[segment]
            denominator = float(np.sum(sample_weight))
            if not np.isfinite(denominator) or denominator <= 0.0:
                continue
            center = float(np.sum(segment * sample_weight) / denominator) / FS_HZ
            times.append(pulse * pulse_spacing_s + center)
            values.append(np.sum(sample_weight * decoded[segment]) / denominator)
            envelopes.append(float(np.sum(sample_weight * power[segment]) / denominator))
            weights.append(denominator * max(float(pulse_snr[pulse]), 1e-3))
            pulse_indices.append(pulse)
    return {
        "absolute_time_s": np.asarray(times, dtype=float),
        "data": np.asarray(values, dtype=np.complex128),
        "envelope": np.asarray(envelopes, dtype=float),
        "weight": np.asarray(weights, dtype=float),
        "pulse_index": np.asarray(pulse_indices, dtype=int),
    }


def rms_baud_amplitudes(raw_pulses: np.ndarray, pulse_templates: np.ndarray) -> np.ndarray:
    """Estimate one pulse amplitude prior from decoded baud-voltage RMS."""
    bauds = baud_measurements(
        raw_pulses,
        pulse_templates,
        np.zeros(len(raw_pulses)),
        PAIR_SPACING_S,
        np.ones(len(raw_pulses)),
    )
    amplitudes = np.full(len(raw_pulses), np.nan)
    for pulse in range(len(raw_pulses)):
        select = bauds["pulse_index"] == pulse
        values = bauds["data"][select]
        amplitudes[pulse] = np.sqrt(np.nanmean(np.abs(values) ** 2))
    return amplitudes


def common_bin_pulse_pair_doppler(
    raw_pulses: np.ndarray,
    pulse_templates: np.ndarray,
    pulse_spacing_s: float,
    frequency_prior_hz: float,
    wavelength_m: float,
    pulse_snr: np.ndarray | None = None,
    zero_pad_factor: int = 16,
) -> dict:
    """Measure pulse-pair Doppler phase at one shared baud-domain FFT bin."""
    raw_pulses = np.asarray(raw_pulses, dtype=np.complex128)
    pulse_templates = np.asarray(pulse_templates, dtype=np.complex128)
    if raw_pulses.ndim != 2 or raw_pulses.shape[0] != 2:
        raise ValueError("raw_pulses must have shape (2, n_fast)")
    if pulse_templates.shape != raw_pulses.shape:
        raise ValueError("pulse_templates must match raw_pulses")
    if pulse_snr is None:
        pulse_snr = np.ones(2, dtype=float)
    if int(zero_pad_factor) < 4:
        raise ValueError("zero_pad_factor must be at least four")

    bauds = baud_measurements(
        raw_pulses,
        pulse_templates,
        np.zeros(2),
        float(pulse_spacing_s),
        np.asarray(pulse_snr, dtype=float),
    )
    pulse_measurements = []
    baud_count = []
    local_times = []
    for pulse in range(2):
        select = bauds["pulse_index"] == pulse
        use = (
            select
            & np.isfinite(bauds["data"].real)
            & np.isfinite(bauds["data"].imag)
            & np.isfinite(bauds["envelope"])
            & np.isfinite(bauds["weight"])
            & (bauds["envelope"] > 0.0)
            & (bauds["weight"] > 0.0)
        )
        if np.count_nonzero(use) < 4:
            raise RuntimeError("too few finite baud averages for pulse-pair Doppler")
        pulse_measurements.append(use)
        baud_count.append(int(np.count_nonzero(use)))
        local_times.append(bauds["absolute_time_s"][use] - pulse * pulse_spacing_s)

    all_local_time = np.unique(np.concatenate(local_times))
    differences = np.diff(all_local_time)
    differences = differences[np.isfinite(differences) & (differences > 0.0)]
    if len(differences) == 0:
        raise RuntimeError("baud times do not define a Doppler grid")
    baud_spacing_s = float(np.median(differences))
    n_baud = max(baud_count)
    nfft = int(2 ** np.ceil(np.log2(max(int(zero_pad_factor) * n_baud, 2))))
    frequency_step_hz = 1.0 / (nfft * baud_spacing_s)
    common_frequency_hz = (
        np.rint(float(frequency_prior_hz) / frequency_step_hz) * frequency_step_hz
    )

    responses = []
    response_norm = []
    for pulse, use in enumerate(pulse_measurements):
        local_time = bauds["absolute_time_s"][use] - pulse * pulse_spacing_s
        basis = bauds["envelope"][use] * np.exp(
            1j * 2.0 * np.pi * common_frequency_hz * local_time
        )
        weight = bauds["weight"][use]
        numerator = np.sum(weight * bauds["data"][use] * np.conj(basis))
        denominator = np.sum(weight * np.abs(basis) ** 2)
        response = numerator / max(float(denominator), np.finfo(float).tiny)
        responses.append(response)
        response_norm.append(
            abs(numerator)
            / max(
                float(np.sum(weight * np.abs(bauds["data"][use]) * np.abs(basis))),
                np.finfo(float).tiny,
            )
        )
    responses = np.asarray(responses, dtype=np.complex128)
    phase_rad = float(np.angle(responses[1] * np.conj(responses[0])))
    wrapped_frequency_hz = phase_rad / (2.0 * np.pi * float(pulse_spacing_s))
    frequency_ambiguity_hz = 1.0 / float(pulse_spacing_s)
    resolved_frequency_hz = wrapped_frequency_hz + np.rint(
        (float(frequency_prior_hz) - wrapped_frequency_hz) / frequency_ambiguity_hz
    ) * frequency_ambiguity_hz
    return {
        "phase_rad": phase_rad,
        "wrapped_frequency_hz": wrapped_frequency_hz,
        "frequency_ambiguity_hz": frequency_ambiguity_hz,
        "wrapped_velocity_mps": wrapped_frequency_hz * float(wavelength_m) / 2.0,
        "velocity_ambiguity_mps": frequency_ambiguity_hz * float(wavelength_m) / 2.0,
        "resolved_frequency_hz": resolved_frequency_hz,
        "resolved_velocity_mps": resolved_frequency_hz * float(wavelength_m) / 2.0,
        "common_frequency_hz": common_frequency_hz,
        "frequency_step_hz": frequency_step_hz,
        "nfft": nfft,
        "baud_count": np.asarray(baud_count, dtype=int),
        "response": responses,
        "coherence": float(min(response_norm)),
    }


def fit_three_pulse_envelope(
    raw_pulses: np.ndarray,
    pulse_templates: np.ndarray,
    pulse_spacing_s: float,
    frequency_guess_hz: float,
    acceleration_guess_mps2: float,
    wavelength_m: float,
    acceleration_half_width_mps2: float = 2.0e4,
    frequency_half_width_hz: float = 1.5e4,
    pulse_snr: np.ndarray | None = None,
    matched_filter_amplitudes: np.ndarray | None = None,
    fixed_acceleration_mps2: float | None = None,
    acceleration_phase_rad: float | None = None,
    acceleration_phase_std_rad: float | None = None,
    initial_parameters: np.ndarray | None = None,
    max_nfev: int = 2000,
) -> dict:
    """Jointly fit three amplitudes, common phase, frequency, and acceleration."""
    raw_pulses = np.asarray(raw_pulses, dtype=np.complex128)
    pulse_templates = np.asarray(pulse_templates, dtype=np.complex128)
    if raw_pulses.ndim != 2 or raw_pulses.shape[0] != 3:
        raise ValueError("raw_pulses must have shape (3, n_fast)")
    if pulse_templates.shape != raw_pulses.shape:
        raise ValueError("pulse_templates must match raw_pulses")

    if pulse_snr is None:
        pulse_snr = np.ones(3, dtype=float)
    pulse_snr = np.asarray(pulse_snr, dtype=float)
    if pulse_snr.shape != (3,):
        raise ValueError("pulse_snr must contain one value for each pulse")
    if matched_filter_amplitudes is not None:
        matched_filter_amplitudes = np.asarray(matched_filter_amplitudes, dtype=float)
        if matched_filter_amplitudes.shape != (3,):
            raise ValueError("matched_filter_amplitudes must contain one value for each pulse")

    initial_bauds = baud_measurements(
        raw_pulses, pulse_templates, np.zeros(3), pulse_spacing_s, pulse_snr
    )
    absolute_time = initial_bauds["absolute_time_s"]
    data = initial_bauds["data"]
    envelope = initial_bauds["envelope"]
    weight = initial_bauds["weight"]
    pulse_index = initial_bauds["pulse_index"]
    use = (
        np.isfinite(data.real)
        & np.isfinite(data.imag)
        & np.isfinite(envelope)
        & np.isfinite(weight)
        & (envelope > 0.0)
        & (weight > 0.0)
    )
    if np.count_nonzero(use) < 8:
        raise RuntimeError("too few finite baud-averaged samples in the three-pulse envelope")

    time_center_s = float(np.sum(weight[use] * absolute_time[use]) / np.sum(weight[use]))
    time_s = absolute_time - time_center_s
    scale = float(np.nanmax(np.abs(data[use])))
    data = data / max(scale, 1e-30)
    envelope = envelope / max(float(np.nanmax(envelope[use])), 1e-30)
    weight = weight / max(float(np.nanmax(weight[use])), 1e-30)
    sqrt_weight = np.sqrt(weight[use])

    if initial_parameters is None:
        trial = np.asarray(
            [0.0, 0.0, 0.0, 0.0, frequency_guess_hz, acceleration_guess_mps2]
        )
        unit_model = complex_envelope_model(trial, time_s, envelope, pulse_index, wavelength_m)
        aligned_data = data
        projections = []
        for pulse in range(3):
            pulse_use = use & (pulse_index == pulse)
            numerator = np.sum(
                weight[pulse_use]
                * aligned_data[pulse_use]
                * np.conj(unit_model[pulse_use])
            )
            denominator = np.sum(weight[pulse_use] * np.abs(unit_model[pulse_use]) ** 2)
            projections.append(numerator / max(float(denominator), 1e-30))
        common_phase = float(np.angle(np.sum(projections)))
        projected_amplitudes = np.maximum(np.abs(projections), 1e-6)
        if matched_filter_amplitudes is not None and np.all(matched_filter_amplitudes > 0.0):
            geometric_mean = float(np.prod(projected_amplitudes) ** (1.0 / 3.0))
            matched_geometric_mean = float(np.prod(matched_filter_amplitudes) ** (1.0 / 3.0))
            projected_amplitudes = (
                geometric_mean * matched_filter_amplitudes / matched_geometric_mean
            )
        initial_parameters = np.asarray(
            [
                np.log(projected_amplitudes[0]),
                np.log(projected_amplitudes[1]),
                np.log(projected_amplitudes[2]),
                common_phase,
                frequency_guess_hz,
                acceleration_guess_mps2,
            ],
            dtype=float,
        )
    else:
        initial_parameters = np.asarray(initial_parameters, dtype=float).copy()
    starting_parameters = initial_parameters.copy()

    def residual_full(parameters: np.ndarray) -> np.ndarray:
        model = complex_envelope_model(parameters, time_s, envelope, pulse_index, wavelength_m)
        residual = sqrt_weight * (data[use] - model[use])
        components = [residual.real, residual.imag]
        if acceleration_phase_rad is not None:
            if acceleration_phase_std_rad is None or acceleration_phase_std_rad <= 0.0:
                raise ValueError(
                    "acceleration_phase_std_rad must be positive with a phase measurement"
                )
            predicted_phase = (
                4.0
                * np.pi
                * parameters[5]
                * float(pulse_spacing_s) ** 2
                / float(wavelength_m)
            )
            wrapped_residual = np.angle(
                np.exp(1j * (float(acceleration_phase_rad) - predicted_phase))
            )
            components.append(
                np.asarray([wrapped_residual / float(acceleration_phase_std_rad)])
            )
        return np.concatenate(components)

    lower = np.asarray(
        [
            max(-12.0, initial_parameters[0] + np.log(0.9)),
            max(-12.0, initial_parameters[1] + np.log(0.9)),
            max(-12.0, initial_parameters[2] + np.log(0.9)),
            -np.inf,
            frequency_guess_hz - frequency_half_width_hz,
            acceleration_guess_mps2 - acceleration_half_width_mps2,
        ]
    )
    upper = np.asarray(
        [
            min(5.0, initial_parameters[0] + np.log(1.1)),
            min(5.0, initial_parameters[1] + np.log(1.1)),
            min(5.0, initial_parameters[2] + np.log(1.1)),
            np.inf,
            frequency_guess_hz + frequency_half_width_hz,
            acceleration_guess_mps2 + acceleration_half_width_mps2,
        ]
    )
    if fixed_acceleration_mps2 is None:
        fit = least_squares(
            residual_full,
            np.clip(initial_parameters, lower + 1e-9, upper - 1e-9),
            bounds=(lower, upper),
            x_scale=np.asarray([1.0, 1.0, 1.0, 1.0, 2.0e3, 5.0e3]),
            max_nfev=int(max_nfev),
        )
        parameters = fit.x
        jacobian = fit.jac
    else:
        fixed_acceleration_mps2 = float(fixed_acceleration_mps2)

        def residual_nuisance(nuisance: np.ndarray) -> np.ndarray:
            parameters = np.insert(nuisance, 5, fixed_acceleration_mps2)
            return residual_full(parameters)

        nuisance_lower = np.delete(lower, 5)
        nuisance_upper = np.delete(upper, 5)
        nuisance_initial = np.clip(
            np.delete(initial_parameters, 5), nuisance_lower + 1e-9, nuisance_upper - 1e-9
        )
        fit = least_squares(
            residual_nuisance,
            nuisance_initial,
            bounds=(nuisance_lower, nuisance_upper),
            x_scale=np.asarray([1.0, 1.0, 1.0, 1.0, 2.0e3]),
            max_nfev=min(int(max_nfev), 1200),
        )
        parameters = np.insert(fit.x, 5, fixed_acceleration_mps2)
        jacobian = fit.jac

    final_data = data
    model = complex_envelope_model(parameters, time_s, envelope, pulse_index, wavelength_m)
    residual_vector = residual_full(parameters)
    weighted_sse = float(np.sum(residual_vector**2))
    degrees_of_freedom = max(len(residual_vector) - len(fit.x), 1)
    covariance = np.full((len(fit.x), len(fit.x)), np.nan)
    parameter_noise_influence = np.full((len(fit.x), len(residual_vector)), np.nan)
    try:
        inverse_information = np.linalg.inv(jacobian.T @ jacobian)
        residual_variance = weighted_sse / degrees_of_freedom
        covariance = inverse_information * residual_variance
        parameter_noise_influence = (
            inverse_information @ jacobian.T * np.sqrt(residual_variance)
        )
    except np.linalg.LinAlgError:
        pass
    velocity_mps = float(parameters[4] * wavelength_m / 2.0)
    velocity_acceleration_covariance = np.asarray(
        [
            [
                covariance[4, 4] * (wavelength_m / 2.0) ** 2,
                covariance[4, 5] * wavelength_m / 2.0,
            ],
            [
                covariance[5, 4] * wavelength_m / 2.0,
                covariance[5, 5],
            ],
        ]
    )
    velocity_acceleration_noise_influence = np.vstack(
        (
            parameter_noise_influence[4] * wavelength_m / 2.0,
            parameter_noise_influence[5],
        )
    )
    residual_pulse_index = np.r_[
        pulse_index[use],
        pulse_index[use],
        np.full(len(residual_vector) - 2 * np.count_nonzero(use), -1, dtype=int),
    ]
    return {
        "parameters": parameters,
        "starting_parameters": starting_parameters,
        "parameter_covariance": covariance,
        "velocity_mps": velocity_mps,
        "acceleration_mps2": float(parameters[5]),
        "velocity_acceleration_covariance": velocity_acceleration_covariance,
        "velocity_acceleration_noise_influence": velocity_acceleration_noise_influence,
        "residual_pulse_index": residual_pulse_index,
        "model": model,
        "data": final_data,
        "envelope": envelope,
        "weight": weight,
        "use": use,
        "time_s": time_s,
        "absolute_time_s": absolute_time,
        "pulse_index": pulse_index,
        "time_center_s": time_center_s,
        "scale": scale,
        "weighted_sse": weighted_sse,
        "degrees_of_freedom": degrees_of_freedom,
        "success": bool(fit.success),
        "message": str(fit.message),
    }


def fit_three_pulse_acceleration_aliases(
    raw_pulses: np.ndarray,
    pulse_templates: np.ndarray,
    pulse_spacing_s: float,
    frequency_guess_hz: float,
    acceleration_phase_guess_mps2: float,
    wavelength_m: float,
    *,
    acceleration_limits_mps2: tuple[float, float] = (-1.5e5, 1.5e5),
    acceleration_phase_rad: float | None = None,
    acceleration_phase_std_rad: float | None = None,
    pulse_snr: np.ndarray | None = None,
    matched_filter_amplitudes: np.ndarray | None = None,
    frequency_half_width_hz: float | None = None,
) -> dict:
    """Search the FFT-selected velocity alias and every phase-compatible acceleration.

    Three equally spaced pulses repeat their quadratic-phase information when
    acceleration changes by ``wavelength / (2 * spacing**2)``.  The established
    pulse-to-pulse phase progression pins narrow solutions in every acceleration
    alias cell.  The established FFT Doppler chooses the frequency alias, so the
    local frequency search is limited to half of its ``1 / pulse_spacing``
    ambiguity interval.  The complex-envelope objective then distinguishes the
    remaining acceleration cells.
    """
    lower_limit, upper_limit = map(float, acceleration_limits_mps2)
    if not lower_limit < upper_limit:
        raise ValueError("acceleration_limits_mps2 must be increasing")
    alias_spacing = float(wavelength_m) / (2.0 * float(pulse_spacing_s) ** 2)
    alias_half_width = 0.499 * alias_spacing
    frequency_alias_spacing_hz = 1.0 / float(pulse_spacing_s)
    if frequency_half_width_hz is None:
        frequency_half_width_hz = 0.499 * frequency_alias_spacing_hz
    minimum_alias = int(
        np.ceil((lower_limit - acceleration_phase_guess_mps2) / alias_spacing)
    )
    maximum_alias = int(
        np.floor((upper_limit - acceleration_phase_guess_mps2) / alias_spacing)
    )
    alias_number = np.arange(minimum_alias, maximum_alias + 1, dtype=int)
    alias_center = acceleration_phase_guess_mps2 + alias_number * alias_spacing
    fits = []
    for center in alias_center:
        local_half_width = min(
            alias_half_width,
            center - lower_limit if center > lower_limit else alias_half_width,
            upper_limit - center if center < upper_limit else alias_half_width,
        )
        local_half_width = max(float(local_half_width), 1.0)
        fits.append(
            fit_three_pulse_envelope(
                raw_pulses,
                pulse_templates,
                pulse_spacing_s,
                frequency_guess_hz,
                float(center),
                wavelength_m,
                acceleration_half_width_mps2=local_half_width,
                frequency_half_width_hz=frequency_half_width_hz,
                pulse_snr=pulse_snr,
                matched_filter_amplitudes=matched_filter_amplitudes,
                acceleration_phase_rad=acceleration_phase_rad,
                acceleration_phase_std_rad=acceleration_phase_std_rad,
            )
        )
    weighted_sse = np.asarray([fit["weighted_sse"] for fit in fits], dtype=float)
    selected_index = int(np.nanargmin(weighted_sse))
    residual_variance = weighted_sse[selected_index] / max(
        int(fits[selected_index]["degrees_of_freedom"]), 1
    )
    delta_chi2 = (weighted_sse - weighted_sse[selected_index]) / max(
        float(residual_variance), np.finfo(float).tiny
    )
    frequency_std_hz = np.asarray(
        [np.sqrt(fit["parameter_covariance"][4, 4]) for fit in fits], dtype=float
    )
    acceleration_std_mps2 = np.asarray(
        [np.sqrt(fit["parameter_covariance"][5, 5]) for fit in fits], dtype=float
    )
    return {
        "selected_fit": fits[selected_index],
        "selected_index": selected_index,
        "alias_number": alias_number,
        "alias_center_mps2": alias_center,
        "alias_spacing_mps2": alias_spacing,
        "fit_acceleration_mps2": np.asarray(
            [fit["parameters"][5] for fit in fits], dtype=float
        ),
        "fit_frequency_hz": np.asarray(
            [fit["parameters"][4] for fit in fits], dtype=float
        ),
        "fit_velocity_mps": np.asarray(
            [fit["parameters"][4] * wavelength_m / 2.0 for fit in fits], dtype=float
        ),
        "fit_velocity_std_mps": frequency_std_hz * wavelength_m / 2.0,
        "fit_acceleration_std_mps2": acceleration_std_mps2,
        "frequency_alias_spacing_hz": frequency_alias_spacing_hz,
        "velocity_alias_spacing_mps": wavelength_m / (2.0 * pulse_spacing_s),
        "frequency_half_width_hz": float(frequency_half_width_hz),
        "weighted_sse": weighted_sse,
        "delta_chi2": delta_chi2,
        "residual_variance": residual_variance,
        "fits": fits,
    }


def select_eight_ms_triplet(decoded: dict, pairs: dict) -> tuple[int, int, int, int]:
    delta_t = np.asarray(pairs["delta_t_s"], dtype=float)
    valid = np.isfinite(delta_t) & (np.abs(delta_t - PAIR_SPACING_S) <= PAIR_SPACING_TOLERANCE_S)
    if not np.any(valid):
        raise RuntimeError("no same-beam 8-ms pulse pair is available")
    pair_score = np.asarray(pairs["coherence"], dtype=float) / np.maximum(
        np.asarray(pairs["phase_doppler_std_mps"], dtype=float), 1e-6
    )
    pair_score[~valid] = -np.inf
    best = None
    for first_row in np.flatnonzero(valid):
        middle = int(pairs["current_index"][first_row])
        candidates = np.flatnonzero(valid & (np.asarray(pairs["previous_index"], dtype=int) == middle))
        for second_row in candidates:
            score = min(float(pair_score[first_row]), float(pair_score[second_row]))
            if best is None or score > best[0]:
                best = (score, int(first_row), int(second_row))
    if best is None:
        raise RuntimeError("no three-pulse same-beam 8-ms sequence is available")
    _, first_row, second_row = best
    previous = int(pairs["previous_index"][first_row])
    middle = int(pairs["current_index"][first_row])
    current = int(pairs["current_index"][second_row])
    response = decoded["response"]
    channel_score = (
        np.abs(response[middle] * np.conj(response[previous]))
        + np.abs(response[current] * np.conj(response[middle]))
    )
    channel = int(np.nanargmax(channel_score))
    return previous, middle, current, channel


def acceleration_profile(
    fit: dict,
    raw_pulses: np.ndarray,
    pulse_templates: np.ndarray,
    pulse_spacing_s: float,
    frequency_guess_hz: float,
    acceleration_guess_mps2: float,
    wavelength_m: float,
    points: int,
    half_width_mps2: float,
    pulse_snr: np.ndarray,
    matched_filter_amplitudes: np.ndarray,
) -> dict:
    acceleration = np.linspace(
        acceleration_guess_mps2 - half_width_mps2,
        acceleration_guess_mps2 + half_width_mps2,
        int(points),
    )
    objective = np.full(len(acceleration), np.nan)
    frequency = np.full(len(acceleration), np.nan)
    initial = fit["parameters"].copy()
    for index, value in enumerate(acceleration):
        profile_initial = initial.copy()
        profile_initial[:3] = fit["starting_parameters"][:3]
        fixed = fit_three_pulse_envelope(
            raw_pulses,
            pulse_templates,
            pulse_spacing_s,
            frequency_guess_hz,
            acceleration_guess_mps2,
            wavelength_m,
            acceleration_half_width_mps2=half_width_mps2,
            pulse_snr=pulse_snr,
            matched_filter_amplitudes=matched_filter_amplitudes,
            fixed_acceleration_mps2=value,
            initial_parameters=profile_initial,
        )
        objective[index] = fixed["weighted_sse"]
        frequency[index] = fixed["parameters"][4]
        initial = fixed["parameters"]
    residual_variance = fit["weighted_sse"] / fit["degrees_of_freedom"]
    delta_chi2 = (objective - np.nanmin(objective)) / max(float(residual_variance), 1e-30)
    return {
        "acceleration_mps2": acceleration,
        "frequency_hz": frequency,
        "weighted_sse": objective,
        "delta_chi2": delta_chi2,
    }


def plot_result(
    sample_idx: int,
    fit: dict,
    conditioned_fit: dict,
    profile: dict,
    metadata: dict,
    output: Path,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.5), sharey=True, constrained_layout=True)
    pulse_index = fit["pulse_index"]
    use = fit["use"]
    local_time_us = 1e6 * (
        fit["absolute_time_s"] - pulse_index * metadata["pulse_spacing_s"]
    )
    for column, pulse in enumerate((0, 1, 2)):
        ax = axes[column]
        ax.set_title((r"Pulse $n$", r"Pulse $n+1$", r"Pulse $n+2$")[pulse])
        for component, color, label in (
            (np.real, "C0", "real"),
            (np.imag, "C1", "imaginary"),
        ):
            select = use & (pulse_index == pulse)
            ax.scatter(
                local_time_us[select],
                component(fit["data"][select]),
                s=18,
                facecolor="none",
                edgecolor=color,
                alpha=0.8,
                label=f"{label} measurement",
                zorder=3,
            )
            order = np.argsort(local_time_us[select])
            ax.plot(
                local_time_us[select][order],
                component(fit["model"][select])[order],
                color=color,
                lw=1.5,
                label=f"{label} model",
            )
        ax.set_xlabel("Fast time within pulse (us)")
        ax.grid(alpha=0.2, lw=0.5)
        ax.legend(frameon=False, loc="best")
        if pulse == 2:
            ax.text(
                0.03,
                0.03,
                rf"$a_{{\rm fit}}={fit['parameters'][5] / 1e3:.2f}$ km s$^{{-2}}$"
                + "\n"
                + rf"$v_0={fit['parameters'][4] * pc.wavelength / 2e3:.3f}$ km s$^{{-1}}$"
                + "\n"
                + rf"$a_{{\rm phase,0}}={metadata['acceleration_guess_mps2'] / 1e3:.2f}$ km s$^{{-2}}$",
                transform=ax.transAxes,
                ha="left",
                va="bottom",
            )
    axes[0].set_ylabel("Normalized baud-averaged complex voltage")

    timestamp = dt.datetime.fromtimestamp(sample_idx / 1e6, tz=dt.timezone.utc).isoformat(timespec="milliseconds")
    fig.suptitle(
        f"Three-pulse complex-envelope fit, {timestamp}; "
        f"beam {metadata['beam_id']}, channel {metadata['channel']}, "
        f"separation {1e3 * metadata['pulse_spacing_s']:.3f} ms"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=190)
    plt.close(fig)


def write_result(
    path: Path,
    sample_idx: int,
    fit: dict,
    conditioned_fit: dict,
    profile: dict,
    metadata: dict,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as handle:
        handle.attrs["schema"] = "pansy.three_pulse_complex_envelope.v1"
        handle.attrs["sample_idx"] = int(sample_idx)
        handle.attrs["model"] = (
            "separate positive pulse amplitudes; shared phase; continuous quadratic phase; "
            "range gates fixed to established matched-filter estimates"
        )
        for name, value in metadata.items():
            handle.attrs[name] = value
        handle.create_dataset("parameters", data=fit["parameters"])
        handle.create_dataset("starting_parameters", data=fit["starting_parameters"])
        handle["parameters"].attrs["names"] = np.asarray(
            [
                "log_amplitude_0",
                "log_amplitude_1",
                "log_amplitude_2",
                "phase0_rad",
                "frequency_hz",
                "acceleration_mps2",
            ],
            dtype=h5py.string_dtype("utf-8"),
        )
        handle.create_dataset("parameter_covariance", data=fit["parameter_covariance"])
        handle.create_dataset("velocity_mps", data=fit["velocity_mps"])
        handle.create_dataset("acceleration_mps2", data=fit["acceleration_mps2"])
        handle.create_dataset(
            "velocity_acceleration_covariance",
            data=fit["velocity_acceleration_covariance"],
        )
        for name in ("data", "model", "envelope", "weight", "use", "time_s", "pulse_index"):
            handle.create_dataset(name, data=fit[name])
        conditioned = handle.create_group("conditioned_on_beat_phase_acceleration")
        conditioned.create_dataset("parameters", data=conditioned_fit["parameters"])
        conditioned.create_dataset("parameter_covariance", data=conditioned_fit["parameter_covariance"])
        conditioned.attrs["fixed_acceleration_mps2"] = float(conditioned_fit["parameters"][5])
        profile_group = handle.create_group("acceleration_profile")
        for name, value in profile.items():
            profile_group.create_dataset(name, data=value)


def main() -> int:
    from interferometer_alias_diagnostics import load_cut
    from inter_pulse_phase_deceleration import (
        baud_averaged_beat_pairs,
        decoded_pulse_responses,
        diagnostic_path,
        fractional_segment,
        load_selected,
        precise_matched_filter_estimates,
        robust_line,
    )

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-idx", type=int, required=True)
    parser.add_argument("--base", type=Path, default=Path("/mnt/data/juha/pansy"))
    parser.add_argument("--output-dir", type=Path, default=Path("test_plots/three_pulse_complex_envelope"))
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    parser.add_argument("--profile-points", type=int, default=161)
    parser.add_argument("--acceleration-half-width-mps2", type=float, default=2.0e4)
    parser.add_argument(
        "--initial-fit-h5",
        type=Path,
        help="Established high-resolution range/Doppler product used to initialize the fit",
    )
    parser.add_argument(
        "--prior-profile-h5",
        type=Path,
        help="Established phase-aware profile providing the local triplet acceleration seed",
    )
    args = parser.parse_args()

    cut = load_cut(args.base / "metadata/cut", args.sample_idx)
    selected_path = diagnostic_path(args.base, args.sample_idx)
    if selected_path.exists():
        import fit_best_alias_physics_models as physics

        hypothesis = load_selected(selected_path)
        precise = precise_matched_filter_estimates(cut, hypothesis, args.snr_threshold)
        decoded = decoded_pulse_responses(
            cut,
            hypothesis,
            args.snr_threshold,
            precise_range_km=precise["range_km"],
            precise_doppler_mps=precise["doppler_mps"],
        )
        seed_source = "trajectory model"
    else:
        initial_fit_h5 = args.initial_fit_h5
        if initial_fit_h5 is None:
            initial_fit_h5 = (
                Path("data/server_initial_fits")
                / f"highres_fft_i2_p16_{args.sample_idx}.h5"
            )
        if not initial_fit_h5.exists():
            raise FileNotFoundError(
                "the selected-fit diagnostic is absent and no established "
                f"high-resolution initialization exists at {initial_fit_h5}"
            )
        with h5py.File(initial_fit_h5, "r") as initial_handle:
            raw_idx = np.asarray(initial_handle["raw_idx"], dtype=int)
            precise_range_km = np.asarray(initial_handle["range_km"], dtype=float)
            precise_doppler_mps = np.asarray(initial_handle["doppler_mps"], dtype=float)
        raw_tx_idx = np.asarray(cut["tx_idx"], dtype=float)
        hypothesis = {
            "selected_range_time_component": -1,
            "t_rel_s": (raw_tx_idx[raw_idx] - raw_tx_idx[raw_idx[0]]) / FS_HZ,
        }
        decoded = decoded_pulse_responses(
            cut,
            hypothesis,
            args.snr_threshold,
            precise_range_km=precise_range_km,
            precise_doppler_mps=precise_doppler_mps,
        )
        seed_source = f"established high-resolution product {initial_fit_h5}"
    pairs = baud_averaged_beat_pairs(decoded, same_beam=True)
    pair_weight = np.asarray(pairs["coherence"], dtype=float) / np.maximum(
        np.asarray(pairs["phase_doppler_std_mps"], dtype=float), 1.0
    ) ** 2
    phase_fit = robust_line(pairs["time_s"], pairs["phase_doppler_mps"], pair_weight)
    global_acceleration_guess = float(phase_fit["acceleration_mps2"])
    previous, middle, current, channel = select_eight_ms_triplet(decoded, pairs)
    pulse_indices = np.asarray([previous, middle, current], dtype=int)
    pulse_deltas_s = np.diff(decoded["tx_idx"][pulse_indices]) / FS_HZ
    if np.max(np.abs(pulse_deltas_s - PAIR_SPACING_S)) > PAIR_SPACING_TOLERANCE_S:
        raise RuntimeError("selected three-pulse sequence is not uniformly spaced by 8 ms")
    pulse_spacing_s = float(np.mean(pulse_deltas_s))
    acceleration_guess = global_acceleration_guess
    acceleration_seed_source = "event-wide robust beat-phase line"
    prior_profile_h5 = args.prior_profile_h5
    if prior_profile_h5 is None:
        prior_profile_h5 = (
            Path("data/server_initial_fits")
            / f"mass_profile_phase_aware_{args.sample_idx}.h5"
        )
    if prior_profile_h5.exists():
        with h5py.File(prior_profile_h5, "r") as prior_handle:
            phase_group = prior_handle["phase_acceleration"]
            support = np.asarray(phase_group["support_observation_indices"], dtype=int)
            target_support = np.asarray([previous, middle, middle, current], dtype=int)
            matches = np.flatnonzero(np.all(support == target_support[None, :], axis=1))
            if len(matches):
                acceleration_guess = float(
                    phase_group["best_model_radial_acceleration_mps2"][int(matches[0])]
                )
                acceleration_seed_source = (
                    f"local phase-aware production model {prior_profile_h5}"
                )
    raw_indices = np.asarray(decoded["raw_idx"][pulse_indices], dtype=int)
    n_fast = decoded["z_tx"].shape[1]
    raw_pulses = np.stack(
        [
            fractional_segment(
                decoded["z_rx"][raw_index, channel],
                decoded["range_gate"][observation],
                n_fast,
            )
            for raw_index, observation in zip(raw_indices, pulse_indices)
        ]
    )
    pulse_templates = decoded["z_tx"][raw_indices].astype(np.complex128)
    triplet_center_time_s = float(decoded["tx_idx"][middle] / FS_HZ)
    if "physics_model" in hypothesis:
        trajectory_doppler = physics.predicted_doppler(
            hypothesis["physics_model"], hypothesis["physics_velocity"]
        ) * 1e3
        velocity_guess = float(
            np.interp(triplet_center_time_s, decoded["tx_idx"] / FS_HZ, trajectory_doppler)
        )
    else:
        velocity_guess = float(
            decoded["coarse_doppler_mps"][middle]
        )
    frequency_guess = 2.0 * velocity_guess / pc.wavelength
    pulse_snr = np.asarray(decoded["snr"][pulse_indices], dtype=float)
    matched_filter_amplitudes = rms_baud_amplitudes(raw_pulses, pulse_templates)
    fit = fit_three_pulse_envelope(
        raw_pulses,
        pulse_templates,
        pulse_spacing_s,
        frequency_guess,
        acceleration_guess,
        pc.wavelength,
        acceleration_half_width_mps2=args.acceleration_half_width_mps2,
        pulse_snr=pulse_snr,
        matched_filter_amplitudes=matched_filter_amplitudes,
    )
    conditioned_fit = fit_three_pulse_envelope(
        raw_pulses,
        pulse_templates,
        pulse_spacing_s,
        frequency_guess,
        acceleration_guess,
        pc.wavelength,
        acceleration_half_width_mps2=args.acceleration_half_width_mps2,
        pulse_snr=pulse_snr,
        matched_filter_amplitudes=matched_filter_amplitudes,
        fixed_acceleration_mps2=acceleration_guess,
        initial_parameters=fit["starting_parameters"],
    )
    profile = acceleration_profile(
        fit,
        raw_pulses,
        pulse_templates,
        pulse_spacing_s,
        frequency_guess,
        acceleration_guess,
        pc.wavelength,
        args.profile_points,
        args.acceleration_half_width_mps2,
        pulse_snr,
        matched_filter_amplitudes,
    )
    metadata = {
        "previous_index": previous,
        "middle_index": middle,
        "current_index": current,
        "beam_id": int(decoded["beam_id"][current]),
        "channel": channel,
        "pulse_spacing_s": pulse_spacing_s,
        "frequency_guess_hz": frequency_guess,
        "velocity_guess_mps": velocity_guess,
        "acceleration_guess_mps2": acceleration_guess,
        "global_acceleration_guess_mps2": global_acceleration_guess,
        "acceleration_seed_source": acceleration_seed_source,
        "frequency_seed_source": seed_source,
        "range_gate_0_samples": float(decoded["range_gate"][previous]),
        "range_gate_1_samples": float(decoded["range_gate"][middle]),
        "range_gate_2_samples": float(decoded["range_gate"][current]),
        "matched_filter_amplitude_0": float(matched_filter_amplitudes[0]),
        "matched_filter_amplitude_1": float(matched_filter_amplitudes[1]),
        "matched_filter_amplitude_2": float(matched_filter_amplitudes[2]),
        "pulse_snr_0": float(pulse_snr[0]),
        "pulse_snr_1": float(pulse_snr[1]),
        "pulse_snr_2": float(pulse_snr[2]),
    }
    covariance = np.asarray(fit["parameter_covariance"], dtype=float)
    metadata["frequency_acceleration_correlation"] = float(
        covariance[4, 5] / np.sqrt(covariance[4, 4] * covariance[5, 5])
    )
    metadata["conditioned_frequency_hz"] = float(conditioned_fit["parameters"][4])
    stem = f"three_pulse_complex_envelope_{args.sample_idx}"
    plot_result(
        args.sample_idx,
        fit,
        conditioned_fit,
        profile,
        metadata,
        args.output_dir / f"{stem}.png",
    )
    write_result(
        args.output_dir / f"{stem}.h5",
        args.sample_idx,
        fit,
        conditioned_fit,
        profile,
        metadata,
    )
    print(
        f"sample={args.sample_idx} triplet=({previous},{middle},{current}) beam={metadata['beam_id']} channel={channel} "
        f"dt_ms={1e3 * pulse_spacing_s:.6f} f_guess_hz={frequency_guess:.6f} "
        f"a_guess_mps2={acceleration_guess:.6f} f_fit_hz={fit['parameters'][4]:.6f} "
        f"a_fit_mps2={fit['parameters'][5]:.6f} "
        f"f_conditioned_hz={conditioned_fit['parameters'][4]:.6f} "
        f"amplitude_ratios=(1,{np.exp(fit['parameters'][1] - fit['parameters'][0]):.6f},"
        f"{np.exp(fit['parameters'][2] - fit['parameters'][0]):.6f})",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
