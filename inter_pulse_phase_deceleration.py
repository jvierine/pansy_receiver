#!/usr/bin/env python3
"""Estimate meteor deceleration from phase-coherent decoded pulse pairs."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_best_alias_physics_models as physics
import pansy_config as pc
from interferometer_alias_diagnostics import amp_scale, load_cut, recompute_cut_observables
from plot_interferometric_disambiguation import split_observations_by_range_time, subset_pulse_observations
from range_interpolation_study import DRG_KM, PowerOnlyRangeDopplerSearch, quadratic_peak_2d


FS_HZ = 1e6


def load_selected(path: Path) -> dict:
    with h5py.File(path, "r") as handle:
        label = handle.attrs["selected_hypothesis"]
        label = label.decode() if isinstance(label, bytes) else str(label)
        group = handle["hypotheses"][label]
        out = {
            "label": label,
            "selected_range_time_component": int(handle.attrs.get("selected_range_time_component", -1)),
            "t_rel_s": np.asarray(group["t_rel_s"], dtype=float),
            "range_km": np.asarray(group["range_km"], dtype=float),
            "doppler_mps": np.asarray(group["doppler_mps"], dtype=float),
            "physics_model": np.asarray(group["physics_ceplecha_model"], dtype=float),
            "physics_velocity": np.asarray(group["physics_ceplecha_velocity_km_s"], dtype=float),
        }
    order = np.argsort(out["t_rel_s"])
    for name, values in list(out.items()):
        if isinstance(values, np.ndarray) and len(values) == len(order):
            out[name] = values[order]
    return out


def diagnostic_measurement_clock(cut: dict, hyp: dict, snr_threshold: float) -> dict:
    observations = recompute_cut_observables(cut, interp=1)
    good = np.asarray(observations["snr"], dtype=float) > snr_threshold
    clock = {
        name: values[good] if isinstance(values, np.ndarray) and len(values) == len(good) else values
        for name, values in observations.items()
    }
    component = hyp["selected_range_time_component"]
    if component >= 0:
        segments = split_observations_by_range_time(clock, min_points=3)
        if component < len(segments):
            clock = subset_pulse_observations(clock, segments[component])
    return clock


def diagnostic_to_raw_pulses(cut: dict, clock: dict, t_rel_s: np.ndarray) -> np.ndarray:
    raw_tx = np.asarray(cut["tx_idx"], dtype=float)
    absolute_tx = float(np.asarray(clock["tx_idx"])[0]) + np.asarray(t_rel_s) * FS_HZ
    raw_index = np.asarray([np.argmin(np.abs(raw_tx - value)) for value in absolute_tx], dtype=int)
    error_s = (raw_tx[raw_index] - absolute_tx) / FS_HZ
    if np.nanmax(np.abs(error_s)) > 0.002:
        raise RuntimeError(f"diagnostic/raw pulse clock mismatch: {np.nanmax(np.abs(error_s)):.6f} s")
    return raw_index


def precise_matched_filter_estimates(cut: dict, hyp: dict, snr_threshold: float) -> dict:
    clock = diagnostic_measurement_clock(cut, hyp, snr_threshold)
    raw_idx = diagnostic_to_raw_pulses(cut, clock, hyp["t_rel_s"])
    z_rx = np.asarray(cut["zrx_echoes_re"], np.float32) + 1j * np.asarray(cut["zrx_echoes_im"], np.float32)
    z_tx = np.asarray(cut["ztx_pulses_re"], np.float32) + 1j * np.asarray(cut["ztx_pulses_im"], np.float32)
    z_rx *= amp_scale()[None, :, None]
    delays = np.asarray(cut["delays"], dtype=float)
    search = PowerOnlyRangeDopplerSearch(
        txlen=z_tx.shape[1], echolen=z_rx.shape[2], n_channels=z_rx.shape[1], interp=2, fft_pad=16
    )
    ranges = np.full(len(raw_idx), np.nan)
    dopplers = np.full(len(raw_idx), np.nan)
    for obs_i, pulse in enumerate(raw_idx):
        power = search.mf(z_tx[pulse], z_rx[pulse])
        range_profile = np.nanmax(power, axis=1)
        range_index = int(np.nanargmax(range_profile))
        doppler_index = int(np.nanargmax(power[range_index]))
        dr, dd = quadratic_peak_2d(np.log(np.maximum(power, np.nanmedian(power))), range_index, doppler_index)
        ranges[obs_i] = (delays[pulse] + (range_index + dr) / 2.0) * DRG_KM
        dopplers[obs_i] = search.dopv[doppler_index] + dd * (search.dopv[1] - search.dopv[0])
    return {"raw_idx": raw_idx, "range_km": ranges, "doppler_mps": dopplers}


def fractional_segment(values: np.ndarray, start: float, length: int) -> np.ndarray:
    """Linearly interpolate a complex fast-time segment at a fractional delay."""
    sample = float(start) + np.arange(length, dtype=np.float64)
    grid = np.arange(values.shape[-1], dtype=np.float64)
    if sample[0] < 0.0 or sample[-1] > grid[-1]:
        return np.full(length, np.nan + 1j * np.nan, dtype=np.complex128)
    return np.interp(sample, grid, values.real) + 1j * np.interp(sample, grid, values.imag)


def decoded_pulse_responses(
    cut: dict,
    hyp: dict,
    snr_threshold: float,
    precise_range_km: np.ndarray,
    precise_doppler_mps: np.ndarray,
) -> dict:
    """Range-align and decode retained meteor pulses while preserving phase."""
    clock = diagnostic_measurement_clock(cut, hyp, snr_threshold)
    raw_idx = diagnostic_to_raw_pulses(cut, clock, hyp["t_rel_s"])
    absolute_tx = float(np.asarray(clock["tx_idx"])[0]) + np.asarray(hyp["t_rel_s"]) * FS_HZ
    clock_tx = np.asarray(clock["tx_idx"], dtype=float)
    clock_idx = np.asarray([np.argmin(np.abs(clock_tx - value)) for value in absolute_tx], dtype=int)
    z_rx = np.asarray(cut["zrx_echoes_re"], np.float32) + 1j * np.asarray(
        cut["zrx_echoes_im"], np.float32
    )
    z_tx = np.asarray(cut["ztx_pulses_re"], np.float32) + 1j * np.asarray(
        cut["ztx_pulses_im"], np.float32
    )
    z_rx *= amp_scale()[None, :, None]
    delays = np.asarray(cut["delays"], dtype=np.float64)
    tx_idx = np.asarray(cut["tx_idx"], dtype=np.float64)
    beam_id = np.asarray(cut["beam_id"], dtype=np.int64)
    range_km = np.asarray(precise_range_km, dtype=np.float64)
    coarse = np.asarray(precise_doppler_mps, dtype=np.float64)
    n_obs = len(raw_idx)
    n_ch = z_rx.shape[1]
    n_fast = z_tx.shape[1]
    fast_time_s = np.arange(n_fast, dtype=np.float64) / FS_HZ
    response = np.full((n_obs, n_ch), np.nan + 1j * np.nan, dtype=np.complex128)
    range_gate = range_km / DRG_KM - delays[raw_idx]
    decoded = np.full((n_obs, n_ch, n_fast), np.nan + 1j * np.nan, dtype=np.complex64)
    derotated = np.full_like(decoded, np.nan + 1j * np.nan)

    for obs_i, pulse in enumerate(raw_idx):
        tx = z_tx[pulse].astype(np.complex128)
        envelope = np.abs(tx) ** 2
        use = envelope > 0.05 * np.nanmax(envelope)
        if not np.any(use):
            continue
        doppler_hz = 2.0 * coarse[obs_i] / pc.wavelength
        derotation = np.exp(-1j * 2.0 * np.pi * doppler_hz * fast_time_s)
        for channel in range(n_ch):
            echo = fractional_segment(z_rx[pulse, channel], range_gate[obs_i], n_fast)
            dec = echo * np.conj(tx)
            derot = dec * derotation
            decoded[obs_i, channel] = dec
            derotated[obs_i, channel] = derot
            response[obs_i, channel] = np.sum(derot[use]) / np.sum(use)

    return {
        "raw_idx": raw_idx,
        "tx_idx": tx_idx[raw_idx],
        "beam_id": beam_id[raw_idx],
        "range_gate": range_gate,
        "range_km": range_km,
        "coarse_doppler_mps": coarse,
        "snr": np.asarray(clock["snr"], dtype=float)[clock_idx],
        "xc": np.asarray(clock["xc"])[clock_idx],
        "response": response,
        "decoded": decoded,
        "derotated": derotated,
        "z_rx": z_rx,
        "z_tx": z_tx,
        "fast_time_s": fast_time_s,
    }


def decoded_waveform_pairs(decoded: dict, same_beam: bool, search_hz: float = 2000.0) -> dict:
    """Fit the residual fast-time sinusoid in d[n+1] conj(d[n])."""
    tx_s = decoded["tx_idx"] / FS_HZ
    beam = decoded["beam_id"]
    wave = decoded["decoded"].astype(np.complex128)
    previous_by_beam: dict[int, int] = {}
    frequencies = np.linspace(-float(search_hz), float(search_hz), 4001)
    fast_time = decoded["fast_time_s"]
    basis = np.exp(-1j * 2.0 * np.pi * frequencies[:, None] * fast_time[None, :])
    rows = []
    for cur in range(len(tx_s)):
        if same_beam:
            prev = previous_by_beam.get(int(beam[cur]))
            previous_by_beam[int(beam[cur])] = cur
        else:
            prev = cur - 1 if cur else None
        if prev is None:
            continue
        delta_t = tx_s[cur] - tx_s[prev]
        if not np.isfinite(delta_t) or delta_t <= 0.0 or delta_t > 0.025:
            continue
        channel_product = wave[cur] * np.conj(wave[prev])
        amplitude = np.abs(channel_product)
        finite = np.isfinite(channel_product.real) & np.isfinite(channel_product.imag)
        channel_product = np.where(finite, channel_product, 0.0j)
        amplitude = np.where(finite, amplitude, 0.0)
        product = np.sum(channel_product, axis=0)
        weight = np.sum(amplitude, axis=0)
        if not np.any(weight > 0.0):
            continue
        threshold = 0.10 * np.nanmax(weight)
        use = np.isfinite(product.real) & np.isfinite(product.imag) & (weight > threshold)
        if np.count_nonzero(use) < 12:
            continue
        phasor = product[use] / np.maximum(np.abs(product[use]), 1e-30)
        sample_weight = weight[use] / np.nanmax(weight[use])
        spectrum = basis[:, use] @ (sample_weight * phasor)
        power = np.abs(spectrum) ** 2
        peak = int(np.nanargmax(power))
        frequency_hz = float(frequencies[peak])
        if 0 < peak < len(frequencies) - 1:
            ym, y0, yp = np.log(np.maximum(power[peak - 1 : peak + 2], 1e-300))
            denom = ym - 2.0 * y0 + yp
            if np.isfinite(denom) and abs(denom) > 1e-15:
                frequency_hz += float(0.5 * (ym - yp) / denom) * (frequencies[1] - frequencies[0])
        model_rotation = np.exp(-1j * 2.0 * np.pi * frequency_hz * fast_time[use])
        coefficient = np.sum(sample_weight * phasor * model_rotation) / np.sum(sample_weight)
        residual_phase = np.angle(phasor * model_rotation * np.exp(-1j * np.angle(coefficient)))
        phase_rms = float(np.sqrt(np.sum(sample_weight * residual_phase**2) / np.sum(sample_weight)))
        centered_time = fast_time[use] - np.sum(sample_weight * fast_time[use]) / np.sum(sample_weight)
        information = np.sum(sample_weight * centered_time**2)
        frequency_std_hz = phase_rms / (2.0 * np.pi * np.sqrt(max(information, 1e-30)))
        delta_velocity = frequency_hz * pc.wavelength / 2.0
        delta_velocity_std = frequency_std_hz * pc.wavelength / 2.0
        rows.append(
            (
                prev,
                cur,
                0.5 * (tx_s[prev] + tx_s[cur]),
                delta_t,
                int(beam[prev]),
                int(beam[cur]),
                frequency_hz,
                frequency_std_hz,
                delta_velocity,
                delta_velocity_std,
                delta_velocity / delta_t,
                delta_velocity_std / delta_t,
                phase_rms,
                float(abs(coefficient)),
            )
        )
    names = (
        "previous_index",
        "current_index",
        "time_s",
        "delta_t_s",
        "previous_beam_id",
        "current_beam_id",
        "delta_frequency_hz",
        "delta_frequency_std_hz",
        "delta_velocity_mps",
        "delta_velocity_std_mps",
        "acceleration_mps2",
        "acceleration_std_mps2",
        "phase_residual_rms_rad",
        "coherence",
    )
    if not rows:
        return {name: np.asarray([]) for name in names}
    return {name: np.asarray(values) for name, values in zip(names, zip(*rows))}


def baud_averaged_beat_pairs(decoded: dict, same_beam: bool = True) -> dict:
    """Average the nearly constant decoded beat phasor and recover pair Doppler."""
    tx_s = decoded["tx_idx"] / FS_HZ
    beam = decoded["beam_id"]
    wave = decoded["decoded"].astype(np.complex128)
    raw_idx = decoded["raw_idx"]
    coarse = decoded["coarse_doppler_mps"]
    previous_by_beam: dict[int, int] = {}
    rows = []
    for cur in range(len(tx_s)):
        if same_beam:
            prev = previous_by_beam.get(int(beam[cur]))
            previous_by_beam[int(beam[cur])] = cur
        else:
            prev = cur - 1 if cur else None
        if prev is None:
            continue
        delta_t = tx_s[cur] - tx_s[prev]
        if not np.isfinite(delta_t) or delta_t <= 0.0 or delta_t > 0.025:
            continue
        channel_product = wave[cur] * np.conj(wave[prev])
        finite = np.isfinite(channel_product.real) & np.isfinite(channel_product.imag)
        product = np.sum(np.where(finite, channel_product, 0.0j), axis=0)
        tx_prev = decoded["z_tx"][int(raw_idx[prev])]
        tx_cur = decoded["z_tx"][int(raw_idx[cur])]
        joint_weight = np.abs(tx_prev) ** 2 * np.abs(tx_cur) ** 2
        segments = envelope_segments(joint_weight)
        _centers, baud_values = weighted_baud_values(product, joint_weight, segments)
        baud_weight = np.asarray([np.sum(joint_weight[segment]) for segment in segments], dtype=float)
        good = np.isfinite(baud_values.real) & np.isfinite(baud_values.imag) & (baud_weight > 0.0)
        if np.count_nonzero(good) < 4:
            continue
        baud_values = baud_values[good]
        baud_weight = baud_weight[good]
        baud_weight /= np.nanmax(baud_weight)
        beat = np.sum(baud_weight * baud_values) / np.sum(baud_weight)
        if not np.isfinite(beat.real) or not np.isfinite(beat.imag) or abs(beat) == 0.0:
            continue
        phase = float(np.angle(beat))
        phase_residual = np.angle(baud_values * np.exp(-1j * phase))
        phase_rms = float(np.sqrt(np.sum(baud_weight * phase_residual**2) / np.sum(baud_weight)))
        effective_n = float(np.sum(baud_weight) ** 2 / np.sum(baud_weight**2))
        phase_std = phase_rms / np.sqrt(max(effective_n, 1.0))
        coherence = float(abs(np.sum(baud_weight * baud_values)) / np.sum(baud_weight * np.abs(baud_values)))
        wrapped_velocity = phase * pc.wavelength / (4.0 * np.pi * delta_t)
        ambiguity = pc.wavelength / (2.0 * delta_t)
        coarse_pair = float(0.5 * (coarse[prev] + coarse[cur]))
        velocity = wrapped_velocity + np.rint((coarse_pair - wrapped_velocity) / ambiguity) * ambiguity
        velocity_std = phase_std * pc.wavelength / (4.0 * np.pi * delta_t)
        rows.append(
            (
                prev,
                cur,
                0.5 * (tx_s[prev] + tx_s[cur]),
                delta_t,
                int(beam[cur]),
                beat.real,
                beat.imag,
                phase,
                phase_std,
                coherence,
                effective_n,
                ambiguity,
                coarse_pair,
                velocity,
                velocity_std,
            )
        )
    names = (
        "previous_index",
        "current_index",
        "time_s",
        "delta_t_s",
        "beam_id",
        "beat_real",
        "beat_imag",
        "beat_phase_rad",
        "beat_phase_std_rad",
        "coherence",
        "effective_bauds",
        "ambiguity_mps",
        "coarse_doppler_mps",
        "phase_doppler_mps",
        "phase_doppler_std_mps",
    )
    if not rows:
        return {name: np.asarray([]) for name in names}
    return {name: np.asarray(values) for name, values in zip(names, zip(*rows))}


def pair_responses(decoded: dict, same_beam: bool) -> dict:
    """Cross-multiply consecutive decoded responses and recover Doppler branch."""
    tx_s = decoded["tx_idx"] / FS_HZ
    beam = decoded["beam_id"]
    response = decoded["response"]
    coarse = decoded["coarse_doppler_mps"]
    prev_for_beam: dict[int, int] = {}
    rows = []
    for cur in range(len(tx_s)):
        if same_beam:
            prev = prev_for_beam.get(int(beam[cur]))
            prev_for_beam[int(beam[cur])] = cur
        else:
            prev = cur - 1 if cur else None
        if prev is None:
            continue
        delta_t = tx_s[cur] - tx_s[prev]
        if not np.isfinite(delta_t) or delta_t <= 0.0 or delta_t > 0.025:
            continue
        channel_cross = response[cur] * np.conj(response[prev])
        finite = np.isfinite(channel_cross.real) & np.isfinite(channel_cross.imag)
        if not np.any(finite):
            continue
        cross = np.sum(channel_cross[finite])
        norm = np.sum(np.abs(channel_cross[finite]))
        if norm <= 0.0 or abs(cross) <= 0.0:
            continue
        phase = float(np.angle(cross))
        wrapped_velocity = phase * pc.wavelength / (4.0 * np.pi * delta_t)
        ambiguity = pc.wavelength / (2.0 * delta_t)
        coarse_pair = float(0.5 * (coarse[prev] + coarse[cur]))
        velocity = wrapped_velocity + np.rint((coarse_pair - wrapped_velocity) / ambiguity) * ambiguity
        coarse_phase = 4.0 * np.pi * coarse_pair * delta_t / pc.wavelength
        residual_phase = float(np.angle(np.exp(1j * (phase - coarse_phase))))
        rows.append(
            (
                prev,
                cur,
                0.5 * (tx_s[prev] + tx_s[cur]),
                delta_t,
                int(beam[cur]),
                phase,
                residual_phase,
                abs(cross),
                abs(cross) / norm,
                ambiguity,
                coarse_pair,
                velocity,
            )
        )
    names = (
        "previous_index",
        "current_index",
        "time_s",
        "delta_t_s",
        "beam_id",
        "phase_rad",
        "residual_phase_rad",
        "pair_amplitude",
        "channel_coherence",
        "ambiguity_mps",
        "coarse_doppler_mps",
        "phase_doppler_mps",
    )
    if not rows:
        return {name: np.asarray([]) for name in names}
    columns = zip(*rows)
    return {name: np.asarray(values) for name, values in zip(names, columns)}


def robust_line(time_s: np.ndarray, velocity_mps: np.ndarray, weight: np.ndarray) -> dict:
    """Fit velocity intercept and acceleration with simple robust clipping."""
    time_s = np.asarray(time_s, dtype=float)
    velocity_mps = np.asarray(velocity_mps, dtype=float)
    weight = np.asarray(weight, dtype=float)
    good = np.isfinite(time_s) & np.isfinite(velocity_mps) & np.isfinite(weight) & (weight > 0)
    center = float(np.nanmedian(time_s[good]))
    x = time_s - center
    keep = good.copy()
    coef = np.asarray([np.nan, np.nan])
    for _ in range(5):
        if np.count_nonzero(keep) < 5:
            break
        coef = np.polyfit(x[keep], velocity_mps[keep], 1, w=np.sqrt(weight[keep]))
        residual = velocity_mps - np.polyval(coef, x)
        median = np.nanmedian(residual[keep])
        sigma = 1.4826 * np.nanmedian(np.abs(residual[keep] - median))
        if not np.isfinite(sigma) or sigma <= 0.0:
            break
        keep = good & (np.abs(residual - median) < 4.0 * sigma)
    model = np.polyval(coef, x)
    rms = float(np.sqrt(np.nanmean((velocity_mps[keep] - model[keep]) ** 2)))
    return {
        "center_time_s": center,
        "intercept_mps": float(coef[1]),
        "acceleration_mps2": float(coef[0]),
        "model_mps": model,
        "keep": keep,
        "rms_mps": rms,
    }


def analyze_event(cut: dict, hyp: dict, snr_threshold: float, precise: dict | None = None) -> dict:
    if precise is None:
        precise = precise_matched_filter_estimates(cut, hyp, snr_threshold)
    decoded = decoded_pulse_responses(
        cut,
        hyp,
        snr_threshold,
        precise_range_km=precise["range_km"],
        precise_doppler_mps=precise["doppler_mps"],
    )
    same = pair_responses(decoded, same_beam=True)
    adjacent = pair_responses(decoded, same_beam=False)
    waveform_same = decoded_waveform_pairs(decoded, same_beam=True)
    waveform_adjacent = decoded_waveform_pairs(decoded, same_beam=False)
    beat_pairs = baud_averaged_beat_pairs(decoded, same_beam=True)
    if len(same["time_s"]) < 5:
        raise RuntimeError("too few same-beam decoded pulse pairs")
    weight = same["pair_amplitude"] * np.maximum(same["channel_coherence"], 1e-3) ** 2
    fit = robust_line(same["time_s"], same["phase_doppler_mps"], weight)
    t_obs_s = decoded["tx_idx"] / FS_HZ
    model_doppler = physics.predicted_doppler(hyp["physics_model"], hyp["physics_velocity"]) * 1e3
    model_at_pair = np.interp(same["time_s"], t_obs_s, model_doppler)
    model_fit = robust_line(same["time_s"], model_at_pair, np.ones(len(model_at_pair)))
    beat_weight = beat_pairs["coherence"] / np.maximum(beat_pairs["phase_doppler_std_mps"], 1.0) ** 2
    beat_fit = robust_line(beat_pairs["time_s"], beat_pairs["phase_doppler_mps"], beat_weight)
    beat_model = np.interp(beat_pairs["time_s"], t_obs_s, model_doppler)
    return {
        "decoded": decoded,
        "same_beam_pairs": same,
        "adjacent_pairs": adjacent,
        "waveform_same_beam_pairs": waveform_same,
        "waveform_adjacent_pairs": waveform_adjacent,
        "baud_averaged_beat_pairs": beat_pairs,
        "baud_beat_fit": beat_fit,
        "baud_beat_model_doppler_mps": beat_model,
        "fit": fit,
        "model_doppler_mps": model_doppler,
        "model_at_pair_mps": model_at_pair,
        "model_fit": model_fit,
    }


def best_stage_pair(result: dict) -> tuple[int, int]:
    pairs = result["waveform_same_beam_pairs"]
    score = pairs["coherence"] / np.maximum(pairs["delta_frequency_std_hz"], 1e-12)
    index = int(np.nanargmax(score))
    return int(pairs["previous_index"][index]), int(pairs["current_index"][index])


def normalized(values: np.ndarray) -> np.ndarray:
    scale = np.nanmax(np.abs(values))
    return values / scale if np.isfinite(scale) and scale > 0 else values


def envelope_segments(weight: np.ndarray, threshold: float = 0.05) -> list[np.ndarray]:
    normalized_weight = np.asarray(weight, dtype=float) / max(float(np.nanmax(weight)), 1e-30)
    indices = np.flatnonzero(normalized_weight > threshold)
    if len(indices) == 0:
        return []
    breaks = np.flatnonzero(np.diff(indices) > 1) + 1
    return [segment for segment in np.split(indices, breaks) if len(segment) >= 2]


def weighted_baud_values(values: np.ndarray, weight: np.ndarray, segments: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    centers = []
    averages = []
    for segment in segments:
        baud_weight = np.asarray(weight[segment], dtype=float)
        denominator = np.sum(baud_weight)
        if not np.isfinite(denominator) or denominator <= 0.0:
            continue
        centers.append(float(np.sum(segment * baud_weight) / denominator))
        averages.append(np.sum(baud_weight * values[segment]) / denominator)
    return np.asarray(centers), np.asarray(averages)


def plot_baud_averages(sample_idx: int, result: dict, output: Path) -> None:
    """Show one SNR-weighted complex value and phase for each transmit baud."""
    decoded = result["decoded"]
    prev, cur = best_stage_pair(result)
    channel_cross = decoded["response"][cur] * np.conj(decoded["response"][prev])
    channel = int(np.nanargmax(np.abs(channel_cross)))
    raw_prev = int(decoded["raw_idx"][prev])
    raw_cur = int(decoded["raw_idx"][cur])
    first = decoded["decoded"][prev, channel].astype(np.complex128)
    second = decoded["decoded"][cur, channel].astype(np.complex128)
    product = second * np.conj(first)
    weight_first = np.abs(decoded["z_tx"][raw_prev]) ** 2
    weight_second = np.abs(decoded["z_tx"][raw_cur]) ** 2
    joint_weight = weight_first * weight_second
    segments = envelope_segments(joint_weight)
    rows = []
    for values, weight, label in (
        (first, weight_first, r"$d_n$"),
        (second, weight_second, r"$d_{n+1}$"),
        (product, joint_weight, r"$q_n=d_{n+1}d_n^*$"),
    ):
        center, average = weighted_baud_values(values, weight, segments)
        rows.append((center, average, label))

    fig, axes = plt.subplots(3, 1, figsize=(9.0, 8.4), sharex=True, constrained_layout=True)
    for ax, (center, average, label) in zip(axes, rows):
        scaled = normalized(average)
        ax.plot(center, scaled.real, "o-", ms=5, lw=1.0, color="C0", label="real")
        ax.plot(center, scaled.imag, "o-", ms=5, lw=1.0, color="C1", label="imaginary")
        ax.set_ylabel(label + "\nnormalized complex value")
        ax.grid(alpha=0.2, lw=0.5)
        ax.legend(frameon=False, loc="upper left", ncol=2)
        phase_ax = ax.twinx()
        phase = np.angle(average)
        phase_ax.plot(center, phase, "o-", ms=4, lw=1.0, color="C3", label="phase")
        phase_ax.set_ylabel("Phase (rad)", color="C3")
        phase_ax.set_ylim(-np.pi, np.pi)
        phase_ax.tick_params(axis="y", colors="C3")
        phase_ax.legend(frameon=False, loc="upper right")
    axes[-1].set_xlabel("SNR-weighted baud center (fast-time sample)")
    beam = int(decoded["beam_id"][cur])
    delta_t_ms = (decoded["tx_idx"][cur] - decoded["tx_idx"][prev]) / 1e3
    timestamp = dt.datetime.fromtimestamp(sample_idx / 1e6, tz=dt.timezone.utc).isoformat(timespec="milliseconds")
    fig.suptitle(
        f"SNR-weighted decoded baud values, {timestamp}; channel {channel}, beam {beam}, "
        f"pulse separation {delta_t_ms:.1f} ms"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=190)
    plt.close(fig)


def plot_baud_beat_deceleration(sample_idx: int, result: dict, output: Path) -> None:
    pairs = result["baud_averaged_beat_pairs"]
    fit = result["baud_beat_fit"]
    zero = float(result["decoded"]["tx_idx"][0]) / FS_HZ
    time_s = pairs["time_s"] - zero
    keep = np.asarray(fit["keep"], dtype=bool)
    fig, axes = plt.subplots(3, 1, figsize=(9.0, 8.5), sharex=True, constrained_layout=True)
    ax = axes[0]
    scatter = ax.scatter(time_s, pairs["beat_phase_rad"], c=pairs["coherence"], s=14, cmap="viridis")
    ax.errorbar(time_s, pairs["beat_phase_rad"], yerr=pairs["beat_phase_std_rad"], fmt="none", ecolor="0.75", lw=0.5)
    ax.set_ylabel(r"$\arg Q_n$ (rad)")
    ax.set_ylim(-np.pi, np.pi)
    fig.colorbar(scatter, ax=ax, label="Baud coherence")

    ax = axes[1]
    ax.errorbar(
        time_s,
        pairs["phase_doppler_mps"] / 1e3,
        yerr=pairs["phase_doppler_std_mps"] / 1e3,
        fmt=".",
        ms=3,
        color="black",
        ecolor="0.8",
        label="SNR-weighted beat phase",
    )
    ax.plot(time_s, result["baud_beat_model_doppler_mps"] / 1e3, color="C0", lw=1.1, label="trajectory model")
    ax.plot(time_s[keep], fit["model_mps"][keep] / 1e3, color="C3", lw=1.2, label="robust linear fit")
    ax.set_ylabel("Radial velocity (km/s)")
    ax.legend(frameon=False)

    ax = axes[2]
    residual = pairs["phase_doppler_mps"] - result["baud_beat_model_doppler_mps"]
    ax.errorbar(
        time_s,
        residual,
        yerr=pairs["phase_doppler_std_mps"],
        fmt=".",
        ms=3,
        color="black",
        ecolor="0.8",
    )
    ax.axhline(0.0, color="C0", lw=1.0)
    ax.set_ylabel("Velocity residual (m/s)")
    ax.set_xlabel("Time (s)")
    ax.text(
        0.02,
        0.96,
        rf"phase acceleration {fit['acceleration_mps2']:.0f} m s$^{{-2}}$"
        + "\n"
        + rf"model {result['model_fit']['acceleration_mps2']:.0f} m s$^{{-2}}$"
        + "\n"
        + rf"velocity RMS {fit['rms_mps']:.0f} m s$^{{-1}}$",
        transform=ax.transAxes,
        va="top",
    )
    for ax in axes:
        ax.grid(alpha=0.2, lw=0.5)
    timestamp = dt.datetime.fromtimestamp(sample_idx / 1e6, tz=dt.timezone.utc).isoformat(timespec="milliseconds")
    fig.suptitle(f"Deceleration from SNR-weighted mean decoded beat phasor: {timestamp}")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=190)
    plt.close(fig)


def plot_decoded_waveforms(sample_idx: int, result: dict, output: Path) -> None:
    """Plot the two decoded echoes and their direct sample-wise product."""
    decoded = result["decoded"]
    pairs = result["waveform_same_beam_pairs"]
    prev, cur = best_stage_pair(result)
    channel_cross = decoded["response"][cur] * np.conj(decoded["response"][prev])
    channel = int(np.nanargmax(np.abs(channel_cross)))
    first = decoded["decoded"][prev, channel].astype(np.complex128)
    second = decoded["decoded"][cur, channel].astype(np.complex128)
    product = second * np.conj(first)
    sample = np.arange(len(first))
    row = np.flatnonzero((pairs["previous_index"] == prev) & (pairs["current_index"] == cur))[0]
    delta_frequency = float(pairs["delta_frequency_hz"][row])
    first_n = normalized(first)
    second_n = normalized(second)
    product_n = normalized(product)

    fig, axes = plt.subplots(3, 1, figsize=(9.0, 8.2), sharex=True, constrained_layout=True)
    for ax, values, label in (
        (axes[0], first_n, r"$d_n[k]$"),
        (axes[1], second_n, r"$d_{n+1}[k]$"),
        (axes[2], product_n, r"$q_n[k]=d_{n+1}[k]d_n^*[k]$"),
    ):
        ax.plot(sample, values.real, color="C0", lw=1.0, label="real")
        ax.plot(sample, values.imag, color="C1", lw=1.0, label="imaginary")
        ax.plot(sample, np.abs(values), color="0.25", lw=0.8, alpha=0.7, label="magnitude")
        ax.set_ylabel(label + "\nnormalized voltage")
        ax.grid(alpha=0.2, lw=0.5)
        ax.legend(frameon=False, ncol=3, loc="upper right")
    good = np.isfinite(product.real) & (np.abs(product) > 0.10 * np.nanmax(np.abs(product)))
    phase_ax = axes[2].twinx()
    phase = np.unwrap(np.angle(product[good]))
    fast_time = decoded["fast_time_s"][good]
    phase_model = 2.0 * np.pi * delta_frequency * fast_time
    phase_model += np.nanmedian(phase - phase_model)
    phase_ax.scatter(sample[good], phase, s=9, color="C3", alpha=0.55, label="unwrapped phase")
    phase_ax.plot(sample[good], phase_model, color="C3", lw=1.3, label=rf"fit $\Delta f={delta_frequency:.1f}$ Hz")
    phase_ax.set_ylabel("Product phase (rad)", color="C3")
    phase_ax.legend(frameon=False, loc="lower right")
    axes[2].set_xlabel("Fast-time sample at 1 MHz")
    beam = int(decoded["beam_id"][cur])
    delta_t_ms = (decoded["tx_idx"][cur] - decoded["tx_idx"][prev]) / 1e3
    timestamp = dt.datetime.fromtimestamp(sample_idx / 1e6, tz=dt.timezone.utc).isoformat(timespec="milliseconds")
    fig.suptitle(
        f"Decoded echo cancellation, {timestamp}; channel {channel}, beam {beam}, "
        f"pulse separation {delta_t_ms:.1f} ms"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=190)
    plt.close(fig)


def plot_event(sample_idx: int, result: dict, output: Path) -> None:
    decoded = result["decoded"]
    pairs = result["same_beam_pairs"]
    prev, cur = best_stage_pair(result)
    channel_cross = decoded["response"][cur] * np.conj(decoded["response"][prev])
    channel = int(np.nanargmax(np.abs(channel_cross)))
    raw_prev = int(decoded["raw_idx"][prev])
    raw_cur = int(decoded["raw_idx"][cur])
    gate_prev = decoded["range_gate"][prev]
    gate_cur = decoded["range_gate"][cur]
    n_fast = decoded["z_tx"].shape[1]
    echo_prev = fractional_segment(decoded["z_rx"][raw_prev, channel], gate_prev, n_fast)
    echo_cur = fractional_segment(decoded["z_rx"][raw_cur, channel], gate_cur, n_fast)
    tx_prev = decoded["z_tx"][raw_prev]
    tx_cur = decoded["z_tx"][raw_cur]
    dec_prev = decoded["decoded"][prev, channel]
    dec_cur = decoded["decoded"][cur, channel]
    der_prev = decoded["derotated"][prev, channel]
    der_cur = decoded["derotated"][cur, channel]
    product = dec_cur * np.conj(dec_prev)
    sample = np.arange(n_fast)
    time_zero = decoded["tx_idx"][0] / FS_HZ
    waveform_pairs = result["waveform_same_beam_pairs"]
    adjacent_waveform_pairs = result["waveform_adjacent_pairs"]
    t_pair = pairs["time_s"] - time_zero
    t_wave = waveform_pairs["time_s"] - time_zero
    t_adjacent = adjacent_waveform_pairs["time_s"] - time_zero
    t_obs = decoded["tx_idx"] / FS_HZ - time_zero

    fig, axes = plt.subplots(3, 2, figsize=(12.0, 10.0), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(sample, normalized(echo_prev).real, color="C0", lw=0.9, label="pulse n")
    ax.plot(sample, normalized(echo_cur).real, color="C1", lw=0.9, label="pulse n+1, same beam")
    ax.set_title(f"Range-aligned raw voltage, channel {channel}")
    ax.set_ylabel("Normalized real voltage")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(sample, np.angle(tx_prev), color="C0", lw=0.8, label="TX reference n")
    ax.plot(sample, np.angle(tx_cur), color="C1", lw=0.8, label="TX reference n+1")
    ax.set_title("Measured transmit-reference phase")
    ax.set_ylabel("Phase (rad)")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(sample, normalized(dec_prev).real, color="C0", lw=0.9, label="decoded n")
    ax.plot(sample, normalized(dec_cur).real, color="C1", lw=0.9, label="decoded n+1")
    ax.set_title(r"Code decoding: $e_n[k]x_n^*[k]$")
    ax.set_ylabel("Normalized real voltage")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.plot(sample, normalized(der_prev).real, color="C0", lw=0.9, label="Doppler-removed n")
    ax.plot(sample, normalized(der_cur).real, color="C1", lw=0.9, label="Doppler-removed n+1")
    ax.set_title("After within-pulse coarse-Doppler removal")
    ax.set_ylabel("Normalized real voltage")
    ax.legend(frameon=False)

    ax = axes[2, 0]
    product_good = np.isfinite(product.real) & (np.abs(product) > 0.05 * np.nanmax(np.abs(product)))
    product_phase = np.unwrap(np.angle(product[product_good]))
    stage_row = np.flatnonzero(
        (waveform_pairs["previous_index"] == prev) & (waveform_pairs["current_index"] == cur)
    )[0]
    delta_frequency_hz = waveform_pairs["delta_frequency_hz"][stage_row]
    stage_time = decoded["fast_time_s"][product_good]
    phase_model = 2.0 * np.pi * delta_frequency_hz * stage_time
    phase_model += np.nanmedian(product_phase - phase_model)
    ax.scatter(sample[product_good], product_phase, c=np.abs(product[product_good]), s=12, cmap="plasma")
    ax.plot(sample[product_good], phase_model, color="black", lw=1.2, label=rf"fit: $\Delta f={delta_frequency_hz:.1f}$ Hz")
    ax.set_title(r"Decoded product; its fast-time slope measures $f_{n+1}-f_n$")
    ax.set_ylabel("Product phase (rad)")
    ax.legend(frameon=False)

    ax = axes[2, 1]
    ax.errorbar(
        t_adjacent,
        adjacent_waveform_pairs["acceleration_mps2"],
        yerr=adjacent_waveform_pairs["acceleration_std_mps2"],
        fmt=".",
        ms=3,
        color="0.65",
        ecolor="0.85",
        label="literal adjacent echoes",
    )
    ax.errorbar(
        t_wave,
        waveform_pairs["acceleration_mps2"],
        yerr=waveform_pairs["acceleration_std_mps2"],
        fmt=".",
        ms=4,
        color="black",
        ecolor="0.75",
        label="next echo in same beam",
    )
    ax.axhline(result["model_fit"]["acceleration_mps2"], color="C0", lw=1.2, label="trajectory-model radial acceleration")
    ax.set_title("Acceleration from residual fast-time sinusoid")
    ax.set_ylabel(r"Radial acceleration (m s$^{-2}$)")
    ax.legend(frameon=False, fontsize=8)

    for ax in axes.ravel():
        ax.set_xlabel("Fast-time sample" if ax is not axes[2, 1] else "Time (s)")
        ax.grid(alpha=0.2, lw=0.5)
    timestamp = dt.datetime.fromtimestamp(sample_idx / 1e6, tz=dt.timezone.utc).isoformat(timespec="milliseconds")
    fig.suptitle(f"Inter-pulse decoded-phase deceleration test: {timestamp}")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def write_h5(sample_idx: int, result: dict, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        handle.attrs["sample_idx"] = sample_idx
        handle.attrs["phase_acceleration_mps2"] = result["fit"]["acceleration_mps2"]
        handle.attrs["trajectory_model_radial_acceleration_mps2"] = result["model_fit"]["acceleration_mps2"]
        handle.attrs["phase_velocity_rms_mps"] = result["fit"]["rms_mps"]
        handle.attrs["phase_convention"] = "current decoded response times conjugate previous decoded response"
        for group_name in (
            "same_beam_pairs",
            "adjacent_pairs",
            "waveform_same_beam_pairs",
            "waveform_adjacent_pairs",
            "baud_averaged_beat_pairs",
        ):
            group = handle.create_group(group_name)
            for name, values in result[group_name].items():
                group.create_dataset(name, data=values)
        measurement = handle.create_group("measurement")
        for name in ("raw_idx", "tx_idx", "beam_id", "range_gate", "range_km", "coarse_doppler_mps", "snr", "xc", "response"):
            measurement.create_dataset(name, data=result["decoded"][name])
        measurement.create_dataset("trajectory_model_doppler_mps", data=result["model_doppler_mps"])


def diagnostic_path(base: Path, sample_idx: int) -> Path:
    day = dt.datetime.fromtimestamp(sample_idx / 1e6, tz=dt.timezone.utc).strftime("%Y-%m-%d")
    return base / "events" / day / f"pansy_disambiguation_diagnostics_{sample_idx}.h5"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-idx", type=int, action="append", required=True)
    parser.add_argument("--base", type=Path, default=Path("/mnt/data/juha/pansy"))
    parser.add_argument("--output-dir", type=Path, default=Path("test_plots/inter_pulse_phase"))
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    parser.add_argument("--precise-cache-dir", type=Path)
    args = parser.parse_args()
    for sample_idx in args.sample_idx:
        diag = diagnostic_path(args.base, sample_idx)
        hyp = load_selected(diag)
        cut = load_cut(args.base / "metadata/cut", sample_idx)
        precise = None
        if args.precise_cache_dir is not None:
            cache_path = args.precise_cache_dir / f"highres_fft_i2_p16_{sample_idx}.h5"
            if cache_path.exists():
                with h5py.File(cache_path, "r") as cache:
                    precise = {
                        "raw_idx": np.asarray(cache["raw_idx"], dtype=int),
                        "range_km": np.asarray(cache["range_km"], dtype=float),
                        "doppler_mps": np.asarray(cache["doppler_mps"], dtype=float),
                    }
        result = analyze_event(cut, hyp, args.snr_threshold, precise=precise)
        stem = f"inter_pulse_phase_{sample_idx}"
        plot_event(sample_idx, result, args.output_dir / f"{stem}.png")
        plot_decoded_waveforms(sample_idx, result, args.output_dir / f"{stem}_waveforms.png")
        plot_baud_averages(sample_idx, result, args.output_dir / f"{stem}_baud_averages.png")
        plot_baud_beat_deceleration(sample_idx, result, args.output_dir / f"{stem}_baud_beat_deceleration.png")
        write_h5(sample_idx, result, args.output_dir / f"{stem}.h5")
        print(
            sample_idx,
            f"n_pairs={len(result['same_beam_pairs']['time_s'])}",
            f"phase_acceleration_mps2={result['fit']['acceleration_mps2']:.6g}",
            f"model_acceleration_mps2={result['model_fit']['acceleration_mps2']:.6g}",
            f"phase_velocity_rms_mps={result['fit']['rms_mps']:.6g}",
            f"waveform_pairs={len(result['waveform_same_beam_pairs']['time_s'])}",
            f"median_waveform_acceleration_mps2={np.nanmedian(result['waveform_same_beam_pairs']['acceleration_mps2']):.6g}",
            f"median_waveform_acceleration_std_mps2={np.nanmedian(result['waveform_same_beam_pairs']['acceleration_std_mps2']):.6g}",
            f"baud_beat_acceleration_mps2={result['baud_beat_fit']['acceleration_mps2']:.6g}",
            f"baud_beat_velocity_rms_mps={result['baud_beat_fit']['rms_mps']:.6g}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
