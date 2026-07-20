"""Calibrated coherent receiver-module addition using a fitted arrival direction."""

from __future__ import annotations

import numpy as np

import pansy_config as pc
import pansy_interferometry as pint


def estimate_event_module_voltage_gain(response: np.ndarray) -> np.ndarray:
    """Estimate residual module voltage gains from matched-filter responses.

    A robust two-way log-amplitude decomposition removes the common meteor
    amplitude at each pulse and leaves one constant gain per receiver module.
    """
    amplitude = np.abs(np.asarray(response, dtype=np.complex128))
    valid = np.isfinite(amplitude) & (amplitude > 0.0)
    log_amplitude = np.full(amplitude.shape, np.nan, dtype=float)
    log_amplitude[valid] = np.log(amplitude[valid])
    pulse_level = np.nanmedian(log_amplitude, axis=1, keepdims=True)
    log_gain = np.nanmedian(log_amplitude - pulse_level, axis=0)
    finite = np.isfinite(log_gain)
    if not np.any(finite):
        return np.ones(amplitude.shape[1], dtype=float)
    log_gain[~finite] = np.nanmedian(log_gain[finite])
    log_gain -= np.mean(log_gain)
    return np.exp(log_gain)


def enu_line_of_sight_to_arrival_uvw(position_enu: np.ndarray) -> np.ndarray:
    """Convert an ENU radar-to-target vector to the interferometer arrival convention.

    The catalogue trajectory uses positive Up, while ``pansy_interferometry`` uses
    a wave-arrival vector with negative ``w`` for targets above the horizon.
    """
    position = np.asarray(position_enu, dtype=float)
    if position.shape[-1] != 3:
        raise ValueError("position_enu must have a final axis of length three")
    norm = np.linalg.norm(position, axis=-1, keepdims=True)
    if np.any(~np.isfinite(norm)) or np.any(norm <= 0.0):
        raise ValueError("position_enu must contain finite, nonzero vectors")
    arrival = position / norm
    arrival = arrival.copy()
    arrival[..., 2] *= -1.0
    return arrival


def aoa_steering_vector(
    arrival_uvw: np.ndarray,
    beam_id: int,
    *,
    phasecal: np.ndarray | None = None,
    antpos_m: np.ndarray | None = None,
    wavelength_m: float = pc.wavelength,
) -> np.ndarray:
    """Return module voltage multipliers for the catalogue AoA convention."""
    arrival = np.asarray(arrival_uvw, dtype=float)
    if arrival.shape != (3,):
        raise ValueError("arrival_uvw must have shape (3,)")
    if not np.all(np.isfinite(arrival)):
        raise ValueError("arrival_uvw must be finite")
    norm = float(np.linalg.norm(arrival))
    if not np.isclose(norm, 1.0, rtol=0.0, atol=1e-6):
        raise ValueError("arrival_uvw must be a unit vector")
    if wavelength_m <= 0.0:
        raise ValueError("wavelength_m must be positive")

    positions = np.asarray(pint.get_antpos() if antpos_m is None else antpos_m, dtype=float)
    calibrations = np.asarray(
        pint.get_phasecal() if phasecal is None else phasecal, dtype=float
    )
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("antpos_m must have shape (n_modules, 3)")
    if calibrations.ndim == 2:
        if not 0 <= int(beam_id) < calibrations.shape[0]:
            raise ValueError(f"beam_id {beam_id} is outside the phase calibration")
        calibrations = calibrations[int(beam_id)]
    if calibrations.shape != (positions.shape[0],):
        raise ValueError("phasecal must contain one phase per receiver module")

    geometric_phase = (2.0 * np.pi / float(wavelength_m)) * (positions @ arrival)
    return np.exp(1j * (calibrations + geometric_phase))


def coherent_add_modules(
    voltage: np.ndarray,
    arrival_uvw: np.ndarray,
    beam_id: int,
    *,
    module_axis: int = 0,
    module_voltage_gain: np.ndarray | None = None,
    module_noise_variance: np.ndarray | None = None,
    normalization: str = "sum",
    phasecal: np.ndarray | None = None,
    antpos_m: np.ndarray | None = None,
    wavelength_m: float = pc.wavelength,
) -> np.ndarray:
    """Gain-correct, AoA-steer, and coherently add receiver-module voltages.

    ``normalization='sum'`` preserves the physical coherent voltage sum and gives
    an ideal power-SNR gain equal to the number of equal-noise modules.  The
    ``weighted_mean`` option is useful when a scale-invariant downstream fit uses
    explicitly supplied per-module noise variances.
    """
    values = np.asarray(voltage, dtype=np.complex128)
    values = np.moveaxis(values, module_axis, 0)
    steering = aoa_steering_vector(
        arrival_uvw,
        beam_id,
        phasecal=phasecal,
        antpos_m=antpos_m,
        wavelength_m=wavelength_m,
    )
    if values.shape[0] != len(steering):
        raise ValueError("the module axis does not match the receiver geometry")

    if module_voltage_gain is None:
        gain = np.ones(len(steering), dtype=float)
    else:
        gain = np.asarray(module_voltage_gain, dtype=float)
        if gain.shape != steering.shape or np.any(~np.isfinite(gain)) or np.any(gain <= 0.0):
            raise ValueError("module_voltage_gain must be finite and positive per module")
    calibrated = values / gain.reshape((-1,) + (1,) * (values.ndim - 1))

    if module_noise_variance is None:
        precision = np.ones(len(steering), dtype=float)
    else:
        variance = np.asarray(module_noise_variance, dtype=float)
        if variance.shape != steering.shape or np.any(~np.isfinite(variance)) or np.any(variance <= 0.0):
            raise ValueError("module_noise_variance must be finite and positive per module")
        precision = gain**2 / variance
    multipliers = precision * steering
    result = np.sum(
        calibrated * multipliers.reshape((-1,) + (1,) * (values.ndim - 1)), axis=0
    )
    if normalization == "sum":
        return result
    if normalization == "weighted_mean":
        return result / np.sum(precision)
    raise ValueError("normalization must be 'sum' or 'weighted_mean'")


def coherent_power_gain(steered_module_values: np.ndarray, module_axis: int = -1) -> np.ndarray:
    """Return ``|sum z_m|^2 / sum |z_m|^2`` for already-steered values."""
    values = np.asarray(steered_module_values, dtype=np.complex128)
    numerator = np.abs(np.sum(values, axis=module_axis)) ** 2
    denominator = np.sum(np.abs(values) ** 2, axis=module_axis)
    return numerator / np.maximum(denominator, np.finfo(float).tiny)
