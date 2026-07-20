"""Quality metrics for the fixed 50-event three-pulse fit benchmark."""

from __future__ import annotations

import numpy as np


QUALITY_SCALES = {
    "up_rms_km": 0.056,
    "ns_rms_km": 0.200,
    "ew_rms_km": 0.250,
    "doppler_rms_mps": 129.0,
    "acceleration_rms_km_s2": 3.0,
}


def rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.inf
    return float(np.sqrt(np.mean(values**2)))


def variance_weighted_score(metrics: dict[str, float]) -> float:
    """Return the sum of squared RMS values normalized by benchmark variance."""
    return float(
        sum((float(metrics[name]) / scale) ** 2 for name, scale in QUALITY_SCALES.items())
    )


def event_fit_quality(
    measured_position_km: np.ndarray,
    model_position_km: np.ndarray,
    position_keep: np.ndarray,
    triplet_time_s: np.ndarray,
    measured_doppler_mps: np.ndarray,
    measured_acceleration_mps2: np.ndarray,
    velocity_keep: np.ndarray,
    acceleration_keep: np.ndarray,
    trajectory_time_s: np.ndarray,
    model_doppler_mps: np.ndarray,
    model_acceleration_mps2: np.ndarray,
) -> dict[str, float]:
    position_keep = np.asarray(position_keep, dtype=bool)
    position_residual = (
        np.asarray(measured_position_km, dtype=float)
        - np.asarray(model_position_km, dtype=float)
    )[position_keep]
    position_rms = [rms(position_residual[:, component]) for component in range(3)]

    triplet_time_s = np.asarray(triplet_time_s, dtype=float)
    velocity_residual = np.asarray(measured_doppler_mps, dtype=float) - np.interp(
        triplet_time_s,
        np.asarray(trajectory_time_s, dtype=float),
        np.asarray(model_doppler_mps, dtype=float),
    )
    acceleration_residual = np.asarray(measured_acceleration_mps2, dtype=float) - np.interp(
        triplet_time_s,
        np.asarray(trajectory_time_s, dtype=float),
        np.asarray(model_acceleration_mps2, dtype=float),
    )
    metrics = {
        "up_rms_km": position_rms[2],
        "ns_rms_km": position_rms[1],
        "ew_rms_km": position_rms[0],
        "doppler_rms_mps": rms(velocity_residual[np.asarray(velocity_keep, dtype=bool)]),
        "acceleration_rms_km_s2": rms(
            acceleration_residual[np.asarray(acceleration_keep, dtype=bool)]
        )
        / 1e3,
    }
    metrics["variance_weighted_score"] = variance_weighted_score(metrics)
    return metrics
