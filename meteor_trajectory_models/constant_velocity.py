"""Constant-velocity meteor trajectory model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class ConstantVelocityResult:
    """Integrated constant-velocity trajectory state."""

    time_s: np.ndarray
    position_m: np.ndarray
    velocity_mps: np.ndarray

    @property
    def speed_mps(self) -> np.ndarray:
        return np.linalg.norm(self.velocity_mps, axis=1)


def _as_vector3(values: Iterable[float], name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {arr.shape}")
    return arr


def integrate_constant_velocity(
    position0_m: Iterable[float],
    velocity0_mps: Iterable[float],
    *,
    t_span_s: tuple[float, float] = (0.0, 0.5),
    sample_dt_s: float = 1e-3,
) -> ConstantVelocityResult:
    """Evaluate a straight-line constant-velocity trajectory."""

    position0 = _as_vector3(position0_m, "position0_m")
    velocity0 = _as_vector3(velocity0_mps, "velocity0_mps")
    if sample_dt_s <= 0.0:
        raise ValueError("sample_dt_s must be positive")
    t0, t1 = map(float, t_span_s)
    n_samples = max(2, int(np.floor(abs(t1 - t0) / sample_dt_s)) + 1)
    time_s = np.linspace(t0, t1, n_samples)
    dt = (time_s - t0)[:, None]
    position_m = position0[None, :] + dt * velocity0[None, :]
    velocity_mps = np.repeat(velocity0[None, :], n_samples, axis=0)
    return ConstantVelocityResult(time_s=time_s, position_m=position_m, velocity_mps=velocity_mps)
