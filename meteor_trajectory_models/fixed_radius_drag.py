"""Fixed-radius atmospheric-drag meteor trajectory model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np
from scipy.integrate import solve_ivp

from .ceplecha import HeightFunction, local_up_height


DensityFunction = Callable[[float], float]


@dataclass(frozen=True)
class FixedRadiusDragResult:
    """Integrated fixed-radius drag trajectory state."""

    time_s: np.ndarray
    position_m: np.ndarray
    velocity_mps: np.ndarray
    radius_m: float
    success: bool
    message: str

    @property
    def speed_mps(self) -> np.ndarray:
        return np.linalg.norm(self.velocity_mps, axis=1)


def _as_vector3(values: Iterable[float], name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {arr.shape}")
    return arr


def integrate_fixed_radius_drag(
    position0_m: Iterable[float],
    velocity0_mps: Iterable[float],
    radius_m: float,
    atmosphere_density: DensityFunction,
    *,
    meteoroid_density_kg_m3: float = 1000.0,
    t_span_s: tuple[float, float] = (0.0, 0.5),
    sample_dt_s: float = 1e-3,
    height_function: HeightFunction = local_up_height,
    rtol: float = 1e-8,
    atol: float = 1e-10,
) -> FixedRadiusDragResult:
    """Integrate a fixed-radius meteoroid under atmospheric drag.

    This is the ``msis_trajectory``-type model used when the measurements
    justify deceleration but not an additional ablation/shrinking-radius
    parameter.  The atmospheric density is supplied by the caller, so projects
    can use MSIS, a tabulated profile, or a test profile without coupling this
    package to one atmosphere library.
    """

    position0 = _as_vector3(position0_m, "position0_m")
    velocity0 = _as_vector3(velocity0_mps, "velocity0_mps")
    radius = float(radius_m)
    if radius <= 0.0:
        raise ValueError("radius_m must be positive")
    if meteoroid_density_kg_m3 <= 0.0:
        raise ValueError("meteoroid_density_kg_m3 must be positive")
    if sample_dt_s <= 0.0:
        raise ValueError("sample_dt_s must be positive")

    y0 = np.concatenate([position0, velocity0])
    t0, t1 = map(float, t_span_s)
    n_samples = max(2, int(np.floor(abs(t1 - t0) / sample_dt_s)) + 1)
    t_eval = np.linspace(t0, t1, n_samples)

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        pos = y[:3]
        vel = y[3:6]
        speed = float(np.linalg.norm(vel))
        if speed <= 0.0:
            return np.zeros(6, dtype=np.float64)
        height = float(height_function(pos))
        rho_air = max(float(atmosphere_density(height)), 0.0)
        accel = -(3.0 / 4.0) * rho_air / (meteoroid_density_kg_m3 * radius) * speed * vel
        return np.concatenate([vel, accel])

    solution = solve_ivp(rhs, (t0, t1), y0, t_eval=t_eval, rtol=rtol, atol=atol)
    state = solution.y.T
    return FixedRadiusDragResult(
        time_s=solution.t,
        position_m=state[:, :3],
        velocity_mps=state[:, 3:6],
        radius_m=radius,
        success=bool(solution.success),
        message=str(solution.message),
    )
