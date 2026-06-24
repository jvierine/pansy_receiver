"""Ceplecha-style meteoroid drag and ablation trajectory model.

This module is based on the coupled drag/ablation prototype in
``/Users/jvi019/src/meteor/alpha_beta_est.py``.  The old implementation used a
first-order Euler loop.  This package version keeps the same physical model but
uses :func:`scipy.integrate.solve_ivp` so catalogue projects can share one
well-defined implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np
from scipy.integrate import solve_ivp


DensityFunction = Callable[[float], float]
HeightFunction = Callable[[np.ndarray], float]


@dataclass(frozen=True)
class CeplechaResult:
    """Integrated trajectory state."""

    time_s: np.ndarray
    position_m: np.ndarray
    velocity_mps: np.ndarray
    radius_m: np.ndarray
    mass_kg: np.ndarray
    success: bool
    message: str

    @property
    def speed_mps(self) -> np.ndarray:
        return np.linalg.norm(self.velocity_mps, axis=1)


def spherical_mass(radius_m: float, density_kg_m3: float) -> float:
    """Return the mass of a compact sphere."""

    return (4.0 / 3.0) * np.pi * density_kg_m3 * radius_m**3


def spherical_radius(mass_kg: np.ndarray | float, density_kg_m3: float) -> np.ndarray:
    """Return compact-sphere radius from mass and density."""

    mass = np.asarray(mass_kg, dtype=np.float64)
    return np.maximum(mass / ((4.0 / 3.0) * np.pi * density_kg_m3), 0.0) ** (1.0 / 3.0)


def local_up_height(position_m: np.ndarray) -> float:
    """Default height function for local ENU coordinates."""

    return float(position_m[2])


def _as_vector3(values: Iterable[float], name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {arr.shape}")
    return arr


def integrate_ceplecha(
    position0_m: Iterable[float],
    velocity0_mps: Iterable[float],
    radius0_m: float,
    atmosphere_density: DensityFunction,
    *,
    meteoroid_density_kg_m3: float = 1000.0,
    ablation_sigma_kg_j: float = 1e-8,
    t_span_s: tuple[float, float] = (0.0, 0.5),
    sample_dt_s: float = 1e-3,
    height_function: HeightFunction = local_up_height,
    rtol: float = 1e-8,
    atol: float = 1e-10,
) -> CeplechaResult:
    """Integrate a spherical meteoroid with drag and ablation.

    Parameters
    ----------
    position0_m, velocity0_mps
        Initial Cartesian state.  By default the third coordinate is treated as
        height above the local tangent plane, matching the historical MAARSY
        prototype.
    radius0_m
        Initial meteoroid radius in meters.
    atmosphere_density
        Callable returning atmospheric mass density in kg m^-3 from height in
        meters.
    meteoroid_density_kg_m3
        Bulk meteoroid density used to convert mass to radius.
    ablation_sigma_kg_j
        Mass-loss coefficient in kg J^-1.  The historical Sanya/MAARSY prototype
        commonly used 1e-8 kg J^-1 as a fixed value.
    t_span_s
        Integration interval.
    sample_dt_s
        Output sampling interval.
    height_function
        Callable mapping position vector to height in meters.

    Notes
    -----
    The model equations use the molecular free-flow drag convention

    ``dv/dt = -(3/4) rho_a(h) / (rho_m r) |v| v``

    and

    ``dm/dt = -(1/2) rho_a(h) pi r^2 sigma |v|^3``.

    The acceleration expression is the vector form of the original scalar-speed
    update, with the direction allowed to follow the velocity vector.
    """

    position0 = _as_vector3(position0_m, "position0_m")
    velocity0 = _as_vector3(velocity0_mps, "velocity0_mps")
    radius0 = float(radius0_m)
    if radius0 <= 0.0:
        raise ValueError("radius0_m must be positive")
    if meteoroid_density_kg_m3 <= 0.0:
        raise ValueError("meteoroid_density_kg_m3 must be positive")
    if sample_dt_s <= 0.0:
        raise ValueError("sample_dt_s must be positive")

    mass0 = spherical_mass(radius0, meteoroid_density_kg_m3)
    y0 = np.concatenate([position0, velocity0, [mass0]])
    t0, t1 = map(float, t_span_s)
    n_samples = max(2, int(np.floor(abs(t1 - t0) / sample_dt_s)) + 1)
    t_eval = np.linspace(t0, t1, n_samples)

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        pos = y[:3]
        vel = y[3:6]
        mass = max(float(y[6]), 0.0)
        speed = float(np.linalg.norm(vel))
        if mass <= 0.0 or speed <= 0.0:
            return np.zeros(7, dtype=np.float64)

        radius = float(spherical_radius(mass, meteoroid_density_kg_m3))
        height = float(height_function(pos))
        rho_air = max(float(atmosphere_density(height)), 0.0)

        accel = -(3.0 / 4.0) * rho_air / (meteoroid_density_kg_m3 * radius) * speed * vel
        dm_dt = -0.5 * rho_air * np.pi * radius**2 * ablation_sigma_kg_j * speed**3
        return np.concatenate([vel, accel, [dm_dt]])

    def mass_depleted(_t: float, y: np.ndarray) -> float:
        return float(y[6])

    mass_depleted.terminal = True
    mass_depleted.direction = -1.0

    solution = solve_ivp(
        rhs,
        (t0, t1),
        y0,
        t_eval=t_eval,
        events=mass_depleted,
        rtol=rtol,
        atol=atol,
    )

    state = solution.y.T
    mass = np.maximum(state[:, 6], 0.0)
    radius = spherical_radius(mass, meteoroid_density_kg_m3)
    return CeplechaResult(
        time_s=solution.t,
        position_m=state[:, :3],
        velocity_mps=state[:, 3:6],
        radius_m=radius,
        mass_kg=mass,
        success=bool(solution.success),
        message=str(solution.message),
    )
