"""Spherical-harmonic filtering utilities for HEALPix radiant maps."""

from __future__ import annotations

import healpy as hp
import numpy as np
from matplotlib.colors import LogNorm, SymLogNorm


def required_lmax(
    nside: int,
    zonal_pass: float,
    zonal_taper: float,
    meridional_pass: float,
    meridional_taper: float,
) -> int:
    """Return the highest degree at which the separable taper can be nonzero."""
    supported = np.ceil(zonal_pass + zonal_taper + meridional_pass + meridional_taper)
    return min(3 * int(nside) - 1, max(0, int(supported)))


def raised_cosine_lowpass_weights(frequency: np.ndarray, pass_bin: float, taper_width: float) -> np.ndarray:
    """Return raised-cosine low-pass weights evaluated at nonnegative frequencies."""
    frequency = np.asarray(frequency, dtype=np.float64)
    pass_bin = max(float(pass_bin), 0.0)
    taper_width = max(float(taper_width), 1e-6)
    stop_bin = pass_bin + taper_width
    weights = np.ones(frequency.shape, dtype=np.float64)
    weights[frequency >= stop_bin] = 0.0
    transition = (frequency > pass_bin) & (frequency < stop_bin)
    if np.any(transition):
        x = (frequency[transition] - pass_bin) / taper_width
        weights[transition] = 0.5 * (1.0 + np.cos(np.pi * x))
    return weights


def spherical_harmonic_weights(
    lmax: int,
    zonal_pass: float,
    zonal_taper: float,
    meridional_pass: float,
    meridional_taper: float,
) -> np.ndarray:
    """Filter packed ``a_lm`` coefficients by order and latitudinal degree."""
    ell, emm = hp.Alm.getlm(int(lmax))
    zonal = raised_cosine_lowpass_weights(emm, zonal_pass, zonal_taper)
    meridional = raised_cosine_lowpass_weights(ell - emm, meridional_pass, meridional_taper)
    return zonal * meridional


def spherical_harmonic_split(
    raw_count: np.ndarray,
    nside: int,
    zonal_pass: float,
    zonal_taper: float,
    meridional_pass: float,
    meridional_taper: float,
    *,
    lmax: int | None = None,
    iterations: int = 3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Split a RING-ordered HEALPix count map into low- and high-frequency parts."""
    nside = int(nside)
    values = np.nan_to_num(np.asarray(raw_count, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    if values.shape != (hp.nside2npix(nside),):
        raise ValueError(f"expected a HEALPix map with {hp.nside2npix(nside)} pixels")
    if lmax is None:
        lmax = required_lmax(
            nside,
            zonal_pass,
            zonal_taper,
            meridional_pass,
            meridional_taper,
        )
    lmax = int(lmax)
    alm = hp.map2alm(values, lmax=lmax, iter=int(iterations), pol=False)
    weights = spherical_harmonic_weights(
        lmax,
        zonal_pass,
        zonal_taper,
        meridional_pass,
        meridional_taper,
    )
    lowpass = hp.alm2map(alm * weights, nside=nside, lmax=lmax, pol=False)
    highpass = values - lowpass
    return lowpass, highpass, weights


def iterative_positive_spherical_harmonic_split(
    raw_count: np.ndarray,
    nside: int,
    zonal_pass: float,
    zonal_taper: float,
    meridional_pass: float,
    meridional_taper: float,
    *,
    lmax: int | None = None,
    iterations: int = 12,
    positive_fraction: float = 0.5,
    transform_iterations: int = 3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Isolate positive compact excess while iteratively refitting the background."""
    values = np.nan_to_num(np.asarray(raw_count, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    if iterations < 1:
        raise ValueError("iterations must be at least one")
    if not 0.0 < positive_fraction <= 1.0:
        raise ValueError("positive_fraction must be in (0, 1]")

    working = values.copy()
    history = np.empty((int(iterations), 3), dtype=np.float64)
    weights = None
    for iteration in range(int(iterations)):
        smooth, residual, weights = spherical_harmonic_split(
            working,
            nside,
            zonal_pass,
            zonal_taper,
            meridional_pass,
            meridional_taper,
            lmax=lmax,
            iterations=transform_iterations,
        )
        positive = np.clip(residual, 0.0, None)
        extracted = positive_fraction * positive
        working = np.clip(working - extracted, 0.0, None)
        history[iteration] = (
            float(np.sum(positive)),
            float(np.sum(extracted)),
            float(np.max(positive, initial=0.0)),
        )

    background, _, weights = spherical_harmonic_split(
        working,
        nside,
        zonal_pass,
        zonal_taper,
        meridional_pass,
        meridional_taper,
        lmax=lmax,
        iterations=transform_iterations,
    )
    background = np.clip(background, 0.0, None)
    positive_highpass = np.clip(values - background, 0.0, None)
    return background, positive_highpass, weights, history


def positive_log_norm(values: np.ndarray, percentile: float = 99.8) -> LogNorm:
    positive = np.asarray(values, dtype=np.float64)
    positive = positive[np.isfinite(positive) & (positive > 0.0)]
    vmax = max(2.0, float(np.nanpercentile(positive, percentile))) if positive.size else 2.0
    return LogNorm(vmin=1.0, vmax=vmax)


def signed_residual_norm(values: np.ndarray, percentile: float = 99.5) -> SymLogNorm:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    scale = float(np.nanpercentile(np.abs(finite), percentile)) if finite.size else 1.0
    scale = max(scale, 1.0)
    return SymLogNorm(linthresh=1.0, linscale=0.8, vmin=-scale, vmax=scale, base=10)
