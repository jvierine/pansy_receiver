"""Radius-grid refinement helpers for dynamic-mass likelihood profiles."""

from __future__ import annotations

import numpy as np


def adaptive_profile_radii(
    radius_um: np.ndarray,
    delta_chi2: np.ndarray,
    *,
    spacing_dex: float = 0.02,
    selection_delta_chi2: float = 25.0,
    interval_threshold: float = 3.841458820694124,
) -> np.ndarray:
    """Return extra log-spaced radii near profile minima and crossings."""
    radius_um = np.asarray(radius_um, dtype=float)
    delta_chi2 = np.asarray(delta_chi2, dtype=float)
    if len(radius_um) < 3 or spacing_dex <= 0.0:
        return np.empty(0, dtype=float)
    order = np.argsort(radius_um)
    radius_um = radius_um[order]
    delta_chi2 = delta_chi2[order]
    log_radius = np.log10(radius_um)
    finite = np.isfinite(delta_chi2)
    if not np.any(finite):
        return np.empty(0, dtype=float)

    best = int(np.nanargmin(delta_chi2))
    intervals: list[tuple[float, float]] = []

    # Refine the connected low-chi-square basin. For an open tail, stop a few
    # coarse cells beyond the minimum; the distant flat tail remains coarse.
    selected = finite & (delta_chi2 <= selection_delta_chi2)
    left = best
    while left > 0 and selected[left - 1]:
        left -= 1
    right = best
    while right + 1 < len(radius_um) and selected[right + 1]:
        right += 1
    basin_left = max(0, left - 1)
    basin_right = min(len(radius_um) - 1, right + 1)
    if left == 0:
        basin_left = max(0, best - 2)
    if right == len(radius_um) - 1:
        basin_right = min(len(radius_um) - 1, best + 2)
    intervals.append((log_radius[basin_left], log_radius[basin_right]))

    # Preserve disconnected secondary minima and resolve both the 95% and
    # diagnostic-selection threshold crossings.
    for index in range(1, len(radius_um) - 1):
        if (
            finite[index - 1 : index + 2].all()
            and delta_chi2[index] <= selection_delta_chi2
            and delta_chi2[index] <= delta_chi2[index - 1]
            and delta_chi2[index] <= delta_chi2[index + 1]
            and (
                delta_chi2[index] < delta_chi2[index - 1]
                or delta_chi2[index] < delta_chi2[index + 1]
            )
        ):
            intervals.append((log_radius[index - 1], log_radius[index + 1]))
    for threshold in (interval_threshold, selection_delta_chi2):
        for index in range(len(radius_um) - 1):
            if not (finite[index] and finite[index + 1]):
                continue
            if (delta_chi2[index] - threshold) * (
                delta_chi2[index + 1] - threshold
            ) <= 0.0:
                intervals.append((log_radius[index], log_radius[index + 1]))

    new_log_radius = []
    for lower, upper in intervals:
        count = max(2, int(np.ceil((upper - lower) / spacing_dex)) + 1)
        new_log_radius.extend(np.linspace(lower, upper, count))
    candidates = np.unique(np.power(10.0, np.asarray(new_log_radius, dtype=float)))
    if len(candidates) > 1:
        distinct = np.r_[True, np.diff(np.log10(candidates)) > 1e-8]
        candidates = candidates[distinct]
    duplicate = np.any(
        np.isclose(candidates[:, None], radius_um[None, :], rtol=1e-10, atol=0.0),
        axis=1,
    )
    return candidates[~duplicate]
