#!/usr/bin/env python3
"""Focused tests for adaptive fixed-radius profile sampling."""

import numpy as np

from mass_profile_grid import adaptive_profile_radii


def test_refines_closed_likelihood_basin_and_crossings():
    radius = np.logspace(0.0, 4.0, 41)
    delta = ((np.log10(radius) - 2.0) / 0.12) ** 2
    extra = adaptive_profile_radii(radius, delta, spacing_dex=0.02)
    combined = np.sort(np.r_[radius, extra])
    relevant = combined[(combined >= 50.0) & (combined <= 200.0)]
    assert len(extra) > 0
    assert np.max(np.diff(np.log10(relevant))) <= 0.0200001


def test_does_not_densify_an_entire_open_tail():
    radius = np.logspace(0.0, 4.0, 41)
    delta = np.maximum(0.0, 4.0 - 5.0 * (np.log10(radius) - 2.0))
    extra = adaptive_profile_radii(radius, delta, spacing_dex=0.02)
    assert len(extra) < 40
    assert np.count_nonzero(extra > 1e3) == 0


def test_overlapping_refinement_intervals_do_not_duplicate_radii():
    radius = np.logspace(0.0, 4.0, 41)
    delta = ((np.log10(radius) - 2.0) / 0.2) ** 2
    extra = adaptive_profile_radii(radius, delta, spacing_dex=0.02)
    assert np.min(np.diff(np.log10(extra))) > 1e-8


if __name__ == "__main__":
    test_refines_closed_likelihood_basin_and_crossings()
    test_does_not_densify_an_entire_open_tail()
    test_overlapping_refinement_intervals_do_not_duplicate_radii()
    print("adaptive profile tests passed")
