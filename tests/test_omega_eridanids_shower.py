import numpy as np

from plot_omega_eridanids_shower import (
    ang2pix_ring,
    circular_mean_std_deg,
    ecliptic_to_equatorial_deg,
)


def test_ang2pix_ring_matches_radview_reference_pixels():
    # Avoid points exactly on HEALPix boundaries, where implementations may
    # choose either adjacent pixel because of floating-point rounding.
    lon = np.asarray([-73.5, -70.2, 1.0, 181.0, 45.0, 275.0])
    lat = np.asarray([-49.53125, -48.6, 1.0, -1.0, 20.0, -20.0])
    expected = np.asarray([10861, 10754, 5952, 6272, 3920, 8225])

    np.testing.assert_array_equal(ang2pix_ring(32, lon, lat), expected)


def test_circular_statistics_cross_zero_degrees():
    mean, std = circular_mean_std_deg(np.asarray([359.0, 0.0, 1.0]))

    assert np.isclose(mean, 0.0, atol=1e-12)
    assert np.isclose(std, 1.0)


def test_ecliptic_to_equatorial_preserves_vernal_equinox():
    ra, dec = ecliptic_to_equatorial_deg(np.asarray([0.0]), np.asarray([0.0]))

    np.testing.assert_allclose(ra, 0.0, atol=1e-12)
    np.testing.assert_allclose(dec, 0.0, atol=1e-12)
