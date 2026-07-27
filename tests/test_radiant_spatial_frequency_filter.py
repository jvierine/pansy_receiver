import healpy as hp
import numpy as np

from radiant_spatial_frequency_filter import (
    iterative_positive_spherical_harmonic_split,
    spherical_harmonic_weights,
)


def test_spherical_harmonic_weights_are_bounded_and_preserve_dc():
    weights = spherical_harmonic_weights(23, 1.0, 5.0, 3.0, 7.0)
    assert np.all((weights >= 0.0) & (weights <= 1.0))
    assert weights[hp.Alm.getidx(23, 0, 0)] == 1.0


def test_iterative_positive_split_returns_positive_compact_excess():
    nside = 8
    npix = hp.nside2npix(nside)
    raw = np.full(npix, 10.0)
    source_pixel = hp.ang2pix(nside, 45.0, -20.0, lonlat=True)
    raw[source_pixel] += 100.0

    low, high, weights, history = iterative_positive_spherical_harmonic_split(
        raw,
        nside,
        zonal_pass=1.0,
        zonal_taper=4.0,
        meridional_pass=2.0,
        meridional_taper=5.0,
        iterations=4,
        positive_fraction=0.5,
        transform_iterations=1,
    )

    assert low.shape == raw.shape
    assert high.shape == raw.shape
    assert weights.ndim == 1
    assert history.shape == (4, 3)
    assert np.all(np.isfinite(low))
    assert np.all(high >= 0.0)
    assert np.argmax(high) == source_pixel
    assert high[source_pixel] > 50.0
    assert np.all(history[:, 0] >= 0.0)
    assert np.allclose(history[:, 1], 0.5 * history[:, 0])
