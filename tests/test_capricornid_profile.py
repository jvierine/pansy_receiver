import numpy as np

from plot_capricornid_conjugate_stream import ecliptic_to_equatorial_deg


def test_ecliptic_to_equatorial_cardinal_directions():
    ra_deg, dec_deg = ecliptic_to_equatorial_deg(
        np.asarray([0.0, 90.0]),
        np.asarray([0.0, 0.0]),
    )

    np.testing.assert_allclose(ra_deg, [0.0, 90.0], atol=1e-8)
    np.testing.assert_allclose(dec_deg, [0.0, 23.4392911], atol=1e-8)
