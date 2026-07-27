import numpy as np

from compare_ser_radiant_inclination_geometry import (
    HIGH_INCLINATION_DEG,
    LOW_INCLINATION_DEG,
    circular_mean_deg,
    geocentric_radiants,
    matched_inclination_clouds,
)


def test_matched_clouds_change_only_mean_inclination():
    elements = np.asarray(
        [
            [3.4, 0.7, 90.0, 290.0, 315.0, 45.0, 0.9],
            [3.5, 0.72, 92.0, 291.0, 320.0, 40.0, 0.95],
        ]
    )
    clouds = matched_inclination_clouds(elements)
    assert np.allclose(clouds[0, :, [0, 1, 3, 4, 5, 6]], clouds[1, :, [0, 1, 3, 4, 5, 6]])
    assert np.isclose(np.mean(clouds[0, :, 2]), HIGH_INCLINATION_DEG)
    assert np.isclose(np.mean(clouds[1, :, 2]), LOW_INCLINATION_DEG)


def test_circular_mean_handles_longitude_wrap():
    assert np.isclose(circular_mean_deg([359.0, 1.0]), 0.0, atol=1e-10)


def test_radiant_projection_is_finite():
    elements = np.asarray([[3.44, 0.712, 91.22, 290.12, 317.12, 42.45, 0.899]])
    radiant = geocentric_radiants(elements, np.asarray([109.5]))
    assert radiant.shape == (1, 2)
    assert np.all(np.isfinite(radiant))
    assert -180.0 <= radiant[0, 0] <= 180.0
    assert -90.0 <= radiant[0, 1] <= 90.0
