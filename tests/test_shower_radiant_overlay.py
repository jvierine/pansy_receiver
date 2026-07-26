from types import SimpleNamespace

import numpy as np

from shower_radiant_overlay import shower_healpy_points


def test_shower_healpy_points_match_realtime_longitude_convention():
    showers = [
        SimpleNamespace(
            radiant_solar_ecliptic_lon_deg=0.0,
            radiant_ecliptic_lat_deg=12.0,
            code="HEL",
            name="Helion test",
        ),
        SimpleNamespace(
            radiant_solar_ecliptic_lon_deg=270.0,
            radiant_ecliptic_lat_deg=-8.0,
            code="APX",
            name="Apex test",
        ),
        SimpleNamespace(
            radiant_solar_ecliptic_lon_deg=180.0,
            radiant_ecliptic_lat_deg=3.0,
            code="",
            name="Antihelion test",
        ),
    ]

    lon, lat, labels = shower_healpy_points(showers, longitude_offset_deg=90.0)

    np.testing.assert_allclose(lon, [90.0, 0.0, -90.0])
    np.testing.assert_allclose(lat, [12.0, -8.0, 3.0])
    assert labels == ["HEL", "APX", "Antih"]
