from __future__ import annotations

import h5py
import numpy as np

from export_pansy_catalogue import (
    _init_level2,
    _init_level3,
    append_level2,
    append_level3,
    validate_month,
)
from orbit_metadata_table import EVENT_DTYPE, PATH_DTYPE


def test_monthly_level2_level3_join_and_units(tmp_path):
    events = np.zeros(2, dtype=EVENT_DTYPE)
    events["sample_idx"] = [1_700_000_000_100_000, 1_700_000_001_200_000]
    events["v_g_km_s"] = [42.0, 43.0]
    events["radiant_ra_deg"] = [10.0, 11.0]
    events["radiant_dec_deg"] = [-20.0, -21.0]
    events["radiant_ecliptic_lon_deg"] = [30.0, 31.0]
    events["radiant_ecliptic_lat_deg"] = [-5.0, -6.0]
    events["radiant_sun_ecliptic_lon_deg"] = [100.0, 101.0]
    events["radiant_sun_ecliptic_lat_deg"] = [-5.0, -6.0]
    events["initial_state_gcrs_m_mps"][:, 3] = [40_000.0, 41_000.0]
    events["fit_parameter_covariance"][:, :6, :6] = np.eye(6)
    events["kepler"] = np.arange(14, dtype=np.float32).reshape(2, 7)
    events["kepler_std"] = 0.1
    events["kepler_covariance"] = np.eye(7)
    events["selected_hypothesis"] = b"k0"

    paths = np.zeros(3, dtype=PATH_DTYPE)
    paths["sample_idx"] = [events["sample_idx"][0], events["sample_idx"][0], events["sample_idx"][1]]
    paths["t_rel_s"] = [0.0, 0.008, 0.0]
    paths["position_enu_km"] = [[3.0, 4.0, 0.0], [0.0, 0.0, 6.0], [0.0, 8.0, 15.0]]
    paths["snr"] = [10.0, 100.0, 1.0]
    paths["doppler_mps"] = [1.0, 2.0, 3.0]
    paths["beam_id"] = [0, 0, 1]
    paths["selection_keep"] = [True, False, True]

    level2_path = tmp_path / "level2.h5"
    level3_path = tmp_path / "level3.h5"
    with _init_level2(level2_path, "1", "test") as level2:
        append_level2(level2, events, paths)
    with _init_level3(level3_path, "1", "test") as level3:
        append_level3(level3, events)

    assert validate_month(level2_path, level3_path) == {"events": 2, "measurements": 3}
    with h5py.File(level2_path, "r") as h:
        np.testing.assert_array_equal(h["events/measurement_count"][()], [2, 1])
        np.testing.assert_allclose(h["measurements/range_km"][()], [5.0, 6.0, 17.0])
        np.testing.assert_allclose(h["measurements/snr_db"][()], [10.0, 20.0, 0.0])
        np.testing.assert_array_equal(
            h["measurements/time_sample_idx"][()],
            [events["sample_idx"][0], events["sample_idx"][0] + 8000, events["sample_idx"][1]],
        )
    with h5py.File(level3_path, "r") as h:
        np.testing.assert_allclose(h["events/solar_longitude_deg"][()], [280.0, 281.0])
        np.testing.assert_allclose(h["events/sun_centered_ecliptic_longitude_deg"][()], [290.0, 290.0])
        np.testing.assert_allclose(h["events/semi_major_axis_AU"][()], [0.0, 7.0])
        np.testing.assert_allclose(h["events/v0_km_s"][()], [0.0, 0.0])
        assert h["events/geocentric_state_gcrs_m_mps"].shape == (2, 6)
