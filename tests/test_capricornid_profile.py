import numpy as np
import pytest

from plot_capricornid_conjugate_stream import (
    PASSAGES,
    ecliptic_to_equatorial_deg,
    radiant_color_values,
    validate_passage_sources,
)


def test_ecliptic_to_equatorial_cardinal_directions():
    ra_deg, dec_deg = ecliptic_to_equatorial_deg(
        np.asarray([0.0, 90.0]),
        np.asarray([0.0, 0.0]),
    )

    np.testing.assert_allclose(ra_deg, [0.0, 90.0], atol=1e-8)
    np.testing.assert_allclose(dec_deg, [0.0, 23.4392911], atol=1e-8)


def test_radiant_color_values_use_documented_kepler_columns():
    rows = np.zeros(
        2,
        dtype=[("vg", "f8"), ("kepler", "f8", (7,))],
    )
    rows["vg"] = [21.0, 25.0]
    rows["kepler"] = np.asarray(
        [
            [2.5, 0.74, 7.6, 132.0, 264.0, 275.0, 0.63],
            [3.0, 0.78, 9.1, 313.0, 266.0, 273.0, 0.66],
        ]
    )

    np.testing.assert_allclose(radiant_color_values(rows, "vg"), [21.0, 25.0])
    np.testing.assert_allclose(radiant_color_values(rows, "a"), [2.5, 3.0])
    np.testing.assert_allclose(radiant_color_values(rows, "inclination"), [7.6, 9.1])
    np.testing.assert_allclose(radiant_color_values(rows, "omega"), [264.0, 266.0])
    np.testing.assert_allclose(radiant_color_values(rows, "q"), [0.63, 0.66])


def test_passage_source_validation_rejects_missing_pansy_catalogue():
    rows = np.zeros(3, dtype=[("dataset", "U8")])
    rows["dataset"] = "MAARSY"

    with pytest.raises(RuntimeError, match="missing PANSY"):
        validate_passage_sources(rows, PASSAGES[1])


def test_passage_source_validation_accepts_combined_catalogue():
    rows = np.zeros(3, dtype=[("dataset", "U8")])
    rows["dataset"] = ["PANSY", "MAARSY", "PANSY"]

    validate_passage_sources(rows, PASSAGES[1])
