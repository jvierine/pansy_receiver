from types import SimpleNamespace

import numpy as np

from plot_paper_radiant_results import (
    LATITUDE_BINS,
    LONGITUDE_BINS,
    corrected_flux_histogram,
    fit_zenith_exponent,
    observing_time_samples,
    radiant_healpix_histogram,
    right_side_contour_label_positions,
)
from radiant_visibility import (
    altitude_from_radec_deg,
    ecliptic_lonlat_to_radec_deg,
    radiant_exposure_hours_grid,
    radiant_exposure_hours_points,
    wrap360,
)


def test_observing_time_samples_preserves_partial_bin_duration():
    intervals = np.asarray([[100.0, 1000.0], [1100.0, 1700.0]])
    epoch, hours = observing_time_samples(intervals, cadence_seconds=900.0)
    np.testing.assert_allclose(epoch, [450.0, 1350.0])
    np.testing.assert_allclose(hours, [800.0 / 3600.0, 700.0 / 3600.0])


def test_exposure_grid_matches_direct_altitude_calculation():
    epoch = np.asarray([1_750_000_000.0])
    sun = np.asarray([123.0])
    plot_lon = np.asarray([-120.0, 15.0])
    beta = np.asarray([-25.0, 20.0])
    _, _, exposure = radiant_exposure_hours_grid(epoch, sun, [0.25], plot_lon, beta)

    for iy, beta_deg in enumerate(beta):
        for ix, plot_deg in enumerate(plot_lon):
            lambda_minus_sun = wrap360(-90.0 - plot_deg)
            ra, dec = ecliptic_lonlat_to_radec_deg(wrap360(sun[0] + lambda_minus_sun), beta_deg)
            altitude = altitude_from_radec_deg(ra, dec, epoch[0])
            expected = 0.25 if altitude > 0.0 else 0.0
            assert exposure[iy, ix] == expected


def test_exposure_points_matches_paired_grid_samples():
    epoch = np.asarray([1_750_000_000.0])
    sun = np.asarray([123.0])
    plot_lon = np.asarray([-120.0, 15.0])
    beta = np.asarray([-25.0, 20.0])
    _, _, grid = radiant_exposure_hours_grid(epoch, sun, [0.25], plot_lon, beta)
    points = radiant_exposure_hours_points(epoch, sun, [0.25], plot_lon, beta)
    np.testing.assert_allclose(points, [grid[0, 0], grid[1, 1]])


def test_healpix_histogram_preserves_total_count():
    rows = np.zeros(
        4,
        dtype=[("lambda_minus_sun_deg", "f8"), ("radiant_beta_ecliptic_deg", "f8")],
    )
    rows["lambda_minus_sun_deg"] = [0.0, 0.0, 180.0, 270.0]
    rows["radiant_beta_ecliptic_deg"] = [0.0, 0.0, 30.0, -45.0]
    count = radiant_healpix_histogram(rows, nside=8)
    assert count.shape == (12 * 8**2,)
    assert np.sum(count) == len(rows)


def test_apex_symmetry_fit_recovers_known_exponent():
    dtype = np.dtype(
        [
            ("plot_longitude_deg", "f8"),
            ("radiant_beta_ecliptic_deg", "f8"),
            ("speed_km_s", "f8"),
        ]
    )
    rows = np.zeros(3, dtype=dtype)
    rows["plot_longitude_deg"] = 0.0
    rows["radiant_beta_ecliptic_deg"] = [20.0, 20.0, -20.0]
    rows["speed_km_s"] = 73.0
    cos_z = np.asarray([1.0, 1.0, 0.5])
    exposure = np.ones((LATITUDE_BINS, LONGITUDE_BINS), dtype=np.float64)

    alpha, objective, north_flux, south_flux = fit_zenith_exponent(
        rows,
        cos_z,
        exposure,
        min_cos_z=0.1,
        alpha_bounds=(0.0, 3.0),
    )
    assert abs(alpha - 1.0) < 2e-3
    assert objective < 1e-8
    np.testing.assert_allclose(north_flux, south_flux, rtol=2e-3)
    _, _, flux = corrected_flux_histogram(rows, cos_z, exposure, 0.1, alpha)
    assert np.sum(flux) > 0.0


def test_contour_labels_prefer_dark_right_side():
    xcenters = np.asarray([-120.0, 60.0, 110.0, 145.0])
    ycenters = np.asarray([-20.0, 20.0])
    raw_hist = np.asarray([[0.0, 100.0, 0.0, 20.0], [0.0, 100.0, 50.0, 20.0]])
    contour_set = SimpleNamespace(
        allsegs=[
            [
                np.deg2rad(
                    np.asarray(
                        [
                            [-120.0, -20.0],
                            [60.0, -20.0],
                            [110.0, -20.0],
                            [145.0, -20.0],
                        ]
                    )
                )
            ]
        ]
    )

    positions = right_side_contour_label_positions(contour_set, raw_hist, xcenters, ycenters)

    np.testing.assert_allclose(np.rad2deg(positions[0]), [110.0, -20.0])
