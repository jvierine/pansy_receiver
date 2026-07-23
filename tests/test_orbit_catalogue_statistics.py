import sys
from pathlib import Path

import h5py
import numpy as np


PANSY_RECEIVER_ROOT = Path(__file__).resolve().parents[1]
if str(PANSY_RECEIVER_ROOT) not in sys.path:
    sys.path.insert(0, str(PANSY_RECEIVER_ROOT))


def test_solar_longitude_count_density_uses_bin_width():
    from plot_orbit_catalogue_statistics import count_density_from_counts

    counts_by_year = np.asarray([[10, 20], [4, 8]], dtype=np.int32)
    edges_deg = np.asarray([0.0, 2.0, 6.0], dtype=np.float32)

    density_by_year, all_density = count_density_from_counts(counts_by_year, edges_deg)

    np.testing.assert_allclose(density_by_year, [[5.0, 5.0], [2.0, 2.0]])
    np.testing.assert_allclose(all_density, [7.0, 7.0])


def test_solar_longitude_count_rate_uses_measurement_hours():
    from plot_orbit_catalogue_statistics import count_rate_from_counts_and_exposure

    counts_by_year = np.asarray([[10.0, 8.0], [6.0, 4.0]], dtype=np.float32)
    all_counts = np.asarray([16.0, 12.0], dtype=np.float32)
    hours_by_year = np.asarray([[5.0, 4.0], [10.0, 8.0]], dtype=np.float32)
    all_hours = np.asarray([15.0, 12.0], dtype=np.float32)

    rate_by_year, all_rate = count_rate_from_counts_and_exposure(
        counts_by_year,
        all_counts,
        hours_by_year,
        all_hours,
    )

    np.testing.assert_allclose(rate_by_year[0, 0], 2.0)
    assert np.isnan(rate_by_year[0, 1])
    np.testing.assert_allclose(rate_by_year[1], [0.6, 0.5])
    np.testing.assert_allclose(all_rate, [16.0 / 15.0, 1.0])


def test_solar_longitude_exposure_is_clipped_to_catalogue_time_span(tmp_path):
    from plot_orbit_catalogue_statistics import mesomode_exposure_by_solar_longitude

    base = 1_750_000_000.0
    sidecar = tmp_path / "mesomode_intervals.h5"
    with h5py.File(sidecar, "w") as h:
        intervals = h.create_group("intervals")
        intervals.create_dataset("t0_unix", data=[base, base + 200.0])
        intervals.create_dataset("t1_unix", data=[base + 100.0, base + 400.0])

    by_year, total = mesomode_exposure_by_solar_longitude(
        sidecar,
        edges=np.asarray([0.0, 360.0]),
        years=np.asarray([2025]),
        start_unix_s=base + 50.0,
        stop_unix_s=base + 300.0,
    )

    np.testing.assert_allclose(by_year, [[150.0 / 3600.0]], rtol=1e-6)
    np.testing.assert_allclose(total, [150.0 / 3600.0], rtol=1e-6)


def test_height_velocity_histogram_spans_60_to_160_km():
    from plot_orbit_catalogue_statistics import histogram_height_velocity

    counts, height_edges, speed_edges = histogram_height_velocity(
        np.asarray([59.9, 60.1, 159.9, 160.1]),
        np.asarray([20.0, 20.0, 20.0, 20.0]),
    )

    assert height_edges[0] == 60.0
    assert height_edges[-1] == 160.0
    assert speed_edges[0] == 0.0
    assert speed_edges[-1] == 80.0
    assert counts.sum() == 2


def test_height_velocity_histogram_can_use_custom_height_range():
    from plot_orbit_catalogue_statistics import histogram_height_velocity

    counts, height_edges, _speed_edges = histogram_height_velocity(
        np.asarray([49.9, 50.1, 179.9, 180.1]),
        np.asarray([20.0, 20.0, 20.0, 20.0]),
        height_min_km=50.0,
        height_max_km=180.0,
    )

    assert height_edges[0] == 50.0
    assert height_edges[-1] == 180.0
    assert counts.sum() == 2


def test_initial_detection_speed_uses_fitted_local_velocity():
    import orbit_metadata_table as omt
    from plot_orbit_catalogue_statistics import initial_detection_speed_km_s

    events = np.zeros(2, dtype=omt.EVENT_DTYPE)
    events["fit_parameters"][0, 3:6] = [3_000.0, 4_000.0, 12_000.0]
    events["fit_parameters"][1, 3:6] = [np.nan, 0.0, 0.0]

    speed = initial_detection_speed_km_s(events)

    np.testing.assert_allclose(speed[0], 13.0)
    assert np.isnan(speed[1])


def test_measurement_histogram_input_excludes_robust_fit_outliers():
    import orbit_metadata_table as omt
    from plot_orbit_catalogue_statistics import measurement_height_velocity_arrays

    sample_idx = 1_740_000_000_000_000
    events = np.zeros(1, dtype=omt.EVENT_DTYPE)
    events["sample_idx"] = sample_idx
    events["initial_detection_height_km"] = 100.0
    events["v_g_km_s"] = 40.0
    events["radiant_sun_ecliptic_lon_deg"] = 120.0
    events["orbit_solution_type"] = b"dasst_winning_alias"
    events["combined_score"] = 0.5
    events["frac_e_gt_1"] = 0.0
    events["n_uncertainty_samples"] = 3

    paths = np.zeros(2, dtype=omt.PATH_DTYPE)
    paths["sample_idx"] = sample_idx
    paths["position_enu_km"][:, 2] = [100.0, 40.0]
    paths["selection_keep"] = [True, False]

    result = measurement_height_velocity_arrays(
        events,
        paths,
        max_radiant_sigma_deg=None,
        min_sample_idx=sample_idx - 1,
        max_sample_idx=sample_idx + 1,
        max_combined_score=1.5,
        max_frac_e_gt_1=0.5,
        min_uncertainty_samples=3,
        max_initial_state_position_sigma_m=1000.0,
        max_initial_state_radiant_angle_sigma_deg=3.0,
    )

    np.testing.assert_array_equal(result["measurement_height_km"], [100.0])


def test_selected_closest_tx_alias_sample_indices_requires_selected_alias_to_be_tx_nearest():
    import orbit_metadata_table as omt
    from plot_orbit_catalogue_statistics import selected_closest_tx_alias_sample_indices

    events = np.zeros(2, dtype=omt.EVENT_DTYPE)
    events["sample_idx"] = [10, 20]
    events["selected_hypothesis"] = [b"H02", b"H01"]

    aliases = np.zeros(4, dtype=omt.ALIAS_DTYPE)
    aliases["sample_idx"] = [10, 10, 20, 20]
    aliases["hypothesis_label"] = [b"H01", b"H02", b"H01", b"H02"]
    aliases["tx_beam_snr_weighted_mean_dc"] = [0.1, 0.02, 0.4, 0.03]

    result = selected_closest_tx_alias_sample_indices(events, aliases)

    np.testing.assert_array_equal(result, [10])


def test_measurement_histogram_input_excludes_points_outside_plot_domain():
    import orbit_metadata_table as omt
    from plot_orbit_catalogue_statistics import measurement_height_velocity_arrays

    sample_idx = 1_740_000_000_000_000
    events = np.zeros(1, dtype=omt.EVENT_DTYPE)
    events["sample_idx"] = sample_idx
    events["initial_detection_height_km"] = 100.0
    events["v_g_km_s"] = 40.0
    events["radiant_sun_ecliptic_lon_deg"] = 120.0
    events["orbit_solution_type"] = b"dasst_winning_alias"
    events["combined_score"] = 0.5
    events["frac_e_gt_1"] = 0.0
    events["n_uncertainty_samples"] = 3

    paths = np.zeros(3, dtype=omt.PATH_DTYPE)
    paths["sample_idx"] = sample_idx
    paths["position_enu_km"][:, 2] = [59.0, 100.0, 161.0]
    paths["selection_keep"] = True

    result = measurement_height_velocity_arrays(
        events,
        paths,
        max_radiant_sigma_deg=None,
        min_sample_idx=sample_idx - 1,
        max_sample_idx=sample_idx + 1,
        max_combined_score=1.5,
        max_frac_e_gt_1=0.5,
        min_uncertainty_samples=3,
        max_initial_state_position_sigma_m=1000.0,
        max_initial_state_radiant_angle_sigma_deg=3.0,
    )

    np.testing.assert_array_equal(result["measurement_height_km"], [100.0])


def test_path_arc_length_uses_selected_points_along_best_fit_direction():
    import orbit_metadata_table as omt
    from plot_orbit_catalogue_statistics import path_arc_lengths_km, sample_indices_with_min_arc_length

    paths = np.zeros(5, dtype=omt.PATH_DTYPE)
    paths["sample_idx"] = [10, 10, 10, 20, 20]
    paths["position_enu_km"] = [
        [0.0, 0.0, 100.0],
        [9.0, 0.1, 100.0],
        [30.0, 0.0, 100.0],
        [0.0, 0.0, 100.0],
        [10.0, 0.0, 100.0],
    ]
    paths["selection_keep"] = [True, True, False, True, True]

    sample, arc = path_arc_lengths_km(paths)

    np.testing.assert_array_equal(sample, [10, 20])
    np.testing.assert_allclose(arc, [9.0005, 10.0], rtol=1e-4)
    np.testing.assert_array_equal(sample_indices_with_min_arc_length(paths, 9.5), [20])


def test_measurement_histogram_input_requires_min_arc_length_when_requested():
    import orbit_metadata_table as omt
    from plot_orbit_catalogue_statistics import measurement_height_velocity_arrays, sample_indices_with_min_arc_length

    events = np.zeros(2, dtype=omt.EVENT_DTYPE)
    events["sample_idx"] = [10, 20]
    events["initial_detection_height_km"] = 100.0
    events["v_g_km_s"] = 40.0
    events["radiant_sun_ecliptic_lon_deg"] = 120.0
    events["orbit_solution_type"] = b"dasst_winning_alias"
    events["combined_score"] = 0.5
    events["frac_e_gt_1"] = 0.0
    events["n_uncertainty_samples"] = 3

    paths = np.zeros(4, dtype=omt.PATH_DTYPE)
    paths["sample_idx"] = [10, 10, 20, 20]
    paths["position_enu_km"] = [
        [0.0, 0.0, 100.0],
        [10.0, 0.0, 100.0],
        [0.0, 0.0, 100.0],
        [20.0, 0.0, 100.0],
    ]
    paths["selection_keep"] = True
    long_sample = sample_indices_with_min_arc_length(paths, 15.0)

    result = measurement_height_velocity_arrays(
        events,
        paths,
        max_radiant_sigma_deg=None,
        min_sample_idx=0,
        max_sample_idx=30,
        max_combined_score=1.5,
        max_frac_e_gt_1=0.5,
        min_uncertainty_samples=3,
        max_initial_state_position_sigma_m=1000.0,
        max_initial_state_radiant_angle_sigma_deg=3.0,
        min_arc_sample_idx=long_sample,
    )

    np.testing.assert_array_equal(result["measurement_event_sample_idx"], [20, 20])
    np.testing.assert_array_equal(result["measurement_height_km"], [100.0, 100.0])


def test_measurement_histogram_input_excludes_points_far_from_tx_beam_center():
    import orbit_metadata_table as omt
    import pansy_gain as pgain
    from plot_orbit_catalogue_statistics import measurement_height_velocity_arrays

    sample_idx = 1_740_000_000_000_000
    events = np.zeros(1, dtype=omt.EVENT_DTYPE)
    events["sample_idx"] = sample_idx
    events["initial_detection_height_km"] = 100.0
    events["v_g_km_s"] = 40.0
    events["radiant_sun_ecliptic_lon_deg"] = 120.0
    events["orbit_solution_type"] = b"dasst_winning_alias"
    events["combined_score"] = 0.5
    events["frac_e_gt_1"] = 0.0
    events["n_uncertainty_samples"] = 3

    beam_vec = np.asarray(pgain.tx_beam_unit_vectors()[0], dtype=np.float64)
    beam_vec[2] *= -1.0
    far_vec = beam_vec.copy()
    far_vec[:2] += np.asarray([np.sin(np.deg2rad(20.0)), 0.0])
    far_vec /= np.linalg.norm(far_vec)

    paths = np.zeros(2, dtype=omt.PATH_DTYPE)
    paths["sample_idx"] = sample_idx
    paths["position_enu_km"] = np.vstack([100.0 * beam_vec, 100.0 * far_vec])
    paths["beam_id"] = 0
    paths["selection_keep"] = True

    result = measurement_height_velocity_arrays(
        events,
        paths,
        max_radiant_sigma_deg=None,
        min_sample_idx=sample_idx - 1,
        max_sample_idx=sample_idx + 1,
        max_combined_score=1.5,
        max_frac_e_gt_1=0.5,
        min_uncertainty_samples=3,
        max_initial_state_position_sigma_m=1000.0,
        max_initial_state_radiant_angle_sigma_deg=3.0,
        max_tx_beam_angle_deg=10.0,
    )

    np.testing.assert_allclose(result["measurement_height_km"], [100.0 * beam_vec[2]], rtol=1e-6)


def test_catalogue_arrays_can_require_initial_point_near_tx_beam_center():
    import orbit_metadata_table as omt
    import pansy_gain as pgain
    from plot_orbit_catalogue_statistics import catalogue_arrays, sample_indices_with_max_initial_tx_beam_angle

    events = np.zeros(2, dtype=omt.EVENT_DTYPE)
    events["sample_idx"] = [10, 20]
    events["initial_detection_height_km"] = 100.0
    events["v_g_km_s"] = 40.0
    events["radiant_sun_ecliptic_lon_deg"] = 120.0

    beam_vec = np.asarray(pgain.tx_beam_unit_vectors()[0], dtype=np.float64)
    beam_vec[2] *= -1.0
    far_vec = beam_vec.copy()
    far_vec[:2] += np.asarray([np.sin(np.deg2rad(20.0)), 0.0])
    far_vec /= np.linalg.norm(far_vec)

    paths = np.zeros(4, dtype=omt.PATH_DTYPE)
    paths["sample_idx"] = [10, 10, 20, 20]
    paths["t_rel_s"] = [0.0, 1.0, 0.0, 1.0]
    paths["position_enu_km"] = np.vstack([100.0 * beam_vec, 100.0 * far_vec, 100.0 * far_vec, 100.0 * beam_vec])
    paths["beam_id"] = 0
    paths["selection_keep"] = True

    initial_tx_sample_idx = sample_indices_with_max_initial_tx_beam_angle(paths, 10.0)
    result = catalogue_arrays(
        events,
        max_radiant_sigma_deg=None,
        min_sample_idx=0,
        max_sample_idx=30,
        initial_tx_sample_idx=initial_tx_sample_idx,
    )

    np.testing.assert_array_equal(initial_tx_sample_idx, [10])
    np.testing.assert_array_equal(result["sample_idx"], [10])


def test_low_height_diagnostics_identifies_suspicious_selected_measurements():
    import orbit_metadata_table as omt
    from plot_orbit_catalogue_statistics import measurement_low_height_diagnostics

    sample_idx = 1_740_000_000_000_000
    events = np.zeros(1, dtype=omt.EVENT_DTYPE)
    events["sample_idx"] = sample_idx
    events["initial_detection_height_km"] = 100.0
    events["v_g_km_s"] = 40.0
    events["radiant_sun_ecliptic_lon_deg"] = 120.0
    events["orbit_solution_type"] = b"dasst_winning_alias"
    events["combined_score"] = 0.5
    events["frac_e_gt_1"] = 0.0
    events["n_uncertainty_samples"] = 3

    paths = np.zeros(4, dtype=omt.PATH_DTYPE)
    paths["sample_idx"] = sample_idx
    paths["position_enu_km"][:, 2] = [42.0, 80.0, 100.0, 140.0]
    paths["selection_keep"] = True

    result = measurement_low_height_diagnostics(
        events,
        paths,
        max_radiant_sigma_deg=None,
        min_sample_idx=sample_idx - 1,
        max_sample_idx=sample_idx + 1,
        max_combined_score=1.5,
        max_frac_e_gt_1=0.5,
        min_uncertainty_samples=3,
        max_initial_state_position_sigma_m=1000.0,
        max_initial_state_radiant_angle_sigma_deg=3.0,
        low_height_threshold_km=60.0,
    )

    assert len(result) == 1
    assert result["sample_idx"][0] == sample_idx
    assert result["n_selected_measurements"][0] == 4
    assert result["n_below_threshold"][0] == 1
    assert result["min_height_km"][0] == 42.0
