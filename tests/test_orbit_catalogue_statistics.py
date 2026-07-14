import sys
from pathlib import Path

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


def test_height_velocity_histogram_spans_20_to_180_km():
    from plot_orbit_catalogue_statistics import histogram_height_velocity

    counts, height_edges, speed_edges = histogram_height_velocity(
        np.asarray([19.9, 20.1, 179.9, 180.1]),
        np.asarray([20.0, 20.0, 20.0, 20.0]),
    )

    assert height_edges[0] == 20.0
    assert height_edges[-1] == 180.0
    assert speed_edges[0] == 0.0
    assert speed_edges[-1] == 80.0
    assert counts.sum() == 2


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
    paths["position_enu_km"][:, 2] = [19.0, 100.0, 181.0]
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
