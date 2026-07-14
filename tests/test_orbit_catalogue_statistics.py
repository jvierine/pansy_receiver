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
