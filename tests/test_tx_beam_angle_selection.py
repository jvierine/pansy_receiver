import numpy as np


def _track(hypothesis_id, angle_deg, distance_dc, reduced_chi2):
    return {
        "hypothesis_id": hypothesis_id,
        "reason": "kept",
        "unique_pulses": 12,
        "measurement_range_span_km": 3.0,
        "fixed_velocity_reduced_chi2": reduced_chi2,
        "fixed_velocity_params": np.asarray(
            [0.0, 0.0, 100e3, 10e3, 0.0, -10e3], dtype=np.float64
        ),
        "tx_beam_weighted_mean_deg": angle_deg,
        "tx_beam_snr_weighted_mean_dc": distance_dc,
        "line_reduced_chi2": 0.5,
        "line_length_adjusted_reduced_chi2": 0.5,
        "ballistic_pulse_coverage_reject": False,
        "common_pulse_completion_reject": False,
        "linearity_reject": False,
        "descent_reject": False,
        "low_detection_altitude_reject": False,
        "ballistic_low_start_altitude_reject": False,
    }


def test_far_provisional_winner_is_replaced_by_tx_nearest_hypothesis():
    import plot_interferometric_disambiguation as disamb

    far_good_fit = _track(1, angle_deg=30.0, distance_dc=0.50, reduced_chi2=0.5)
    nearest_bad_fit = _track(2, angle_deg=3.0, distance_dc=0.05, reduced_chi2=5.0)

    disamb.score_combined_hypotheses([far_good_fit, nearest_bad_fit])

    assert nearest_bad_fit["combined_rank"] == 0
    assert nearest_bad_fit["combined_tx_angle_override_selected"]
    assert nearest_bad_fit["combined_tx_angle_override_triggered"]
    assert nearest_bad_fit["combined_reject"]
    assert far_good_fit["combined_rank"] == 1
    assert far_good_fit["tx_beam_angle_reject"]


def test_nearest_hypothesis_over_limit_is_ranked_first_but_rejected():
    import plot_interferometric_disambiguation as disamb

    far_good_fit = _track(1, angle_deg=35.0, distance_dc=0.60, reduced_chi2=0.5)
    nearest_bad_fit = _track(2, angle_deg=20.0, distance_dc=0.30, reduced_chi2=5.0)

    disamb.score_combined_hypotheses([far_good_fit, nearest_bad_fit])

    assert nearest_bad_fit["combined_rank"] == 0
    assert nearest_bad_fit["combined_tx_angle_override_selected"]
    assert nearest_bad_fit["tx_beam_angle_reject"]
    assert nearest_bad_fit["combined_reject"]


def test_valid_tx_nearest_hypothesis_remains_selected_without_override():
    import plot_interferometric_disambiguation as disamb

    farther = _track(1, angle_deg=12.0, distance_dc=0.20, reduced_chi2=0.5)
    nearest = _track(2, angle_deg=3.0, distance_dc=0.05, reduced_chi2=0.7)

    disamb.score_combined_hypotheses([farther, nearest])

    assert nearest["combined_rank"] == 0
    assert not nearest["combined_tx_angle_override_triggered"]
    assert not nearest["combined_reject"]
