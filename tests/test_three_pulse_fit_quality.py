import numpy as np

from three_pulse_fit_quality import QUALITY_SCALES
from three_pulse_fit_quality import event_fit_quality
from three_pulse_fit_quality import variance_weighted_score


def test_reference_standard_deviations_have_score_five():
    assert variance_weighted_score(QUALITY_SCALES) == 5.0


def test_event_fit_quality_uses_requested_units_and_masks():
    metrics = event_fit_quality(
        measured_position_km=np.asarray([[0.25, 0.2, 0.056], [99.0, 99.0, 99.0]]),
        model_position_km=np.zeros((2, 3)),
        position_keep=np.asarray([True, False]),
        triplet_time_s=np.asarray([0.0, 1.0]),
        measured_doppler_mps=np.asarray([129.0, 999.0]),
        measured_acceleration_mps2=np.asarray([3000.0, 999000.0]),
        velocity_keep=np.asarray([True, False]),
        acceleration_keep=np.asarray([True, False]),
        trajectory_time_s=np.asarray([0.0, 1.0]),
        model_doppler_mps=np.zeros(2),
        model_acceleration_mps2=np.zeros(2),
    )
    assert metrics["ew_rms_km"] == 0.25
    assert metrics["ns_rms_km"] == 0.2
    assert metrics["up_rms_km"] == 0.056
    assert metrics["doppler_rms_mps"] == 129.0
    assert metrics["acceleration_rms_km_s2"] == 3.0
    assert metrics["variance_weighted_score"] == 5.0
