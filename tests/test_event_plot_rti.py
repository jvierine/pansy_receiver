import numpy as np

from fit_full_event_three_pulse_complex_envelope import uniform_ipp_rti


def test_uniform_ipp_rti_leaves_missing_pulses_as_nan_columns():
    ipp_samples = 8000
    tx_idx = np.array([1_000_000, 1_008_000, 1_024_000], dtype=np.int64)
    rti = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])

    uniform_tx_idx, uniform_rti = uniform_ipp_rti(tx_idx, rti, ipp_samples)

    np.testing.assert_array_equal(
        uniform_tx_idx,
        np.array([1_000_000, 1_008_000, 1_016_000, 1_024_000]),
    )
    np.testing.assert_array_equal(uniform_rti[[0, 1, 3]], rti)
    assert np.isnan(uniform_rti[2]).all()

