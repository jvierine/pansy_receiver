import numpy as np

from fit_full_event_three_pulse_complex_envelope import global_alias_scores


def test_global_alias_score_balances_local_and_trajectory_chi2():
    aliases = {
        "delta_chi2": np.asarray([0.0, 2.0, 12.0]),
        "fit_velocity_mps": np.asarray([400.0, 0.0, 0.0]),
        "fit_acceleration_mps2": np.asarray([0.0, 0.0, 9_000.0]),
        "fit_success": np.asarray([True, True, True]),
    }

    score, local, trajectory = global_alias_scores(
        aliases,
        model_velocity_mps=0.0,
        model_acceleration_mps2=0.0,
        velocity_sigma_mps=100.0,
        acceleration_sigma_mps2=3_000.0,
    )

    # The local best is four Doppler sigmas off the global trajectory, so the
    # mildly worse local alias is the correct joint choice.
    assert np.argmin(local) == 0
    assert np.argmin(score) == 1
    np.testing.assert_allclose(trajectory, [16.0, 0.0, 9.0])

