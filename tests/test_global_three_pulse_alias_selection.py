import numpy as np

from fit_full_event_three_pulse_complex_envelope import (
    global_alias_scores,
    select_global_alias,
)


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


def test_acceleration_branch_is_profiled_over_velocity_aliases():
    aliases = {
        # Candidate 0 has the best absolute velocity, but the wrong acceleration
        # branch.  Candidate 1 has the right acceleration branch despite a bad
        # velocity alias.  Candidate 2 is that same acceleration branch with the
        # best available velocity alias, and should be selected second.
        "alias_number": np.asarray([-1, 0, 0]),
        "velocity_alias_number": np.asarray([0, -1, 0]),
        "delta_chi2": np.asarray([0.0, 0.5, 2.0]),
        "fit_velocity_mps": np.asarray([0.0, 400.0, 30.0]),
        "fit_acceleration_mps2": np.asarray([-50_000.0, 0.0, 0.0]),
        "fit_success": np.asarray([True, True, True]),
    }

    selected, acceleration_alias, *_ = select_global_alias(
        aliases,
        model_velocity_mps=0.0,
        model_acceleration_mps2=0.0,
        velocity_sigma_mps=100.0,
        acceleration_sigma_mps2=3_000.0,
    )

    assert acceleration_alias == 0
    assert selected == 2
