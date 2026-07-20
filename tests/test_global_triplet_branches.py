import numpy as np

from fit_full_event_three_pulse_complex_envelope import (
    smooth_global_branch_indices,
)


def test_smooth_branch_selection_allows_an_alias_boundary_crossing():
    candidates = [
        {
            "velocity_mps": np.zeros(3),
            "acceleration_mps2": np.asarray([-50e3, 0.0, 50e3]),
        }
        for _ in range(5)
    ]
    selected = smooth_global_branch_indices(
        candidates,
        predicted_velocity_mps=np.zeros(5),
        predicted_acceleration_mps2=np.asarray([5e3, 15e3, 25e3, 35e3, 45e3]),
        velocity_sigma_mps=1e3,
        acceleration_sigma_mps2=10e3,
    )

    np.testing.assert_array_equal(selected, np.asarray([1, 1, 1, 2, 2]))

