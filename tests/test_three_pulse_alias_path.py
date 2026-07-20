import numpy as np

from fit_full_event_three_pulse_complex_envelope import (
    select_global_doppler_alias_path,
)


def test_global_alias_path_removes_independent_sideband_switches():
    model = np.zeros(6)
    base = np.asarray([200.0, -200.0, 200.0, -200.0, 200.0, -200.0])
    aliases = np.asarray([-1, 0, 1])
    rows = [
        {
            "fit_velocity_mps": value + 400.0 * aliases,
            "fit_velocity_std_mps": np.full(3, 5.0),
            "delta_chi2": np.zeros(3),
        }
        for value in base
    ]

    path = select_global_doppler_alias_path(
        rows,
        np.arange(len(rows), dtype=float) * 0.008,
        model,
        transition_sigma_mps=20.0,
        anchor_sigma_mps=400.0,
    )
    selected = np.asarray(
        [rows[index]["fit_velocity_mps"][state] for index, state in enumerate(path["selected_index"])]
    )

    np.testing.assert_allclose(np.diff(selected), 0.0)
    assert np.all(np.abs(selected) == 200.0)
