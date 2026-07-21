import numpy as np

from fit_full_event_three_pulse_complex_envelope import (
    fft_locked_acceleration_limits,
    select_local_alias,
)


def test_fft_lock_spans_exactly_one_acceleration_alias():
    spacing = 49_800.0
    lower, upper = fft_locked_acceleration_limits(8_000.0, spacing)

    assert lower > 8_000.0 - 0.5 * spacing
    assert upper < 8_000.0 + 0.5 * spacing
    np.testing.assert_allclose((lower + upper) / 2.0, 8_000.0)


def test_local_chi2_selects_velocity_within_locked_acceleration_branch():
    aliases = {
        "alias_number": np.asarray([0, 0, 0]),
        "weighted_sse": np.asarray([12.0, 2.0, 9.0]),
        "fit_success": np.asarray([True, True, True]),
    }

    assert select_local_alias(aliases) == 1


def test_local_selection_rejects_multiple_acceleration_branches():
    aliases = {
        "alias_number": np.asarray([-1, 0]),
        "weighted_sse": np.asarray([1.0, 2.0]),
        "fit_success": np.asarray([True, True]),
    }

    with np.testing.assert_raises(RuntimeError):
        select_local_alias(aliases)
