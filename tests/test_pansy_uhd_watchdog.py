import importlib.util
from pathlib import Path

import numpy as np


WATCHDOG_PATH = Path(__file__).resolve().parents[1] / "receiver" / "scripts" / "pansy_uhd_watchdog.py"
SPEC = importlib.util.spec_from_file_location("pansy_uhd_watchdog", WATCHDOG_PATH)
watchdog = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(watchdog)


def test_phase_vector_accepts_small_channel_offsets():
    reference = np.deg2rad([0.0, -5.4, -6.3, -25.7, -7.4, -4.4, -1.4, -174.0])
    measured = reference + np.deg2rad([0.0, -4.8, 1.7, -1.6, -3.4, -3.5, 6.4, 0.8])
    status = watchdog.classify_phase_vector(
        measured,
        np.ones(8),
        reference,
        threshold_deg=30.0,
        min_valid_channels=7,
    )
    assert status["ok"] is True
    np.testing.assert_allclose(status["max_abs_deg"], 6.4, atol=1e-10)


def test_phase_vector_rejects_half_cycle_usrp_offset():
    reference = np.deg2rad([0.0, -5.4, -6.3, -25.7, -7.4, -4.4, -1.4, -174.0])
    measured = reference.copy()
    measured[2] += np.deg2rad(-169.0)
    measured[3] += np.deg2rad(169.0)
    status = watchdog.classify_phase_vector(
        measured,
        np.ones(8),
        reference,
        threshold_deg=30.0,
        min_valid_channels=7,
    )
    assert status["ok"] is False
    assert status["max_abs_deg"] > 160.0
