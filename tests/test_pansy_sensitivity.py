import sys
from pathlib import Path

import numpy as np


PANSY_RECEIVER_ROOT = Path(__file__).resolve().parents[1]
if str(PANSY_RECEIVER_ROOT) not in sys.path:
    sys.path.insert(0, str(PANSY_RECEIVER_ROOT))


def test_single_pulse_minimum_detectable_rcs():
    from pansy_sensitivity import minimum_detectable_rcs_m2

    transmit_gain = 10.0 ** (37.0 / 10.0)
    rcs_m2 = minimum_detectable_rcs_m2(
        snr_linear=7.0,
        range_m=100e3,
        system_temperature_k=4550.0,
        pulse_length_s=128e-6,
        peak_power_w=500e3,
        transmit_gain_linear=transmit_gain,
        receive_gain_linear=transmit_gain / 54.0,
        frequency_hz=47e6,
    )

    np.testing.assert_allclose(rcs_m2, 7.204277073777118e-5, rtol=1e-12)
    np.testing.assert_allclose(10.0 * np.log10(rcs_m2), -41.42409592729447, rtol=1e-12)
