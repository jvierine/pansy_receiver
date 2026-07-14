#!/usr/bin/env python3
"""Calculate the PANSY single-pulse minimum detectable radar cross section."""

import numpy as np
from scipy.constants import c, k


def minimum_detectable_rcs_m2(
    *,
    snr_linear: float,
    range_m: float,
    system_temperature_k: float,
    pulse_length_s: float,
    peak_power_w: float,
    transmit_gain_linear: float,
    receive_gain_linear: float,
    frequency_hz: float,
) -> float:
    wavelength_m = c / frequency_hz
    receive_effective_area_m2 = receive_gain_linear * wavelength_m**2 / (4.0 * np.pi)
    noise_bandwidth_hz = 1.0 / pulse_length_s
    return float(
        snr_linear
        * (4.0 * np.pi) ** 2
        * range_m**4
        * k
        * system_temperature_k
        * noise_bandwidth_hz
        / (peak_power_w * transmit_gain_linear * receive_effective_area_m2)
    )


def main() -> None:
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
    print(f"minimum_detectable_rcs_m2 {rcs_m2:.9g}")
    print(f"minimum_detectable_rcs_dbsm {10.0 * np.log10(rcs_m2):.6f}")
    print(f"snr_threshold_db {10.0 * np.log10(7.0):.6f}")


if __name__ == "__main__":
    main()
