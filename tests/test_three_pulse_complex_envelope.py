import numpy as np

from aod_complex_doppler_fit import event_module_voltage_gain
from fit_three_pulse_complex_envelope import (
    baud_measurements,
    complex_envelope_model,
    fit_three_pulse_envelope,
)


def test_event_module_voltage_gain_removes_common_meteor_amplitude():
    pulse_amplitude = np.geomspace(0.2, 8.0, 31)
    module_gain = np.asarray([0.7, 0.9, 1.0, 1.1, 1.3, 1.7, 2.2])
    phase = np.exp(1j * np.linspace(-1.0, 1.0, len(module_gain)))
    response = pulse_amplitude[:, None] * module_gain[None, :] * phase[None, :]

    estimated = event_module_voltage_gain(response)
    expected = module_gain / np.exp(np.mean(np.log(module_gain)))

    np.testing.assert_allclose(estimated, expected, rtol=1e-12, atol=1e-12)


def test_three_pulse_fit_recovers_frequency_and_acceleration_from_good_initial_branch():
    rng = np.random.default_rng(42)
    wavelength_m = 6.378562936170213
    n_fast = 132
    spacing_s = 0.008
    local_time = np.arange(n_fast) / 1e6
    absolute_time = np.concatenate(
        (local_time, spacing_s + local_time, 2.0 * spacing_s + local_time)
    )
    template = np.zeros(n_fast, dtype=complex)
    for baud, start in enumerate(range(2, 122, 10)):
        shape = np.sin(np.linspace(0.15, np.pi - 0.15, 7))
        template[start : start + 7] = shape * np.exp(1j * rng.uniform(-np.pi, np.pi))
    templates = np.stack((template, template, template))
    dummy = baud_measurements(templates, templates, np.zeros(3), spacing_s, np.ones(3))
    center = np.sum(dummy["weight"] * dummy["absolute_time_s"]) / np.sum(dummy["weight"])
    truth = np.asarray([np.log(0.8), np.log(1.25), np.log(0.95), 0.7, -6800.0, 8500.0])
    frequency_rate = 2.0 * truth[5] / wavelength_m
    phase = truth[3] + 2.0 * np.pi * (
        truth[4] * (absolute_time - center) + 0.5 * frequency_rate * (absolute_time - center) ** 2
    )
    raw = np.empty((3, n_fast), dtype=complex)
    for pulse in (0, 1, 2):
        samples = slice(pulse * n_fast, (pulse + 1) * n_fast)
        raw[pulse] = np.exp(truth[pulse]) * template * np.exp(1j * phase[samples])
    noise = 0.001 * (rng.normal(size=raw.shape) + 1j * rng.normal(size=raw.shape))
    raw += noise

    fit = fit_three_pulse_envelope(
        raw,
        templates,
        spacing_s,
        frequency_guess_hz=-6798.0,
        acceleration_guess_mps2=8400.0,
        wavelength_m=wavelength_m,
        acceleration_half_width_mps2=5000.0,
        frequency_half_width_hz=1000.0,
        pulse_snr=np.asarray([20.0, 30.0, 25.0]),
        matched_filter_amplitudes=np.asarray([0.8, 1.25, 0.95]),
    )

    assert fit["success"]
    assert abs(fit["parameters"][4] - truth[4]) < 5.0
    assert abs(fit["parameters"][5] - truth[5]) < 500.0
    assert abs(np.exp(fit["parameters"][1] - fit["parameters"][0]) - 1.25 / 0.8) < 0.12
    propagated = (
        fit["velocity_acceleration_noise_influence"]
        @ fit["velocity_acceleration_noise_influence"].T
    )
    np.testing.assert_allclose(
        propagated,
        fit["velocity_acceleration_covariance"],
        rtol=1e-10,
        atol=1e-10,
    )
    assert fit["residual_pulse_index"].shape == (
        fit["velocity_acceleration_noise_influence"].shape[1],
    )


def test_independent_pulse_amplitudes_do_not_introduce_independent_phase():
    source = complex_envelope_model(
        np.asarray([0.0, np.log(2.0), np.log(3.0), 0.3, 1000.0, 500.0]),
        np.asarray([-0.008, 0.0, 0.008]),
        np.ones(3, dtype=float),
        np.asarray([0, 1, 2]),
        wavelength_m=6.0,
    )
    assert np.isclose(abs(source[1]) / abs(source[0]), 2.0)
