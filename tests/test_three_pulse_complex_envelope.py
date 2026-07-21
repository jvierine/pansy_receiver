import numpy as np

from aod_complex_doppler_fit import event_module_voltage_gain
from fit_three_pulse_complex_envelope import (
    FS_HZ,
    common_bin_pulse_pair_doppler,
    baud_measurements,
    complex_envelope_model,
    fit_three_pulse_acceleration_aliases,
    fit_three_pulse_envelope,
    gapped_two_pulse_fft_doppler,
    independent_pulse_fft_doppler,
)


def test_common_bin_pulse_pair_recovers_frequency_alias():
    sample_count = 80
    template = np.zeros(sample_count, dtype=np.complex128)
    for start in (5, 15, 25, 35, 45, 55):
        template[start : start + 3] = 1.0
    pulse_spacing_s = 0.008
    frequency_hz = 12345.0
    fast_time_s = np.arange(sample_count) / FS_HZ
    raw_pulses = np.stack(
        [
            template
            * np.exp(
                1j
                * 2.0
                * np.pi
                * frequency_hz
                * (pulse * pulse_spacing_s + fast_time_s)
            )
            for pulse in range(2)
        ]
    )
    result = common_bin_pulse_pair_doppler(
        raw_pulses,
        np.stack((template, template)),
        pulse_spacing_s,
        frequency_hz + 20.0,
        0.6,
        zero_pad_factor=16,
    )
    assert result["nfft"] >= 16 * np.max(result["baud_count"])
    assert result["frequency_step_hz"] < 1e4
    assert abs(result["resolved_frequency_hz"] - frequency_hz) < 1e-6
    assert 0.0 < result["coherence"] <= 1.0


def test_gapped_two_pulse_fft_uses_native_separation():
    sample_count = 80
    template = np.zeros(sample_count, dtype=np.complex128)
    for start in (5, 15, 25, 35, 45, 55):
        template[start : start + 3] = 1.0
    pulse_spacing_s = 0.008
    frequency_hz = 12345.0
    fast_time_s = np.arange(sample_count) / FS_HZ
    raw_pulses = np.stack(
        [
            template
            * np.exp(
                1j
                * 2.0
                * np.pi
                * frequency_hz
                * (pulse * pulse_spacing_s + fast_time_s)
            )
            for pulse in range(2)
        ]
    )
    result = gapped_two_pulse_fft_doppler(
        raw_pulses,
        np.stack((template, template)),
        pulse_spacing_s,
        frequency_hz + 20.0,
        0.6,
        zero_pad_factor=16,
    )
    assert result["timeline_length"] > 0.9 * pulse_spacing_s * result["baud_sample_rate_hz"]
    assert result["occupied_samples"] == 12
    assert result["physical_aperture_s"] > pulse_spacing_s
    assert abs(result["resolved_frequency_hz"] - frequency_hz) < 10.0


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


def synthetic_three_pulse_voltage(seed=73):
    rng = np.random.default_rng(seed)
    wavelength_m = 6.378562936170213
    n_fast = 132
    spacing_s = 0.008
    local_time = np.arange(n_fast) / FS_HZ
    absolute_time = np.concatenate(
        (local_time, spacing_s + local_time, 2.0 * spacing_s + local_time)
    )
    template = np.zeros(n_fast, dtype=complex)
    for start in range(2, 122, 10):
        shape = np.sin(np.linspace(0.15, np.pi - 0.15, 7))
        template[start : start + 7] = shape * np.exp(1j * rng.uniform(-np.pi, np.pi))
    templates = np.stack((template, template, template))
    dummy = baud_measurements(templates, templates, np.zeros(3), spacing_s, np.ones(3))
    center = np.sum(dummy["weight"] * dummy["absolute_time_s"]) / np.sum(dummy["weight"])
    truth = np.asarray([np.log(0.8), np.log(1.25), np.log(0.95), 0.7, -6800.0, 8500.0])
    frequency_rate = 2.0 * truth[5] / wavelength_m
    phase = truth[3] + 2.0 * np.pi * (
        truth[4] * (absolute_time - center)
        + 0.5 * frequency_rate * (absolute_time - center) ** 2
    )
    raw = np.empty((3, n_fast), dtype=complex)
    for pulse in (0, 1, 2):
        samples = slice(pulse * n_fast, (pulse + 1) * n_fast)
        raw[pulse] = np.exp(truth[pulse]) * template * np.exp(1j * phase[samples])
    raw += 2e-4 * (rng.normal(size=raw.shape) + 1j * rng.normal(size=raw.shape))
    return raw, templates, spacing_s, wavelength_m, truth


def test_independent_pulse_fft_mean_seeds_three_pulse_frequency():
    raw, templates, spacing_s, wavelength_m, truth = synthetic_three_pulse_voltage()
    result = independent_pulse_fft_doppler(
        raw,
        templates,
        spacing_s,
        wavelength_m,
        pulse_snr=np.asarray([20.0, 30.0, 25.0]),
        zero_pad_factor=16,
    )

    assert result["pulse_frequency_hz"].shape == (3,)
    assert np.all(result["nfft"] >= 16 * 12)
    assert abs(result["mean_frequency_hz"] - truth[4]) < 300.0


def test_joint_velocity_acceleration_alias_search_recovers_wrong_seed_branches():
    raw, templates, spacing_s, wavelength_m, truth = synthetic_three_pulse_voltage()
    velocity_alias_hz = 1.0 / spacing_s
    acceleration_alias_mps2 = wavelength_m / (2.0 * spacing_s**2)
    result = fit_three_pulse_acceleration_aliases(
        raw,
        templates,
        spacing_s,
        frequency_guess_hz=truth[4] + velocity_alias_hz,
        acceleration_phase_guess_mps2=truth[5] + acceleration_alias_mps2,
        wavelength_m=wavelength_m,
        acceleration_limits_mps2=(-1.0e5, 1.0e5),
        pulse_snr=np.asarray([20.0, 30.0, 25.0]),
        matched_filter_amplitudes=np.exp(truth[:3]),
        velocity_alias_half_count=1,
        maximum_reseed_passes=3,
    )

    fit = result["selected_fit"]
    assert fit["success"]
    assert result["velocity_alias_number"][result["selected_index"]] == -1
    assert abs(fit["parameters"][4] - truth[4]) < 5.0
    assert abs(fit["parameters"][5] - truth[5]) < 500.0
    assert result["reseed_passes"] >= 2
    candidate_count = len(result["weighted_sse"])
    for name in (
        "alias_number",
        "velocity_alias_number",
        "fit_velocity_mps",
        "fit_acceleration_mps2",
        "fit_velocity_variance_mps2",
        "fit_velocity_acceleration_covariance_mps3_s2",
        "fit_acceleration_variance_mps4_s4",
        "fit_success",
        "degrees_of_freedom",
        "delta_chi2",
    ):
        assert len(result[name]) == candidate_count
    assert np.all(np.isfinite(result["weighted_sse"]))
