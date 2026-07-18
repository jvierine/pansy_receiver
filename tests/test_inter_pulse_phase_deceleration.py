import numpy as np

from inter_pulse_phase_deceleration import decoded_waveform_pairs, pair_responses, robust_line


def test_pair_phase_recovers_velocity_branch_and_acceleration():
    wavelength = 6.378562936170213
    n = 80
    dt = 0.0064
    time_s = np.arange(n) * dt
    velocity = -22000.0 + 1200.0 * time_s
    phase = np.zeros(n)
    phase[1:] = np.cumsum(4.0 * np.pi * velocity[1:] * dt / wavelength)
    response = np.exp(1j * phase)[:, None] * np.ones((1, 7))
    decoded = {
        "tx_idx": time_s * 1e6,
        "beam_id": np.zeros(n, dtype=int),
        "response": response,
        "coarse_doppler_mps": velocity + 40.0,
    }
    pairs = pair_responses(decoded, same_beam=True)
    fit = robust_line(pairs["time_s"], pairs["phase_doppler_mps"], np.ones(n - 1))
    assert np.max(np.abs(pairs["phase_doppler_mps"] - velocity[1:])) < 1e-6
    assert abs(fit["acceleration_mps2"] - 1200.0) < 1e-6


def test_waveform_product_cancels_common_doppler_and_recovers_delta_frequency():
    n_fast = 132
    fast_time = np.arange(n_fast) / 1e6
    delta_frequency_hz = 73.0
    common_frequency_hz = -6800.0
    previous = np.exp(1j * 2.0 * np.pi * common_frequency_hz * fast_time)
    current = np.exp(1j * 2.0 * np.pi * (common_frequency_hz + delta_frequency_hz) * fast_time)
    wave = np.stack([previous, current])[:, None, :] * np.ones((1, 7, 1))
    decoded = {
        "tx_idx": np.asarray([0.0, 6400.0]),
        "beam_id": np.asarray([0, 0]),
        "decoded": wave,
        "fast_time_s": fast_time,
    }
    pairs = decoded_waveform_pairs(decoded, same_beam=True)
    assert len(pairs["delta_frequency_hz"]) == 1
    assert abs(pairs["delta_frequency_hz"][0] - delta_frequency_hz) < 0.2
