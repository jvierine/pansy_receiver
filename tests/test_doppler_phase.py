import numpy as np
import pytest

from doppler_phase import (
    circular_phase_residual,
    fit_wrapped_phase_time,
    phase_fit_doppler_mps,
    pulse_pair_phase_samples,
    pulse_phase_diagnostics,
)


def test_pulse_pair_phase_samples_uses_same_beam_five_ipp_pairs():
    wavelength_m = 4.0
    tx_idx = np.arange(10, dtype=np.float64) * 1_600.0
    beam_id = np.arange(len(tx_idx), dtype=np.int64) % 5
    phase = np.arange(len(tx_idx), dtype=np.float64) * 0.4
    response = np.exp(1j * phase)[:, None]
    coarse = np.zeros(len(tx_idx), dtype=np.float64)

    zpp, phase_dt_s, _coarse_phase, pair_snr = pulse_pair_phase_samples(
        tx_idx,
        beam_id,
        coarse,
        response,
        wavelength_m,
        snr=np.arange(1, len(tx_idx) + 1, dtype=np.float64),
    )

    assert np.all(np.isnan(zpp[:5].real))
    assert np.allclose(phase_dt_s[5:], 0.008)
    assert np.angle(zpp[5:]) == pytest.approx(np.full(5, 2.0))
    assert pair_snr[5:].tolist() == [1.0, 2.0, 3.0, 4.0, 5.0]


def test_fit_wrapped_phase_time_fits_circular_phase_without_branch_correction():
    tx_idx = np.arange(12, dtype=np.float64) * 8_000.0
    beam_id = np.zeros(len(tx_idx), dtype=np.int64)
    true_phase = 2.8 - 35.0 * (tx_idx / 1e6 - tx_idx[0] / 1e6)
    zpp = np.exp(1j * true_phase).astype(np.complex64)
    zpp[0] = np.nan + 1j * np.nan

    fit_phase, residual, slope, intercept = fit_wrapped_phase_time(tx_idx, beam_id, zpp)

    good = np.isfinite(residual)
    assert np.sqrt(np.mean(residual[good] ** 2)) < 1e-5
    assert slope[good] == pytest.approx(-35.0, abs=1e-4)
    assert np.isfinite(intercept[good]).all()
    assert circular_phase_residual(np.angle(zpp[good]), fit_phase[good]) == pytest.approx(residual[good])


def test_phase_fit_doppler_and_diagnostics():
    fit_phase = np.array([np.nan, np.pi / 2.0])
    dt = np.array([np.nan, 0.008])
    doppler = phase_fit_doppler_mps(fit_phase, dt, wavelength_m=4.0)
    assert np.isnan(doppler[0])
    assert doppler[1] == pytest.approx(62.5)

    diag = pulse_phase_diagnostics(np.array([0.0, 1.0]), np.array([0.1, -0.2]))
    assert diag["n_phase"] == 2
    assert diag["phase_residual_rms_rad"] == pytest.approx(np.sqrt(0.025))
