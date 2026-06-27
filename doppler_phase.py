import numpy as np
import scipy.optimize as so


def circular_phase_residual(observed_phase_rad, model_phase_rad):
    """Circular phase residual, matching the MAARSY zpp fitting convention."""
    return np.angle(np.exp(1j * observed_phase_rad) * np.exp(-1j * model_phase_rad))


def pulse_pair_phase_samples(tx_idx, beam_id, coarse_doppler_mps, pulse_response, wavelength_m, snr=None):
    """Build wrapped same-beam pulse-to-pulse phase samples.

    Each valid output sample corresponds to the phase difference between the
    current pulse and the previous pulse with the same beam ID.  No phase
    phase branch correction is attempted.
    """
    tx_idx = np.asarray(tx_idx, dtype=np.float64)
    beam_id = np.asarray(beam_id, dtype=np.int64)
    coarse_doppler_mps = np.asarray(coarse_doppler_mps, dtype=np.float64)
    pulse_response = np.asarray(pulse_response, dtype=np.complex64)
    if snr is None:
        snr = np.ones(len(tx_idx), dtype=np.float64)
    else:
        snr = np.asarray(snr, dtype=np.float64)

    zpp = np.full(len(tx_idx), np.nan + 1j * np.nan, dtype=np.complex64)
    phase_dt_s = np.full(len(tx_idx), np.nan, dtype=np.float64)
    coarse_phase_rad = np.full(len(tx_idx), np.nan, dtype=np.float64)
    pair_snr = np.full(len(tx_idx), np.nan, dtype=np.float64)

    tx_time_s = tx_idx / 1e6
    for beam in np.unique(beam_id):
        idx = np.where(beam_id == beam)[0]
        if len(idx) < 2:
            continue
        idx = idx[np.argsort(tx_time_s[idx])]
        for prev_i, cur_i in zip(idx[:-1], idx[1:]):
            dt_s = tx_time_s[cur_i] - tx_time_s[prev_i]
            if not np.isfinite(dt_s) or dt_s <= 0.0:
                continue
            cross = np.sum(pulse_response[cur_i] * np.conj(pulse_response[prev_i]))
            if not np.isfinite(cross.real) or not np.isfinite(cross.imag) or np.abs(cross) == 0.0:
                continue

            zpp[cur_i] = cross / np.abs(cross)
            phase_dt_s[cur_i] = dt_s
            coarse_hz = 2.0 * coarse_doppler_mps[cur_i] / wavelength_m
            coarse_phase_rad[cur_i] = np.angle(np.exp(1j * 2.0 * np.pi * coarse_hz * dt_s))
            pair_snr[cur_i] = min(snr[prev_i], snr[cur_i])

    return zpp, phase_dt_s, coarse_phase_rad, pair_snr


def fit_wrapped_phase_time(tx_idx, beam_id, zpp, snr=None):
    """Fit wrapped pulse-to-pulse phase versus time with circular residuals."""
    tx_idx = np.asarray(tx_idx, dtype=np.float64)
    beam_id = np.asarray(beam_id, dtype=np.int64)
    zpp = np.asarray(zpp, dtype=np.complex64)
    if snr is None:
        snr = np.ones(len(tx_idx), dtype=np.float64)
    else:
        snr = np.asarray(snr, dtype=np.float64)

    fit_phase_rad = np.full(len(tx_idx), np.nan, dtype=np.float64)
    fit_residual_rad = np.full(len(tx_idx), np.nan, dtype=np.float64)
    phase_slope_rad_s = np.full(len(tx_idx), np.nan, dtype=np.float64)
    phase_intercept_rad = np.full(len(tx_idx), np.nan, dtype=np.float64)
    tx_time_s = tx_idx / 1e6

    for beam in np.unique(beam_id):
        idx = np.where(beam_id == beam)[0]
        finite = idx[np.isfinite(zpp[idx].real) & np.isfinite(zpp[idx].imag)]
        if len(finite) < 2:
            continue

        t = tx_time_s[finite] - tx_time_s[finite[0]]
        observed = np.angle(zpp[finite])
        weight = np.sqrt(np.maximum(snr[finite], 1e-6))
        intercept0 = np.angle(np.sum(weight * np.exp(1j * observed)))

        def residual(par):
            model = par[0] + par[1] * t
            return weight * circular_phase_residual(observed, model)

        opt = so.least_squares(residual, np.array([intercept0, 0.0]), loss="linear")
        model_phase = opt.x[0] + opt.x[1] * t
        fit_phase_rad[finite] = np.angle(np.exp(1j * model_phase))
        fit_residual_rad[finite] = circular_phase_residual(observed, model_phase)
        phase_intercept_rad[finite] = np.angle(np.exp(1j * opt.x[0]))
        phase_slope_rad_s[finite] = opt.x[1]

    return fit_phase_rad, fit_residual_rad, phase_slope_rad_s, phase_intercept_rad


def phase_fit_doppler_mps(fit_phase_rad, phase_dt_s, wavelength_m):
    """Convert wrapped pulse-pair phase to aliased Doppler velocity."""
    fit_phase_rad = np.asarray(fit_phase_rad, dtype=np.float64)
    phase_dt_s = np.asarray(phase_dt_s, dtype=np.float64)
    return fit_phase_rad / (2.0 * np.pi * phase_dt_s) * wavelength_m / 2.0


def pulse_phase_diagnostics(observed_phase_rad, fit_residual_rad):
    good = np.isfinite(observed_phase_rad) & np.isfinite(fit_residual_rad)
    if np.count_nonzero(good) == 0:
        return {
            "n_phase": 0,
            "phase_residual_median_rad": np.nan,
            "phase_residual_mad_rad": np.nan,
            "phase_residual_rms_rad": np.nan,
        }

    resid = fit_residual_rad[good]
    med = float(np.median(resid))
    return {
        "n_phase": int(np.count_nonzero(good)),
        "phase_residual_median_rad": med,
        "phase_residual_mad_rad": float(np.median(np.abs(resid - med))),
        "phase_residual_rms_rad": float(np.sqrt(np.mean(resid**2))),
    }
