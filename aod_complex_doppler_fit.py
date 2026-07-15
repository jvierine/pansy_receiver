#!/usr/bin/env python3
"""Fit Doppler from an AoD-coherently summed decoded meteor echo."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

import estimate_meteor_rcs as emr
from interferometer_alias_diagnostics import RangeDopplerSearch, amp_scale, load_cut
import pansy_config as pc
import pansy_interferometry as pint


def beamform_echo(echo_ch: np.ndarray, uvw: np.ndarray, beam: int) -> tuple[np.ndarray, tuple[int, int]]:
    """Try calibration/geometry sign conventions and keep the strongest one."""
    antpos = pint.get_antpos()
    phasecal = pint.get_phasecal()[beam]
    geom_phase = (2.0 * np.pi / pc.wavelength) * (antpos @ uvw)
    best = None
    for s_cal in (-1, 1):
        for s_geom in (-1, 1):
            weights = np.exp(1j * (s_cal * phasecal + s_geom * geom_phase))
            summed = np.sum(echo_ch * weights[:, None], axis=0) / len(weights)
            amp = float(np.nanmax(np.abs(summed)))
            if best is None or amp > best[0]:
                best = (amp, summed, (s_cal, s_geom))
    assert best is not None
    return best[1], best[2]


def fit_weighted_complex_sinusoid(decoded: np.ndarray, tx: np.ndarray, coarse_doppler_mps: float):
    """Fit decoded ~= C |tx|^2 exp(i 2 pi f t), weighted by TX envelope."""
    ns = np.arange(len(tx), dtype=np.float64)
    t = ns / 1e6
    env2 = np.abs(tx) ** 2
    env2 = env2 / max(float(np.nanmax(env2)), 1e-30)
    weights = env2.copy()
    weights[env2 <= 0.05] = 0.0

    coarse_hz = coarse_doppler_mps * 2.0 * pc.freq / sc.c
    freqs = np.linspace(coarse_hz - 12000.0, coarse_hz + 12000.0, 2401)
    chi = np.empty(len(freqs), dtype=np.float64)
    coef = np.empty(len(freqs), dtype=np.complex128)
    for i, freq in enumerate(freqs):
        model_basis = env2 * np.exp(1j * 2.0 * np.pi * freq * t)
        den = np.sum(weights * np.abs(model_basis) ** 2)
        c = np.sum(weights * np.conj(model_basis) * decoded) / max(float(den), 1e-30)
        resid = decoded - c * model_basis
        chi[i] = np.sum(weights * np.abs(resid) ** 2) / max(float(np.sum(weights)), 1e-30)
        coef[i] = c
    best = int(np.argmin(chi))
    freq = float(freqs[best])
    doppler = freq * sc.c / (2.0 * pc.freq)
    model = coef[best] * env2 * np.exp(1j * 2.0 * np.pi * freq * t)
    return {
        "sample": ns,
        "time_s": t,
        "env2": env2,
        "weights": weights,
        "freq_hz": freq,
        "doppler_mps": float(doppler),
        "coef": coef[best],
        "model": model,
        "resid": decoded - model,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-idx", type=int, required=True)
    parser.add_argument("--cut-dir", type=Path, default=Path("/mnt/data/juha/pansy/metadata/cut"))
    parser.add_argument("--orbit-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cut = load_cut(args.cut_dir, args.sample_idx)
    z_rx = np.asarray(cut["zrx_echoes_re"], dtype=np.complex64) + 1j * np.asarray(cut["zrx_echoes_im"], dtype=np.complex64)
    z_tx = np.asarray(cut["ztx_pulses_re"], dtype=np.complex64) + 1j * np.asarray(cut["ztx_pulses_im"], dtype=np.complex64)
    z_rx = z_rx * amp_scale()[None, :, None]
    tx_idx = np.asarray(cut["tx_idx"], dtype=np.float64)
    beam_id = np.asarray(cut["beam_id"], dtype=np.int64)

    events = [event for event in emr.read_orbit_metadata_file(args.orbit_file) if int(event["sample_idx"]) == args.sample_idx]
    if not events:
        raise RuntimeError(f"sample {args.sample_idx} not found in {args.orbit_file}")
    event = events[0]
    path_t = np.asarray(event["t_rel_s"], dtype=np.float64)
    path_uvw = np.asarray(event["uvw"], dtype=np.float64)
    t_rel = (tx_idx - args.sample_idx) / 1e6
    uvw = np.full((len(tx_idx), 3), np.nan)
    for i, tt in enumerate(t_rel):
        uvw[i] = path_uvw[int(np.nanargmin(np.abs(path_t - tt)))]

    rds = RangeDopplerSearch(txlen=z_tx.shape[1], echolen=z_rx.shape[2], interp=1, n_channels=z_rx.shape[1])
    coarse = []
    for pulse in range(z_rx.shape[0]):
        mf, pprof, peak_dopv, noise_floor, _xc, rgmax = rds.mf(z_tx[pulse], z_rx[pulse])
        doppler = float(rds.dopv[int(np.argmax(mf[rgmax]))])
        snr = (float(pprof[rgmax]) - float(noise_floor)) / max(float(noise_floor), 1e-30)
        coarse.append((snr, int(rgmax), doppler))
    coarse = np.asarray(coarse, dtype=[("snr", "f8"), ("range_bin", "i4"), ("doppler_mps", "f8")])

    valid = np.all(np.isfinite(uvw), axis=1) & (beam_id >= 0) & (beam_id < 5)
    pulse = int(next(i for i in np.argsort(coarse["snr"])[::-1] if valid[i]))
    beam = int(beam_id[pulse])
    range_bin = int(coarse["range_bin"][pulse])
    summed, signs = beamform_echo(z_rx[pulse].astype(np.complex128), uvw[pulse], beam)
    tx = z_tx[pulse].astype(np.complex128)
    decoded = summed[range_bin : range_bin + len(tx)] * np.conj(tx)
    fit = fit_weighted_complex_sinusoid(decoded, tx, float(coarse["doppler_mps"][pulse]))

    sample = fit["sample"]
    model = fit["model"]
    env_mask = fit["env2"] > 0.05
    data_norm = decoded[env_mask] / np.maximum(fit["env2"][env_mask], 1e-6)
    model_norm = fit["coef"] * np.exp(1j * 2.0 * np.pi * fit["freq_hz"] * fit["time_s"][env_mask])
    phase_data = np.unwrap(np.angle(data_norm))
    phase_model = np.unwrap(np.angle(model_norm))
    phase_data -= np.nanmedian(phase_data - phase_model)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(4, 1, figsize=(8.2, 9.0), sharex=True, constrained_layout=True)
    axes[0].plot(sample, np.abs(tx) / max(float(np.nanmax(np.abs(tx))), 1e-30), color="0.2", label="|TX pulse|")
    axes[0].plot(sample, fit["env2"], color="tab:blue", label="fit weight envelope")
    axes[0].set_ylabel("normalized")
    axes[0].legend(frameon=False, loc="upper right")
    axes[0].set_title(
        f"AoD coherent sum + weighted complex sinusoid, sample {args.sample_idx}, "
        f"pulse {pulse}, beam {beam}, SNR {10*np.log10(max(coarse['snr'][pulse], 1e-30)):.1f} dB"
    )

    axes[1].plot(sample, decoded.real, ".", ms=3, color="tab:blue", alpha=0.75, label="decoded real")
    axes[1].plot(sample, model.real, "-", lw=1.4, color="navy", label="model real")
    axes[1].plot(sample, decoded.imag, ".", ms=3, color="tab:orange", alpha=0.75, label="decoded imag")
    axes[1].plot(sample, model.imag, "-", lw=1.4, color="darkorange", label="model imag")
    axes[1].set_ylabel("decoded amplitude")
    axes[1].legend(frameon=False, ncol=2, loc="upper right")

    axes[2].plot(sample[env_mask], phase_data, ".", ms=4, color="tab:green", label="decoded / envelope")
    axes[2].plot(sample[env_mask], phase_model, "-", lw=1.4, color="black", label="sinusoid phase")
    axes[2].set_ylabel("unwrapped phase (rad)")
    axes[2].legend(frameon=False, loc="upper left")

    axes[3].plot(sample, np.abs(fit["resid"]), ".", ms=3, color="0.3", label="|residual|")
    axes[3].plot(sample, np.abs(decoded), "-", lw=1.0, color="0.75", label="|decoded|")
    axes[3].set_ylabel("amplitude")
    axes[3].set_xlabel("sample within decoded TX pulse")
    axes[3].legend(frameon=False, loc="upper right")
    for ax in axes:
        ax.grid(True, color="0.9", lw=0.6)
    fig.text(
        0.01,
        0.01,
        f"coarse Doppler {coarse['doppler_mps'][pulse]:.1f} m/s; "
        f"sinusoid Doppler {fit['doppler_mps']:.1f} m/s; "
        f"f={fit['freq_hz']:.1f} Hz; signs cal={signs[0]}, geom={signs[1]}",
        fontsize=9,
    )
    fig.savefig(args.output, dpi=170)
    print(args.output)
    print(f"pulse {pulse} beam {beam} range_bin {range_bin}")
    print(f"coarse_doppler_mps {coarse['doppler_mps'][pulse]:.6g}")
    print(f"fit_doppler_mps {fit['doppler_mps']:.6g}")
    print(f"fit_frequency_hz {fit['freq_hz']:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
