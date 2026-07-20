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
from pansy_coherent import (
    coherent_add_modules,
    estimate_event_module_voltage_gain,
)


def receiver_noise_variance(z_rx: np.ndarray) -> np.ndarray:
    """Robust complex-noise variance per pulse and receiver module."""
    differences = np.diff(np.asarray(z_rx, dtype=np.complex128), axis=-1) / np.sqrt(2.0)
    center = np.median(differences.real, axis=-1) + 1j * np.median(
        differences.imag, axis=-1
    )
    variance = np.median(np.abs(differences - center[..., None]) ** 2, axis=-1) / np.log(2.0)
    fallback = np.nanmedian(variance, axis=0)
    bad = ~np.isfinite(variance) | (variance <= 0.0)
    variance[bad] = np.broadcast_to(fallback, variance.shape)[bad]
    return variance


def event_receiver_noise_variance(
    z_rx: np.ndarray,
    raw_indices: np.ndarray,
    range_gates: np.ndarray,
    decoded_length: int,
) -> np.ndarray:
    """Estimate module noise while excluding every decoded echo interval."""
    z_rx = np.asarray(z_rx, dtype=np.complex128)
    variance = np.full(z_rx.shape[:2], np.nan)
    by_raw: dict[int, list[int]] = {}
    for observation, raw_index in enumerate(raw_indices):
        by_raw.setdefault(int(raw_index), []).append(observation)
    for raw_index in range(z_rx.shape[0]):
        exclude = np.zeros(z_rx.shape[2], dtype=bool)
        for observation in by_raw.get(raw_index, []):
            start = int(np.floor(range_gates[observation])) - 6
            stop = int(np.ceil(range_gates[observation])) + decoded_length + 6
            exclude[max(start, 0) : min(stop, len(exclude))] = True
        for channel in range(z_rx.shape[1]):
            values = z_rx[raw_index, channel, ~exclude]
            center = np.median(values.real) + 1j * np.median(values.imag)
            estimate = np.median(np.abs(values - center) ** 2) / np.log(2.0)
            if not np.isfinite(estimate) or estimate <= 0.0:
                differences = np.diff(z_rx[raw_index, channel]) / np.sqrt(2.0)
                center = np.median(differences.real) + 1j * np.median(differences.imag)
                estimate = np.median(np.abs(differences - center) ** 2) / np.log(2.0)
            variance[raw_index, channel] = estimate
    fallback = np.nanmedian(variance, axis=0)
    bad = ~np.isfinite(variance) | (variance <= 0.0)
    variance[bad] = np.broadcast_to(fallback, variance.shape)[bad]
    return variance


def event_module_voltage_gain(response: np.ndarray) -> np.ndarray:
    """Estimate residual per-module voltage gain from one event.

    The input voltages have already received the static mesocal correction.  A
    robust two-way log-amplitude decomposition removes the common meteor
    amplitude at each pulse and leaves one constant residual gain per module.
    """
    return estimate_event_module_voltage_gain(response)


def fractional_delay_fft(x: np.ndarray, delay_samples: float) -> np.ndarray:
    """Apply a fractional sample delay using an FFT phase ramp."""
    if abs(delay_samples) < 1e-12:
        return x
    freq = np.fft.fftfreq(len(x))
    phase = np.exp(-1j * 2.0 * np.pi * freq * delay_samples)
    return np.fft.ifft(np.fft.fft(x) * phase)


def fit_weighted_complex_sinusoid(decoded: np.ndarray, tx: np.ndarray, coarse_doppler_mps: float, search_hz: float):
    """Fit decoded ~= C |tx|^2 exp(i 2 pi f t), weighted by TX envelope."""
    ns = np.arange(len(tx), dtype=np.float64)
    t = ns / 1e6
    env2 = np.abs(tx) ** 2
    env2 = env2 / max(float(np.nanmax(env2)), 1e-30)
    weights = env2.copy()
    weights[env2 <= 0.05] = 0.0

    coarse_hz = coarse_doppler_mps * 2.0 * pc.freq / sc.c
    freqs = np.linspace(coarse_hz - search_hz, coarse_hz + search_hz, 4001)
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
        "chi": float(chi[best]),
    }


def envelope_segments(env2: np.ndarray, threshold: float = 0.05) -> list[np.ndarray]:
    """Contiguous above-threshold TX-envelope blobs."""
    above = np.asarray(env2) > threshold
    idx = np.flatnonzero(above)
    if len(idx) == 0:
        return []
    breaks = np.where(np.diff(idx) > 1)[0] + 1
    return [seg for seg in np.split(idx, breaks) if len(seg) >= 2]


def fit_baud_complex_sinusoid(decoded: np.ndarray, tx: np.ndarray, coarse_doppler_mps: float, search_hz: float):
    """Estimate one matched-filter complex amplitude per TX baud blob, then fit Doppler."""
    ns = np.arange(len(tx), dtype=np.float64)
    t = ns / 1e6
    env2 = np.abs(tx) ** 2
    env2 = env2 / max(float(np.nanmax(env2)), 1e-30)
    segs = envelope_segments(env2)
    if len(segs) < 4:
        raise RuntimeError("too few TX envelope blobs for baud Doppler fit")
    baud_t = []
    baud_z = []
    baud_w = []
    for seg in segs:
        weight = env2[seg]
        den = np.sum(weight)
        baud_t.append(float(np.sum(t[seg] * weight) / max(float(den), 1e-30)))
        baud_z.append(np.sum(weight * decoded[seg]) / max(float(den), 1e-30))
        baud_w.append(float(den))
    baud_t = np.asarray(baud_t)
    baud_z = np.asarray(baud_z)
    baud_w = np.asarray(baud_w)
    baud_w = baud_w / max(float(np.nanmax(baud_w)), 1e-30)

    coarse_hz = coarse_doppler_mps * 2.0 * pc.freq / sc.c
    freqs = np.linspace(coarse_hz - search_hz, coarse_hz + search_hz, 4001)
    chi = np.empty(len(freqs), dtype=np.float64)
    coef = np.empty(len(freqs), dtype=np.complex128)
    for i, freq in enumerate(freqs):
        basis = np.exp(1j * 2.0 * np.pi * freq * baud_t)
        den = np.sum(baud_w * np.abs(basis) ** 2)
        c = np.sum(baud_w * np.conj(basis) * baud_z) / max(float(den), 1e-30)
        resid = baud_z - c * basis
        chi[i] = np.sum(baud_w * np.abs(resid) ** 2) / max(float(np.sum(baud_w)), 1e-30)
        coef[i] = c
    best = int(np.argmin(chi))
    freq = float(freqs[best])
    sample_model = coef[best] * env2 * np.exp(1j * 2.0 * np.pi * freq * t)
    baud_model = coef[best] * np.exp(1j * 2.0 * np.pi * freq * baud_t)
    return {
        "sample": ns,
        "time_s": t,
        "env2": env2,
        "freq_hz": freq,
        "doppler_mps": float(freq * sc.c / (2.0 * pc.freq)),
        "coef": coef[best],
        "model": sample_model,
        "resid": decoded - sample_model,
        "chi": float(chi[best]),
        "baud_t": baud_t,
        "baud_sample": baud_t * 1e6,
        "baud_z": baud_z,
        "baud_w": baud_w,
        "baud_model": baud_model,
        "n_bauds": int(len(baud_t)),
    }


def baud_fit_snr_db(fit: dict) -> float:
    """Model-to-residual SNR for the baud-integrated complex amplitudes."""
    signal = np.sum(fit["baud_w"] * np.abs(fit["baud_model"]) ** 2)
    resid = np.sum(fit["baud_w"] * np.abs(fit["baud_z"] - fit["baud_model"]) ** 2)
    return float(10.0 * np.log10(max(float(signal), 1e-300) / max(float(resid), 1e-300)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-idx", type=int, required=True)
    parser.add_argument("--cut-dir", type=Path, default=Path("/mnt/data/juha/pansy/metadata/cut"))
    parser.add_argument("--orbit-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--range-search-samples", type=int, default=8)
    parser.add_argument("--fractional-delay-step", type=float, default=0.05)
    parser.add_argument("--frequency-search-hz", type=float, default=20000.0)
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
    tx = z_tx[pulse].astype(np.complex128)
    best = None
    frac_delays = np.arange(-0.5, 0.5001, args.fractional_delay_step)
    summed_i = coherent_add_modules(
        z_rx[pulse].astype(np.complex128), uvw[pulse], beam
    )
    for _ in range(1):
        for range_bin_i in range(
            int(coarse["range_bin"][pulse]) - args.range_search_samples,
            int(coarse["range_bin"][pulse]) + args.range_search_samples + 1,
        ):
            if range_bin_i < 0 or range_bin_i + len(tx) > len(summed_i):
                continue
            for frac_delay in frac_delays:
                shifted_i = fractional_delay_fft(summed_i, frac_delay)
                decoded_i = shifted_i[range_bin_i : range_bin_i + len(tx)] * np.conj(tx)
                fit_i = fit_baud_complex_sinusoid(
                    decoded_i,
                    tx,
                    float(coarse["doppler_mps"][pulse]),
                    search_hz=float(args.frequency_search_hz),
                )
                snr_db_i = baud_fit_snr_db(fit_i)
                candidate = {
                    "score": snr_db_i,
                    "range_bin": range_bin_i,
                    "frac_delay": float(frac_delay),
                    "summed": summed_i,
                    "decoded": decoded_i,
                    "fit": fit_i,
                }
                if best is None or candidate["score"] > best["score"]:
                    best = candidate
    if best is None:
        raise RuntimeError("no valid range candidates")
    range_bin = best["range_bin"]
    frac_delay = best["frac_delay"]
    decoded = best["decoded"]
    fit = best["fit"]

    sample = fit["sample"]
    model = fit["model"]
    baud_phase = np.unwrap(np.angle(fit["baud_z"]))
    baud_model_phase = np.unwrap(np.angle(fit["baud_model"]))
    baud_phase -= np.nanmedian(baud_phase - baud_model_phase)

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
    sign_text = (
        "catalogue AoA steering\n"
        f"baud-fit SNR {best['score']:.1f} dB @ {range_bin}{frac_delay:+.2f}"
    )
    axes[0].text(
        0.01,
        0.96,
        sign_text,
        transform=axes[0].transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 2.5},
    )

    axes[1].plot(sample, decoded.real, ".", ms=3, color="tab:blue", alpha=0.75, label="decoded real")
    axes[1].plot(sample, model.real, "-", lw=1.4, color="navy", label="model real")
    axes[1].plot(sample, decoded.imag, ".", ms=3, color="tab:orange", alpha=0.75, label="decoded imag")
    axes[1].plot(sample, model.imag, "-", lw=1.4, color="darkorange", label="model imag")
    axes[1].set_ylabel("decoded amplitude")
    axes[1].legend(frameon=False, ncol=2, loc="upper right")

    axes[2].scatter(
        fit["baud_sample"],
        baud_phase,
        s=16.0 + 34.0 * fit["baud_w"],
        color="tab:purple",
        label="baud-integrated amplitude phase",
        zorder=3,
    )
    axes[2].plot(fit["baud_sample"], baud_model_phase, "-", lw=1.5, color="black", label="baud sinusoid phase")
    axes[2].set_ylabel("unwrapped phase (rad)")
    axes[2].legend(frameon=False, loc="upper left")

    axes[3].plot(fit["baud_sample"], np.abs(fit["baud_z"]), "o", ms=5, color="tab:blue", label="|baud amplitude|")
    axes[3].plot(
        fit["baud_sample"],
        np.abs(fit["baud_z"] - fit["baud_model"]),
        "o",
        ms=5,
        color="0.3",
        label="|baud residual|",
    )
    axes[3].set_ylabel("baud amplitude")
    axes[3].set_xlabel("sample within decoded TX pulse")
    axes[3].legend(frameon=False, loc="upper right")
    for ax in axes:
        ax.grid(True, color="0.9", lw=0.6)
    fig.text(
        0.01,
        0.01,
        f"coarse Doppler {coarse['doppler_mps'][pulse]:.1f} m/s; "
        f"sinusoid Doppler {fit['doppler_mps']:.1f} m/s; "
        f"f={fit['freq_hz']:.1f} Hz; range {range_bin}{frac_delay:+.2f}; "
        f"catalogue AoA steering; baud-fit SNR={best['score']:.1f} dB; chi={fit['chi']:.3g}",
        fontsize=9,
    )
    fig.savefig(args.output, dpi=170)
    print(args.output)
    print(f"pulse {pulse} beam {beam} range_bin {range_bin} fractional_delay {frac_delay:.6g}")
    print(f"coarse_doppler_mps {coarse['doppler_mps'][pulse]:.6g}")
    print(f"fit_doppler_mps {fit['doppler_mps']:.6g}")
    print(f"fit_frequency_hz {fit['freq_hz']:.6g}")
    print("steering_convention catalogue_arrival_uvw")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
