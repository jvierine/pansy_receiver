#!/usr/bin/env python3
"""Standalone short-lag within-pulse Doppler test for one PANSY meteor cut."""

from __future__ import annotations

import argparse
from pathlib import Path

import digital_rf as drf
import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as c
import scipy.optimize as opt

import pansy_config as pc
import stuffr
from doppler_phase import circular_phase_residual
from interferometer_alias_diagnostics import amp_scale, load_cut, recompute_cut_observables


DRG_KM = c.c / 2.0 / 1e6 / 1e3


def discover_sample(cut_dir: Path) -> int:
    dm = drf.DigitalMetadataReader(str(cut_dir))
    bounds = dm.get_bounds()
    records = dm.read(bounds[1] - 24 * 3600 * 1_000_000, bounds[1])
    if not records:
        records = dm.read(bounds[0], bounds[1])
    if not records:
        raise RuntimeError(f"no cut metadata records found in {cut_dir}")
    return int(sorted(records.keys())[-1])


def circular_lag_slope(lag_s, lag_product, omega0_rad_s):
    good = np.isfinite(lag_product.real) & np.isfinite(lag_product.imag) & (np.abs(lag_product) > 0.0)
    if np.count_nonzero(good) < 2:
        return np.nan, np.nan, np.nan
    lag_s = np.asarray(lag_s, dtype=np.float64)[good]
    lag_product = np.asarray(lag_product, dtype=np.complex128)[good]
    weight = np.sqrt(np.maximum(np.abs(lag_product), 1e-12))
    phase = np.angle(lag_product)

    def residual(par):
        return weight * circular_phase_residual(phase, par[0] * lag_s)

    result = opt.least_squares(residual, np.asarray([omega0_rad_s], dtype=np.float64), loss="linear")
    omega = float(result.x[0])
    res = circular_phase_residual(phase, omega * lag_s)
    rms = float(np.sqrt(np.mean(res**2)))
    return omega, rms, int(np.count_nonzero(good))


def moving_complex_average(values, weights, half_width):
    values = np.asarray(values, dtype=np.complex128)
    weights = np.asarray(weights, dtype=np.float64)
    out = np.full(values.shape, np.nan + 1j * np.nan, dtype=np.complex128)
    for i in range(values.shape[0]):
        lo = max(0, i - half_width)
        hi = min(values.shape[0], i + half_width + 1)
        lag_good = np.isfinite(values[lo:hi].real) & np.isfinite(values[lo:hi].imag)
        pulse_good = np.any(lag_good, axis=1) & np.isfinite(weights[lo:hi])
        if np.count_nonzero(pulse_good) == 0:
            continue
        ww = np.maximum(weights[lo:hi][pulse_good], 1e-6)
        vv = values[lo:hi][pulse_good].copy()
        vv[~lag_good[pulse_good]] = np.nan + 1j * np.nan
        num = np.nansum(ww[:, None] * vv, axis=0)
        den = np.nansum(ww[:, None] * np.isfinite(vv.real), axis=0)
        good_den = den > 0.0
        out[i, good_den] = num[good_den] / den[good_den]
    return out


def compute_short_lag_case(cut, lags=(4, 8, 13), neighbor_half_width=2, tx_conjugate=False):
    obs = recompute_cut_observables(cut)
    z_rx = np.asarray(cut["zrx_echoes_re"], dtype=np.complex64) + 1j * np.asarray(cut["zrx_echoes_im"], dtype=np.complex64)
    z_tx = np.asarray(cut["ztx_pulses_re"], dtype=np.complex64) + 1j * np.asarray(cut["ztx_pulses_im"], dtype=np.complex64)
    scale = amp_scale()
    for i in range(z_rx.shape[1]):
        z_rx[:, i, :] *= scale[i]

    lags = np.asarray(lags, dtype=np.int64)
    delays = np.asarray(cut["delays"], dtype=np.float64)
    range_gate = np.rint(np.asarray(obs["range_km"], dtype=np.float64) / DRG_KM - delays).astype(np.int64)
    txlen = z_tx.shape[1]
    n_pulse = z_rx.shape[0]
    n_lag = len(lags)
    lag_product = np.full((n_pulse, n_lag), np.nan + 1j * np.nan, dtype=np.complex128)

    for pi in range(n_pulse):
        rg = int(range_gate[pi])
        if rg < 0 or rg + txlen > z_rx.shape[2]:
            continue
        tx = np.conj(z_tx[pi]) if tx_conjugate else z_tx[pi]
        for li, lag in enumerate(lags):
            if lag <= 0 or lag >= txlen:
                continue
            acc = 0.0j
            for chi in range(z_rx.shape[1]):
                base = z_rx[pi, chi, rg : rg + txlen] * tx
                acc += np.sum(base[lag:] * np.conj(base[:-lag]))
            lag_product[pi, li] = acc

    lag_product_avg = moving_complex_average(lag_product, obs["snr"], neighbor_half_width)
    lag_s = lags.astype(np.float64) / 1e6
    coarse_omega = 2.0 * np.pi * (2.0 * np.asarray(obs["doppler_mps"], dtype=np.float64) / pc.wavelength)

    omega = np.full(n_pulse, np.nan, dtype=np.float64)
    omega_rms = np.full(n_pulse, np.nan, dtype=np.float64)
    omega_n = np.zeros(n_pulse, dtype=np.int64)
    omega_avg = np.full(n_pulse, np.nan, dtype=np.float64)
    omega_avg_rms = np.full(n_pulse, np.nan, dtype=np.float64)
    omega_avg_n = np.zeros(n_pulse, dtype=np.int64)
    for pi in range(n_pulse):
        omega[pi], omega_rms[pi], omega_n[pi] = circular_lag_slope(lag_s, lag_product[pi], coarse_omega[pi])
        omega_avg[pi], omega_avg_rms[pi], omega_avg_n[pi] = circular_lag_slope(lag_s, lag_product_avg[pi], coarse_omega[pi])

    short_lag_doppler_mps = omega * pc.wavelength / (4.0 * np.pi)
    short_lag_avg_doppler_mps = omega_avg * pc.wavelength / (4.0 * np.pi)

    err_raw = short_lag_avg_doppler_mps - obs["doppler_mps"]
    err_flip = -short_lag_avg_doppler_mps - obs["doppler_mps"]
    use = np.isfinite(err_raw) & np.isfinite(err_flip) & (obs["snr"] > 9.0)
    sign = 1.0
    if np.count_nonzero(use) and np.nanmedian(np.abs(err_flip[use])) < np.nanmedian(np.abs(err_raw[use])):
        sign = -1.0

    return {
        "tx_idx": np.asarray(obs["tx_idx"], dtype=np.float64),
        "beam_id": np.asarray(obs["beam_id"], dtype=np.int64),
        "range_km": np.asarray(obs["range_km"], dtype=np.float64),
        "range_gate": range_gate,
        "snr": np.asarray(obs["snr"], dtype=np.float64),
        "coarse_doppler_mps": np.asarray(obs["doppler_mps"], dtype=np.float64),
        "lags_samples": lags,
        "lags_seconds": lag_s,
        "lag_product": lag_product,
        "lag_product_avg": lag_product_avg,
        "lag_phase_rad": np.angle(lag_product),
        "lag_phase_avg_rad": np.angle(lag_product_avg),
        "short_lag_doppler_mps": short_lag_doppler_mps,
        "short_lag_residual_rms_rad": omega_rms,
        "short_lag_n": omega_n,
        "short_lag_avg_doppler_mps": short_lag_avg_doppler_mps,
        "short_lag_avg_doppler_sign_aligned_mps": sign * short_lag_avg_doppler_mps,
        "short_lag_avg_residual_rms_rad": omega_avg_rms,
        "short_lag_avg_n": omega_avg_n,
        "sign_alignment": np.asarray(sign, dtype=np.float64),
        "neighbor_half_width": np.asarray(neighbor_half_width, dtype=np.int64),
        "tx_conjugate": np.asarray(int(tx_conjugate), dtype=np.int64),
    }


def write_hdf5(result, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        handle.attrs["schema_version"] = "pansy_short_lag_doppler_test_v1"
        handle.attrs["source_program"] = "short_lag_doppler_test_case.py"
        handle.attrs["wavelength_m"] = pc.wavelength
        handle.attrs["sample_rate_hz"] = 1e6
        for key, value in result.items():
            handle.create_dataset(key, data=value)


def plot_case(result, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    t_rel = result["tx_idx"] / 1e6 - result["tx_idx"][0] / 1e6
    snr_db = 10.0 * np.log10(np.maximum(result["snr"], 1e-6))
    good = np.isfinite(result["short_lag_avg_doppler_sign_aligned_mps"])

    fig, axes = plt.subplots(4, 1, figsize=(9.0, 9.5), sharex=True, constrained_layout=True)
    ax = axes[0]
    sc = ax.scatter(t_rel, result["coarse_doppler_mps"] / 1e3, c=snr_db, s=18, cmap="viridis", label="FFT peak")
    ax.plot(t_rel[good], result["short_lag_avg_doppler_sign_aligned_mps"][good] / 1e3, ".", ms=4.0, color="tab:red", label="short-lag avg")
    ax.set_ylabel("Doppler (km/s)")
    ax.legend(fontsize=8)
    fig.colorbar(sc, ax=ax, label="SNR (dB)")

    ax = axes[1]
    delta = result["short_lag_avg_doppler_sign_aligned_mps"] - result["coarse_doppler_mps"]
    ax.axhline(0.0, color="0.25", lw=0.8)
    ax.plot(t_rel[good], delta[good], ".", color="tab:purple")
    ax.set_ylabel("short-lag - FFT (m/s)")

    ax = axes[2]
    for li, lag in enumerate(result["lags_samples"]):
        ax.plot(t_rel, result["lag_phase_avg_rad"][:, li], ".", ms=3.0, label=f"lag {int(lag)}")
    ax.set_ylabel("Lag phase (rad)")
    ax.set_ylim(-np.pi, np.pi)
    ax.legend(fontsize=8, ncol=len(result["lags_samples"]))

    ax = axes[3]
    ax.plot(t_rel, result["short_lag_avg_residual_rms_rad"], ".", color="tab:orange")
    ax.set_xlabel("Time from first cut pulse (s)")
    ax.set_ylabel("Lag slope RMS (rad)")

    fig.suptitle(
        "Short-lag within-pulse Doppler test "
        f"(lags {','.join(str(int(x)) for x in result['lags_samples'])} samples)"
    )
    fig.savefig(output, dpi=220)
    plt.close(fig)


def write_latex_memo(sample_idx: int, result, figure_path: Path, h5_path: Path, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    use = np.isfinite(result["short_lag_avg_doppler_sign_aligned_mps"]) & (result["snr"] > 9.0)
    delta = result["short_lag_avg_doppler_sign_aligned_mps"] - result["coarse_doppler_mps"]
    rms = float(np.sqrt(np.nanmean(delta[use] ** 2))) if np.count_nonzero(use) else np.nan
    mad = float(1.4826 * np.nanmedian(np.abs(delta[use] - np.nanmedian(delta[use])))) if np.count_nonzero(use) else np.nan
    text = rf"""\documentclass[11pt]{{article}}
\usepackage[margin=1in]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{siunitx}}

\newif\ifshowscriptprovenance
\showscriptprovenancetrue
\newcommand{{\scriptprovenance}}[1]{{\ifshowscriptprovenance\par\small\emph{{Script: #1}}\fi}}

\begin{{document}}
\section*{{Short-lag within-pulse Doppler test}}

Cut sample {sample_idx} begins at {stuffr.unix2datestr(sample_idx / 1e6)} UTC.
The estimator uses lag products at {", ".join(str(int(x)) for x in result["lags_samples"])}
samples, then averages lag products over neighboring pulses before solving for
the wrapped phase slope versus lag time. Data products are stored in
\texttt{{{h5_path.name}}}.

\begin{{figure}}[h]
  \centering
  \includegraphics[width=\linewidth]{{{figure_path.name}}}
  \caption{{Short-lag Doppler comparison against the ordinary within-pulse FFT
  peak estimator. \scriptprovenance{{short_lag_doppler_test_case.py --sample-idx {sample_idx}}}}}
\end{{figure}}

\begin{{tabular}}{{lr}}
\hline
SNR-selected samples & {int(np.count_nonzero(use))} \\
Short-lag minus FFT RMS & \SI{{{rms:.1f}}}{{m.s^{{-1}}}} \\
Short-lag minus FFT MAD & \SI{{{mad:.1f}}}{{m.s^{{-1}}}} \\
\hline
\end{{tabular}}

\end{{document}}
"""
    output.write_text(text)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cut-dir", type=Path, default=Path("data/metadata/cut"))
    parser.add_argument("--sample-idx", type=int)
    parser.add_argument("--lags", type=int, nargs="+", default=[4, 8, 13])
    parser.add_argument("--neighbor-half-width", type=int, default=2)
    parser.add_argument("--tx-conjugate", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("memos/short_lag_doppler_test"))
    args = parser.parse_args()

    sample_idx = args.sample_idx if args.sample_idx is not None else discover_sample(args.cut_dir)
    cut = load_cut(args.cut_dir, sample_idx)
    result = compute_short_lag_case(
        cut,
        lags=args.lags,
        neighbor_half_width=args.neighbor_half_width,
        tx_conjugate=args.tx_conjugate,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fig_path = args.output_dir / f"short_lag_doppler_{sample_idx}.png"
    h5_path = args.output_dir / f"short_lag_doppler_{sample_idx}.h5"
    tex_path = args.output_dir / f"short_lag_doppler_{sample_idx}.tex"
    plot_case(result, fig_path)
    write_hdf5(result, h5_path)
    write_latex_memo(sample_idx, result, fig_path, h5_path, tex_path)

    use = np.isfinite(result["short_lag_avg_doppler_sign_aligned_mps"]) & (result["snr"] > 9.0)
    delta = result["short_lag_avg_doppler_sign_aligned_mps"] - result["coarse_doppler_mps"]
    rms = float(np.sqrt(np.nanmean(delta[use] ** 2))) if np.count_nonzero(use) else np.nan
    mad = float(1.4826 * np.nanmedian(np.abs(delta[use] - np.nanmedian(delta[use])))) if np.count_nonzero(use) else np.nan
    print(f"sample_idx {sample_idx} {stuffr.unix2datestr(sample_idx / 1e6)}")
    print(f"lags_samples {' '.join(str(int(x)) for x in result['lags_samples'])}")
    print(f"neighbor_half_width {args.neighbor_half_width}")
    print(f"sign_alignment {float(result['sign_alignment']):.0f}")
    print(f"snr_selected {int(np.count_nonzero(use))}")
    print(f"short_lag_minus_fft_rms_mps {rms:.3f}")
    print(f"short_lag_minus_fft_mad_mps {mad:.3f}")
    print(f"wrote {fig_path}")
    print(f"wrote {h5_path}")
    print(f"wrote {tex_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
