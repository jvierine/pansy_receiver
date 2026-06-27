#!/usr/bin/env python3
"""Test sinc-shaped sub-bin Doppler fits on one PANSY meteor cut."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

from interferometer_alias_diagnostics import amp_scale, load_cut
from range_interpolation_study import (
    DRG_KM,
    PowerOnlyRangeDopplerSearch,
    load_winning_hypothesis,
    quadratic_vertex,
    robust_poly_residual,
)


def dirichlet_power(rel_bins, frac_bin, n_time, fftlen):
    """Squared finite-record FFT response for a tone offset by frac_bin bins."""
    x = (np.asarray(rel_bins, dtype=np.float64) - float(frac_bin)) / float(fftlen)
    den = np.sin(np.pi * x)
    num = np.sin(np.pi * n_time * x)
    out = np.empty_like(x, dtype=np.float64)
    small = np.abs(den) < 1e-12
    out[small] = 1.0
    out[~small] = (num[~small] / (n_time * den[~small])) ** 2.0
    return out


def fit_sinc_doppler(profile, peak_bin, n_time, fftlen, fit_bins=5):
    """Fit background + amplitude * Dirichlet-power lobe around one Doppler peak."""
    profile = np.asarray(profile, dtype=np.float64)
    lo = max(0, int(peak_bin) - int(fit_bins))
    hi = min(len(profile), int(peak_bin) + int(fit_bins) + 1)
    bins = np.arange(lo, hi, dtype=np.int64)
    rel = bins - int(peak_bin)
    y = np.maximum(profile[bins], 1e-30)
    floor0 = max(float(np.nanmedian(profile)), 1e-30)
    amp0 = max(float(profile[int(peak_bin)] - floor0), floor0)

    def model(par):
        frac, log_amp, log_floor = par
        return np.exp(log_floor) + np.exp(log_amp) * dirichlet_power(rel, frac, n_time, fftlen)

    def residual(par):
        return np.log(y) - np.log(np.maximum(model(par), 1e-30))

    result = opt.least_squares(
        residual,
        np.asarray([0.0, np.log(amp0), np.log(floor0)], dtype=np.float64),
        bounds=([-0.75, np.log(1e-30), np.log(1e-30)], [0.75, np.log(1e40), np.log(1e40)]),
        loss="soft_l1",
        f_scale=0.25,
        max_nfev=200,
    )
    frac = float(result.x[0])
    fit_rms = float(np.sqrt(np.mean(residual(result.x) ** 2.0)))
    return {
        "frac_bin": frac,
        "amplitude": float(np.exp(result.x[1])),
        "floor": float(np.exp(result.x[2])),
        "fit_rms_logpower": fit_rms,
        "success": bool(result.success and np.isfinite(frac) and abs(frac) < 0.749),
        "bins": bins,
        "model_power": model(result.x),
    }


def range_gate_highres_doppler(tx, echo, range_bin, interp, fdec, fftlen, oversample):
    txi = np.repeat(tx, interp)
    fftlen_hi = int(fftlen) * int(oversample)
    n_channels = echo.shape[0]
    txlen = len(txi)
    z_rx = np.zeros((n_channels, echo.shape[1] * interp), dtype=np.complex64)
    idx = np.arange(int(txlen / fdec), dtype=np.int64)
    idx2 = idx * fdec
    power = np.zeros(fftlen_hi, dtype=np.float64)
    for chi in range(n_channels):
        z_rx[chi, :] = np.repeat(echo[chi, :], interp)
        z = z_rx[chi, range_bin : range_bin + txlen] * txi
        zd = np.zeros(len(idx), dtype=np.complex64)
        for fi in range(fdec):
            zd += z[idx2 + fi]
        zf = np.fft.fftshift(np.fft.fft(zd, fftlen_hi))
        power += zf.real**2.0 + zf.imag**2.0
    return power


def estimate_dopplers(cut, pulses, interp=1, fft_pad=1, fit_bins=5, oversample=16):
    z_rx = np.asarray(cut["zrx_echoes_re"], dtype=np.complex64) + 1j * np.asarray(cut["zrx_echoes_im"], dtype=np.complex64)
    z_tx = np.asarray(cut["ztx_pulses_re"], dtype=np.complex64) + 1j * np.asarray(cut["ztx_pulses_im"], dtype=np.complex64)
    scale = amp_scale()
    for i in range(z_rx.shape[1]):
        z_rx[:, i, :] *= scale[i]

    delays = np.asarray(cut["delays"], dtype=np.int64)
    rds = PowerOnlyRangeDopplerSearch(
        txlen=z_tx.shape[1],
        echolen=z_rx.shape[2],
        n_channels=z_rx.shape[1],
        interp=interp,
        fft_pad=fft_pad,
    )
    dop_step = float(rds.dopv[1] - rds.dopv[0])
    n_time = int(rds.txlen / rds.fdec)

    out = {
        "pulse": [],
        "range_km": [],
        "snr": [],
        "doppler_bin_mps": [],
        "doppler_parabolic_mps": [],
        "doppler_sinc_mps": [],
        "doppler_sinc_interp_mps": [],
        "parabolic_frac_bin": [],
        "sinc_frac_bin": [],
        "sinc_interp_frac_bin": [],
        "sinc_fit_rms_logpower": [],
        "sinc_success": [],
        "peak_bin": [],
        "peak_range_bin": [],
    }
    example = None
    for pulse in pulses:
        power = rds.mf(z_tx[pulse, :], z_rx[pulse, :, :])
        noise = float(np.nanmedian(power))
        pprof = np.nanmax(power, axis=1)
        ri = int(np.nanargmax(pprof))
        profile = power[ri, :]
        di = int(np.nanargmax(profile))
        log_profile = np.log(np.maximum(profile, max(noise, 1e-30)))
        frac_parabola = 0.0
        if 0 < di < len(profile) - 1:
            frac_parabola = quadratic_vertex(log_profile[di - 1], log_profile[di], log_profile[di + 1])
        sinc = fit_sinc_doppler(profile, di, n_time=n_time, fftlen=rds.fftlen, fit_bins=fit_bins)
        highres = range_gate_highres_doppler(
            z_tx[pulse, :],
            z_rx[pulse, :, :],
            ri,
            interp=interp,
            fdec=rds.fdec,
            fftlen=rds.fftlen,
            oversample=oversample,
        )
        center_hi = int(di * oversample)
        hi_lo = max(0, center_hi - oversample)
        hi_hi = min(len(highres), center_hi + oversample + 1)
        di_hi = int(hi_lo + np.nanargmax(highres[hi_lo:hi_hi]))
        frac_interp = (di_hi - center_hi) / float(oversample)

        out["pulse"].append(int(pulse))
        out["range_km"].append((delays[pulse] + ri / interp) * DRG_KM)
        out["snr"].append((pprof[ri] - noise) / max(noise, 1e-30))
        out["doppler_bin_mps"].append(rds.dopv[di])
        out["doppler_parabolic_mps"].append(rds.dopv[di] + frac_parabola * dop_step)
        out["doppler_sinc_mps"].append(rds.dopv[di] + sinc["frac_bin"] * dop_step)
        out["doppler_sinc_interp_mps"].append(rds.dopv[di] + frac_interp * dop_step)
        out["parabolic_frac_bin"].append(frac_parabola)
        out["sinc_frac_bin"].append(sinc["frac_bin"])
        out["sinc_interp_frac_bin"].append(frac_interp)
        out["sinc_fit_rms_logpower"].append(sinc["fit_rms_logpower"])
        out["sinc_success"].append(sinc["success"])
        out["peak_bin"].append(di)
        out["peak_range_bin"].append(ri)
        if example is None or out["snr"][-1] > example["snr"]:
            example = {
                "pulse": int(pulse),
                "snr": float(out["snr"][-1]),
                "rel_bins": sinc["bins"] - di,
                "power": profile[sinc["bins"]],
                "model_power": sinc["model_power"],
            }
    return {key: np.asarray(val) for key, val in out.items()}, example, dop_step, n_time


def metric_rows(t, keep, estimates):
    rows = []
    for key, label in [
        ("doppler_bin_mps", "FFT bin max"),
        ("doppler_parabolic_mps", "3-point log parabola"),
        ("doppler_sinc_mps", "sinc lobe fit"),
        ("doppler_sinc_interp_mps", "zero-padded sinc interp"),
    ]:
        _pred, res, use, rms, mad = robust_poly_residual(t, estimates[key], degree=2, mask=keep)
        rows.append(
            {
                "key": key,
                "label": label,
                "doppler_rms_mps": rms,
                "doppler_mad_mps": mad,
                "n_used": int(np.count_nonzero(use)),
                "residual_mps": res,
                "use": use,
            }
        )
    rows.sort(key=lambda row: row["doppler_rms_mps"])
    return rows


def write_h5(path, sample_idx, hyp, estimates, rows, dop_step, n_time, oversample):
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h:
        h.attrs["source_program"] = "sinc_doppler_peak_fit_test.py"
        h.attrs["schema_version"] = "sinc_doppler_peak_fit_test_v1"
        h.attrs["sample_idx"] = int(sample_idx)
        h.attrs["selected_hypothesis"] = str(hyp["label"])
        h.attrs["doppler_bin_spacing_mps"] = float(abs(dop_step))
        h.attrs["matched_filter_time_samples"] = int(n_time)
        h.attrs["sinc_interp_oversample"] = int(oversample)
        data = h.create_group("pulse_estimates")
        for key, value in estimates.items():
            data.create_dataset(key, data=value)
        metrics = h.create_group("metrics")
        string_dtype = h5py.string_dtype(encoding="utf-8")
        metrics.create_dataset("label", data=np.asarray([r["label"] for r in rows], dtype=string_dtype))
        metrics.create_dataset("key", data=np.asarray([r["key"] for r in rows], dtype=string_dtype))
        metrics.create_dataset("doppler_rms_mps", data=np.asarray([r["doppler_rms_mps"] for r in rows]))
        metrics.create_dataset("doppler_mad_mps", data=np.asarray([r["doppler_mad_mps"] for r in rows]))
        metrics.create_dataset("n_used", data=np.asarray([r["n_used"] for r in rows], dtype=np.int64))


def plot_results(path, t, keep, estimates, rows, example, sample_idx, hyp_label):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0), constrained_layout=True)

    ax = axes[0, 0]
    for row in rows:
        ax.plot(t, row["residual_mps"], ".", ms=4, label=f"{row['label']}: RMS {row['doppler_rms_mps']:.0f} m/s")
    ax.axhline(0.0, color="0.35", lw=0.8)
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Doppler residual after quadratic trend (m/s)")
    ax.set_title("Pulse-wise Doppler residual metric")
    residual_scale = []
    for row in rows:
        vals = np.asarray(row["residual_mps"], dtype=np.float64)[np.asarray(row["use"], dtype=bool)]
        residual_scale.extend(vals[np.isfinite(vals)])
    if residual_scale:
        lim = max(1500.0, 1.5 * float(np.nanpercentile(np.abs(residual_scale), 95.0)))
        ax.set_ylim(-lim, lim)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    labels = [row["label"] for row in rows]
    rms = [row["doppler_rms_mps"] for row in rows]
    ax.bar(np.arange(len(rows)), rms, color=["tab:blue", "tab:green", "tab:orange", "tab:purple"][: len(rows)])
    ax.set_xticks(np.arange(len(rows)))
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("Doppler residual RMS (m/s)")
    ax.set_title("Lower is smoother against the trajectory proxy")
    ax.grid(True, axis="y", alpha=0.25)

    ax = axes[1, 0]
    ax.plot(t, estimates["doppler_bin_mps"] / 1e3, ".", ms=4, label="FFT bin max")
    ax.plot(t, estimates["doppler_sinc_mps"] / 1e3, ".", ms=4, label="sinc lobe fit")
    ax.plot(t, estimates["doppler_sinc_interp_mps"] / 1e3, ".", ms=4, label="zero-padded sinc interp")
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Doppler velocity (km/s)")
    ax.set_title("Doppler estimates")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    if example is not None:
        ax.plot(example["rel_bins"], example["power"], "o", label=f"pulse {example['pulse']}")
        ax.plot(example["rel_bins"], example["model_power"], "-", label="sinc lobe model")
        ax.set_yscale("log")
    ax.set_xlabel("Doppler bin offset from peak")
    ax.set_ylabel("Matched-filter power")
    ax.set_title("Example local Doppler peak fit")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)

    fig.suptitle(f"Sinc Doppler peak fit test, sample {sample_idx}, {hyp_label}")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_latex_table(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("% Generated by sinc_doppler_peak_fit_test.py\n")
        f.write("\\begin{tabular}{lrrr}\n")
        f.write("\\hline\n")
        f.write("Estimator & Doppler RMS (m s$^{-1}$) & Doppler MAD (m s$^{-1}$) & pulses\\\\\n")
        f.write("\\hline\n")
        for row in rows:
            f.write(f"{row['label']} & {row['doppler_rms_mps']:.1f} & {row['doppler_mad_mps']:.1f} & {row['n_used']}\\\\\n")
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")


def main():
    parser = argparse.ArgumentParser(description="Test sinc-shaped sub-bin Doppler fits on one meteor cut.")
    parser.add_argument("--sample-idx", type=int, default=1746489746272806)
    parser.add_argument("--cut-dir", type=Path, default=Path("data/metadata/cut"))
    parser.add_argument("--diagnostics", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=Path("test_plots/sinc_doppler"))
    parser.add_argument("--interp", type=int, default=1)
    parser.add_argument("--fft-pad", type=int, default=1)
    parser.add_argument("--fit-bins", type=int, default=5)
    parser.add_argument("--oversample", type=int, default=16)
    args = parser.parse_args()

    sample_idx = int(args.sample_idx)
    diag = args.diagnostics or Path(f"test_plots/pansy_disambiguation_diagnostics_{sample_idx}.h5")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cut = load_cut(args.cut_dir, sample_idx)
    hyp = load_winning_hypothesis(diag)
    pulses = np.asarray(hyp["pulse"], dtype=np.int64)
    t = np.asarray(hyp["t_rel_s"], dtype=np.float64)
    keep = np.asarray(hyp["ballistic_keep"], dtype=bool)
    estimates, example, dop_step, n_time = estimate_dopplers(
        cut, pulses, interp=args.interp, fft_pad=args.fft_pad, fit_bins=args.fit_bins, oversample=args.oversample
    )
    estimates["t_rel_s"] = t
    estimates["ballistic_keep"] = keep.astype(np.int8)
    rows = metric_rows(t, keep, estimates)

    h5_path = args.output_dir / f"sinc_doppler_peak_fit_{sample_idx}.h5"
    figure_path = args.output_dir / f"sinc_doppler_peak_fit_{sample_idx}.png"
    tex_path = args.output_dir / f"sinc_doppler_peak_fit_metrics_{sample_idx}.tex"
    write_h5(h5_path, sample_idx, hyp, estimates, rows, dop_step, n_time, args.oversample)
    plot_results(figure_path, t, keep, estimates, rows, example, sample_idx, hyp["label"])
    write_latex_table(tex_path, rows)

    print(f"selected_hypothesis {hyp['label']}")
    print(f"n_pulses {len(pulses)} n_metric {int(np.count_nonzero(keep))}")
    print(f"doppler_bin_spacing_mps {abs(dop_step):.3f}")
    for row in rows:
        print(
            f"{row['key']} doppler_rms_mps {row['doppler_rms_mps']:.3f} "
            f"doppler_mad_mps {row['doppler_mad_mps']:.3f} n_used {row['n_used']}"
        )
    print(f"h5 {h5_path}")
    print(f"figure {figure_path}")
    print(f"tex {tex_path}")


if __name__ == "__main__":
    main()
