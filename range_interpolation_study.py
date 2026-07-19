#!/usr/bin/env python3
"""Range/Doppler peak interpolation study for one PANSY meteor cut."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as c
import scipy.signal as signal

import pansy_config as pc
from interferometer_alias_diagnostics import amp_scale, load_cut


DRG_KM = c.c / 2.0 / 1e6 / 1e3


METHOD_NAMES = {
    "nearest_range_km": "bin maximum",
    "range_gauss_km": "range log-quadratic",
    "peak2d_range_km": "2-D log-quadratic",
}


def quadratic_vertex(y_minus, y0, y_plus):
    """Return sub-bin vertex offset for three samples of a parabola."""
    denom = y_minus - 2.0 * y0 + y_plus
    if not np.isfinite(denom) or abs(denom) < 1e-12:
        return 0.0
    delta = 0.5 * (y_minus - y_plus) / denom
    if not np.isfinite(delta) or abs(delta) > 1.0:
        return 0.0
    return float(delta)


def quadratic_peak_2d(log_power, ri, di):
    """Fit a local 2-D quadratic to a 3x3 log-power patch."""
    nr, nd = log_power.shape
    if ri <= 0 or ri >= nr - 1 or di <= 0 or di >= nd - 1:
        return 0.0, 0.0
    xs = []
    ys = []
    for dr in (-1, 0, 1):
        for dd in (-1, 0, 1):
            x = float(dr)
            y = float(dd)
            xs.append([1.0, x, y, x * x, x * y, y * y])
            ys.append(float(log_power[ri + dr, di + dd]))
    coeff, *_ = np.linalg.lstsq(np.asarray(xs), np.asarray(ys), rcond=None)
    _a, bx, by, cxx, cxy, cyy = coeff
    hess = np.asarray([[2.0 * cxx, cxy], [cxy, 2.0 * cyy]], dtype=np.float64)
    grad = np.asarray([bx, by], dtype=np.float64)
    try:
        delta = -np.linalg.solve(hess, grad)
    except np.linalg.LinAlgError:
        return 0.0, 0.0
    if not np.all(np.isfinite(delta)) or np.any(np.abs(delta) > 1.0):
        return 0.0, 0.0
    return float(delta[0]), float(delta[1])


class PowerOnlyRangeDopplerSearch:
    """Range-Doppler matched-filter power, without interferometric cross spectra."""

    def __init__(self, txlen, echolen, n_channels, interp=1, fft_pad=1):
        self.interp = int(interp)
        self.txlen = int(txlen) * self.interp
        self.echolen = int(echolen) * self.interp
        self.n_channels = int(n_channels)
        self.fdec = 8 * self.interp
        self.fftlen = 512
        if int(self.txlen * self.interp / self.fdec) > self.fftlen:
            self.fftlen = int(self.txlen * self.interp / self.fdec)
        self.fftlen *= int(fft_pad)
        self.n_rg = self.echolen - self.txlen
        self.rg = np.arange(self.n_rg, dtype=np.int64)
        self.idx = np.arange(self.txlen, dtype=np.int64)
        self.idx_mat = self.rg[:, None] + self.idx[None, :]
        self.z_rx = np.zeros((self.n_channels, self.echolen), dtype=np.complex64)
        self.fvec = np.fft.fftshift(np.fft.fftfreq(self.fftlen, d=self.fdec / (self.interp * 1e6)))
        self.dopv = self.fvec * c.c / 2.0 / pc.freq

    def decim(self, values):
        new_width = int(self.txlen / self.fdec)
        out = np.zeros((self.n_rg, new_width), dtype=np.complex64)
        idx = np.arange(new_width, dtype=np.int64)
        idx2 = idx * self.fdec
        for i in range(self.fdec):
            out[:, idx] += values[:, idx2 + i]
        return out

    def mf(self, tx, echo):
        txi = np.repeat(tx, self.interp)
        for chi in range(self.n_channels):
            self.z_rx[chi, :] = np.repeat(echo[chi, :], self.interp)
        power = np.zeros((self.n_rg, self.fftlen), dtype=np.float32)
        for chi in range(self.n_channels):
            z = self.z_rx[chi, self.idx_mat] * txi[None, :]
            zd = self.decim(z)
            zf = np.fft.fftshift(np.fft.fft(zd, self.fftlen, axis=1), axes=1)
            power += zf.real**2 + zf.imag**2
        return power


class FFTFractionalRangeDopplerSearch:
    """FFT-bandlimited fractional range sampling of the complex ambiguity function."""

    def __init__(self, txlen, echolen, n_channels, range_oversample=4, fft_pad=1):
        self.txlen = int(txlen)
        self.echolen = int(echolen)
        self.n_channels = int(n_channels)
        self.range_oversample = int(range_oversample)
        if self.range_oversample < 1:
            raise ValueError("range_oversample must be positive")
        self.fdec = 8
        self.fftlen = 512 * int(fft_pad)
        self.n_rg_native = self.echolen - self.txlen
        self.n_rg = self.n_rg_native * self.range_oversample
        self.rg = np.arange(self.n_rg_native, dtype=np.int64)
        self.idx = np.arange(self.txlen, dtype=np.int64)
        self.idx_mat = self.rg[:, None] + self.idx[None, :]
        self.fvec = np.fft.fftshift(
            np.fft.fftfreq(self.fftlen, d=self.fdec / 1e6)
        )
        self.dopv = self.fvec * c.c / 2.0 / pc.freq

    def decim(self, values):
        width = self.txlen // self.fdec
        trimmed = values[:, : width * self.fdec]
        return trimmed.reshape(self.n_rg_native, width, self.fdec).sum(axis=2)

    def mf(self, tx, echo):
        power = np.zeros((self.n_rg, self.fftlen), dtype=np.float32)
        for channel in range(self.n_channels):
            product = echo[channel, self.idx_mat] * tx[None, :]
            doppler_spectrum = np.fft.fftshift(
                np.fft.fft(self.decim(product), self.fftlen, axis=1), axes=1
            )
            fractional_ambiguity = signal.resample(
                doppler_spectrum, self.n_rg, axis=0
            )
            power += (
                fractional_ambiguity.real**2 + fractional_ambiguity.imag**2
            ).astype(np.float32)
        return power


def robust_poly_residual(t, y, degree=2, mask=None, n_iter=5):
    t = np.asarray(t, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    good = np.isfinite(t) & np.isfinite(y)
    if mask is not None:
        good &= np.asarray(mask, dtype=bool)
    x = t - np.nanmean(t[good])
    keep = good.copy()
    for _ in range(n_iter):
        if np.count_nonzero(keep) < degree + 2:
            break
        coeff = np.polyfit(x[keep], y[keep], degree)
        pred = np.polyval(coeff, x)
        res = y - pred
        med = np.nanmedian(res[keep])
        sig = 1.4826 * np.nanmedian(np.abs(res[keep] - med))
        if not np.isfinite(sig) or sig <= 0.0:
            break
        keep = good & (np.abs(res - med) < 4.0 * sig)
    coeff = np.polyfit(x[keep], y[keep], degree)
    pred = np.polyval(coeff, x)
    res = y - pred
    use = keep & np.isfinite(res)
    rms = float(np.sqrt(np.nanmean(res[use] ** 2)))
    mad = float(1.4826 * np.nanmedian(np.abs(res[use] - np.nanmedian(res[use]))))
    return pred, res, use, rms, mad


def estimate_for_config(cut, pulses, interp, fft_pad):
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
    out = {
        "nearest_range_km": [],
        "range_gauss_km": [],
        "peak2d_range_km": [],
        "nearest_doppler_mps": [],
        "peak2d_doppler_mps": [],
        "snr": [],
    }
    for pulse in pulses:
        mf = rds.mf(z_tx[pulse, :], z_rx[pulse, :, :])
        noise = float(np.nanmedian(mf))
        pprof = np.nanmax(mf, axis=1)
        ri = int(np.nanargmax(pprof))
        di = int(np.nanargmax(mf[ri, :]))
        log_prof = np.log(np.maximum(pprof, max(noise, 1e-12)))
        dr_gauss = 0.0
        if 0 < ri < len(pprof) - 1:
            dr_gauss = quadratic_vertex(log_prof[ri - 1], log_prof[ri], log_prof[ri + 1])
        log_mf = np.log(np.maximum(mf, max(noise, 1e-12)))
        dr_2d, dd_2d = quadratic_peak_2d(log_mf, ri, di)
        out["nearest_range_km"].append((delays[pulse] + ri / interp) * DRG_KM)
        out["range_gauss_km"].append((delays[pulse] + (ri + dr_gauss) / interp) * DRG_KM)
        out["peak2d_range_km"].append((delays[pulse] + (ri + dr_2d) / interp) * DRG_KM)
        out["nearest_doppler_mps"].append(rds.dopv[di])
        dop_step = float(rds.dopv[1] - rds.dopv[0]) if len(rds.dopv) > 1 else 0.0
        out["peak2d_doppler_mps"].append(rds.dopv[di] + dd_2d * dop_step)
        out["snr"].append((pprof[ri] - noise) / max(noise, 1e-12))
    return {key: np.asarray(val, dtype=np.float64) for key, val in out.items()}


def load_winning_hypothesis(diag_h5):
    with h5py.File(diag_h5, "r") as h:
        label = h.attrs.get("selected_hypothesis", "")
        if not label:
            scored = []
            for name, grp in h["hypotheses"].items():
                if "combined_rank" in grp.attrs:
                    scored.append((int(grp.attrs["combined_rank"]), name))
            label = sorted(scored)[0][1]
        grp = h["hypotheses"][label]
        return {
            "label": str(label),
            "pulse": grp["pulse"][()],
            "t_rel_s": grp["t_rel_s"][()],
            "range_km": grp["range_km"][()],
            "doppler_mps": grp["doppler_mps"][()],
            "snr": grp["snr"][()],
            "ballistic_keep": grp["ballistic_keep"][()] if "ballistic_keep" in grp else np.ones(len(grp["pulse"]), dtype=bool),
        }


def write_latex_table(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("% Generated by range_interpolation_study.py\n")
        f.write("\\begin{tabular}{lrrrrr}\n")
        f.write("\\hline\n")
        f.write("Estimator & grid (m) & Doppler step (m s$^{-1}$) & range RMS (m) & range MAD (m) & Doppler RMS (m s$^{-1}$)\\\\\n")
        f.write("\\hline\n")
        for row in rows:
            f.write(
                f"{row['label']} & {row['range_grid_m']:.1f} & {row['doppler_step_mps']:.1f} & "
                f"{row['range_rms_m']:.1f} & {row['range_mad_m']:.1f} & {row['doppler_rms_mps']:.1f}\\\\\n"
            )
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")


def symmetric_residual_ylim(residuals, mask):
    values = np.asarray(residuals, dtype=np.float64)[np.asarray(mask, dtype=bool)]
    values = values[np.isfinite(values)]
    if values.size == 0:
        return (-1.0, 1.0)
    lim = float(np.nanpercentile(np.abs(values), 95.0) * 1.25)
    lim = max(lim, 250.0)
    return (-lim, lim)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-idx", type=int, default=1746489746272806)
    parser.add_argument("--cut-dir", type=Path, default=Path("data/metadata/cut"))
    parser.add_argument("--diagnostics", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=Path("test_plots/range_interpolation"))
    args = parser.parse_args()

    sample_idx = int(args.sample_idx)
    diag = args.diagnostics or Path(f"test_plots/pansy_disambiguation_diagnostics_{sample_idx}.h5")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cut = load_cut(args.cut_dir, sample_idx)
    hyp = load_winning_hypothesis(diag)
    pulses = np.asarray(hyp["pulse"], dtype=np.int64)
    t = np.asarray(hyp["t_rel_s"], dtype=np.float64)
    keep = np.asarray(hyp["ballistic_keep"], dtype=bool)

    configs = [
        (1, 1),
        (2, 1),
        (4, 1),
        (8, 1),
        (4, 2),
        (4, 4),
    ]
    estimates = {}
    rows = []
    for interp, fft_pad in configs:
        est = estimate_for_config(cut, pulses, interp=interp, fft_pad=fft_pad)
        estimates[(interp, fft_pad)] = est
        rds_tmp = PowerOnlyRangeDopplerSearch(
            txlen=np.asarray(cut["ztx_pulses_re"]).shape[1],
            echolen=np.asarray(cut["zrx_echoes_re"]).shape[2],
            n_channels=np.asarray(cut["zrx_echoes_re"]).shape[1],
            interp=interp,
            fft_pad=fft_pad,
        )
        dop_step = float(abs(rds_tmp.dopv[1] - rds_tmp.dopv[0]))
        for key, label in [
            ("nearest_range_km", "bin max"),
            ("range_gauss_km", "range Gaussian"),
            ("peak2d_range_km", "2-D Gaussian"),
        ]:
            pred_r, res_r, use_r, rms_r, mad_r = robust_poly_residual(t, est[key] * 1e3, degree=2, mask=keep)
            dop_key = "peak2d_doppler_mps" if key == "peak2d_range_km" else "nearest_doppler_mps"
            _pred_d, _res_d, _use_d, rms_d, _mad_d = robust_poly_residual(t, est[dop_key], degree=2, mask=keep)
            rows.append(
                {
                    "interp": interp,
                    "fft_pad": fft_pad,
                    "method": key,
                    "label": f"I{interp} P{fft_pad} {label}",
                    "range_grid_m": DRG_KM * 1e3 / interp,
                    "doppler_step_mps": dop_step,
                    "range_rms_m": rms_r,
                    "range_mad_m": mad_r,
                    "doppler_rms_mps": rms_d,
                    "n_used": int(np.count_nonzero(use_r)),
                }
            )

    rows.sort(key=lambda r: (r["range_rms_m"], r["doppler_rms_mps"]))
    csv_path = args.output_dir / f"range_interpolation_metrics_{sample_idx}.csv"
    with csv_path.open("w") as f:
        cols = ["label", "interp", "fft_pad", "method", "range_grid_m", "doppler_step_mps", "range_rms_m", "range_mad_m", "doppler_rms_mps", "n_used"]
        f.write(",".join(cols) + "\n")
        for row in rows:
            f.write(",".join(str(row[c]) for c in cols) + "\n")
    write_latex_table(args.output_dir / f"range_interpolation_metrics_{sample_idx}.tex", rows[:12])

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.0), constrained_layout=True)
    show = [
        ((1, 1), "nearest_range_km", "I1 bin maximum"),
        ((4, 1), "nearest_range_km", "I4 bin maximum"),
        ((4, 1), "range_gauss_km", "I4 range log-quadratic"),
        ((8, 1), "peak2d_range_km", "I8 2-D log-quadratic"),
    ]
    residual_ylim = None
    for cfg, key, label in show:
        est = estimates[cfg]
        pred, res, use, rms, mad = robust_poly_residual(t, est[key] * 1e3, degree=2, mask=keep)
        axes[0, 0].plot(t, res, ".", ms=4, label=f"{label}: RMS {rms:.0f} m")
        if cfg == (1, 1) and key == "nearest_range_km":
            residual_ylim = symmetric_residual_ylim(res, use)
    axes[0, 0].axhline(0.0, color="0.3", lw=0.8)
    axes[0, 0].set_xlabel("Time since cut start (s)")
    axes[0, 0].set_ylabel("Range residual after quadratic trend (m)")
    if residual_ylim is not None:
        axes[0, 0].set_ylim(*residual_ylim)
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(True, alpha=0.25)

    best = rows[0]
    base = [r for r in rows if r["interp"] == 1 and r["fft_pad"] == 1 and r["method"] == "nearest_range_km"][0]
    axes[0, 1].bar(["current", "best"], [base["range_rms_m"], best["range_rms_m"]], color=["0.55", "tab:blue"])
    axes[0, 1].set_ylabel("Range residual RMS (m)")
    axes[0, 1].set_title(f"Best: {best['label']}")
    axes[0, 1].grid(True, axis="y", alpha=0.25)

    for method, marker in [("nearest_range_km", "o"), ("range_gauss_km", "s"), ("peak2d_range_km", "^")]:
        subset = [r for r in rows if r["fft_pad"] == 1 and r["method"] == method]
        subset.sort(key=lambda r: r["interp"])
        axes[1, 0].plot(
            [r["range_grid_m"] for r in subset],
            [r["range_rms_m"] for r in subset],
            marker=marker,
            label=METHOD_NAMES[method],
        )
    axes[1, 0].invert_xaxis()
    axes[1, 0].set_xlabel("Range grid spacing (m)")
    axes[1, 0].set_ylabel("Range residual RMS (m)")
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].grid(True, alpha=0.25)

    subset = [r for r in rows if r["interp"] == 4 and r["method"] == "peak2d_range_km"]
    subset.sort(key=lambda r: r["fft_pad"])
    axes[1, 1].plot([r["doppler_step_mps"] for r in subset], [r["doppler_rms_mps"] for r in subset], "o-", color="tab:purple")
    axes[1, 1].invert_xaxis()
    axes[1, 1].set_xlabel("Doppler grid spacing (m/s)")
    axes[1, 1].set_ylabel("Doppler residual RMS (m/s)")
    axes[1, 1].grid(True, alpha=0.25)
    fig.suptitle(f"PANSY range/Doppler interpolation study, sample {sample_idx}, {hyp['label']}")
    fig.savefig(args.output_dir / f"range_interpolation_summary_{sample_idx}.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.5, 4.8), constrained_layout=True)
    cfg = (8, 1)
    est = estimates[cfg]
    pred, res, use, rms, mad = robust_poly_residual(t, est["peak2d_range_km"] * 1e3, degree=2, mask=keep)
    ax.plot(t, est["peak2d_range_km"], ".", ms=5, label="I8 2-D Gaussian range")
    ax.plot(t, pred / 1e3, "-", color="black", lw=1.3, label="quadratic trend used for residual metric")
    ax.set_xlabel("Time since cut start (s)")
    ax.set_ylabel("Range (km)")
    ax.set_title(f"Best high-resolution range track: RMS {rms:.1f} m, MAD {mad:.1f} m")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(args.output_dir / f"range_interpolation_track_{sample_idx}.png", dpi=220)
    plt.close(fig)

    print(f"selected_hypothesis {hyp['label']}")
    print(f"n_pulses {len(pulses)} n_metric {int(np.count_nonzero(keep))}")
    print(f"csv {csv_path}")
    print(f"best {best['label']} range_rms_m {best['range_rms_m']:.3f} doppler_rms_mps {best['doppler_rms_mps']:.3f}")
    print(f"current {base['label']} range_rms_m {base['range_rms_m']:.3f} doppler_rms_mps {base['doppler_rms_mps']:.3f}")


if __name__ == "__main__":
    main()
