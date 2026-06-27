#!/usr/bin/env python3
"""Standalone wrapped pulse-to-pulse Doppler phase test for one PANSY cut."""

from __future__ import annotations

import argparse
from pathlib import Path

import digital_rf as drf
import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as c
import scipy.fftpack as fp

import pansy_config as pc
import stuffr
from doppler_phase import (
    fit_wrapped_phase_time,
    phase_fit_doppler_mps,
    pulse_pair_phase_samples,
    pulse_phase_diagnostics,
)
from interferometer_alias_diagnostics import amp_scale, load_cut


class ResponseRangeDopplerSearch:
    """Matched-filter bank with complex channel response at the selected peak."""

    def __init__(self, txlen, echolen, n_channels, interp=1):
        self.txlen = int(txlen * interp)
        self.fdec = int(8 * interp)
        self.fftlen = 512
        if int(self.txlen * interp / self.fdec) > self.fftlen:
            self.fftlen = int(self.txlen * interp / self.fdec)
        self.echolen = int(echolen * interp)
        self.n_channels = int(n_channels)
        self.interp = int(interp)
        self.n_rg = self.echolen - self.txlen
        self.rg = np.arange(self.n_rg, dtype=np.int64)
        self.idx = np.arange(self.txlen, dtype=np.int64)
        self.idx_mat = np.zeros((self.n_rg, self.txlen), dtype=np.int64)
        for ri in range(self.n_rg):
            self.idx_mat[ri, :] = self.idx + self.rg[ri]
        self.z_rx = np.zeros((self.n_channels, self.echolen), dtype=np.complex64)
        self.rangev = self.rg * c.c / 2.0 / 1e6 / 1e3
        self.fvec = np.fft.fftshift(np.fft.fftfreq(self.fftlen, d=self.fdec / (interp * 1e6)))
        self.dopv = self.fvec * c.c / 2.0 / pc.freq

    def decim(self, z):
        new_width = int(self.txlen / self.fdec)
        z2 = np.zeros((self.n_rg, new_width), dtype=np.complex64)
        idx = np.arange(new_width, dtype=np.int64)
        idx2 = np.arange(new_width, dtype=np.int64) * self.fdec
        for i in range(self.fdec):
            z2[:, idx] += z[:, idx2 + i]
        return z2

    def mf(self, tx, echo):
        txi = np.repeat(tx, self.interp)
        for chi in range(self.n_channels):
            self.z_rx[chi, :] = np.repeat(echo[chi, :], self.interp)

        mf = np.zeros((self.n_rg, self.fftlen), dtype=np.float32)
        zf_channels = []
        for chi in range(self.n_channels):
            z = self.z_rx[chi, self.idx_mat] * txi[None, :]
            zd = self.decim(z)
            zf = np.fft.fftshift(fp.fft(zd, self.fftlen, axis=1), axes=1)
            zf_channels.append(zf)
            mf += zf.real**2.0 + zf.imag**2.0

        noise_floor = float(np.median(mf))
        pprof = np.max(mf, axis=1)
        dop_idx = np.argmax(mf, axis=1)
        max_rg = int(np.argmax(pprof))
        max_dop_idx = int(dop_idx[max_rg])
        response = np.asarray([zf[max_rg, max_dop_idx] for zf in zf_channels], dtype=np.complex64)
        return {
            "pprof": pprof,
            "noise_floor": noise_floor,
            "range_gate": max_rg,
            "doppler_mps": float(self.dopv[max_dop_idx]),
            "response": response,
        }


def discover_sample(cut_dir: Path) -> int:
    dm = drf.DigitalMetadataReader(str(cut_dir))
    bounds = dm.get_bounds()
    records = dm.read(bounds[1] - 24 * 3600 * 1_000_000, bounds[1])
    if not records:
        records = dm.read(bounds[0], bounds[1])
    if not records:
        raise RuntimeError(f"no cut metadata records found in {cut_dir}")
    return int(sorted(records.keys())[-1])


def compute_pulse_phase_case(cut, interp=1):
    z_rx = np.asarray(cut["zrx_echoes_re"], dtype=np.complex64) + 1j * np.asarray(cut["zrx_echoes_im"], dtype=np.complex64)
    z_tx = np.asarray(cut["ztx_pulses_re"], dtype=np.complex64) + 1j * np.asarray(cut["ztx_pulses_im"], dtype=np.complex64)
    scale = amp_scale()
    for i in range(z_rx.shape[1]):
        z_rx[:, i, :] *= scale[i]

    tx_idx = np.asarray(cut["tx_idx"], dtype=np.float64)
    beam_id = np.asarray(cut["beam_id"], dtype=np.int64)
    delays = np.asarray(cut["delays"], dtype=np.int64)
    rds = ResponseRangeDopplerSearch(
        txlen=z_tx.shape[1],
        echolen=z_rx.shape[2],
        n_channels=z_rx.shape[1],
        interp=interp,
    )

    coarse = []
    ranges = []
    snr = []
    response = []
    drg_km = c.c / 1e6 / 2.0 / 1e3
    for pulse in range(z_rx.shape[0]):
        out = rds.mf(z_tx[pulse, :], z_rx[pulse, :, :])
        coarse.append(out["doppler_mps"])
        ranges.append((delays[pulse] + out["range_gate"] / interp) * drg_km)
        snr.append((out["pprof"][out["range_gate"]] - out["noise_floor"]) / out["noise_floor"])
        response.append(out["response"])

    coarse = np.asarray(coarse, dtype=np.float64)
    snr = np.asarray(snr, dtype=np.float64)
    response = np.asarray(response, dtype=np.complex64)
    zpp, phase_dt, coarse_phase, pair_snr = pulse_pair_phase_samples(
        tx_idx,
        beam_id,
        coarse,
        response,
        pc.wavelength,
        snr,
    )
    fit_phase, fit_residual, fit_phase_slope, fit_phase_intercept = fit_wrapped_phase_time(
        tx_idx,
        beam_id,
        zpp,
        pair_snr,
    )
    fit_doppler = phase_fit_doppler_mps(fit_phase, phase_dt, pc.wavelength)

    return {
        "tx_idx": tx_idx,
        "beam_id": beam_id,
        "range_km": np.asarray(ranges, dtype=np.float64),
        "snr": snr,
        "pair_snr": pair_snr,
        "coarse_doppler_mps": coarse,
        "zpp": zpp,
        "pulse_pair_phase_rad": np.angle(zpp),
        "phase_dt_s": phase_dt,
        "coarse_phase_rad": coarse_phase,
        "phase_fit_rad": fit_phase,
        "phase_fit_residual_rad": fit_residual,
        "phase_fit_slope_rad_s": fit_phase_slope,
        "phase_fit_intercept_rad": fit_phase_intercept,
        "phase_fit_doppler_mps": fit_doppler,
        "diagnostics": pulse_phase_diagnostics(np.angle(zpp), fit_residual),
    }


def write_hdf5_summary(result, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        for key, value in result.items():
            if key == "diagnostics":
                group = handle.create_group("diagnostics")
                for diag_key, diag_value in value.items():
                    group.attrs[diag_key] = diag_value
            else:
                handle.create_dataset(key, data=value)


def plot_case(result, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    t_rel = result["tx_idx"] / 1e6 - result["tx_idx"][0] / 1e6
    good = np.isfinite(result["phase_fit_doppler_mps"])
    snr_db = 10.0 * np.log10(np.maximum(result["snr"], 1e-6))

    fig, axes = plt.subplots(3, 1, figsize=(8.5, 8.2), sharex=True, constrained_layout=True)
    ax = axes[0]
    sc = ax.scatter(
        t_rel,
        result["coarse_doppler_mps"] / 1e3,
        c=snr_db,
        s=22,
        cmap="viridis",
        label="within-pulse matched filter",
    )
    ax.scatter(
        t_rel[good],
        result["phase_fit_doppler_mps"][good] / 1e3,
        marker="x",
        s=28,
        color="tab:red",
        label="wrapped phase fit",
    )
    ax.set_ylabel("Doppler velocity (km/s)")
    ax.legend(fontsize=8)
    fig.colorbar(sc, ax=ax, label="SNR (dB)")

    ax = axes[1]
    delta = result["phase_fit_doppler_mps"] - result["coarse_doppler_mps"]
    ax.axhline(0.0, color="0.25", lw=0.8)
    ax.plot(t_rel[good], delta[good], ".", color="tab:purple")
    ax.set_ylabel("Fit - coarse (m/s)")

    ax = axes[2]
    for bid in np.unique(result["beam_id"]):
        use = result["beam_id"] == bid
        ax.plot(t_rel[use], result["range_km"][use], ".", label=f"beam {int(bid)}")
    ax.set_xlabel("Time from first cut pulse (s)")
    ax.set_ylabel("Range (km)")
    ax.legend(fontsize=7, ncol=5)

    fig.suptitle("Pulse-to-pulse Doppler test")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_phase_progression(result, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    t_rel = result["tx_idx"] / 1e6 - result["tx_idx"][0] / 1e6
    phase_cycles = result["pulse_pair_phase_rad"] / (2.0 * np.pi)
    fit_cycles = result["phase_fit_rad"] / (2.0 * np.pi)
    residual_cycles = result["phase_fit_residual_rad"] / (2.0 * np.pi)
    snr_db = 10.0 * np.log10(np.maximum(result["snr"], 1e-6))

    fig, axes = plt.subplots(3, 1, figsize=(8.5, 8.2), sharex=True, constrained_layout=True)
    ax = axes[0]
    fit_label_added = False
    for bid in np.unique(result["beam_id"]):
        use = (result["beam_id"] == bid) & np.isfinite(phase_cycles)
        if not np.any(use):
            continue
        sc = ax.scatter(t_rel[use], phase_cycles[use], s=24, label=f"beam {int(bid)}")
        fit_use = use & np.isfinite(fit_cycles)
        color = sc.get_facecolors()[0] if len(sc.get_facecolors()) else None
        ax.plot(
            t_rel[fit_use],
            fit_cycles[fit_use],
            lw=1.5,
            color=color,
            label="SNR-weighted circular fit" if not fit_label_added else None,
        )
        fit_label_added = True
    ax.set_ylabel("Pulse-pair phase (cycles)")
    ax.set_ylim(-0.55, 0.55)
    ax.legend(fontsize=7, ncol=5)

    ax = axes[1]
    good = np.isfinite(residual_cycles)
    sc = ax.scatter(t_rel[good], residual_cycles[good], c=snr_db[good], s=24, cmap="viridis")
    ax.axhline(0.0, color="0.25", lw=0.8)
    ax.set_ylabel("Circular residual (cycles)")
    fig.colorbar(sc, ax=ax, label="SNR (dB)")

    ax = axes[2]
    good_fit = np.isfinite(result["phase_fit_doppler_mps"])
    ax.plot(t_rel, result["coarse_doppler_mps"] / 1e3, ".", color="0.5", label="within-pulse")
    ax.plot(t_rel[good_fit], result["phase_fit_doppler_mps"][good_fit] / 1e3, ".", color="tab:red", label="phase fit")
    ax.set_xlabel("Time from first cut pulse (s)")
    ax.set_ylabel("Aliased phase Doppler (km/s)")
    ax.legend(fontsize=8)

    fig.suptitle("Wrapped pulse-to-pulse phase fit")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def write_latex_memo(sample_idx: int, result, figure_path: Path, phase_figure_path: Path, h5_path: Path, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    diag = result["diagnostics"]
    fig_rel = figure_path.name
    phase_fig_rel = phase_figure_path.name
    h5_rel = h5_path.name
    text = rf"""\documentclass[11pt]{{article}}
\usepackage[margin=1in]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{siunitx}}
\usepackage{{hyperref}}

\newif\ifshowscriptprovenance
\showscriptprovenancetrue
\newcommand{{\scriptprovenance}}[1]{{\ifshowscriptprovenance\par\small\emph{{Script: #1}}\fi}}

\begin{{document}}
\section*{{Pulse-to-pulse Doppler phase test}}

Cut sample {sample_idx} begins at {stuffr.unix2datestr(sample_idx / 1e6)} UTC.  The
test compares the ordinary within-pulse matched-filter Doppler estimate with a
wrapped pulse-to-pulse phase fit. Same-beam pulses are paired so fixed beam
phase offsets cancel. The fit is performed directly on wrapped complex phase
samples using the circular residual convention from \texttt{{msis\_trajectory.py}}.

\begin{{figure}}[h]
  \centering
  \includegraphics[width=\linewidth]{{{fig_rel}}}
  \caption{{Doppler comparison for one meteor cut.
  \scriptprovenance{{pulse_phase_doppler_test_case.py --sample-idx {sample_idx}}}}}
\end{{figure}}

\begin{{figure}}[h]
  \centering
  \includegraphics[width=\linewidth]{{{phase_fig_rel}}}
  \caption{{Same-beam pulse-to-pulse phase samples and an SNR-weighted circular
  phase fit. The residual is computed as
  \texttt{{angle(exp(1j*angle(zpp))*exp(-1j*angle(model\_zpp)))}}.
  \scriptprovenance{{pulse_phase_doppler_test_case.py --sample-idx {sample_idx}}}}}
\end{{figure}}

\begin{{tabular}}{{lr}}
\hline
Pulses with phase sample & {diag["n_phase"]} \\
Median circular residual & \SI{{{diag["phase_residual_median_rad"]:.3f}}}{{rad}} \\
MAD circular residual & \SI{{{diag["phase_residual_mad_rad"]:.3f}}}{{rad}} \\
RMS circular residual & \SI{{{diag["phase_residual_rms_rad"]:.3f}}}{{rad}} \\
\hline
\end{{tabular}}

Per-pulse values are written to \texttt{{{h5_rel}}}.

\end{{document}}
"""
    output.write_text(text)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cut-dir", type=Path, default=Path("data/metadata/cut"))
    parser.add_argument("--sample-idx", type=int)
    parser.add_argument("--interp", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=Path("memos/pulse_phase_doppler_test"))
    args = parser.parse_args()

    sample_idx = args.sample_idx if args.sample_idx is not None else discover_sample(args.cut_dir)
    cut = load_cut(args.cut_dir, sample_idx)
    result = compute_pulse_phase_case(cut, interp=args.interp)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure_path = args.output_dir / f"pulse_phase_doppler_{sample_idx}.png"
    phase_figure_path = args.output_dir / f"pulse_phase_progression_{sample_idx}.png"
    h5_path = args.output_dir / f"pulse_phase_doppler_{sample_idx}.h5"
    tex_path = args.output_dir / f"pulse_phase_doppler_{sample_idx}.tex"

    plot_case(result, figure_path)
    plot_phase_progression(result, phase_figure_path)
    write_hdf5_summary(result, h5_path)
    write_latex_memo(sample_idx, result, figure_path, phase_figure_path, h5_path, tex_path)

    diag = result["diagnostics"]
    print(f"sample_idx {sample_idx} {stuffr.unix2datestr(sample_idx / 1e6)}")
    print(f"wrote {figure_path}")
    print(f"wrote {phase_figure_path}")
    print(f"wrote {h5_path}")
    print(f"wrote {tex_path}")
    print(
        "phase fit: n={n_phase} median_residual={phase_residual_median_rad:.3f} rad "
        "mad_residual={phase_residual_mad_rad:.3f} rad rms_residual={phase_residual_rms_rad:.3f} rad".format(**diag)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
