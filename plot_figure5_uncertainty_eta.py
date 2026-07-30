#!/usr/bin/env python3
"""Build paper Figure 5 from uncertainty and eta-Aquariids HDF5 products."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
import numpy as np


ETA_SC_LON_DEG = 293.3
ETA_BETA_DEG = 7.9


def histogram_density(values: np.ndarray, bins: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = values[np.isfinite(values) & (values > 0.0)]
    counts, edges = np.histogram(values, bins=bins)
    widths = np.diff(np.log10(edges))
    density = counts / np.maximum(np.sum(counts) * widths, np.finfo(float).tiny)
    return density, edges


def load_uncertainty(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        angle = np.asarray(h5["initial_state_radiant_angle_sigma_deg"], dtype=np.float64)
        frac = np.asarray(h5["fractional_initial_state_velocity_sigma"], dtype=np.float64)
    return angle, frac


def load_eta_histogram(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5:
        xedges = np.asarray(h5["x_edges_sun_centered_lon_deg"], dtype=np.float64)
        yedges = np.asarray(h5["beta_edges_deg"], dtype=np.float64)
        count = np.asarray(h5["count"], dtype=np.float64)
    return xedges, yedges, count


def plot_figure(
    uncertainty_h5: Path,
    eta_h5: Path,
    output_png: Path,
    output_pdf: Path | None,
) -> None:
    angle, frac = load_uncertainty(uncertainty_h5)
    xedges, yedges, count = load_eta_histogram(eta_h5)

    with mpl.rc_context(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
        }
    ):
        fig = plt.figure(figsize=(4.8, 4.45), constrained_layout=False)
        gs = fig.add_gridspec(
            2,
            2,
            width_ratios=(1.0, 0.045),
            height_ratios=(1.0, 0.78),
            left=0.15,
            right=0.89,
            bottom=0.12,
            top=0.88,
            hspace=0.34,
            wspace=0.10,
        )

        ax0 = fig.add_subplot(gs[0, :])
        bins = np.geomspace(0.15, 60.0, 70)
        density, edges = histogram_density(angle, bins)
        peak_i = int(np.nanargmax(density))
        peak_angle = float(np.sqrt(edges[peak_i] * edges[peak_i + 1]))
        ax0.stairs(density, edges, color="#2f5f9f", linewidth=1.8)
        ax0.axvline(peak_angle, color="black", lw=1.4, ls="--")
        ax0.set_xscale("log")
        ax0.set_xlabel("Velocity-vector angular uncertainty (deg)")
        ax0.set_ylabel("Probability density per decade")
        ax0.grid(True, which="both", alpha=0.25)

        frac_good = frac[np.isfinite(frac) & (frac > 0.0)]
        top = ax0.twiny()
        frac_bins = np.geomspace(
            max(1e-5, np.nanpercentile(frac_good, 0.1)),
            max(2e-5, np.nanpercentile(frac_good, 99.8)),
            70,
        )
        frac_density, frac_edges = histogram_density(frac_good, frac_bins)
        frac_peak_i = int(np.nanargmax(frac_density))
        frac_peak = float(np.sqrt(frac_edges[frac_peak_i] * frac_edges[frac_peak_i + 1]))
        top.stairs(frac_density, frac_edges, color="#c45a1a", linewidth=1.5)
        top.set_xscale("log")
        top.set_xlabel(r"Fractional velocity uncertainty $\sigma_v/v_0$")
        top.tick_params(axis="x", colors="#9f4613")
        top.xaxis.label.set_color("#9f4613")
        ax0.text(
            0.97,
            0.95,
            f"N = {len(angle):,}\nangle peak = {peak_angle:.2f} deg\n$\\sigma_v/v_0$ peak = {frac_peak:.3f}",
            transform=ax0.transAxes,
            ha="right",
            va="top",
            fontsize=8,
        )

        ax1 = fig.add_subplot(gs[1, 0])
        cax = fig.add_subplot(gs[1, 1])
        positive = count[count > 0.0]
        vmax = max(1.0, float(np.nanpercentile(positive, 99.5))) if len(positive) else 1.0
        cmap = plt.get_cmap("magma").copy()
        cmap.set_bad("black")
        mesh = ax1.pcolormesh(
            xedges,
            yedges,
            np.ma.masked_where(count <= 0.0, count),
            cmap=cmap,
            norm=LogNorm(vmin=1.0, vmax=vmax),
            shading="flat",
        )
        ax1.scatter(ETA_SC_LON_DEG, ETA_BETA_DEG, marker="+", s=140, color="cyan", linewidth=1.5)
        ax1.set_aspect("equal", adjustable="box")
        ax1.set_xlim(xedges[-1], xedges[0])
        ax1.set_ylim(yedges[0], yedges[-1])
        ax1.set_xlabel(r"Sun-centered ecliptic longitude, $\lambda_g-\lambda_\odot$ (deg)")
        ax1.set_ylabel(r"Ecliptic latitude, $\beta_g$ (deg)")
        ax1.text(
            0.03,
            0.04,
            f"0.1 deg bins; N={int(np.nansum(count))}",
            transform=ax1.transAxes,
            ha="left",
            va="bottom",
            fontsize=8,
        )
        cbar = fig.colorbar(mesh, cax=cax)
        cbar.set_label(r"Raw count per 0.1$^\circ$ bin", fontsize=9)
        cbar.set_ticks([1, 3, 10, 30])
        cbar.ax.yaxis.set_major_formatter(LogFormatterMathtext())
        cbar.ax.tick_params(labelsize=8, pad=1)

        ax0.text(0.02, 0.92, "a)", transform=ax0.transAxes, fontsize=11, fontweight="bold")
        ax1.text(0.02, 0.88, "b)", transform=ax1.transAxes, fontsize=11, fontweight="bold", color="white")

        output_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_png, dpi=260)
        if output_pdf is not None:
            fig.savefig(output_pdf)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uncertainty-h5", type=Path, default=Path("figs/radiant_velocity_uncertainty_actual.h5"))
    parser.add_argument(
        "--eta-h5",
        type=Path,
        default=Path("figs/eta_aquariids_snapshot/eta_aquariids_radiant_full_activity_ls14_100_bin0.1deg.h5"),
    )
    parser.add_argument("--output-png", type=Path, default=Path("figs/paper_figure5_uncertainty_eta.png"))
    parser.add_argument("--output-pdf", type=Path, default=Path("figs/paper_figure5_uncertainty_eta.pdf"))
    args = parser.parse_args()
    plot_figure(args.uncertainty_h5, args.eta_h5, args.output_png, args.output_pdf)
    print(args.output_png)
    print(args.output_pdf)


if __name__ == "__main__":
    main()
