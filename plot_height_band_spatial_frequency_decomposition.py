#!/usr/bin/env python3
"""Append spherical-harmonic radiant components below the height-band figure."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

from healpix_hammer import centered_plot_longitude_deg, render_healpix_hammer
from radiant_spatial_frequency_filter import (
    iterative_positive_spherical_harmonic_split,
    positive_log_norm,
    required_lmax,
)


SIDECAR = Path("figs/paper_radiant_results_current/paper_radiant_results.h5")
TOP_IMAGE = Path("/Users/jvi019/src/pansy_paper/paper_height_band_radiants_top.png")
OUTPUT = Path("/Users/jvi019/src/pansy_paper/paper_height_band_radiants.png")
FILTER_OUTPUT = Path("figs/paper_radiant_results_current/radiant_spherical_harmonic_split.h5")
PLOT_CENTER_LONGITUDE_DEG = -90.0


def add_source_markers(ax) -> None:
    for lon_deg, marker, color, size in (
        (0.0, "o", "#ffd21f", 22),
        (270.0, r"$\otimes$", "black", 45),
        (180.0, "o", "black", 22),
    ):
        ax.scatter(
            np.deg2rad(centered_plot_longitude_deg(lon_deg, PLOT_CENTER_LONGITUDE_DEG)),
            0.0,
            marker=marker,
            s=size,
            color=color,
            edgecolor="black" if marker == "o" else None,
            linewidth=0.3 if marker == "o" else 0.0,
            zorder=12,
        )


def style_hammer(ax):
    tick_pos = np.asarray([-90.0, 0.0, 90.0])
    tick_labels = [r"$0^\circ$", r"$270^\circ$", r"$180^\circ$"]
    ax.set_xticks(np.deg2rad(tick_pos))
    ax.set_xticklabels(tick_labels)
    ax.set_yticks(np.deg2rad([-60.0, -30.0, 0.0, 30.0, 60.0]))
    ax.set_yticklabels([r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$30^\circ$", r"$60^\circ$"])
    ax.grid(True, alpha=0.25, lw=0.45)
    ax.set_xlabel(r"Sun-centered ecliptic longitude, $\lambda_g-\lambda_\odot$")
    ax.set_ylabel(r"Ecliptic latitude, $\beta_g$")


def plot_healpix(ax, values, nside, *, cmap, norm, positive_only=True):
    mesh = render_healpix_hammer(
        ax,
        values,
        nside,
        cmap=cmap,
        norm=norm,
        center_longitude_deg=PLOT_CENTER_LONGITUDE_DEG,
        positive_only=positive_only,
    )
    style_hammer(ax)
    add_source_markers(ax)
    return mesh


def save_filter_product(path: Path, args, raw, lowpass, highpass, weights, history, nside, lmax) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h5:
        h5.attrs["schema"] = "radiant_spherical_harmonic_split_v1"
        h5.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        h5.attrs["source_h5"] = str(args.sidecar)
        h5.attrs["healpix_nside"] = np.int32(nside)
        h5.attrs["healpix_ordering"] = "RING"
        h5.attrs["lmax"] = np.int32(lmax)
        h5.attrs["zonal_frequency"] = "absolute spherical-harmonic order |m|"
        h5.attrs["meridional_frequency"] = "latitudinal degree l-|m|"
        for name in ("zonal_pass", "zonal_taper", "meridional_pass", "meridional_taper"):
            h5.attrs[name] = np.float32(getattr(args, name))
        h5.attrs["positive_iterations"] = np.int32(args.positive_iterations)
        h5.attrs["positive_fraction"] = np.float32(args.positive_fraction)
        h5.attrs["high_frequency_definition"] = "positive part of raw count minus iteratively refitted low-frequency background"
        h5.create_dataset("raw_count", data=np.asarray(raw, dtype=np.float32), compression="gzip", shuffle=True)
        h5.create_dataset("low_frequency_count", data=np.asarray(lowpass, dtype=np.float32), compression="gzip", shuffle=True)
        h5.create_dataset("positive_high_frequency_count", data=np.asarray(highpass, dtype=np.float32), compression="gzip", shuffle=True)
        h5.create_dataset("alm_lowpass_weight", data=np.asarray(weights, dtype=np.float32), compression="gzip", shuffle=True)
        history_dataset = h5.create_dataset(
            "iteration_diagnostics",
            data=np.asarray(history, dtype=np.float32),
            compression="gzip",
            shuffle=True,
        )
        history_dataset.attrs["columns"] = "positive_residual_sum,extracted_sum,positive_residual_max"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sidecar", type=Path, default=SIDECAR)
    parser.add_argument("--top-image", type=Path, default=TOP_IMAGE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--filter-output", type=Path, default=FILTER_OUTPUT)
    parser.add_argument("--zonal-pass", type=float, default=1.0)
    parser.add_argument("--zonal-taper", type=float, default=19.1)
    parser.add_argument("--meridional-pass", type=float, default=6.4)
    parser.add_argument("--meridional-taper", type=float, default=23.9)
    parser.add_argument("--positive-iterations", type=int, default=12)
    parser.add_argument("--positive-fraction", type=float, default=0.5)
    parser.add_argument("--lmax", type=int)
    args = parser.parse_args()

    top_image = mpimg.imread(args.top_image)
    with h5py.File(args.sidecar, "r") as h5:
        raw = np.asarray(h5["healpix_raw_count"], dtype=np.float64)
        nside = int(h5.attrs["healpix_nside"])
    lmax = int(args.lmax) if args.lmax is not None else required_lmax(
        nside,
        args.zonal_pass,
        args.zonal_taper,
        args.meridional_pass,
        args.meridional_taper,
    )
    lowpass, highpass, weights, history = iterative_positive_spherical_harmonic_split(
        raw,
        nside,
        args.zonal_pass,
        args.zonal_taper,
        args.meridional_pass,
        args.meridional_taper,
        lmax=lmax,
        iterations=args.positive_iterations,
        positive_fraction=args.positive_fraction,
    )
    save_filter_product(args.filter_output, args, raw, lowpass, highpass, weights, history, nside, lmax)

    fig = plt.figure(figsize=(12.0, 8.4), constrained_layout=False)
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.08, top=0.98, wspace=0.13, hspace=0.18)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.54, 1.0])
    ax_top = fig.add_subplot(gs[0, :])
    ax_low = fig.add_subplot(gs[1, 0], projection="hammer")
    ax_high = fig.add_subplot(gs[1, 1], projection="hammer")

    ax_top.imshow(top_image)
    ax_top.set_axis_off()

    mesh_low = plot_healpix(
        ax_low,
        np.where(lowpass > 0.0, lowpass, np.nan),
        nside,
        cmap="magma",
        norm=positive_log_norm(lowpass),
    )
    ax_low.set_title("Low spherical-harmonic frequencies", fontsize=11)
    cb_low = fig.colorbar(mesh_low, ax=ax_low, orientation="horizontal", pad=0.08, fraction=0.05)
    cb_low.set_label("Model count per HEALPix pixel")

    mesh_high = plot_healpix(
        ax_high,
        highpass,
        nside,
        cmap="magma",
        norm=positive_log_norm(highpass),
    )
    ax_high.set_title("Positive high-frequency excess", fontsize=11)
    cb_high = fig.colorbar(mesh_high, ax=ax_high, orientation="horizontal", pad=0.08, fraction=0.05)
    cb_high.set_label("Positive excess per HEALPix pixel")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=240)
    plt.close(fig)
    print(
        f"{args.output}\n{args.filter_output}\n"
        f"zonal |m|={args.zonal_pass:.2f}+{args.zonal_taper:.2f}; "
        f"meridional l-|m|={args.meridional_pass:.2f}+{args.meridional_taper:.2f}; "
        f"positive extraction={args.positive_fraction:.2f} x {args.positive_iterations}; lmax={lmax}"
    )


if __name__ == "__main__":
    main()
