#!/usr/bin/env python3
"""Interactively tune the two Figure 7 HEALPix color-scale maxima."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.widgets import Button, Slider

from healpix_hammer import render_healpix_hammer
from plot_paper_radiant_results import (
    FIGURE7_SPEED_VMAX,
    FIGURE7_EXPOSURE_CONTOUR_HOURS,
    FIGURE7_ZENITH_VMAX,
    add_exposure_contours,
    add_source_markers,
    hide_apex_meridian,
    plot_all_radiants,
    style_hammer,
)


DEFAULT_SIDECAR = Path("figs/paper_radiant_results_current/paper_radiant_results.h5")
DEFAULT_OUTPUT = Path("figs/paper_radiant_results_current/paper_radiant_distribution_corrected_tuned.png")


def positive(values: np.ndarray) -> np.ndarray:
    return values[np.isfinite(values) & (values > 0.0)]


def slider_limits(values: np.ndarray, initial: float) -> tuple[float, float, float]:
    data = positive(values)
    low = float(np.nanpercentile(data, 90.0))
    high = float(np.nanmax(data))
    return np.log10(low), np.log10(high), np.log10(np.clip(initial, low, high))


def percentile(values: np.ndarray, value: float) -> float:
    data = np.sort(positive(values))
    return 100.0 * float(np.searchsorted(data, value, side="right")) / len(data)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sidecar", type=Path, default=DEFAULT_SIDECAR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--left-vmax", type=float, default=FIGURE7_ZENITH_VMAX)
    parser.add_argument("--right-vmax", type=float, default=FIGURE7_SPEED_VMAX)
    parser.add_argument("--save-only", action="store_true", help="Render the selected settings without opening the GUI")
    args = parser.parse_args()

    with h5py.File(args.sidecar, "r") as h5:
        raw_hist = np.asarray(h5["raw_count"], dtype=np.float64)
        zenith_rate = np.asarray(h5["healpix_zenith_only_rate_h_inv"], dtype=np.float64)
        speed_rate = np.asarray(h5["healpix_debiased_rate_h_inv"], dtype=np.float64)
        exposure_hours = np.asarray(h5["radiant_exposure_hours"], dtype=np.float64)
        xedges = np.asarray(h5["plot_longitude_edges_deg"], dtype=np.float64)
        yedges = np.asarray(h5["ecliptic_latitude_edges_deg"], dtype=np.float64)
        nside = int(h5.attrs["healpix_nside"])
        alpha = float(h5.attrs["fitted_zenith_exponent_alpha"])

    if args.save_only:
        plot_all_radiants(
            np.empty(0),
            raw_hist,
            zenith_rate,
            speed_rate,
            nside,
            exposure_hours,
            xedges,
            yedges,
            alpha,
            args.output,
            zenith_vmax=args.left_vmax,
            speed_vmax=args.right_vmax,
        )
        print(f"Saved {args.output} with left vmax={args.left_vmax:g}, right vmax={args.right_vmax:g}")
        return

    left_limits = slider_limits(zenith_rate, args.left_vmax)
    right_limits = slider_limits(speed_rate, args.right_vmax)
    left_vmin = float(np.nanpercentile(positive(zenith_rate), 5.0))
    right_vmin = float(np.nanpercentile(positive(speed_rate), 5.0))

    fig = plt.figure(figsize=(12.0, 6.8))
    fig.subplots_adjust(left=0.06, right=0.98, top=0.94, bottom=0.25, wspace=0.15)
    ax0 = fig.add_subplot(121, projection="hammer")
    ax1 = fig.add_subplot(122, projection="hammer")
    norm0 = LogNorm(vmin=left_vmin, vmax=10.0 ** left_limits[2])
    norm1 = LogNorm(vmin=right_vmin, vmax=10.0 ** right_limits[2])
    mesh0 = render_healpix_hammer(ax0, zenith_rate, nside, cmap="magma", norm=norm0)
    mesh1 = render_healpix_hammer(ax1, speed_rate, nside, cmap="magma", norm=norm1)

    for ax in (ax0, ax1):
        style_hammer(ax)
        hide_apex_meridian(ax)
        ax.set_xticklabels([])
        add_source_markers(ax)
        ax.set_xlabel(r"Sun-centered ecliptic longitude, $\lambda-\lambda_\odot$")
        ax.set_ylabel(r"Ecliptic latitude, $\beta$")
    ax0.set_title(rf"Exposure + zenith corrected rate ($\alpha={alpha:.2f}$)", fontsize=10)
    ax1.set_title(r"With additional $1/v_g^3$ speed weight", fontsize=10)
    add_exposure_contours(
        ax0,
        exposure_hours,
        raw_hist,
        xedges,
        yedges,
        contour_levels=FIGURE7_EXPOSURE_CONTOUR_HOURS,
    )
    cb0 = fig.colorbar(mesh0, ax=ax0, orientation="horizontal", pad=0.10, fraction=0.045)
    cb1 = fig.colorbar(mesh1, ax=ax1, orientation="horizontal", pad=0.10, fraction=0.045)
    cb0.set_label(r"Radiant rate (h$^{-1}$ per HEALPix pixel)")
    cb1.set_label(r"Speed-weighted radiant rate (h$^{-1}$ per HEALPix pixel)")

    slider0 = Slider(
        fig.add_axes([0.11, 0.135, 0.34, 0.032]),
        "Left log10(vmax)",
        left_limits[0],
        left_limits[1],
        valinit=left_limits[2],
    )
    slider1 = Slider(
        fig.add_axes([0.58, 0.135, 0.34, 0.032]),
        "Right log10(vmax)",
        right_limits[0],
        right_limits[1],
        valinit=right_limits[2],
    )
    readout = fig.text(0.5, 0.085, "", ha="center", va="center", fontsize=10)
    status = fig.text(0.5, 0.025, f"Save target: {args.output}", ha="center", va="center", fontsize=9)

    def selected() -> tuple[float, float]:
        return 10.0 ** slider0.val, 10.0 ** slider1.val

    def update(_value=None) -> None:
        left_vmax, right_vmax = selected()
        norm0.vmax = left_vmax
        norm1.vmax = right_vmax
        mesh0.changed()
        mesh1.changed()
        cb0.update_normal(mesh0)
        cb1.update_normal(mesh1)
        readout.set_text(
            f"Left vmax={left_vmax:.4g} ({percentile(zenith_rate, left_vmax):.3f} percentile)    "
            f"Right vmax={right_vmax:.4g} ({percentile(speed_rate, right_vmax):.3f} percentile)"
        )
        fig.canvas.draw_idle()

    def save(_event) -> None:
        left_vmax, right_vmax = selected()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        plot_all_radiants(
            np.empty(0),
            raw_hist,
            zenith_rate,
            speed_rate,
            nside,
            exposure_hours,
            xedges,
            yedges,
            alpha,
            args.output,
            zenith_vmax=left_vmax,
            speed_vmax=right_vmax,
        )
        status.set_text(
            f"Saved {args.output} | --figure7-zenith-vmax {left_vmax:.8g} "
            f"--figure7-speed-vmax {right_vmax:.8g}"
        )
        fig.canvas.draw_idle()
        print(status.get_text(), flush=True)

    slider0.on_changed(update)
    slider1.on_changed(update)
    save_button = Button(fig.add_axes([0.42, 0.045, 0.16, 0.04]), "Save paper-style PNG")
    save_button.on_clicked(save)
    update()
    plt.show()


if __name__ == "__main__":
    main()
