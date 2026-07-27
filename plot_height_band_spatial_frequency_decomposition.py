#!/usr/bin/env python3
"""Render the height-band and spherical-harmonic radiant maps as a native 2x2 figure."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from healpix_hammer import centered_plot_longitude_deg, render_healpix_hammer
from radiant_spatial_frequency_filter import (
    iterative_positive_spherical_harmonic_split,
    positive_log_norm,
    required_lmax,
)


SIDECAR = Path("figs/paper_radiant_results_current/paper_radiant_results.h5")
SELECTION_H5 = Path("figs/paper_radiant_results_current/height_band_selection.h5")
OUTPUT = Path("/Users/jvi019/src/pansy_paper/paper_height_band_radiants.png")
FILTER_OUTPUT = Path("figs/paper_radiant_results_current/radiant_spherical_harmonic_split.h5")
PLOT_CENTER_LONGITUDE_DEG = -90.0
UPPER_BAND_COUNT = 165_074
LOWER_BAND_COUNT = 805_800


def add_source_markers(ax) -> None:
    ax.text(
        np.deg2rad(centered_plot_longitude_deg(0.0, PLOT_CENTER_LONGITUDE_DEG)),
        0.0,
        r"$\odot$",
        color="#ffd21f",
        fontsize=8,
        ha="center",
        va="center",
        zorder=12,
    )
    for lon_deg, marker, color, edgecolor, linewidth, size in (
        (270.0, r"$\otimes$", "black", "none", 0.0, 45),
        (180.0, "o", "black", "white", 0.8, 22),
    ):
        ax.scatter(
            np.deg2rad(centered_plot_longitude_deg(lon_deg, PLOT_CENTER_LONGITUDE_DEG)),
            0.0,
            marker=marker,
            s=size,
            color=color,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=12,
        )


def style_hammer(ax) -> None:
    ax.set_xticks(np.deg2rad([-90.0, 0.0, 90.0]))
    ax.set_xticklabels([])
    ax.set_yticks(np.deg2rad([-60.0, -30.0, 0.0, 30.0, 60.0]))
    ax.set_yticklabels(
        [r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$30^\circ$", r"$60^\circ$"]
    )
    ax.grid(True, alpha=0.25, lw=0.45)
    ax.set_xlabel(r"Sun-centered ecliptic longitude, $\lambda_g-\lambda_\odot$")
    ax.set_ylabel(r"Ecliptic latitude, $\beta_g$")


def plot_healpix(ax, values, nside, *, norm):
    mesh = render_healpix_hammer(
        ax,
        values,
        nside,
        cmap="magma",
        norm=norm,
        center_longitude_deg=PLOT_CENTER_LONGITUDE_DEG,
        positive_only=True,
    )
    style_hammer(ax)
    add_source_markers(ax)
    return mesh


def healpix_counts(rows: np.ndarray, mask: np.ndarray, nside: int) -> np.ndarray:
    lon = np.asarray(rows["lambda_minus_sun_deg"], dtype=np.float64)[mask]
    lat = np.asarray(rows["radiant_beta_ecliptic_deg"], dtype=np.float64)[mask]
    good = np.isfinite(lon) & np.isfinite(lat)
    pixel = hp.ang2pix(int(nside), lon[good], lat[good], lonlat=True)
    return np.bincount(pixel, minlength=hp.nside2npix(int(nside))).astype(np.float64)


def derive_height_band_maps(rows: np.ndarray, nside: int):
    height = np.asarray(rows["first_alt_km"], dtype=np.float64)
    speed = np.asarray(rows["speed_km_s"], dtype=np.float64)
    good = np.isfinite(height) & np.isfinite(speed) & (height > 0.0) & (speed > 0.0)
    good_index = np.flatnonzero(good)
    height = height[good]
    speed = speed[good]

    edges = np.arange(0.0, 82.0, 2.0)
    centers = 0.5 * (edges[:-1] + edges[1:])
    median_height = np.asarray(
        [
            np.median(height[(speed >= lower) & (speed < upper)])
            if np.count_nonzero((speed >= lower) & (speed < upper)) >= 100
            else np.nan
            for lower, upper in zip(edges[:-1], edges[1:], strict=True)
        ]
    )
    fit = np.isfinite(median_height)
    ridge_coefficients = np.polyfit(centers[fit], median_height[fit], 2)
    residual = height - np.polyval(ridge_coefficients, speed)

    if len(residual) < UPPER_BAND_COUNT + LOWER_BAND_COUNT:
        raise ValueError("not enough radiant rows to reconstruct both height bands")
    lower_local = np.argpartition(residual, LOWER_BAND_COUNT - 1)[:LOWER_BAND_COUNT]
    upper_local = np.argpartition(
        residual, len(residual) - UPPER_BAND_COUNT
    )[-UPPER_BAND_COUNT:]
    lower_mask = np.zeros(len(rows), dtype=bool)
    upper_mask = np.zeros(len(rows), dtype=bool)
    lower_mask[good_index[lower_local]] = True
    upper_mask[good_index[upper_local]] = True
    if np.any(lower_mask & upper_mask):
        raise RuntimeError("derived upper and lower height bands overlap")

    return (
        healpix_counts(rows, upper_mask, nside),
        healpix_counts(rows, lower_mask, nside),
        np.asarray(rows["sample_idx"][upper_mask], dtype=np.int64),
        np.asarray(rows["sample_idx"][lower_mask], dtype=np.int64),
        ridge_coefficients,
        float(np.min(residual[upper_local])),
        float(np.max(residual[lower_local])),
    )


def load_or_create_height_band_maps(path: Path, rows: np.ndarray, nside: int):
    if path.exists():
        with h5py.File(path, "r") as h5:
            if "upper_healpix_count" in h5 and "lower_healpix_count" in h5:
                return (
                    np.asarray(h5["upper_healpix_count"], dtype=np.float64),
                    np.asarray(h5["lower_healpix_count"], dtype=np.float64),
                )
            if "upper/radiants" in h5 and "lower/radiants" in h5:
                upper_rows = h5["upper/radiants"][()]
                lower_rows = h5["lower/radiants"][()]
                return (
                    healpix_counts(
                        upper_rows, np.ones(len(upper_rows), dtype=bool), nside
                    ),
                    healpix_counts(
                        lower_rows, np.ones(len(lower_rows), dtype=bool), nside
                    ),
                )

    (
        upper_count,
        lower_count,
        upper_sample_idx,
        lower_sample_idx,
        ridge_coefficients,
        upper_threshold,
        lower_threshold,
    ) = derive_height_band_maps(rows, nside)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h5:
        h5.attrs["schema"] = "height_band_healpix_selection_v1"
        h5.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        h5.attrs["healpix_nside"] = np.int32(nside)
        h5.attrs["selection"] = (
            "lower and upper residual tails around the quadratic conditional-median "
            "initial-height ridge, preserving the published band counts"
        )
        h5.attrs["upper_count"] = np.int32(len(upper_sample_idx))
        h5.attrs["lower_count"] = np.int32(len(lower_sample_idx))
        h5.attrs["upper_residual_threshold_km"] = np.float32(upper_threshold)
        h5.attrs["lower_residual_threshold_km"] = np.float32(lower_threshold)
        h5.create_dataset(
            "ridge_coefficients",
            data=np.asarray(ridge_coefficients, dtype=np.float32),
        )
        for name, values in (
            ("upper_sample_idx", upper_sample_idx),
            ("lower_sample_idx", lower_sample_idx),
            ("upper_healpix_count", np.asarray(upper_count, dtype=np.float32)),
            ("lower_healpix_count", np.asarray(lower_count, dtype=np.float32)),
        ):
            h5.create_dataset(name, data=values, compression="gzip", shuffle=True)
    return upper_count, lower_count


def shared_positive_log_norm(*values: np.ndarray) -> LogNorm:
    positive = []
    for value in values:
        value = np.asarray(value, dtype=np.float64)
        positive.append(value[np.isfinite(value) & (value > 0.0)])
    positive = np.concatenate(positive)
    return LogNorm(
        vmin=float(np.min(positive)),
        vmax=float(np.percentile(positive, 99.8)),
    )


def save_filter_product(
    path: Path, args, raw, lowpass, highpass, weights, history, nside, lmax
) -> None:
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
        for name in (
            "zonal_pass",
            "zonal_taper",
            "meridional_pass",
            "meridional_taper",
        ):
            h5.attrs[name] = np.float32(getattr(args, name))
        h5.attrs["positive_iterations"] = np.int32(args.positive_iterations)
        h5.attrs["positive_fraction"] = np.float32(args.positive_fraction)
        h5.attrs["high_frequency_definition"] = (
            "positive part of raw count minus iteratively refitted low-frequency background"
        )
        for name, values in (
            ("raw_count", raw),
            ("low_frequency_count", lowpass),
            ("positive_high_frequency_count", highpass),
            ("alm_lowpass_weight", weights),
        ):
            h5.create_dataset(
                name,
                data=np.asarray(values, dtype=np.float32),
                compression="gzip",
                shuffle=True,
            )
        history_dataset = h5.create_dataset(
            "iteration_diagnostics",
            data=np.asarray(history, dtype=np.float32),
            compression="gzip",
            shuffle=True,
        )
        history_dataset.attrs["columns"] = (
            "positive_residual_sum,extracted_sum,positive_residual_max"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sidecar", type=Path, default=SIDECAR)
    parser.add_argument("--selection-h5", type=Path, default=SELECTION_H5)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--filter-output", type=Path, default=FILTER_OUTPUT)
    parser.add_argument("--zonal-pass", type=float, default=6.1)
    parser.add_argument("--zonal-taper", type=float, default=10.5)
    parser.add_argument("--meridional-pass", type=float, default=32.3)
    parser.add_argument("--meridional-taper", type=float, default=7.4)
    parser.add_argument("--positive-iterations", type=int, default=12)
    parser.add_argument("--positive-fraction", type=float, default=0.5)
    parser.add_argument("--lmax", type=int)
    parser.add_argument(
        "--reuse-filter-output",
        action="store_true",
        help="Reuse the low/high-frequency maps in --filter-output instead of recomputing them",
    )
    args = parser.parse_args()

    with h5py.File(args.sidecar, "r") as h5:
        raw = np.asarray(h5["healpix_raw_count"], dtype=np.float64)
        rows = h5["radiants"][()]
        nside = int(h5.attrs["healpix_nside"])
    upper_count, lower_count = load_or_create_height_band_maps(
        args.selection_h5, rows, nside
    )

    if args.reuse_filter_output:
        with h5py.File(args.filter_output, "r") as h5:
            lowpass = np.asarray(h5["low_frequency_count"], dtype=np.float64)
            highpass = np.asarray(
                h5["positive_high_frequency_count"], dtype=np.float64
            )
            weights = np.asarray(h5["alm_lowpass_weight"], dtype=np.float64)
            history = np.asarray(h5["iteration_diagnostics"], dtype=np.float64)
            lmax = int(h5.attrs["lmax"])
    else:
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
        save_filter_product(
            args.filter_output,
            args,
            raw,
            lowpass,
            highpass,
            weights,
            history,
            nside,
            lmax,
        )

    bin_area_deg2 = hp.nside2pixarea(nside, degrees=True)
    upper_density = upper_count / (np.sum(upper_count) * bin_area_deg2)
    lower_density = lower_count / (np.sum(lower_count) * bin_area_deg2)
    lowpass_density = lowpass / bin_area_deg2
    highpass_density = highpass / bin_area_deg2
    top_norm = shared_positive_log_norm(upper_density, lower_density)

    fig = plt.figure(figsize=(12.0, 8.8), constrained_layout=False)
    fig.subplots_adjust(
        left=0.065,
        right=0.985,
        bottom=0.075,
        top=0.965,
        wspace=0.15,
        hspace=0.38,
    )
    axes = [
        fig.add_subplot(2, 2, index + 1, projection="hammer")
        for index in range(4)
    ]
    panels = (
        (
            axes[0],
            upper_density,
            top_norm,
            f"Upper initial-height band (N={int(np.sum(upper_count)):,})",
            r"Normalized count density (deg$^{-2}$)",
        ),
        (
            axes[1],
            lower_density,
            top_norm,
            f"Lower initial-height band (N={int(np.sum(lower_count)):,})",
            r"Normalized count density (deg$^{-2}$)",
        ),
        (
            axes[2],
            np.where(lowpass_density > 0.0, lowpass_density, np.nan),
            positive_log_norm(lowpass_density),
            "Low spherical-harmonic frequencies",
            r"Model count density (deg$^{-2}$)",
        ),
        (
            axes[3],
            highpass_density,
            positive_log_norm(highpass_density),
            "Positive high-frequency excess",
            r"Positive excess density (deg$^{-2}$)",
        ),
    )
    for panel_label, (ax, values, norm, title, colorbar_label) in zip(
        (r"$\mathbf{a)}$", r"$\mathbf{b)}$", r"$\mathbf{c)}$", r"$\mathbf{d)}$"),
        panels,
        strict=True,
    ):
        mesh = plot_healpix(ax, values, nside, norm=norm)
        ax.set_title(title, fontsize=11)
        ax.text(
            0.105,
            0.79,
            panel_label,
            transform=ax.transAxes,
            color="white",
            fontsize=13,
            ha="left",
            va="top",
            zorder=20,
        )
        colorbar = fig.colorbar(
            mesh,
            ax=ax,
            orientation="horizontal",
            pad=0.12,
            fraction=0.055,
        )
        colorbar.set_label(colorbar_label)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=240)
    plt.close(fig)
    print(
        f"{args.output}\n{args.filter_output}\n{args.selection_h5}\n"
        f"zonal |m|={args.zonal_pass:.2f}+{args.zonal_taper:.2f}; "
        f"meridional l-|m|={args.meridional_pass:.2f}+{args.meridional_taper:.2f}; "
        f"positive extraction={args.positive_fraction:.2f} x "
        f"{args.positive_iterations}; lmax={lmax}"
    )


if __name__ == "__main__":
    main()
