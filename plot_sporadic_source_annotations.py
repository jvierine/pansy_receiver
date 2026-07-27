#!/usr/bin/env python3
"""Plot diagnostic sporadic-source apertures over the raw radiant-count map."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm, SymLogNorm

from healpix_hammer import _pixel_polygons, render_healpix_hammer

SIDECAR = Path("figs/paper_radiant_results_current/paper_radiant_results.h5")
REGION_JSON = Path("figs/sporadic_source_regions_manual.json")
OUTPUT = Path("/Users/j/src/pansy_paper/paper_radiant_distribution_sources.png")
PLOT_CENTER_LONGITUDE_DEG = -90.0


REGION_COLORS = {
    "helion": ("Helion", "#00a6d6"),
    "antihelion": ("Antihelion", "#e69f00"),
    "apex": ("Apex", "#009e73"),
    "narrow_apex": ("Narrow apex", "#d55e00"),
    "northern_toroidal": ("Northern toroidal", "#cc79a7"),
    "southern_toroidal": ("Southern toroidal", "#cc79a7"),
}


def wrap180(deg):
    return (np.asarray(deg, dtype=np.float64) + 180.0) % 360.0 - 180.0


def wrap360(deg):
    return np.asarray(deg, dtype=np.float64) % 360.0


def centered_plot_longitude_deg(lambda_minus_sun_deg):
    signed = wrap180(lambda_minus_sun_deg)
    return -wrap180(signed - PLOT_CENTER_LONGITUDE_DEG)


def load_regions(path: Path = REGION_JSON):
    if path.exists():
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)["regions"]
    return {
        "helion": {"center_lon_deg": 0.0, "half_lon_deg": 30.0, "beta_min_deg": -20.0, "beta_max_deg": 20.0},
        "antihelion": {"center_lon_deg": 180.0, "half_lon_deg": 30.0, "beta_min_deg": -20.0, "beta_max_deg": 20.0},
        "apex": {"center_lon_deg": 270.0, "half_lon_deg": 30.0, "beta_min_deg": -25.0, "beta_max_deg": 25.0},
        "narrow_apex": {"center_lon_deg": 270.0, "half_lon_deg": 15.0, "beta_min_deg": -10.0, "beta_max_deg": 10.0},
        "southern_toroidal": {"center_lon_deg": 270.0, "half_lon_deg": 45.0, "beta_min_deg": -65.0, "beta_max_deg": -35.0},
    }


def style_hammer(ax):
    tick_pos = np.asarray([-90.0, 0.0, 90.0])
    tick_labels = [f"{int(wrap360(PLOT_CENTER_LONGITUDE_DEG - tick))}$^\\circ$" for tick in tick_pos]
    ax.set_xticks(np.deg2rad(tick_pos))
    ax.set_xticklabels(tick_labels)
    ax.set_yticks(np.deg2rad([-60.0, -30.0, 0.0, 30.0, 60.0]))
    ax.set_yticklabels([r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$30^\circ$", r"$60^\circ$"])
    ax.grid(True, alpha=0.25, lw=0.45)


def plot_healpix(ax, values, nside, title, norm):
    mesh = render_healpix_hammer(ax, values, nside, cmap="magma", norm=norm)
    style_hammer(ax)
    ax.set_title(title, fontsize=10)
    return mesh


def render_signed_healpix_hammer(ax, values, nside, *, cmap="coolwarm", norm=None):
    values = np.asarray(values, dtype=np.float64)
    keep = np.isfinite(values)
    pixels = np.flatnonzero(keep)
    polygons, source_indices = _pixel_polygons(int(nside), pixels, PLOT_CENTER_LONGITUDE_DEG)
    collection = PolyCollection(
        polygons,
        array=values[pixels][source_indices],
        cmap=cmap,
        norm=norm,
        edgecolors="none",
        linewidths=0.0,
        antialiased=False,
        rasterized=True,
        zorder=1,
    )
    collection.set_clip_path(ax.patch)
    ax.add_collection(collection)
    ax.set_facecolor("white")
    return collection


def spherical_harmonic_split(raw_count: np.ndarray, lmax: int, lowpass_lmax: int) -> tuple[np.ndarray, np.ndarray]:
    """Return low-pass model and high-pass residual for a full-sky HEALPix map."""
    values = np.asarray(raw_count, dtype=np.float64)
    if values.ndim != 1:
        raise ValueError("raw_count must be a one-dimensional HEALPix map")
    if not np.all(np.isfinite(values)):
        values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    nside = hp.npix2nside(values.size)
    lmax = min(int(lmax), 3 * nside - 1)
    lowpass_lmax = min(int(lowpass_lmax), lmax)
    alm = hp.map2alm(values, lmax=lmax, iter=3, pol=False, use_weights=True)
    alm_low = alm.copy()
    for ell in range(lowpass_lmax + 1, lmax + 1):
        for emm in range(ell + 1):
            alm_low[hp.Alm.getidx(lmax, ell, emm)] = 0.0
    lowpass = hp.alm2map(alm_low, nside=nside, lmax=lmax, pol=False)
    highpass = values - lowpass
    return np.asarray(lowpass, dtype=np.float64), np.asarray(highpass, dtype=np.float64)


def positive_log_norm(values: np.ndarray, floor: float = 1.0) -> LogNorm:
    positive = np.asarray(values, dtype=np.float64)
    positive = positive[np.isfinite(positive) & (positive > 0.0)]
    vmax = max(float(floor) * 2.0, float(np.nanmax(positive))) if positive.size else float(floor) * 2.0
    return LogNorm(vmin=float(floor), vmax=vmax)


def signed_residual_norm(values: np.ndarray) -> SymLogNorm:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    scale = float(np.nanpercentile(np.abs(finite), 99.5)) if finite.size else 1.0
    scale = max(scale, 1.0)
    return SymLogNorm(linthresh=1.0, linscale=0.8, vmin=-scale, vmax=scale, base=10)


def add_panel_label(ax, label: str) -> None:
    ax.text(
        0.015,
        0.965,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        weight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5},
        zorder=30,
    )


def add_source_markers(ax) -> None:
    for lon_deg, marker, color, size in (
        (0.0, "o", "#ffd21f", 22),
        (-90.0, r"$\otimes$", "black", 45),
        (180.0, "o", "black", 22),
    ):
        ax.scatter(
            np.deg2rad(centered_plot_longitude_deg(lon_deg)),
            0.0,
            marker=marker,
            s=size,
            color=color,
            edgecolor="black" if marker == "o" else None,
            linewidth=0.3 if marker == "o" else 0.0,
            zorder=12,
        )


def add_region_box(ax, label, region, color):
    center_lon = region["center_lon_deg"]
    half_width = region["half_lon_deg"]
    beta_min = region["beta_min_deg"]
    beta_max = region["beta_max_deg"]
    lon1 = center_lon - half_width
    lon2 = center_lon + half_width
    lon_edge = np.linspace(lon1, lon2, 160)
    beta_edge = np.linspace(beta_min, beta_max, 80)
    x_top = centered_plot_longitude_deg(lon_edge)
    x_bottom = centered_plot_longitude_deg(lon_edge)
    x_left = centered_plot_longitude_deg(np.full_like(beta_edge, lon1))
    x_right = centered_plot_longitude_deg(np.full_like(beta_edge, lon2))
    ax.plot(np.deg2rad(x_top), np.deg2rad(np.full_like(x_top, beta_max)), color=color, lw=1.2, alpha=0.95)
    ax.plot(np.deg2rad(x_bottom), np.deg2rad(np.full_like(x_bottom, beta_min)), color=color, lw=1.2, alpha=0.95)
    ax.plot(np.deg2rad(x_left), np.deg2rad(beta_edge), color=color, lw=1.2, alpha=0.95)
    ax.plot(np.deg2rad(x_right), np.deg2rad(beta_edge), color=color, lw=1.2, alpha=0.95)
    label_x = centered_plot_longitude_deg(center_lon)
    if label == "Apex":
        label_y = beta_min - 5.0
    elif label == "Narrow apex":
        label_y = beta_max + 5.0
    else:
        label_y = beta_max + 4.0 if beta_max <= 25.0 else beta_max - 4.0
    if beta_max < 0.0:
        label_y = beta_min - 4.0
    font_size = 7.5
    y_offset_points = -0.5 * font_size if label == "Narrow apex" else 0.0
    ax.annotate(
        label,
        xy=(np.deg2rad(label_x), np.deg2rad(label_y)),
        xytext=(0.0, y_offset_points),
        textcoords="offset points",
        color=color,
        ha="center",
        va="center",
        fontsize=font_size,
        weight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.8},
        zorder=20,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sidecar", type=Path, default=SIDECAR)
    parser.add_argument("--regions", type=Path, default=REGION_JSON)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--harmonic-lmax", type=int, default=48, help="Maximum degree used in the spherical harmonic fit")
    parser.add_argument("--lowpass-lmax", type=int, default=8, help="Degree retained in the low-pass component")
    args = parser.parse_args()

    with h5py.File(args.sidecar, "r") as h5:
        raw = np.asarray(h5["healpix_raw_count"], dtype=np.float64)
        nside = int(h5.attrs["healpix_nside"])
        n_rows = len(h5["radiants"])
    regions = load_regions(args.regions)

    lowpass, highpass = spherical_harmonic_split(raw, args.harmonic_lmax, args.lowpass_lmax)
    positive_highpass = np.clip(highpass, 0.0, None)

    raw_norm = positive_log_norm(raw)
    lowpass_norm = positive_log_norm(np.clip(lowpass, 0.0, None))
    highpass_norm = signed_residual_norm(highpass)
    positive_norm = positive_log_norm(positive_highpass)

    fig = plt.figure(figsize=(10.8, 8.0), constrained_layout=False)
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.075, top=0.95, wspace=0.13, hspace=0.42)
    axes = [fig.add_subplot(2, 2, index, projection="hammer") for index in range(1, 5)]
    ax0, ax1, ax2, ax3 = axes

    mesh0 = plot_healpix(ax0, raw, nside, f"Raw counts (N={n_rows:,})", raw_norm)
    mesh1 = render_healpix_hammer(ax1, np.clip(lowpass, 0.0, None), nside, cmap="magma", norm=lowpass_norm)
    style_hammer(ax1)
    ax1.set_title(rf"Low-pass harmonic model ($\ell\leq{args.lowpass_lmax}$)", fontsize=10)
    mesh2 = render_signed_healpix_hammer(ax2, highpass, nside, cmap="coolwarm", norm=highpass_norm)
    style_hammer(ax2)
    ax2.set_title("High-pass residual", fontsize=10)
    mesh3 = render_healpix_hammer(ax3, positive_highpass, nside, cmap="magma", norm=positive_norm)
    style_hammer(ax3)
    ax3.set_title("Positive high-pass enhancements", fontsize=10)

    for label, ax in zip(("a", "b", "c", "d"), axes, strict=True):
        add_panel_label(ax, label)
        ax.set_xlabel(r"$\lambda_g-\lambda_\odot$")
        ax.set_ylabel(r"$\beta_g$")
    add_source_markers(ax0)
    add_source_markers(ax1)
    add_source_markers(ax2)
    add_source_markers(ax3)
    for name, region in regions.items():
        label, color = REGION_COLORS.get(name, (name, "white"))
        add_region_box(ax0, label, region, color)
    cb0 = fig.colorbar(mesh0, ax=ax0, orientation="horizontal", pad=0.08, fraction=0.045)
    cb0.set_label("Count per pixel")
    cb1 = fig.colorbar(mesh1, ax=ax1, orientation="horizontal", pad=0.08, fraction=0.045)
    cb1.set_label("Model count per pixel")
    cb2 = fig.colorbar(mesh2, ax=ax2, orientation="horizontal", pad=0.08, fraction=0.045)
    cb2.set_label("Residual count per pixel")
    cb3 = fig.colorbar(mesh3, ax=ax3, orientation="horizontal", pad=0.08, fraction=0.045)
    cb3.set_label("Positive residual count per pixel")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=240)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
