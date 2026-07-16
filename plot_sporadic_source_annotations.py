#!/usr/bin/env python3
"""Plot diagnostic sporadic-source apertures over the raw radiant-count map."""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


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


def load_regions():
    if REGION_JSON.exists():
        with REGION_JSON.open("r", encoding="utf-8") as fh:
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


def plot_grid(ax, hist, xedges, yedges, title, norm):
    plot_hist = np.asarray(hist, dtype=np.float64).copy()
    plot_hist[~np.isfinite(plot_hist) | (plot_hist <= 0.0)] = float(norm.vmin)
    mesh = ax.pcolormesh(np.deg2rad(xedges), np.deg2rad(yedges), plot_hist, shading="auto", cmap="magma", norm=norm)
    style_hammer(ax)
    ax.set_title(title, fontsize=10)
    return mesh


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
    ax.text(
        np.deg2rad(label_x),
        np.deg2rad(label_y),
        label,
        color=color,
        ha="center",
        va="center",
        fontsize=7.5,
        weight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.8},
        zorder=20,
    )


def main() -> None:
    with h5py.File(SIDECAR, "r") as h5:
        raw = np.asarray(h5["raw_count"], dtype=np.float64)
        xedges = np.asarray(h5["plot_longitude_edges_deg"], dtype=np.float64)
        yedges = np.asarray(h5["ecliptic_latitude_edges_deg"], dtype=np.float64)
        n_rows = len(h5["radiants"])
    regions = load_regions()

    raw_positive = raw[np.isfinite(raw) & (raw > 0.0)]
    raw_norm = LogNorm(vmin=1.0, vmax=max(2.0, float(np.nanmax(raw_positive))))

    fig = plt.figure(figsize=(7.2, 4.8), constrained_layout=True)
    ax0 = fig.add_subplot(111, projection="hammer")
    mesh0 = plot_grid(ax0, raw, xedges, yedges, f"Observed high-quality radiants (N={n_rows:,})", raw_norm)
    ax0.set_xticklabels([])
    add_source_markers(ax0)
    ax0.set_xlabel(r"Sun-centered ecliptic longitude, $\lambda-\lambda_\odot$")
    ax0.set_ylabel(r"Ecliptic latitude, $\beta$")
    for name, region in regions.items():
        label, color = REGION_COLORS.get(name, (name, "white"))
        add_region_box(ax0, label, region, color)
    cb0 = fig.colorbar(mesh0, ax=ax0, orientation="horizontal", pad=0.10, fraction=0.045)
    cb0.set_label("Count per bin")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=240)
    plt.close(fig)
    print(OUTPUT)


if __name__ == "__main__":
    main()
