#!/usr/bin/env python3
"""Animate sun-centered radiant scatter and histogram distributions with IAU shower overlays."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import imageio.v2 as imageio
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

from plot_fitted_radiant_distribution import (
    add_source_markers,
    centered_tick_labels,
    scaled_scatter_area,
)
from shower_radiant_overlay import active_showers, add_shower_overlay_hammer, circular_mean_deg


def read_radiants(path: Path):
    with h5py.File(path, "r") as h:
        arr = h["radiants"][()]
    good = np.isfinite(arr["plot_longitude_deg"])
    good &= np.isfinite(arr["radiant_beta_ecliptic_deg"])
    good &= np.isfinite(arr["speed_km_s"])
    good &= np.isfinite(arr["epoch_unix"])
    return arr[good]


def day_labels(arr):
    days = np.asarray(
        [
            dt.datetime.fromtimestamp(float(t), tz=dt.timezone.utc).strftime("%Y-%m-%d")
            for t in arr["epoch_unix"]
        ]
    )
    return np.unique(days), days


def parse_day(day: str) -> dt.date:
    return dt.datetime.strptime(day, "%Y-%m-%d").date()


def window_label(day: str, window_days: int) -> str:
    if window_days <= 1:
        return day
    half_before = (window_days - 1) // 2
    half_after = window_days - 1 - half_before
    center = parse_day(day)
    start = center - dt.timedelta(days=half_before)
    end = center + dt.timedelta(days=half_after)
    return f"{day}  window {start.isoformat()}..{end.isoformat()}"


def frame_subset(arr, days, day, cumulative: bool, window_days: int):
    if cumulative:
        return arr[days <= day]
    window_days = max(1, int(window_days))
    half_before = (window_days - 1) // 2
    half_after = window_days - 1 - half_before
    center = parse_day(day)
    start = (center - dt.timedelta(days=half_before)).isoformat()
    end = (center + dt.timedelta(days=half_after)).isoformat()
    return arr[(days >= start) & (days <= end)]


def frame_showers(arr, day_rows, shower_catalog: Path | None, peak_tolerance_deg: float):
    rows = day_rows if len(day_rows) else arr
    showers, query = active_showers(
        rows,
        shower_catalog=shower_catalog,
        peak_tolerance_deg=peak_tolerance_deg,
    )
    return showers, float(query)


def setup_hammer():
    fig = plt.figure(figsize=(9.4, 6.2), constrained_layout=True)
    ax = fig.add_subplot(111, projection="hammer")
    tick_pos, tick_labels = centered_tick_labels()
    ax.set_xticks(np.deg2rad(tick_pos))
    ax.set_xticklabels(tick_labels)
    ax.grid(True, alpha=0.42)
    ax.set_xlabel(r"Sun-centered ecliptic longitude, $\lambda-\lambda_\odot$ (apex centered at $-90^\circ$)", labelpad=18)
    ax.set_ylabel(r"Ecliptic latitude, $\beta$")
    return fig, ax


def plot_scatter_frame(rows, day_rows, day_label: str, out: Path, shower_catalog: Path | None, peak_tolerance_deg: float):
    fig, ax = setup_hammer()
    sc = ax.scatter(
        np.deg2rad(rows["plot_longitude_deg"]),
        np.deg2rad(rows["radiant_beta_ecliptic_deg"]),
        c=rows["speed_km_s"],
        s=scaled_scatter_area(len(rows)),
        cmap="turbo",
        vmin=10.0,
        vmax=80.0,
        alpha=0.72,
        linewidths=0.15,
        edgecolors="white",
        zorder=3,
    )
    add_source_markers(ax)
    showers, query = frame_showers(rows, day_rows, shower_catalog, peak_tolerance_deg)
    add_shower_overlay_hammer(ax, showers)
    cb = fig.colorbar(sc, ax=ax, orientation="horizontal", pad=0.14, fraction=0.046)
    cb.set_label("DASST geocentric radiant speed (km/s)")
    ax.set_title(f"{day_label}  N={len(rows)}  shower solar longitude={query:.1f} deg", fontsize=10)
    ax.legend(loc="lower left", fontsize=8, frameon=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_histogram_frame(rows, day_rows, day_label: str, out: Path, shower_catalog: Path | None, peak_tolerance_deg: float, bins: int):
    fig, ax = setup_hammer()
    lon = np.asarray(rows["plot_longitude_deg"], dtype=np.float64)
    lat = np.asarray(rows["radiant_beta_ecliptic_deg"], dtype=np.float64)
    hist, xedges, yedges = np.histogram2d(lon, lat, bins=[bins * 2, bins], range=[[-180, 180], [-90, 90]])
    xx = np.deg2rad(xedges)
    yy = np.deg2rad(yedges)
    plot_hist = hist.T
    plot_count = plot_hist.copy()
    plot_count[~np.isfinite(plot_count) | (plot_count <= 0.0)] = 1.0
    norm = LogNorm(vmin=1, vmax=max(1, float(np.nanmax(plot_count))))
    ax.pcolormesh(xx, yy, plot_count, cmap="plasma", norm=norm, shading="auto")
    add_source_markers(ax)
    showers, query = frame_showers(rows, day_rows, shower_catalog, peak_tolerance_deg)
    add_shower_overlay_hammer(ax, showers, label_color="white")
    mappable = plt.cm.ScalarMappable(norm=norm, cmap="plasma")
    cb = fig.colorbar(mappable, ax=ax, orientation="horizontal", pad=0.14, fraction=0.046)
    cb.set_label("Radiants per bin")
    ax.set_title(f"{day_label}  N={len(rows)}  shower solar longitude={query:.1f} deg", fontsize=10)
    ax.legend(loc="lower left", fontsize=8, frameon=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def make_animation(arr, output_gif: Path, frame_dir: Path, kind: str, shower_catalog: Path | None, peak_tolerance_deg: float, cumulative: bool, window_days: int, bins: int):
    unique_days, days = day_labels(arr)
    frame_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    for day in unique_days:
        rows = frame_subset(arr, days, day, cumulative=cumulative, window_days=window_days)
        day_rows = arr[days == day]
        if len(rows) == 0:
            continue
        out = frame_dir / f"{kind}_{day}.png"
        label = day if cumulative else window_label(day, window_days)
        if kind == "scatter":
            plot_scatter_frame(rows, day_rows, label, out, shower_catalog, peak_tolerance_deg)
        else:
            plot_histogram_frame(rows, day_rows, label, out, shower_catalog, peak_tolerance_deg, bins)
        frames.append(imageio.imread(out))
    output_gif.parent.mkdir(parents=True, exist_ok=True)
    imageio.mimsave(output_gif, frames, duration=0.55, loop=0)
    return len(frames)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5", type=Path, default=Path("server_plots/current_fitted_radiant_distribution.h5"))
    parser.add_argument("--scatter-gif", type=Path, default=Path("server_plots/current_distribution_plots/radiant_scatter_showers.gif"))
    parser.add_argument("--histogram-gif", type=Path, default=Path("server_plots/current_distribution_plots/radiant_histogram_showers.gif"))
    parser.add_argument("--frame-dir", type=Path, default=Path("server_plots/current_distribution_plots/shower_animation_frames"))
    parser.add_argument("--shower-catalog", type=Path, default=None)
    parser.add_argument("--shower-peak-tolerance-deg", type=float, default=5.0)
    parser.add_argument("--bins", type=int, default=144)
    parser.add_argument("--window-days", type=int, default=3, help="Centered day averaging window for non-cumulative frames.")
    parser.add_argument("--cumulative", action="store_true", help="Accumulate all previous days in each frame.")
    args = parser.parse_args()

    arr = read_radiants(args.input_h5)
    order = np.argsort(arr["epoch_unix"])
    arr = arr[order]
    cumulative = args.cumulative
    n_scatter = make_animation(
        arr,
        args.scatter_gif,
        args.frame_dir / "scatter",
        "scatter",
        args.shower_catalog,
        args.shower_peak_tolerance_deg,
        cumulative,
        args.window_days,
        args.bins,
    )
    n_hist = make_animation(
        arr,
        args.histogram_gif,
        args.frame_dir / "histogram",
        "histogram",
        args.shower_catalog,
        args.shower_peak_tolerance_deg,
        cumulative,
        args.window_days,
        args.bins,
    )
    print(f"scatter_frames {n_scatter} {args.scatter_gif}")
    print(f"histogram_frames {n_hist} {args.histogram_gif}")


if __name__ == "__main__":
    main()
