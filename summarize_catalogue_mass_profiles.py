#!/usr/bin/env python3
"""Aggregate and plot a random-subset PANSY mass-profile run."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy.ndimage import gaussian_filter


METEOROID_DENSITY_KG_M3 = 3000.0


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles-dir", type=Path, required=True)
    parser.add_argument(
        "--baseline-profiles-dir",
        type=Path,
        help="Baseline catalogue profiles supplying trajectory quality metadata",
    )
    parser.add_argument("--output-h5", type=Path, required=True)
    parser.add_argument("--output-plot", type=Path, required=True)
    parser.add_argument("--minimum-path-km", type=float, default=15.0)
    parser.add_argument("--minimum-speed-km-s", type=float, default=10.0)
    parser.add_argument("--maximum-speed-km-s", type=float, default=80.0)
    return parser.parse_args()


def scalar(group, name):
    return float(group[name][()])


def load_profiles(profiles_dir: Path, baseline_profiles_dir: Path | None = None):
    rows = []
    lower_statuses = []
    upper_statuses = []
    paths = sorted(profiles_dir.glob("mass_profile_*.h5"))
    for index, path in enumerate(paths, start=1):
        with h5py.File(path, "r") as handle:
            result = handle["result"]
            sample_idx = int(handle.attrs["sample_idx"])
            baseline = None
            if "quality" in handle:
                quality = handle["quality"]
                metadata = handle
            else:
                if baseline_profiles_dir is None:
                    raise RuntimeError(
                        f"{path} has no quality group; --baseline-profiles-dir is required"
                    )
                baseline_path = baseline_profiles_dir / f"mass_profile_{sample_idx}.h5"
                baseline = h5py.File(baseline_path, "r")
                quality = baseline["quality"]
                metadata = baseline

            if "free_best_radius_um" in result:
                best_radius_um = scalar(result, "free_best_radius_um")
                lower_radius_um = scalar(result, "profile_ci95_lower_radius_um")
                upper_radius_um = scalar(result, "profile_ci95_upper_radius_um")
                best_mass_kg = scalar(result, "free_best_mass_kg")
                lower_mass_kg = scalar(result, "profile_ci95_lower_mass_kg")
                upper_mass_kg = scalar(result, "profile_ci95_upper_mass_kg")
                marginal_radius = np.asarray(result["marginal_radius_quantiles_um"], dtype=np.float64)
                marginal_mass = np.asarray(result["marginal_mass_quantiles_kg"], dtype=np.float64)
                lower_status = str(result.attrs["profile_ci95_lower_status"])
                upper_status = str(result.attrs["profile_ci95_upper_status"])
            else:
                best_radius_um = scalar(result, "best_radius_um")
                lower_radius_um = scalar(result, "ci95_lower_radius_um")
                upper_radius_um = scalar(result, "ci95_upper_radius_um")
                best_mass_kg = float(radius_um_to_mass_kg(best_radius_um))
                lower_mass_kg = float(radius_um_to_mass_kg(lower_radius_um))
                upper_mass_kg = float(radius_um_to_mass_kg(upper_radius_um))
                marginal_radius = np.asarray(result["marginal_radius_quantiles_um"], dtype=np.float64)
                marginal_mass = radius_um_to_mass_kg(marginal_radius)
                lower_status = str(result.attrs["ci95_lower_status"])
                upper_status = str(result.attrs["ci95_upper_status"])
            rows.append(
                (
                    sample_idx,
                    float(metadata.attrs["sample_epoch_unix"]),
                    best_radius_um,
                    best_mass_kg,
                    lower_radius_um,
                    upper_radius_um,
                    lower_mass_kg,
                    upper_mass_kg,
                    *marginal_radius,
                    *marginal_mass,
                    scalar(quality, "initial_speed_km_s"),
                    scalar(quality, "fitted_speed_mean_km_s"),
                    scalar(quality, "path_length_km"),
                    scalar(quality, "free_position_3d_rms_km"),
                    scalar(quality, "free_doppler_rms_km_s"),
                    int(quality["n_measurements"][()]),
                )
            )
            lower_statuses.append(lower_status)
            upper_statuses.append(upper_status)
            if baseline is not None:
                baseline.close()
        if index % 10000 == 0:
            print(f"read {index} profiles", flush=True)
    names = (
        "sample_idx",
        "sample_epoch_unix",
        "best_radius_um",
        "best_mass_kg",
        "lower_radius_um",
        "upper_radius_um",
        "lower_mass_kg",
        "upper_mass_kg",
        "marginal_radius_q025_um",
        "marginal_radius_q50_um",
        "marginal_radius_q975_um",
        "marginal_mass_q025_kg",
        "marginal_mass_q50_kg",
        "marginal_mass_q975_kg",
        "initial_speed_km_s",
        "fitted_speed_mean_km_s",
        "path_length_km",
        "position_rms_km",
        "doppler_rms_km_s",
        "n_measurements",
    )
    columns = np.asarray(rows, dtype=np.float64).T
    data = {name: values for name, values in zip(names, columns, strict=True)}
    data["sample_idx"] = data["sample_idx"].astype(np.int64)
    data["n_measurements"] = data["n_measurements"].astype(np.int32)
    data["lower_status"] = np.asarray(lower_statuses)
    data["upper_status"] = np.asarray(upper_statuses)
    return data


def mass_kg_to_radius_um(mass_kg):
    mass_kg = np.maximum(np.asarray(mass_kg, dtype=np.float64), 1e-300)
    return (3.0 * mass_kg / (4.0 * np.pi * METEOROID_DENSITY_KG_M3)) ** (1.0 / 3.0) * 1e6


def radius_um_to_mass_kg(radius_um):
    radius_m = np.maximum(np.asarray(radius_um, dtype=np.float64), 1e-300) * 1e-6
    return (4.0 / 3.0) * np.pi * METEOROID_DENSITY_KG_M3 * radius_m**3


def density_contours(ax, speed, mass_kg, color):
    good = np.isfinite(speed) & np.isfinite(mass_kg) & (mass_kg > 0.0)
    speed = np.asarray(speed[good], dtype=np.float64)
    log_mass = np.log10(np.asarray(mass_kg[good], dtype=np.float64))
    speed_edges = np.linspace(10.0, 80.0, 71)
    mass_edges = np.linspace(-14.0, -2.0, 97)
    histogram, _, _ = np.histogram2d(speed, log_mass, bins=(speed_edges, mass_edges))
    density = gaussian_filter(histogram.T, sigma=(1.4, 1.4), mode="nearest")
    positive = density[density > 0.0]
    if len(positive) == 0:
        return
    ordered = np.sort(positive)[::-1]
    cumulative = np.cumsum(ordered) / np.sum(ordered)
    level_fractions = []
    for enclosed_fraction in (0.95, 0.80, 0.50):
        index = min(int(np.searchsorted(cumulative, enclosed_fraction)), len(ordered) - 1)
        level_fractions.append((float(ordered[index]), enclosed_fraction))
    unique_levels = {}
    for level, fraction in level_fractions:
        unique_levels[level] = min(fraction, unique_levels.get(level, 1.0))
    levels = np.asarray(sorted(unique_levels), dtype=np.float64)
    speed_centers = 0.5 * (speed_edges[:-1] + speed_edges[1:])
    mass_centers = 10.0 ** (0.5 * (mass_edges[:-1] + mass_edges[1:]))
    contours = ax.contour(
        speed_centers,
        mass_centers,
        density,
        levels=levels,
        colors=[color],
        linewidths=np.linspace(0.8, 1.4, len(levels)),
        alpha=0.95,
    )
    labels = {level: f"{100.0 * unique_levels[level]:.0f}%" for level in levels}
    ax.clabel(contours, contours.levels, fmt=labels, inline=True, inline_spacing=2, fontsize=7, colors=[color])


def save_summary(path: Path, data, analysis_mask, args):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["schema"] = "pansy.catalogue_mass_profile_summary.v1"
        handle.attrs["source_profiles_dir"] = str(args.profiles_dir)
        handle.attrs["minimum_path_km"] = float(args.minimum_path_km)
        handle.attrs["minimum_speed_km_s"] = float(args.minimum_speed_km_s)
        handle.attrs["maximum_speed_km_s"] = float(args.maximum_speed_km_s)
        handle.attrs["total_profiles"] = len(analysis_mask)
        handle.attrs["analysis_profiles"] = int(np.count_nonzero(analysis_mask))
        for name, values in data.items():
            if values.dtype.kind in "US":
                handle.create_dataset(name, data=values.astype(object), dtype=h5py.string_dtype("utf-8"))
            else:
                handle.create_dataset(name, data=values, compression="gzip")
        handle.create_dataset("analysis_mask", data=analysis_mask, compression="gzip")
    temporary.replace(path)


def plot_summary(path: Path, data, analysis_mask, minimum_path_km: float):
    speed = data["initial_speed_km_s"]
    mass_limits = (1e-12, 1e-5)

    fig = plt.figure(figsize=(7.4, 4.8), dpi=150)
    grid = fig.add_gridspec(1, 2, width_ratios=(4.6, 1.25), wspace=0.04)
    ax = fig.add_subplot(grid[0, 0])
    ax_mass = fig.add_subplot(grid[0, 1], sharey=ax)
    finite_upper = analysis_mask & np.isfinite(data["upper_mass_kg"])
    density_contours(ax, speed[analysis_mask], data["lower_mass_kg"][analysis_mask], "C0")
    density_contours(ax, speed[finite_upper], data["upper_mass_kg"][finite_upper], "C1")
    ax.set_yscale("log")
    ax.set_ylim(*mass_limits)
    ax.set_xlim(10.0, 80.0)
    ax.set_xlabel(r"Initial fitted speed (km s$^{-1}$)")
    ax.set_ylabel(r"Initial mass $m_0$ (kg)")
    ax.grid(alpha=0.2, which="both", linewidth=0.5)

    mass_bins = np.logspace(np.log10(mass_limits[0]), np.log10(mass_limits[1]), 43)
    ax_mass.hist(
        data["lower_mass_kg"][analysis_mask],
        bins=mass_bins,
        orientation="horizontal",
        histtype="step",
        color="C0",
        linewidth=1.3,
    )
    ax_mass.hist(
        data["upper_mass_kg"][finite_upper],
        bins=mass_bins,
        orientation="horizontal",
        histtype="step",
        color="C1",
        linewidth=1.3,
    )
    ax_mass.set_yscale("log")
    ax_mass.set_ylim(*mass_limits)
    ax_mass.set_xlim(left=0.0)
    ax_mass.xaxis.set_major_locator(MaxNLocator(nbins=3, integer=True, prune="lower"))
    ax_mass.set_xlabel("Count")
    ax_mass.tick_params(axis="y", which="both", left=False, labelleft=False)
    ax_mass.grid(alpha=0.2, which="both", linewidth=0.5)
    radius_axis = ax_mass.secondary_yaxis(
        "right", functions=(mass_kg_to_radius_um, radius_um_to_mass_kg)
    )
    radius_axis.set_yscale("log")
    radius_axis.set_ylabel(r"Initial radius $r_0$ ($\mu$m)")
    ax.legend(
        handles=(
            Line2D([0], [0], color="C0", lw=1.3, label="95% lower bound"),
            Line2D([0], [0], color="C1", lw=1.3, label="Finite 95% upper bound"),
        ),
        loc="lower left",
        frameon=False,
        fontsize=8,
    )
    ax.text(
        0.03,
        0.97,
        rf"Path $\geq{minimum_path_km:g}$ km" + f"\nN = {np.count_nonzero(analysis_mask):,}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    fig.subplots_adjust(left=0.12, right=0.89, bottom=0.14, top=0.98)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def main():
    args = parse_args()
    data = load_profiles(args.profiles_dir, args.baseline_profiles_dir)
    analysis_mask = (
        (data["path_length_km"] >= args.minimum_path_km)
        & (data["initial_speed_km_s"] >= args.minimum_speed_km_s)
        & (data["initial_speed_km_s"] <= args.maximum_speed_km_s)
        & (data["lower_status"] == "bounded")
        & np.isin(data["upper_status"], ["bounded", "open_grid"])
        & np.isfinite(data["best_mass_kg"])
        & np.isfinite(data["lower_mass_kg"])
        & (data["lower_mass_kg"] > 0.0)
    )
    save_summary(args.output_h5, data, analysis_mask, args)
    plot_summary(args.output_plot, data, analysis_mask, args.minimum_path_km)
    print(f"total profiles: {len(analysis_mask)}")
    print(f"analysis profiles: {np.count_nonzero(analysis_mask)}")
    print(f"lower status: {Counter(data['lower_status'])}")
    print(f"upper status: {Counter(data['upper_status'])}")
    print(f"wrote {args.output_h5}")
    print(f"wrote {args.output_plot}")


if __name__ == "__main__":
    main()
