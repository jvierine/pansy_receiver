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
from scipy.ndimage import gaussian_filter


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles-dir", type=Path, required=True)
    parser.add_argument("--output-h5", type=Path, required=True)
    parser.add_argument("--output-plot", type=Path, required=True)
    parser.add_argument("--minimum-path-km", type=float, default=15.0)
    parser.add_argument("--minimum-speed-km-s", type=float, default=10.0)
    parser.add_argument("--maximum-speed-km-s", type=float, default=80.0)
    return parser.parse_args()


def scalar(group, name):
    return float(group[name][()])


def load_profiles(profiles_dir: Path):
    rows = []
    lower_statuses = []
    upper_statuses = []
    paths = sorted(profiles_dir.glob("mass_profile_*.h5"))
    for index, path in enumerate(paths, start=1):
        with h5py.File(path, "r") as handle:
            result = handle["result"]
            quality = handle["quality"]
            rows.append(
                (
                    int(handle.attrs["sample_idx"]),
                    float(handle.attrs["sample_epoch_unix"]),
                    scalar(result, "free_best_radius_um"),
                    scalar(result, "free_best_mass_kg"),
                    scalar(result, "profile_ci95_lower_radius_um"),
                    scalar(result, "profile_ci95_upper_radius_um"),
                    scalar(result, "profile_ci95_lower_mass_kg"),
                    scalar(result, "profile_ci95_upper_mass_kg"),
                    *np.asarray(result["marginal_radius_quantiles_um"], dtype=np.float64),
                    *np.asarray(result["marginal_mass_quantiles_kg"], dtype=np.float64),
                    scalar(quality, "initial_speed_km_s"),
                    scalar(quality, "fitted_speed_mean_km_s"),
                    scalar(quality, "path_length_km"),
                    scalar(quality, "free_position_3d_rms_km"),
                    scalar(quality, "free_doppler_rms_km_s"),
                    int(quality["n_measurements"][()]),
                )
            )
            lower_statuses.append(str(result.attrs["profile_ci95_lower_status"]))
            upper_statuses.append(str(result.attrs["profile_ci95_upper_status"]))
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


def survival_curve(values, mass_grid):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values) | np.isposinf(values)]
    if len(values) == 0:
        return np.full(len(mass_grid), np.nan)
    sorted_values = np.sort(values)
    return (len(sorted_values) - np.searchsorted(sorted_values, mass_grid, side="right")) / len(sorted_values)


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
    levels = []
    for enclosed_fraction in (0.95, 0.80, 0.50):
        index = min(int(np.searchsorted(cumulative, enclosed_fraction)), len(ordered) - 1)
        levels.append(ordered[index])
    levels = np.unique(np.sort(levels))
    speed_centers = 0.5 * (speed_edges[:-1] + speed_edges[1:])
    mass_centers = 10.0 ** (0.5 * (mass_edges[:-1] + mass_edges[1:]))
    ax.contour(
        speed_centers,
        mass_centers,
        density,
        levels=levels,
        colors=[color],
        linewidths=np.linspace(0.8, 1.4, len(levels)),
        alpha=0.95,
    )


def fit_survival_slope(mass_grid, survival, minimum_fraction=0.10, maximum_fraction=0.50):
    keep = (
        np.isfinite(mass_grid)
        & (mass_grid > 0.0)
        & np.isfinite(survival)
        & (survival >= minimum_fraction)
        & (survival <= maximum_fraction)
    )
    if np.count_nonzero(keep) < 5:
        return np.nan, np.nan, np.asarray([], dtype=np.float64)
    coefficients = np.polyfit(np.log10(mass_grid[keep]), np.log10(survival[keep]), 1)
    return float(coefficients[0]), float(coefficients[1]), np.flatnonzero(keep)


def add_slope_marker(ax, slope, intercept, fit_mass):
    if not np.isfinite(slope) or len(fit_mass) < 2:
        return
    x0 = float(np.quantile(fit_mass, 0.58))
    x1 = min(float(np.max(fit_mass)), x0 * 4.0)
    y0 = float(10.0 ** (slope * np.log10(x0) + intercept))
    y1 = float(y0 * (x1 / x0) ** slope)
    ax.plot([x0, x1, x1, x0], [y0, y0, y1, y0], color="C0", lw=0.9)
    ax.text(
        x1 * 1.08,
        np.sqrt(y0 * y1),
        rf"$s={slope:.2f}$",
        color="C0",
        fontsize=8,
        ha="left",
        va="center",
    )


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

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), dpi=150)
    ax = axes[0]
    finite_upper = analysis_mask & np.isfinite(data["upper_mass_kg"])
    density_contours(ax, speed[analysis_mask], data["lower_mass_kg"][analysis_mask], "C0")
    density_contours(ax, speed[finite_upper], data["upper_mass_kg"][finite_upper], "C1")
    ax.set_yscale("log")
    ax.set_xlim(10.0, 80.0)
    ax.set_xlabel(r"Initial fitted speed (km s$^{-1}$)")
    ax.set_ylabel(r"Initial mass $m_0$ (kg)")
    ax.grid(alpha=0.2, which="both", linewidth=0.5)
    ax.legend(
        handles=(
            Line2D([0], [0], color="C0", lw=1.3, label="95% lower bound"),
            Line2D([0], [0], color="C1", lw=1.3, label="Finite 95% upper bound"),
        ),
        loc="lower right",
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

    ax = axes[1]
    valid_lower = data["lower_mass_kg"][analysis_mask]
    valid_upper = data["upper_mass_kg"][analysis_mask]
    positive = valid_lower[np.isfinite(valid_lower) & (valid_lower > 0.0)]
    finite_upper = valid_upper[np.isfinite(valid_upper) & (valid_upper > 0.0)]
    limits = np.r_[positive, finite_upper]
    mass_grid = np.geomspace(np.nanpercentile(limits, 0.2), np.nanpercentile(limits, 99.8), 300)
    lower_survival = survival_curve(valid_lower, mass_grid)
    upper_survival = survival_curve(valid_upper, mass_grid)
    slope, intercept, fit_indices = fit_survival_slope(mass_grid, lower_survival)
    ax.plot(mass_grid, lower_survival, color="C0", lw=1.3, label="95% lower bound")
    ax.plot(mass_grid, upper_survival, color="C1", lw=1.3, label="95% upper bound")
    if len(fit_indices) > 0:
        fit_mass = mass_grid[fit_indices]
        ax.plot(
            fit_mass,
            10.0 ** (slope * np.log10(fit_mass) + intercept),
            color="C0",
            lw=0.9,
            ls=":",
        )
        add_slope_marker(ax, slope, intercept, fit_mass)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 1.0)
    ax.set_xlabel(r"Initial mass threshold $m_0$ (kg)")
    ax.set_ylabel("Cumulative fraction above threshold")
    ax.grid(alpha=0.2, which="both", linewidth=0.5)
    ax.legend(loc="lower left", frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.14, top=0.98, wspace=0.28)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    return slope


def main():
    args = parse_args()
    data = load_profiles(args.profiles_dir)
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
    slope = plot_summary(args.output_plot, data, analysis_mask, args.minimum_path_km)
    with h5py.File(args.output_h5, "r+") as handle:
        handle.attrs["lower_bound_survival_slope"] = slope
        handle.attrs["slope_fit_survival_fraction_min"] = 0.10
        handle.attrs["slope_fit_survival_fraction_max"] = 0.50
    print(f"total profiles: {len(analysis_mask)}")
    print(f"analysis profiles: {np.count_nonzero(analysis_mask)}")
    print(f"lower status: {Counter(data['lower_status'])}")
    print(f"upper status: {Counter(data['upper_status'])}")
    print(f"lower-bound cumulative slope: {slope:.4f}")
    print(f"wrote {args.output_h5}")
    print(f"wrote {args.output_plot}")


if __name__ == "__main__":
    main()
