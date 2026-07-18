#!/usr/bin/env python3
"""Aggregate and plot a random-subset PANSY mass-profile run."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


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
    lower_only = analysis_mask & (data["upper_status"] == "open_grid")
    two_sided = analysis_mask & (data["upper_status"] == "bounded")
    speed = data["initial_speed_km_s"]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), dpi=150)
    ax = axes[0]
    ax.scatter(
        speed[lower_only],
        data["lower_mass_kg"][lower_only],
        s=5,
        color="black",
        alpha=0.14,
        edgecolors="none",
        rasterized=True,
        label="Lower limit; upper interval open",
    )
    ax.scatter(
        speed[two_sided],
        data["best_mass_kg"][two_sided],
        s=6,
        color="C1",
        alpha=0.22,
        edgecolors="none",
        rasterized=True,
        label="Best fit; two-sided interval",
    )
    ax.set_yscale("log")
    ax.set_xlim(10.0, 80.0)
    ax.set_xlabel(r"Initial fitted speed (km s$^{-1}$)")
    ax.set_ylabel(r"Initial mass $m_0$ (kg)")
    ax.grid(alpha=0.2, which="both", linewidth=0.5)
    ax.legend(loc="lower right", frameon=False, fontsize=8)
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
    valid_best = data["best_mass_kg"][analysis_mask]
    valid_upper = data["upper_mass_kg"][analysis_mask]
    positive = valid_lower[np.isfinite(valid_lower) & (valid_lower > 0.0)]
    finite_best = valid_best[np.isfinite(valid_best) & (valid_best > 0.0)]
    finite_upper = valid_upper[np.isfinite(valid_upper) & (valid_upper > 0.0)]
    limits = np.r_[positive, finite_best, finite_upper]
    mass_grid = np.geomspace(np.nanpercentile(limits, 0.2), np.nanpercentile(limits, 99.8), 300)
    lower_survival = survival_curve(valid_lower, mass_grid)
    best_survival = survival_curve(valid_best, mass_grid)
    upper_survival = survival_curve(valid_upper, mass_grid)
    ax.fill_between(mass_grid, lower_survival, upper_survival, color="0.75", alpha=0.55, label="95% interval envelope")
    ax.plot(mass_grid, lower_survival, color="black", lw=1.2, label="Guaranteed above mass")
    ax.plot(mass_grid, best_survival, color="C1", lw=1.2, label="Best fit")
    ax.plot(mass_grid, upper_survival, color="0.35", lw=1.0, ls="--", label="Allowed above mass")
    ax.set_xscale("log")
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(r"Initial mass threshold $m_0$ (kg)")
    ax.set_ylabel("Cumulative fraction above threshold")
    ax.grid(alpha=0.2, which="both", linewidth=0.5)
    ax.legend(loc="lower left", frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.14, top=0.98, wspace=0.28)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


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
    plot_summary(args.output_plot, data, analysis_mask, args.minimum_path_km)
    print(f"total profiles: {len(analysis_mask)}")
    print(f"analysis profiles: {np.count_nonzero(analysis_mask)}")
    print(f"lower status: {Counter(data['lower_status'])}")
    print(f"upper status: {Counter(data['upper_status'])}")
    print(f"wrote {args.output_h5}")
    print(f"wrote {args.output_plot}")


if __name__ == "__main__":
    main()
