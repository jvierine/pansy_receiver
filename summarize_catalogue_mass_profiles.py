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
from scipy.interpolate import UnivariateSpline


METEOROID_DENSITY_KG_M3 = 3000.0
DEFAULT_ZENITH_WEIGHT_ALPHA = 2.095
DEFAULT_MINIMUM_COS_ZENITH = 0.15
DEFAULT_VELOCITY_REFERENCE_KM_S = 73.0
PAPER_AXIS_LABEL_SIZE = 18
PAPER_TICK_LABEL_SIZE = 16
PAPER_ANNOTATION_SIZE = 15
PAPER_CONTOUR_LABEL_SIZE = 14


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
    parser.add_argument(
        "--zenith-weight-alpha",
        type=float,
        default=DEFAULT_ZENITH_WEIGHT_ALPHA,
    )
    parser.add_argument(
        "--minimum-cos-zenith",
        type=float,
        default=DEFAULT_MINIMUM_COS_ZENITH,
    )
    parser.add_argument(
        "--velocity-reference-km-s",
        type=float,
        default=DEFAULT_VELOCITY_REFERENCE_KM_S,
    )
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
            free_parameters = np.asarray(
                metadata["result"]["free_best_parameters7"],
                dtype=np.float64,
            )
            initial_velocity_enu_mps = free_parameters[3:6]
            initial_speed_mps = float(np.linalg.norm(initial_velocity_enu_mps))
            radiant_cos_zenith = (
                -float(initial_velocity_enu_mps[2]) / initial_speed_mps
                if initial_speed_mps > 0.0
                else np.nan
            )
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
                    radiant_cos_zenith,
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
        "radiant_cos_zenith",
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


def density_contours(
    ax,
    speed,
    mass_kg,
    color,
    manual_label_positions=None,
):
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
    manual_label_positions = manual_label_positions or {}
    automatic_levels = []
    for level in contours.levels:
        fraction = unique_levels[level]
        if fraction in manual_label_positions:
            ax.clabel(
                contours,
                levels=[level],
                fmt=labels,
                manual=[manual_label_positions[fraction]],
                inline=True,
                inline_spacing=2,
                fontsize=PAPER_CONTOUR_LABEL_SIZE,
                colors=[color],
            )
        else:
            automatic_levels.append(level)
    if automatic_levels:
        ax.clabel(
            contours,
            levels=automatic_levels,
            fmt=labels,
            inline=True,
            inline_spacing=2,
            fontsize=PAPER_CONTOUR_LABEL_SIZE,
            colors=[color],
        )
    return speed_centers, mass_centers, histogram, density


def plot_smoothed_density_peak(
    ax,
    speed_centers,
    mass_centers,
    histogram,
    density,
    color,
):
    """Plot a smooth conditional-mode ridge through a velocity--mass density."""
    counts_per_speed = np.sum(histogram, axis=1)
    populated = counts_per_speed >= 8
    if np.count_nonzero(populated) < 5:
        return

    peak_indices = np.argmax(density[:, populated], axis=0)
    peak_log_mass = np.log10(mass_centers[peak_indices])
    peak_speed = speed_centers[populated]

    # Fit in log-mass space. The smoothing target allows about 0.15 dex of
    # bin-to-bin scatter while retaining broad changes with entry speed.
    weights = np.sqrt(counts_per_speed[populated])
    weights /= np.median(weights)
    spline = UnivariateSpline(
        peak_speed,
        peak_log_mass,
        w=weights,
        k=min(3, len(peak_speed) - 1),
        s=len(peak_speed) * 0.15**2,
    )
    fit_speed = np.linspace(peak_speed[0], peak_speed[-1], 400)
    ax.plot(
        fit_speed,
        10.0 ** spline(fit_speed),
        color=color,
        linewidth=2.2,
        linestyle="--",
        zorder=5,
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
        handle.attrs["zenith_weight_alpha"] = float(args.zenith_weight_alpha)
        handle.attrs["minimum_cos_zenith"] = float(args.minimum_cos_zenith)
        handle.attrs["velocity_reference_km_s"] = float(
            args.velocity_reference_km_s
        )
        handle.attrs["zenith_weight"] = (
            "max(radiant_cos_zenith, minimum_cos_zenith)^(-zenith_weight_alpha)"
        )
        handle.attrs["velocity_weight"] = (
            "(velocity_reference_km_s / initial_speed_km_s)^3"
        )
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

    # This wide source image is reduced to single-column width in the paper.
    # Use approximately twice the normal Matplotlib type sizes so the final
    # rendered labels remain comparable to the article's 8--9 pt body text.
    fig = plt.figure(figsize=(7.4, 6.0), dpi=200)
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=(4.6, 1.25),
        height_ratios=(1.15, 4.6),
        wspace=0.04,
        hspace=0.04,
    )
    ax = fig.add_subplot(grid[1, 0])
    ax_mass = fig.add_subplot(grid[1, 1], sharey=ax)
    ax_speed = fig.add_subplot(grid[0, 0], sharex=ax)
    finite_upper = analysis_mask & np.isfinite(data["upper_mass_kg"])
    lower_density = density_contours(
        ax,
        speed[analysis_mask],
        data["lower_mass_kg"][analysis_mask],
        "C0",
        manual_label_positions={0.95: (30.0, 5e-10)},
    )
    upper_density = density_contours(
        ax,
        speed[finite_upper],
        data["upper_mass_kg"][finite_upper],
        "C1",
    )
    if lower_density is not None:
        plot_smoothed_density_peak(ax, *lower_density, "C0")
    if upper_density is not None:
        plot_smoothed_density_peak(ax, *upper_density, "C1")
    ax.set_yscale("log")
    ax.set_ylim(*mass_limits)
    ax.set_xlim(10.0, 80.0)
    ax.set_xlabel(r"Initial fitted speed (km s$^{-1}$)", fontsize=PAPER_AXIS_LABEL_SIZE)
    ax.set_ylabel(r"Initial mass $m_0$ (kg)", fontsize=PAPER_AXIS_LABEL_SIZE)
    ax.tick_params(axis="both", which="both", labelsize=PAPER_TICK_LABEL_SIZE)
    ax.grid(alpha=0.2, which="both", linewidth=0.5)

    speed_bins = np.linspace(10.0, 80.0, 36)
    raw_speed_count, _ = np.histogram(speed[analysis_mask], bins=speed_bins)
    debiased_speed_count, _ = np.histogram(
        speed[analysis_mask],
        bins=speed_bins,
        weights=data["debiased_weight"][analysis_mask],
    )
    raw_line = ax_speed.stairs(
        raw_speed_count,
        speed_bins,
        color="0.2",
        linewidth=1.3,
        label="Raw count",
    )
    ax_speed.set_ylabel("Raw count", fontsize=14)
    ax_speed.tick_params(
        axis="both",
        which="both",
        labelsize=13,
        labelbottom=False,
    )
    ax_speed.yaxis.set_major_locator(
        MaxNLocator(nbins=3, integer=True, prune="lower")
    )
    ax_speed.grid(alpha=0.2, linewidth=0.5)
    ax_speed_debiased = ax_speed.twinx()
    debiased_line = ax_speed_debiased.stairs(
        debiased_speed_count,
        speed_bins,
        color="C2",
        linewidth=1.3,
        label=r"$\sum_i w_{z,i}w_{v,i}$",
    )
    ax_speed_debiased.tick_params(
        axis="y",
        colors="C2",
        labelsize=13,
    )
    ax_speed_debiased.yaxis.set_major_locator(
        MaxNLocator(nbins=3, prune="lower")
    )
    ax_speed_debiased.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax_speed_debiased.yaxis.get_offset_text().set_color("C2")
    ax_speed_debiased.yaxis.get_offset_text().set_fontsize(11)
    ax_speed.legend(
        handles=(raw_line, debiased_line),
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        frameon=False,
        fontsize=12,
    )

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
    ax_mass.set_xlabel("Count", fontsize=PAPER_AXIS_LABEL_SIZE)
    ax_mass.tick_params(axis="x", which="both", labelsize=PAPER_TICK_LABEL_SIZE)
    ax_mass.tick_params(axis="y", which="both", left=False, labelleft=False)
    ax_mass.grid(alpha=0.2, which="both", linewidth=0.5)
    radius_axis = ax_mass.secondary_yaxis(
        "right", functions=(mass_kg_to_radius_um, radius_um_to_mass_kg)
    )
    radius_axis.set_yscale("log")
    radius_axis.set_ylabel(
        r"Initial radius $r_0$ ($\mu$m)", fontsize=PAPER_AXIS_LABEL_SIZE
    )
    radius_axis.tick_params(axis="y", which="both", labelsize=PAPER_TICK_LABEL_SIZE)
    ax.legend(
        handles=(
            Line2D(
                [0],
                [0],
                color="C0",
                lw=2.2,
                linestyle="--",
                label="95% lower bound",
            ),
            Line2D(
                [0],
                [0],
                color="C1",
                lw=2.2,
                linestyle="--",
                label="Finite 95% upper bound",
            ),
        ),
        loc="lower left",
        frameon=False,
        fontsize=PAPER_ANNOTATION_SIZE,
    )
    ax.text(
        0.03,
        0.97,
        rf"Path $\geq{minimum_path_km:g}$ km" + f"\nN = {np.count_nonzero(analysis_mask):,}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=PAPER_ANNOTATION_SIZE,
    )
    fig.subplots_adjust(left=0.12, right=0.89, bottom=0.11, top=0.98)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def main():
    args = parse_args()
    data = load_profiles(args.profiles_dir, args.baseline_profiles_dir)
    cos_zenith = np.asarray(data["radiant_cos_zenith"], dtype=np.float64)
    zenith_weight = np.maximum(
        cos_zenith,
        float(args.minimum_cos_zenith),
    ) ** (-float(args.zenith_weight_alpha))
    velocity_weight = (
        float(args.velocity_reference_km_s)
        / np.asarray(data["initial_speed_km_s"], dtype=np.float64)
    ) ** 3
    valid_weight = (
        np.isfinite(zenith_weight)
        & np.isfinite(velocity_weight)
        & (cos_zenith > 0.0)
    )
    zenith_weight[~valid_weight] = np.nan
    velocity_weight[~valid_weight] = np.nan
    data["zenith_weight"] = zenith_weight
    data["velocity_weight"] = velocity_weight
    data["debiased_weight"] = zenith_weight * velocity_weight
    analysis_mask = (
        (data["path_length_km"] >= args.minimum_path_km)
        & (data["initial_speed_km_s"] >= args.minimum_speed_km_s)
        & (data["initial_speed_km_s"] <= args.maximum_speed_km_s)
        & (data["lower_status"] == "bounded")
        & np.isin(data["upper_status"], ["bounded", "open_grid"])
        & np.isfinite(data["best_mass_kg"])
        & np.isfinite(data["lower_mass_kg"])
        & (data["lower_mass_kg"] > 0.0)
        & np.isfinite(data["debiased_weight"])
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
