#!/usr/bin/env python3
"""Compare the radiant compactness of matched SER clouds at two inclinations."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np

from plot_omega_eridanids_shower import load_catalogue, selection_masks


HIGH_INCLINATION_DEG = 91.22
LOW_INCLINATION_DEG = 10.0
MU_SUN_AU3_YR2 = 4.0 * np.pi**2
ELEMENT_COLUMNS = "a_AU,e,i_deg,Omega_deg,omega_deg,nu_deg,q_AU"
RADIANT_COLUMNS = "sun_centered_ecliptic_longitude_deg,ecliptic_latitude_deg"


def wrap180(values):
    return (np.asarray(values, dtype=np.float64) + 180.0) % 360.0 - 180.0


def circular_mean_deg(values) -> float:
    angle = np.deg2rad(np.asarray(values, dtype=np.float64))
    mean = float(np.rad2deg(np.angle(np.mean(np.exp(1j * angle)))) % 360.0)
    return 0.0 if np.isclose(mean, 360.0) else mean


def matched_inclination_clouds(elements: np.ndarray) -> np.ndarray:
    """Preserve the measured cloud while changing only its mean inclination."""
    high = np.asarray(elements, dtype=np.float64).copy()
    high[:, 2] += HIGH_INCLINATION_DEG - np.mean(high[:, 2])
    low = high.copy()
    low[:, 2] += LOW_INCLINATION_DEG - HIGH_INCLINATION_DEG
    return np.stack((high, low), axis=0)


def rotation_matrix(node_deg: float, inc_deg: float, argp_deg: float) -> np.ndarray:
    node, inc, argp = np.deg2rad([node_deg, inc_deg, argp_deg])
    cn, sn = np.cos(node), np.sin(node)
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(argp), np.sin(argp)
    return np.asarray(
        [
            [cn * cw - sn * sw * ci, -cn * sw - sn * cw * ci, sn * si],
            [sn * cw + cn * sw * ci, -sn * sw + cn * cw * ci, -cn * si],
            [sw * si, cw * si, ci],
        ]
    )


def geocentric_radiants(
    elements: np.ndarray,
    solar_longitude_deg: np.ndarray,
) -> np.ndarray:
    """Map heliocentric osculating elements to Sun-centered ecliptic radiants."""
    radiants = np.full((len(elements), 2), np.nan, dtype=np.float64)
    for row, (a, e, inc, node, argp, anomaly, _q) in enumerate(elements):
        p = a * (1.0 - e * e)
        if not np.isfinite(p) or p <= 0.0:
            continue
        anomaly_rad = np.deg2rad(anomaly)
        velocity_perifocal = np.sqrt(MU_SUN_AU3_YR2 / p) * np.asarray(
            [-np.sin(anomaly_rad), e + np.cos(anomaly_rad), 0.0]
        )
        velocity_meteor = rotation_matrix(node, inc, argp) @ velocity_perifocal

        earth_longitude = np.deg2rad((solar_longitude_deg[row] + 180.0) % 360.0)
        velocity_earth = np.sqrt(MU_SUN_AU3_YR2) * np.asarray(
            [-np.sin(earth_longitude), np.cos(earth_longitude), 0.0]
        )
        radiant_vector = -(velocity_meteor - velocity_earth)
        radiant_vector /= np.linalg.norm(radiant_vector)
        longitude = np.rad2deg(np.arctan2(radiant_vector[1], radiant_vector[0])) % 360.0
        latitude = np.rad2deg(np.arcsin(np.clip(radiant_vector[2], -1.0, 1.0)))
        radiants[row] = (
            wrap180(longitude - solar_longitude_deg[row]),
            latitude,
        )
    return radiants


def radiant_offsets(radiants: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    center_lon = circular_mean_deg(radiants[:, 0])
    center_lat = float(np.nanmean(radiants[:, 1]))
    raw = np.column_stack(
        (
            wrap180(radiants[:, 0] - center_lon),
            radiants[:, 1] - center_lat,
        )
    )
    tangent = raw.copy()
    tangent[:, 0] *= np.cos(np.deg2rad(center_lat))
    return raw, tangent


def radiant_metrics(radiants: np.ndarray) -> dict[str, float]:
    raw, tangent = radiant_offsets(radiants)
    covariance = np.cov(tangent, rowvar=False, ddof=1)
    eigenvalues = np.linalg.eigvalsh(covariance)[::-1]
    return {
        "mean_lon_deg": circular_mean_deg(radiants[:, 0]),
        "mean_lat_deg": float(np.nanmean(radiants[:, 1])),
        "sigma_lon_deg": float(np.nanstd(raw[:, 0], ddof=1)),
        "sigma_lat_deg": float(np.nanstd(raw[:, 1], ddof=1)),
        "tangent_major_sigma_deg": float(np.sqrt(eigenvalues[0])),
        "tangent_minor_sigma_deg": float(np.sqrt(eigenvalues[1])),
    }


def add_covariance_ellipse(ax, offsets: np.ndarray, color: str) -> None:
    covariance = np.cov(offsets, rowvar=False, ddof=1)
    values, vectors = np.linalg.eigh(covariance)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    angle = np.rad2deg(np.arctan2(vectors[1, 0], vectors[0, 0]))
    # sqrt(5.991) gives the 95% contour for a two-dimensional Gaussian.
    scale = np.sqrt(5.991)
    ellipse = Ellipse(
        (0.0, 0.0),
        width=2.0 * scale * np.sqrt(values[0]),
        height=2.0 * scale * np.sqrt(values[1]),
        angle=angle,
        fill=False,
        color=color,
        lw=1.5,
    )
    ax.add_patch(ellipse)


def plot_comparison(
    output: Path,
    radiants: np.ndarray,
    metrics: list[dict[str, float]],
) -> None:
    labels = (r"SER geometry, $\langle i\rangle=91.22^\circ$", r"Matched cloud, $\langle i\rangle=10^\circ$")
    colors = ("#1f77b4", "#d95f02")
    offsets = [radiant_offsets(radiants[index])[0] for index in range(2)]
    x_limit = 1.1 * max(np.nanmax(np.abs(value[:, 0])) for value in offsets)
    y_limit = 1.1 * max(np.nanmax(np.abs(value[:, 1])) for value in offsets)

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2), layout="constrained")
    for index in range(2):
        axes[index].scatter(
            offsets[index][:, 0],
            offsets[index][:, 1],
            s=16,
            color=colors[index],
            alpha=0.65,
            linewidths=0,
        )
        add_covariance_ellipse(axes[index], offsets[index], colors[index])
        axes[index].axhline(0.0, color="0.75", lw=0.7)
        axes[index].axvline(0.0, color="0.75", lw=0.7)
        axes[index].set_xlim(-x_limit, x_limit)
        axes[index].set_ylim(-y_limit, y_limit)
        axes[index].set_xlabel(r"$\Delta(\lambda_g-\lambda_\odot)$ (deg)")
        axes[index].set_ylabel(r"$\Delta\beta_g$ (deg)")
        axes[index].set_title(labels[index], fontsize=11)
        axes[index].grid(alpha=0.16)

    metric_names = ("sigma_lon_deg", "sigma_lat_deg", "tangent_major_sigma_deg")
    metric_labels = (
        r"$\sigma_{\lambda'}$",
        r"$\sigma_\beta$",
        r"Tangent-plane $\sigma_{\rm major}$",
    )
    x = np.arange(len(metric_names))
    width = 0.36
    for index in range(2):
        axes[2].bar(
            x + (index - 0.5) * width,
            [metrics[index][name] for name in metric_names],
            width=width,
            color=colors[index],
            label=labels[index],
        )
    axes[2].set_xticks(x, metric_labels)
    axes[2].set_ylabel("Radiant width (deg)")
    axes[2].legend(frameon=False, fontsize=8)
    axes[2].grid(axis="y", alpha=0.2)
    fig.suptitle("Inclination-dependent projection of the same dispersed orbital cloud", fontsize=14)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def save_results(
    output: Path,
    solar_longitude_deg: np.ndarray,
    elements: np.ndarray,
    radiants: np.ndarray,
    metrics: list[dict[str, float]],
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        handle.attrs["description"] = (
            "Static radiant projection of the same SER orbital cloud at mean inclinations 91.22 and 10 deg"
        )
        handle.attrs["element_columns"] = ELEMENT_COLUMNS
        handle.attrs["radiant_columns"] = RADIANT_COLUMNS
        handle.attrs["ensemble_names"] = "SER i=91.22 deg,matched cloud i=10 deg"
        handle.attrs["dynamical_evolution"] = "none"
        handle.create_dataset("solar_longitude_deg", data=solar_longitude_deg.astype(np.float32))
        handle.create_dataset("elements", data=elements.astype(np.float32), compression="gzip", compression_opts=3)
        handle.create_dataset("radiants", data=radiants.astype(np.float32), compression="gzip", compression_opts=3)
        group = handle.create_group("metrics")
        for name in metrics[0]:
            group.create_dataset(
                name,
                data=np.asarray([value[name] for value in metrics], dtype=np.float32),
            )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalogue",
        type=Path,
        default=Path("figs/pansy_maarsy_keplerian_catalogue.h5"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("figs/ser_radiant_inclination_geometry.h5"),
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=Path("figs/ser_radiant_inclination_geometry.png"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = load_catalogue(args.catalogue)
    _shower, core = selection_masks(rows)
    elements = matched_inclination_clouds(rows["kepler"][core])
    solar_longitude = rows["solar"][core]
    radiants = np.stack(
        [geocentric_radiants(cloud, solar_longitude) for cloud in elements],
        axis=0,
    )
    metrics = [radiant_metrics(cloud) for cloud in radiants]
    save_results(args.output, solar_longitude, elements, radiants, metrics)
    plot_comparison(args.plot, radiants, metrics)
    for index, label in enumerate(("SER i=91.22 deg", "matched i=10 deg")):
        print(
            f"{label}: sigma_lon={metrics[index]['sigma_lon_deg']:.3f} deg, "
            f"sigma_lat={metrics[index]['sigma_lat_deg']:.3f} deg, "
            f"major={metrics[index]['tangent_major_sigma_deg']:.3f} deg, "
            f"minor={metrics[index]['tangent_minor_sigma_deg']:.3f} deg"
        )
    print(args.output)
    print(args.plot)


if __name__ == "__main__":
    main()
