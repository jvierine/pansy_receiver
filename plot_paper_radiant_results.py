#!/usr/bin/env python3
"""Generate paper radiant-distribution and shower-summary figures."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from matplotlib.colors import LogNorm


PANSY_LOCATION = EarthLocation(lat=-69.010833 * u.deg, lon=39.599722 * u.deg, height=100.0 * u.m)
PLOT_CENTER_LONGITUDE_DEG = -90.0


@dataclass(frozen=True)
class Shower:
    name: str
    solar_lon_deg: float
    sc_lon_deg: float
    beta_deg: float
    ra_deg: float
    dec_deg: float
    vg_km_s: float
    n: int
    semi_major_axis_au: float
    eccentricity: float
    inclination_deg: float
    note: str


SHOWERS = (
    Shower(
        name=r"$\omega$-Eridanids",
        solar_lon_deg=109.24,
        sc_lon_deg=288.74,
        beta_deg=-48.17,
        ra_deg=52.01,
        dec_deg=-31.36,
        vg_km_s=49.62,
        n=125,
        semi_major_axis_au=6.9259,
        eccentricity=0.7454,
        inclination_deg=92.70,
        note="No close established-shower match",
    ),
    Shower(
        name=r"$\kappa$-Scorpiids",
        solar_lon_deg=1.25,
        sc_lon_deg=240.62,
        beta_deg=-55.12,
        ra_deg=206.83,
        dec_deg=-72.41,
        vg_km_s=39.40,
        n=29,
        semi_major_axis_au=2.3628,
        eccentricity=0.5430,
        inclination_deg=71.77,
        note=r"Nearest broad match: $\lambda$-Centaurids",
    ),
    Shower(
        name=r"$\phi$-Capricornids",
        solar_lon_deg=312.72,
        sc_lon_deg=354.51,
        beta_deg=-8.34,
        ra_deg=312.00,
        dec_deg=-26.51,
        vg_km_s=25.02,
        n=19,
        semi_major_axis_au=2.3073,
        eccentricity=0.7859,
        inclination_deg=7.43,
        note="Possible Daytime Capricornids-Sagittariids substructure",
    ),
)


def wrap180(deg: np.ndarray | float) -> np.ndarray:
    return (np.asarray(deg, dtype=np.float64) + 180.0) % 360.0 - 180.0


def wrap360(deg: np.ndarray | float) -> np.ndarray:
    return np.asarray(deg, dtype=np.float64) % 360.0


def centered_plot_longitude_deg(lambda_minus_sun_deg: np.ndarray | float) -> np.ndarray:
    signed = wrap180(lambda_minus_sun_deg)
    return -wrap180(signed - PLOT_CENTER_LONGITUDE_DEG)


def angular_separation_deg(lon1, lat1, lon2, lat2):
    lon1 = np.deg2rad(lon1)
    lat1 = np.deg2rad(lat1)
    lon2 = np.deg2rad(lon2)
    lat2 = np.deg2rad(lat2)
    s = np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2)
    return np.rad2deg(np.arccos(np.clip(s, -1.0, 1.0)))


def load_filtered_radiants(radiant_h5: Path, statistics_h5: Path) -> np.ndarray:
    with h5py.File(radiant_h5, "r") as h5:
        rows = h5["radiants"][()]
    with h5py.File(statistics_h5, "r") as h5:
        fit_sample_idx = np.asarray(h5["fit_sample_idx"][()], dtype=np.int64)
    keep = np.isin(np.asarray(rows["sample_idx"], dtype=np.int64), fit_sample_idx)
    rows = rows[keep]
    good = (
        np.isfinite(rows["plot_longitude_deg"])
        & np.isfinite(rows["radiant_beta_ecliptic_deg"])
        & np.isfinite(rows["speed_km_s"])
        & (rows["speed_km_s"] > 0.0)
    )
    return rows[good]


def zenith_weight(rows: np.ndarray, min_cos_z: float) -> tuple[np.ndarray, np.ndarray]:
    sky = SkyCoord(
        ra=np.asarray(rows["radiant_ra_gcrs_deg"], dtype=np.float64) * u.deg,
        dec=np.asarray(rows["radiant_dec_gcrs_deg"], dtype=np.float64) * u.deg,
        frame="gcrs",
        obstime=Time(np.asarray(rows["epoch_unix"], dtype=np.float64), format="unix", scale="utc"),
    )
    altaz = sky.transform_to(AltAz(obstime=sky.obstime, location=PANSY_LOCATION))
    cos_z = np.sin(altaz.alt.to_value(u.rad))
    cos_z = np.asarray(cos_z, dtype=np.float64)
    return 1.0 / np.maximum(cos_z, float(min_cos_z)), cos_z


def corrected_weight(rows: np.ndarray, min_cos_z: float) -> tuple[np.ndarray, np.ndarray]:
    wz, cos_z = zenith_weight(rows, min_cos_z=min_cos_z)
    speed = np.asarray(rows["speed_km_s"], dtype=np.float64)
    wv = (73.0 / np.maximum(speed, 1.0)) ** 3
    return wz * wv, cos_z


def radiant_histogram(rows: np.ndarray, weights=None, lon_bins: int = 144, lat_bins: int = 72):
    x = np.asarray(rows["plot_longitude_deg"], dtype=np.float64)
    y = np.asarray(rows["radiant_beta_ecliptic_deg"], dtype=np.float64)
    hist, xedges, yedges = np.histogram2d(
        x,
        y,
        bins=[lon_bins, lat_bins],
        range=[[-180.0, 180.0], [-90.0, 90.0]],
        weights=weights,
    )
    return hist.T, xedges, yedges


def style_hammer(ax):
    tick_pos = np.asarray([-90.0, 0.0, 90.0])
    tick_labels = [f"{int(wrap360(PLOT_CENTER_LONGITUDE_DEG - tick))}$^\\circ$" for tick in tick_pos]
    ax.set_xticks(np.deg2rad(tick_pos))
    ax.set_xticklabels(tick_labels)
    ax.set_yticks(np.deg2rad([-60.0, -30.0, 0.0, 30.0, 60.0]))
    ax.set_yticklabels([r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$30^\circ$", r"$60^\circ$"])
    ax.grid(True, alpha=0.35, lw=0.45)


def plot_hist_panel(ax, rows, weights, title, norm, cmap="magma"):
    hist, xedges, yedges = radiant_histogram(rows, weights=weights)
    x = np.deg2rad(xedges)
    y = np.deg2rad(yedges)
    mesh = ax.pcolormesh(x, y, np.ma.masked_less_equal(hist, 0.0), shading="auto", cmap=cmap, norm=norm)
    style_hammer(ax)
    ax.set_title(title, fontsize=10)
    return mesh


def plot_all_radiants(rows: np.ndarray, weights: np.ndarray, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    raw_hist, *_ = radiant_histogram(rows)
    corr_hist, *_ = radiant_histogram(rows, weights=weights)
    raw_norm = LogNorm(vmin=1.0, vmax=max(2.0, np.nanpercentile(raw_hist[raw_hist > 0], 99.5)))
    corr_norm = LogNorm(vmin=max(1.0, np.nanpercentile(corr_hist[corr_hist > 0], 5)), vmax=np.nanpercentile(corr_hist[corr_hist > 0], 99.5))
    fig = plt.figure(figsize=(10.8, 5.4), constrained_layout=True)
    ax0 = fig.add_subplot(121, projection="hammer")
    ax1 = fig.add_subplot(122, projection="hammer")
    mesh0 = plot_hist_panel(ax0, rows, None, f"Observed high-quality radiants (N={len(rows):,})", raw_norm)
    mesh1 = plot_hist_panel(ax1, rows, weights, "Zenith-angle and velocity corrected", corr_norm)
    for ax in (ax0, ax1):
        ax.set_xlabel(r"Sun-centered ecliptic longitude, $\lambda-\lambda_\odot$")
        ax.set_ylabel(r"Ecliptic latitude, $\beta$")
    cb0 = fig.colorbar(mesh0, ax=ax0, orientation="horizontal", pad=0.10, fraction=0.045)
    cb0.set_label("Count per bin")
    cb1 = fig.colorbar(mesh1, ax=ax1, orientation="horizontal", pad=0.10, fraction=0.045)
    cb1.set_label("Weighted count per bin")
    fig.savefig(out, dpi=240)
    plt.close(fig)


def solar_window_mask(rows: np.ndarray, center_deg: float, half_width_deg: float) -> np.ndarray:
    return np.abs(wrap180(np.asarray(rows["sun_lambda_ecliptic_deg"], dtype=np.float64) - float(center_deg))) <= float(half_width_deg)


def plot_snapshots(rows: np.ndarray, out: Path, half_width_deg: float) -> None:
    fig = plt.figure(figsize=(11.0, 4.4), constrained_layout=True)
    axes = [fig.add_subplot(1, 3, i + 1, projection="hammer") for i in range(3)]
    for ax, shower in zip(axes, SHOWERS, strict=True):
        sub = rows[solar_window_mask(rows, shower.solar_lon_deg, half_width_deg)]
        hist, *_ = radiant_histogram(sub)
        norm = LogNorm(vmin=1.0, vmax=max(2.0, np.nanpercentile(hist[hist > 0], 99.4))) if np.any(hist > 0) else None
        plot_hist_panel(
            ax,
            sub,
            None,
            rf"{shower.name}, $\lambda_\odot={shower.solar_lon_deg:.1f}^\circ\pm{half_width_deg:g}^\circ$",
            norm,
        )
        ax.scatter(
            np.deg2rad(centered_plot_longitude_deg(shower.sc_lon_deg)),
            np.deg2rad(shower.beta_deg),
            marker="+",
            s=170,
            linewidth=1.5,
            color="cyan",
            zorder=10,
        )
        ax.set_xlabel(r"$\lambda-\lambda_\odot$")
    axes[0].set_ylabel(r"$\beta$")
    fig.savefig(out, dpi=240)
    plt.close(fig)


def plot_candidate_showers(rows: np.ndarray, out: Path, half_width_deg: float, radius_deg: float) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.8), constrained_layout=True)
    for ax, shower in zip(axes, SHOWERS, strict=True):
        solar_keep = solar_window_mask(rows, shower.solar_lon_deg, half_width_deg)
        sub = rows[solar_keep]
        x = wrap180(np.asarray(sub["lambda_minus_sun_deg"], dtype=np.float64) - shower.sc_lon_deg)
        y = np.asarray(sub["radiant_beta_ecliptic_deg"], dtype=np.float64)
        near = np.abs(x) <= 24.0
        near &= np.abs(y - shower.beta_deg) <= 18.0
        ax.scatter(x[near], y[near], s=2.0, c="0.55", alpha=0.45, linewidths=0)
        sep = angular_separation_deg(sub["lambda_minus_sun_deg"], sub["radiant_beta_ecliptic_deg"], shower.sc_lon_deg, shower.beta_deg)
        member = sep <= float(radius_deg)
        ax.scatter(x[member], y[member], s=8.0, c="#d62728", alpha=0.9, linewidths=0)
        ax.scatter(0.0, shower.beta_deg, marker="+", s=180, linewidth=1.6, color="black")
        circle = plt.Circle((0.0, shower.beta_deg), radius_deg, color="black", fill=False, lw=0.8, ls="--")
        ax.add_patch(circle)
        ax.set_xlim(-24.0, 24.0)
        ax.set_ylim(shower.beta_deg - 18.0, shower.beta_deg + 18.0)
        ax.set_title(f"{shower.name}\nN={shower.n}, $v_g$={shower.vg_km_s:.1f} km/s", fontsize=10)
        ax.grid(alpha=0.25, lw=0.45)
        ax.set_xlabel(r"$\Delta(\lambda-\lambda_\odot)$ (deg)")
    axes[0].set_ylabel(r"Ecliptic latitude, $\beta$ (deg)")
    fig.savefig(out, dpi=240)
    plt.close(fig)


def write_sidecar(path: Path, rows: np.ndarray, weights: np.ndarray, cos_z: np.ndarray, min_cos_z: float, half_width_deg: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h5:
        h5.attrs["selection"] = "intersection of high-quality radiant table with height-velocity fit_sample_idx gate"
        h5.attrs["zenith_weight"] = "1 / max(cos(zenith_angle), min_cos_z)"
        h5.attrs["velocity_weight"] = "(73 km/s / v_g)^3"
        h5.attrs["min_cos_z"] = float(min_cos_z)
        h5.attrs["snapshot_half_width_solar_longitude_deg"] = float(half_width_deg)
        h5.create_dataset("radiants", data=rows, compression="gzip", compression_opts=3)
        h5.create_dataset("correction_weight", data=weights.astype(np.float32), compression="gzip", compression_opts=3)
        h5.create_dataset("cos_zenith_angle", data=cos_z.astype(np.float32), compression="gzip", compression_opts=3)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radiants-h5", type=Path, required=True)
    parser.add_argument("--statistics-h5", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-cos-z", type=float, default=0.15)
    parser.add_argument("--snapshot-half-width-deg", type=float, default=4.0)
    parser.add_argument("--shower-radius-deg", type=float, default=4.0)
    args = parser.parse_args()

    rows = load_filtered_radiants(args.radiants_h5, args.statistics_h5)
    weights, cos_z = corrected_weight(rows, min_cos_z=args.min_cos_z)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_all_radiants(rows, weights, args.output_dir / "paper_radiant_distribution_corrected.png")
    plot_snapshots(rows, args.output_dir / "paper_radiant_snapshots.png", half_width_deg=args.snapshot_half_width_deg)
    plot_candidate_showers(rows, args.output_dir / "paper_candidate_showers.png", half_width_deg=args.snapshot_half_width_deg, radius_deg=args.shower_radius_deg)
    write_sidecar(args.output_dir / "paper_radiant_results.h5", rows, weights, cos_z, args.min_cos_z, args.snapshot_half_width_deg)
    print(f"paper_radiants {len(rows)}")
    print(args.output_dir)


if __name__ == "__main__":
    main()
