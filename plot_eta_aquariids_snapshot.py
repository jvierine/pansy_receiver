#!/usr/bin/env python3
"""Plot a high-resolution eta-Aquariids radiant snapshot."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np


# IAU MDC shower 31, AdNo 011: Jenniskens (2023).
ETA_SOLAR_LON_DEG = 45.7
ETA_SOLAR_LON_BEGIN_DEG = 14.0
ETA_SOLAR_LON_END_DEG = 100.0
ETA_SC_LON_DEG = 293.3
ETA_BETA_DEG = 7.9
ETA_VG_KM_S = 66.3


def wrap180(deg: np.ndarray | float) -> np.ndarray:
    return (np.asarray(deg, dtype=np.float64) + 180.0) % 360.0 - 180.0


def angular_separation_deg(lon1, lat1, lon2, lat2):
    lon1 = np.deg2rad(lon1)
    lat1 = np.deg2rad(lat1)
    lon2 = np.deg2rad(lon2)
    lat2 = np.deg2rad(lat2)
    dot = np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2)
    return np.rad2deg(np.arccos(np.clip(dot, -1.0, 1.0)))


def load_catalogue(path: Path) -> dict[str, np.ndarray]:
    fields = ("epoch_unix", "solar_longitude_deg", "sun_centered_lon_deg", "ecliptic_lat_deg", "vg_km_s")
    with h5py.File(path, "r") as h5:
        return {name: np.asarray(h5[name], dtype=np.float64) for name in fields}


def utc_days(epoch_unix: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            dt.datetime.fromtimestamp(float(t), tz=dt.timezone.utc).strftime("%Y-%m-%d")
            for t in epoch_unix
        ]
    )


def eta_candidate_mask(data: dict[str, np.ndarray], radiant_radius_deg: float, speed_half_width_kms: float) -> np.ndarray:
    sep = angular_separation_deg(
        data["sun_centered_lon_deg"],
        data["ecliptic_lat_deg"],
        ETA_SC_LON_DEG,
        ETA_BETA_DEG,
    )
    return (
        (np.abs(wrap180(data["solar_longitude_deg"] - ETA_SOLAR_LON_DEG)) <= 12.0)
        & (sep <= float(radiant_radius_deg))
        & (np.abs(data["vg_km_s"] - ETA_VG_KM_S) <= float(speed_half_width_kms))
    )


def peak_day_from_candidates(data: dict[str, np.ndarray], candidate: np.ndarray) -> tuple[str, int]:
    days = utc_days(data["epoch_unix"][candidate])
    if len(days) == 0:
        raise RuntimeError("No eta-Aquariids candidates found with the requested broad selection")
    unique, counts = np.unique(days, return_counts=True)
    order = np.argsort(counts)[::-1]
    return str(unique[order[0]]), int(counts[order[0]])


def plot_snapshot(
    data: dict[str, np.ndarray],
    label: str,
    time_mask: np.ndarray,
    output_png: Path,
    output_pdf: Path,
    output_h5: Path,
    bin_width_deg: float,
    lon_half_width_deg: float,
    beta_min_deg: float,
    beta_max_deg: float,
    show_title: bool,
) -> None:
    lon_offset = wrap180(data["sun_centered_lon_deg"] - ETA_SC_LON_DEG)
    x = ETA_SC_LON_DEG + lon_offset
    y = data["ecliptic_lat_deg"]
    zoom = time_mask & (np.abs(lon_offset) <= lon_half_width_deg) & (y >= beta_min_deg) & (y <= beta_max_deg)
    selected_solar = data["solar_longitude_deg"][time_mask]
    median_solar = float(np.nanmedian(selected_solar)) if len(selected_solar) else np.nan
    min_solar = float(np.nanmin(selected_solar)) if len(selected_solar) else np.nan
    max_solar = float(np.nanmax(selected_solar)) if len(selected_solar) else np.nan

    xedges = np.arange(ETA_SC_LON_DEG - lon_half_width_deg, ETA_SC_LON_DEG + lon_half_width_deg + bin_width_deg, bin_width_deg)
    yedges = np.arange(beta_min_deg, beta_max_deg + bin_width_deg, bin_width_deg)
    hist, _, _ = np.histogram2d(x[zoom], y[zoom], bins=(xedges, yedges))
    hist = hist.T

    output_png.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_h5, "w") as h5:
        h5.attrs["script"] = Path(__file__).name
        h5.attrs["selection_label"] = label
        h5.attrs["histogram_bin_width_deg"] = float(bin_width_deg)
        h5.attrs["eta_solar_longitude_deg"] = ETA_SOLAR_LON_DEG
        h5.attrs["selection_solar_longitude_min_deg"] = min_solar
        h5.attrs["selection_solar_longitude_max_deg"] = max_solar
        h5.attrs["eta_sun_centered_lon_deg"] = ETA_SC_LON_DEG
        h5.attrs["eta_beta_deg"] = ETA_BETA_DEG
        h5.attrs["eta_vg_km_s"] = ETA_VG_KM_S
        h5.attrs["median_solar_longitude_deg"] = median_solar
        h5.attrs["n_zoom_events"] = int(np.count_nonzero(zoom))
        h5.create_dataset("x_edges_sun_centered_lon_deg", data=xedges)
        h5.create_dataset("beta_edges_deg", data=yedges)
        h5.create_dataset("count", data=hist.astype(np.int32))
        h5.create_dataset("zoom_sample_epoch_unix", data=data["epoch_unix"][zoom])
        h5.create_dataset("zoom_sun_centered_lon_deg", data=data["sun_centered_lon_deg"][zoom])
        h5.create_dataset("zoom_beta_deg", data=y[zoom])
        h5.create_dataset("zoom_vg_km_s", data=data["vg_km_s"][zoom])

    positive = hist[hist > 0.0]
    vmax = max(1.0, float(np.nanpercentile(positive, 99.5))) if len(positive) else 1.0
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("black")
    fig, ax = plt.subplots(figsize=(6.7, 4.8), constrained_layout=True)
    mesh = ax.pcolormesh(
        xedges,
        yedges,
        np.ma.masked_where(hist <= 0, hist),
        cmap=cmap,
        norm=LogNorm(vmin=1.0, vmax=vmax),
        shading="flat",
    )
    ax.scatter(ETA_SC_LON_DEG, ETA_BETA_DEG, marker="+", s=120, color="cyan", linewidth=1.4)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(ETA_SC_LON_DEG + lon_half_width_deg, ETA_SC_LON_DEG - lon_half_width_deg)
    ax.set_ylim(beta_min_deg, beta_max_deg)
    ax.set_xlabel(r"Sun-centered ecliptic longitude, $\lambda_g-\lambda_\odot$ (deg)")
    ax.set_ylabel(r"Ecliptic latitude, $\beta_g$ (deg)")
    if show_title:
        ax.set_title(
            rf"$\eta$-Aquariids radiant histogram, {label} ($\lambda_\odot={median_solar:.1f}^\circ$)",
            fontsize=13,
        )
    ax.text(
        0.03,
        0.03,
        f"{bin_width_deg:g} deg bins; N={np.count_nonzero(zoom)}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.08)
    cbar = fig.colorbar(mesh, cax=cax)
    cbar.set_label(rf"Raw count per {bin_width_deg:g}$^\circ$ bin")
    fig.savefig(output_png, dpi=260)
    fig.savefig(output_pdf)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue-h5", type=Path, default=Path("figs/paper_refresh_20260729_current/pansy_keplerian_catalogue.h5"))
    parser.add_argument("--output-dir", type=Path, default=Path("figs/eta_aquariids_snapshot"))
    parser.add_argument("--day", default=None, help="UTC day to plot; default is the data-derived peak candidate day")
    parser.add_argument(
        "--full-activity",
        action="store_true",
        help="Use the full IAU eta-Aquariids activity interval instead of one UTC day.",
    )
    parser.add_argument("--bin-width-deg", type=float, default=0.1)
    parser.add_argument("--lon-half-width-deg", type=float, default=8.0)
    parser.add_argument("--beta-min-deg", type=float, default=ETA_BETA_DEG - 4.0)
    parser.add_argument("--beta-max-deg", type=float, default=ETA_BETA_DEG + 4.0)
    parser.add_argument("--candidate-radius-deg", type=float, default=5.0)
    parser.add_argument("--candidate-speed-half-width-kms", type=float, default=10.0)
    parser.add_argument("--no-title", action="store_true", help="Omit the plot title for use as a paper panel.")
    args = parser.parse_args()

    data = load_catalogue(args.catalogue_h5)
    if args.full_activity:
        time_mask = (data["solar_longitude_deg"] >= ETA_SOLAR_LON_BEGIN_DEG) & (
            data["solar_longitude_deg"] <= ETA_SOLAR_LON_END_DEG
        )
        label = rf"${ETA_SOLAR_LON_BEGIN_DEG:.0f}^\circ \leq \lambda_\odot \leq {ETA_SOLAR_LON_END_DEG:.0f}^\circ$"
        stem = f"eta_aquariids_radiant_full_activity_ls{ETA_SOLAR_LON_BEGIN_DEG:.0f}_{ETA_SOLAR_LON_END_DEG:.0f}_bin{args.bin_width_deg:g}deg"
    else:
        candidate = eta_candidate_mask(data, args.candidate_radius_deg, args.candidate_speed_half_width_kms)
        if args.day is None:
            day, peak_count = peak_day_from_candidates(data, candidate)
            print(f"peak_eta_candidate_day {day} broad_candidate_count {peak_count}")
        else:
            day = args.day
        days = utc_days(data["epoch_unix"])
        time_mask = days == day
        label = f"{day} UTC"
        stem = f"eta_aquariids_radiant_snapshot_{day}_bin{args.bin_width_deg:g}deg"
    output_png = args.output_dir / f"{stem}.png"
    output_pdf = output_png.with_suffix(".pdf")
    output_h5 = output_png.with_suffix(".h5")
    plot_snapshot(
        data,
        label,
        time_mask,
        output_png,
        output_pdf,
        output_h5,
        args.bin_width_deg,
        args.lon_half_width_deg,
        args.beta_min_deg,
        args.beta_max_deg,
        not args.no_title,
    )
    print(output_png)
    print(output_h5)


if __name__ == "__main__":
    main()
