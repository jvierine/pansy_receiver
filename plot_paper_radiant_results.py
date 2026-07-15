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
from scipy.optimize import minimize_scalar

from radiant_visibility import radiant_exposure_hours_grid


PANSY_LOCATION = EarthLocation(lat=-69.010833 * u.deg, lon=39.599722 * u.deg, height=100.0 * u.m)
PLOT_CENTER_LONGITUDE_DEG = -90.0
LONGITUDE_BINS = 144
LATITUDE_BINS = 72
APEX_LONGITUDE_DEG = 270.0
APEX_LONGITUDE_HALF_WIDTH_DEG = 30.0
APEX_BETA_MIN_DEG = 5.0
APEX_BETA_MAX_DEG = 35.0


@dataclass(frozen=True)
class Shower:
    name: str
    solar_lon_deg: float
    snapshot_solar_lon_deg: float | None
    snapshot_half_width_deg: float | None
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
        snapshot_solar_lon_deg=None,
        snapshot_half_width_deg=None,
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
        snapshot_solar_lon_deg=0.0,
        snapshot_half_width_deg=2.0,
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
        snapshot_solar_lon_deg=None,
        snapshot_half_width_deg=None,
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


def load_radiant_epochs(radiant_h5: Path) -> np.ndarray:
    with h5py.File(radiant_h5, "r") as h5:
        return np.asarray(h5["radiants"]["epoch_unix"], dtype=np.float64)


def zenith_cosine(rows: np.ndarray) -> np.ndarray:
    sky = SkyCoord(
        ra=np.asarray(rows["radiant_ra_gcrs_deg"], dtype=np.float64) * u.deg,
        dec=np.asarray(rows["radiant_dec_gcrs_deg"], dtype=np.float64) * u.deg,
        frame="gcrs",
        obstime=Time(np.asarray(rows["epoch_unix"], dtype=np.float64), format="unix", scale="utc"),
    )
    altaz = sky.transform_to(AltAz(obstime=sky.obstime, location=PANSY_LOCATION))
    cos_z = np.sin(altaz.alt.to_value(u.rad))
    return np.asarray(cos_z, dtype=np.float64)


def corrected_weight(rows: np.ndarray, cos_z: np.ndarray, min_cos_z: float, alpha: float) -> np.ndarray:
    speed = np.asarray(rows["speed_km_s"], dtype=np.float64)
    wv = (73.0 / np.maximum(speed, 1.0)) ** 3
    weight = wv / np.maximum(cos_z, float(min_cos_z)) ** float(alpha)
    weight[~np.isfinite(weight) | (cos_z <= 0.0)] = 0.0
    return weight


def radiant_histogram(rows: np.ndarray, weights=None, lon_bins: int = LONGITUDE_BINS, lat_bins: int = LATITUDE_BINS):
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


def merge_intervals(intervals: np.ndarray, max_gap_seconds: float = 0.1) -> np.ndarray:
    intervals = np.asarray(intervals, dtype=np.float64).reshape(-1, 2)
    good = np.all(np.isfinite(intervals), axis=1) & (intervals[:, 1] > intervals[:, 0])
    intervals = intervals[good]
    if len(intervals) == 0:
        return np.empty((0, 2), dtype=np.float64)
    intervals = intervals[np.argsort(intervals[:, 0])]
    merged = [intervals[0].copy()]
    for start, stop in intervals[1:]:
        if start <= merged[-1][1] + float(max_gap_seconds):
            merged[-1][1] = max(merged[-1][1], stop)
        else:
            merged.append(np.asarray([start, stop], dtype=np.float64))
    return np.asarray(merged, dtype=np.float64)


def load_mesomode_intervals(metadata_dir: Path, start_epoch: float, stop_epoch: float) -> np.ndarray:
    import digital_rf

    reader = digital_rf.DigitalMetadataReader(str(metadata_dir))
    start_sample = int(np.floor(float(start_epoch) * 1e6))
    stop_sample = int(np.ceil(float(stop_epoch) * 1e6))
    records = reader.read(start_sample, stop_sample, columns=["start", "end"])
    intervals = []
    for record in records.values():
        start = float(np.asarray(record["start"]).reshape(-1)[0]) / 1e6
        stop = float(np.asarray(record["end"]).reshape(-1)[0]) / 1e6 + 0.032
        start = max(start, float(start_epoch))
        stop = min(stop, float(stop_epoch))
        if stop > start:
            intervals.append((start, stop))
    return merge_intervals(np.asarray(intervals, dtype=np.float64))


def event_occupied_intervals(epoch_unix: np.ndarray, cadence_seconds: float) -> np.ndarray:
    epoch = np.asarray(epoch_unix, dtype=np.float64)
    epoch = epoch[np.isfinite(epoch)]
    if len(epoch) == 0:
        return np.empty((0, 2), dtype=np.float64)
    starts = np.unique(np.floor(epoch / float(cadence_seconds)) * float(cadence_seconds))
    return np.column_stack((starts, starts + float(cadence_seconds)))


def observing_time_samples(intervals: np.ndarray, cadence_seconds: float) -> tuple[np.ndarray, np.ndarray]:
    intervals = merge_intervals(intervals)
    if len(intervals) == 0:
        return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64)
    cadence = float(cadence_seconds)
    first_bin = int(np.floor(np.min(intervals[:, 0]) / cadence))
    last_bin = int(np.floor(np.nextafter(np.max(intervals[:, 1]), -np.inf) / cadence))
    bin_numbers = np.arange(first_bin, last_bin + 1, dtype=np.int64)
    seconds = np.zeros(len(bin_numbers), dtype=np.float64)
    for start, stop in intervals:
        i0 = int(np.floor(start / cadence)) - first_bin
        i1 = int(np.floor(np.nextafter(stop, -np.inf) / cadence)) - first_bin
        for idx in range(i0, i1 + 1):
            bin_start = float(bin_numbers[idx]) * cadence
            seconds[idx] += max(0.0, min(stop, bin_start + cadence) - max(start, bin_start))
    keep = seconds > 0.0
    epoch = (bin_numbers[keep].astype(np.float64) + 0.5) * cadence
    return epoch, seconds[keep] / 3600.0


def interpolate_solar_longitude(rows: np.ndarray, epoch_unix: np.ndarray) -> np.ndarray:
    row_epoch = np.asarray(rows["epoch_unix"], dtype=np.float64)
    row_sun = np.asarray(rows["sun_lambda_ecliptic_deg"], dtype=np.float64)
    good = np.isfinite(row_epoch) & np.isfinite(row_sun)
    order = np.argsort(row_epoch[good])
    row_epoch = row_epoch[good][order]
    row_sun = row_sun[good][order]
    unique_epoch, unique_idx = np.unique(row_epoch, return_index=True)
    unwrapped_sun = np.rad2deg(np.unwrap(np.deg2rad(row_sun[unique_idx])))
    return wrap360(np.interp(np.asarray(epoch_unix, dtype=np.float64), unique_epoch, unwrapped_sun))


def histogram_grid_centers(xedges: np.ndarray, yedges: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])
    xgrid, ygrid = np.meshgrid(xcenters, ycenters)
    physical_lon = wrap360(PLOT_CENTER_LONGITUDE_DEG - xgrid)
    return xcenters, ycenters, physical_lon, ygrid


def corrected_flux_histogram(
    rows: np.ndarray,
    cos_z: np.ndarray,
    exposure_hours: np.ndarray,
    min_cos_z: float,
    alpha: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    weights = corrected_weight(rows, cos_z, min_cos_z=min_cos_z, alpha=alpha)
    weighted_hist, _, _ = radiant_histogram(rows, weights=weights)
    flux = np.zeros_like(weighted_hist, dtype=np.float64)
    np.divide(weighted_hist, exposure_hours, out=flux, where=exposure_hours > 0.0)
    return weights, weighted_hist, flux


def apex_masks(xedges: np.ndarray, yedges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    _, _, physical_lon, beta = histogram_grid_centers(xedges, yedges)
    near_apex = np.abs(wrap180(physical_lon - APEX_LONGITUDE_DEG)) <= APEX_LONGITUDE_HALF_WIDTH_DEG
    north = near_apex & (beta >= APEX_BETA_MIN_DEG) & (beta <= APEX_BETA_MAX_DEG)
    south = near_apex & (beta <= -APEX_BETA_MIN_DEG) & (beta >= -APEX_BETA_MAX_DEG)
    return north, south


def fit_zenith_exponent(
    rows: np.ndarray,
    cos_z: np.ndarray,
    exposure_hours: np.ndarray,
    min_cos_z: float,
    alpha_bounds: tuple[float, float],
) -> tuple[float, float, float, float]:
    _, xedges, yedges = radiant_histogram(rows)
    north, south = apex_masks(xedges, yedges)
    valid_exposure = exposure_hours > 0.0

    def fluxes(alpha: float) -> tuple[float, float]:
        _, _, flux = corrected_flux_histogram(rows, cos_z, exposure_hours, min_cos_z, alpha)
        north_flux = float(np.sum(flux[north & valid_exposure]))
        south_flux = float(np.sum(flux[south & valid_exposure]))
        return north_flux, south_flux

    def objective(alpha: float) -> float:
        north_flux, south_flux = fluxes(alpha)
        if north_flux <= 0.0 or south_flux <= 0.0:
            return np.inf
        return float(np.log(north_flux / south_flux) ** 2)

    result = minimize_scalar(
        objective,
        bounds=(float(alpha_bounds[0]), float(alpha_bounds[1])),
        method="bounded",
        options={"xatol": 1e-3},
    )
    alpha = float(result.x)
    north_flux, south_flux = fluxes(alpha)
    return alpha, float(result.fun), north_flux, south_flux


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
    plot_hist = np.asarray(hist, dtype=np.float64).copy()
    floor = float(norm.vmin) if norm is not None and np.isfinite(norm.vmin) else 1.0
    plot_hist[~np.isfinite(plot_hist) | (plot_hist <= 0.0)] = floor
    mesh = ax.pcolormesh(x, y, plot_hist, shading="auto", cmap=cmap, norm=norm)
    style_hammer(ax)
    ax.set_title(title, fontsize=10)
    return mesh


def plot_grid_panel(ax, hist, xedges, yedges, title, norm, cmap="magma"):
    plot_hist = np.asarray(hist, dtype=np.float64).copy()
    floor = float(norm.vmin) if norm is not None and np.isfinite(norm.vmin) else 1.0
    plot_hist[~np.isfinite(plot_hist) | (plot_hist <= 0.0)] = floor
    mesh = ax.pcolormesh(np.deg2rad(xedges), np.deg2rad(yedges), plot_hist, shading="auto", cmap=cmap, norm=norm)
    style_hammer(ax)
    ax.set_title(title, fontsize=10)
    return mesh


def exposure_contour_levels(exposure_hours: np.ndarray) -> np.ndarray:
    positive = np.asarray(exposure_hours, dtype=np.float64)
    positive = positive[np.isfinite(positive) & (positive > 0.0)]
    if len(positive) == 0:
        return np.asarray([], dtype=np.float64)
    lo = max(1.0, float(np.nanpercentile(positive, 20.0)))
    hi = float(np.nanpercentile(positive, 95.0))
    if hi <= lo:
        return np.asarray([lo], dtype=np.float64)
    targets = np.geomspace(lo, hi, 5)
    magnitude = 10.0 ** np.floor(np.log10(targets))
    levels = np.round(targets / magnitude) * magnitude
    return np.unique(levels[(levels > 0.0) & (levels < np.nanmax(positive))])


def right_side_contour_label_positions(
    contour_set,
    raw_hist: np.ndarray,
    xcenters_deg: np.ndarray,
    ycenters_deg: np.ndarray,
) -> list[tuple[float, float]]:
    """Choose one low-density label location per contour on the map's right side."""
    xcenters_rad = np.deg2rad(np.asarray(xcenters_deg, dtype=np.float64))
    ycenters_rad = np.deg2rad(np.asarray(ycenters_deg, dtype=np.float64))
    density = np.log1p(np.maximum(np.asarray(raw_hist, dtype=np.float64), 0.0))
    density_scale = max(float(np.nanmax(density)), 1.0)
    target_x = np.deg2rad(110.0)
    right_min = np.deg2rad(45.0)
    right_max = np.deg2rad(155.0)
    positions: list[tuple[float, float]] = []

    for segments in contour_set.allsegs:
        candidates = []
        fallback = []
        for segment in segments:
            segment = np.asarray(segment, dtype=np.float64)
            if segment.ndim != 2 or segment.shape[0] == 0:
                continue
            fallback.extend(segment)
            candidates.extend(segment[(segment[:, 0] >= right_min) & (segment[:, 0] <= right_max)])
        points = np.asarray(candidates if candidates else fallback, dtype=np.float64)
        if points.size == 0:
            positions.append((target_x, 0.0))
            continue

        scores = np.empty(len(points), dtype=np.float64)
        for i, (x, y) in enumerate(points):
            ix = int(np.argmin(np.abs(xcenters_rad - x)))
            iy = int(np.argmin(np.abs(ycenters_rad - y)))
            longitude_cost = abs(x - target_x) / np.deg2rad(45.0)
            scores[i] = longitude_cost + 2.5 * density[iy, ix] / density_scale
        positions.append(tuple(points[int(np.argmin(scores))]))
    return positions


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


def add_exposure_contours(
    ax,
    exposure_hours: np.ndarray,
    raw_hist: np.ndarray,
    xedges: np.ndarray,
    yedges: np.ndarray,
) -> None:
    xcenters, ycenters, _, _ = histogram_grid_centers(xedges, yedges)
    levels = exposure_contour_levels(exposure_hours)
    if len(levels) > 0:
        contours = ax.contour(
            np.deg2rad(xcenters),
            np.deg2rad(ycenters),
            exposure_hours,
            levels=levels,
            colors="white",
            linewidths=0.75,
            linestyles="--",
            alpha=0.55,
            zorder=8,
        )
        label_positions = right_side_contour_label_positions(contours, raw_hist, xcenters, ycenters)
        ax.clabel(
            contours,
            levels=contours.levels,
            manual=label_positions,
            fmt=lambda value: f"{value:g} h",
            fontsize=6,
            colors="white",
            inline_spacing=2,
        )
    if np.any(exposure_hours <= 0.0) and np.any(exposure_hours > 0.0):
        zero_boundary = ax.contour(
            np.deg2rad(xcenters),
            np.deg2rad(ycenters),
            (exposure_hours > 0.0).astype(np.float32),
            levels=[0.5],
            colors="white",
            linewidths=0.9,
            linestyles="--",
            alpha=0.72,
            zorder=9,
        )
        zero_position = right_side_contour_label_positions(zero_boundary, raw_hist, xcenters, ycenters)[0]
        ax.clabel(
            zero_boundary,
            levels=[0.5],
            manual=[zero_position],
            fmt={0.5: "0 h"},
            fontsize=6,
            colors="white",
            inline_spacing=2,
        )


def plot_all_radiants(
    rows: np.ndarray,
    raw_hist: np.ndarray,
    flux_hist: np.ndarray,
    exposure_hours: np.ndarray,
    xedges: np.ndarray,
    yedges: np.ndarray,
    alpha: float,
    out: Path,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    raw_norm = LogNorm(vmin=1.0, vmax=max(2.0, np.nanpercentile(raw_hist[raw_hist > 0], 99.5)))
    positive_flux = flux_hist[np.isfinite(flux_hist) & (flux_hist > 0.0)]
    flux_norm = LogNorm(vmin=np.nanpercentile(positive_flux, 5.0), vmax=np.nanpercentile(positive_flux, 99.5))
    fig = plt.figure(figsize=(10.8, 5.4), constrained_layout=True)
    ax0 = fig.add_subplot(121, projection="hammer")
    ax1 = fig.add_subplot(122, projection="hammer")
    mesh0 = plot_grid_panel(ax0, raw_hist, xedges, yedges, f"Observed high-quality radiants (N={len(rows):,})", raw_norm)
    mesh1 = plot_grid_panel(
        ax1,
        flux_hist,
        xedges,
        yedges,
        rf"Exposure-corrected radiant rate ($\alpha={alpha:.2f}$)",
        flux_norm,
    )
    add_exposure_contours(ax0, exposure_hours, raw_hist, xedges, yedges)
    for ax in (ax0, ax1):
        ax.set_xticklabels([])
        add_source_markers(ax)
        ax.set_xlabel(r"Sun-centered ecliptic longitude, $\lambda-\lambda_\odot$")
        ax.set_ylabel(r"Ecliptic latitude, $\beta$")
    cb0 = fig.colorbar(mesh0, ax=ax0, orientation="horizontal", pad=0.10, fraction=0.045)
    cb0.set_label("Count per bin")
    cb1 = fig.colorbar(mesh1, ax=ax1, orientation="horizontal", pad=0.10, fraction=0.045)
    cb1.set_label(r"Debiased radiant rate (h$^{-1}$ per bin)")
    fig.savefig(out, dpi=240)
    plt.close(fig)


def solar_window_mask(rows: np.ndarray, center_deg: float, half_width_deg: float) -> np.ndarray:
    return np.abs(wrap180(np.asarray(rows["sun_lambda_ecliptic_deg"], dtype=np.float64) - float(center_deg))) <= float(half_width_deg)


def solar_longitude_window_mask(solar_longitude_deg: np.ndarray, center_deg: float, half_width_deg: float) -> np.ndarray:
    return np.abs(wrap180(np.asarray(solar_longitude_deg, dtype=np.float64) - float(center_deg))) <= float(half_width_deg)


def plot_snapshots(
    rows: np.ndarray,
    observation_epoch: np.ndarray,
    observation_sun: np.ndarray,
    observation_hours: np.ndarray,
    xedges: np.ndarray,
    yedges: np.ndarray,
    out: Path,
    half_width_deg: float,
) -> None:
    fig = plt.figure(figsize=(11.0, 4.4), constrained_layout=True)
    axes = [fig.add_subplot(1, 3, i + 1, projection="hammer") for i in range(3)]
    xcenters, ycenters, _, _ = histogram_grid_centers(xedges, yedges)
    for ax, shower in zip(axes, SHOWERS, strict=True):
        window_center = shower.snapshot_solar_lon_deg if shower.snapshot_solar_lon_deg is not None else shower.solar_lon_deg
        window_half_width = shower.snapshot_half_width_deg if shower.snapshot_half_width_deg is not None else half_width_deg
        sub = rows[solar_window_mask(rows, window_center, window_half_width)]
        obs_keep = solar_longitude_window_mask(observation_sun, window_center, window_half_width)
        if not np.any(obs_keep):
            raise RuntimeError(f"No observing-time samples found for {shower.name} snapshot")
        _, _, exposure_hours = radiant_exposure_hours_grid(
            observation_epoch[obs_keep],
            observation_sun[obs_keep],
            observation_hours[obs_keep],
            plot_longitude_deg=xcenters,
            beta_deg=ycenters,
        )
        hist, *_ = radiant_histogram(sub)
        norm = LogNorm(vmin=1.0, vmax=max(2.0, np.nanpercentile(hist[hist > 0], 99.4))) if np.any(hist > 0) else None
        plot_hist_panel(
            ax,
            sub,
            None,
            rf"{shower.name}, $\lambda_\odot={window_center:.1f}^\circ\pm{window_half_width:g}^\circ$",
            norm,
        )
        add_exposure_contours(ax, exposure_hours, hist, xedges, yedges)
        add_source_markers(ax)
        ax.scatter(
            np.deg2rad(centered_plot_longitude_deg(shower.sc_lon_deg)),
            np.deg2rad(shower.beta_deg),
            marker="o",
            s=520,
            linewidth=1.6,
            facecolors="none",
            edgecolors="cyan",
            alpha=0.5,
            zorder=10,
        )
        ax.set_xticklabels([])
        ax.set_xlabel(r"$\lambda-\lambda_\odot$")
    axes[0].set_ylabel(r"$\beta$")
    fig.savefig(out, dpi=240)
    plt.close(fig)


def plot_candidate_showers(rows: np.ndarray, out: Path, half_width_deg: float, radius_deg: float) -> None:
    longitude_zoom_deg = 32.0
    latitude_zoom_deg = 24.0
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.8), constrained_layout=True)
    for ax, shower in zip(axes, SHOWERS, strict=True):
        solar_keep = solar_window_mask(rows, shower.solar_lon_deg, half_width_deg)
        sub = rows[solar_keep]
        x = wrap180(np.asarray(sub["lambda_minus_sun_deg"], dtype=np.float64) - shower.sc_lon_deg)
        y = np.asarray(sub["radiant_beta_ecliptic_deg"], dtype=np.float64)
        near = np.abs(x) <= longitude_zoom_deg
        near &= np.abs(y - shower.beta_deg) <= latitude_zoom_deg
        ax.scatter(x[near], y[near], s=8.0, c="0.55", alpha=0.45, linewidths=0)
        sep = angular_separation_deg(sub["lambda_minus_sun_deg"], sub["radiant_beta_ecliptic_deg"], shower.sc_lon_deg, shower.beta_deg)
        member = sep <= float(radius_deg)
        ax.scatter(x[member], y[member], s=32.0, c="#d62728", alpha=0.9, linewidths=0)
        ax.scatter(0.0, shower.beta_deg, marker="+", s=180, linewidth=1.6, color="black")
        circle = plt.Circle((0.0, shower.beta_deg), radius_deg, color="black", fill=False, lw=0.8, ls="--")
        ax.add_patch(circle)
        ax.set_xlim(longitude_zoom_deg, -longitude_zoom_deg)
        ax.set_ylim(shower.beta_deg - latitude_zoom_deg, shower.beta_deg + latitude_zoom_deg)
        ax.set_title(f"{shower.name}\nN={shower.n}, $v_g$={shower.vg_km_s:.1f} km/s", fontsize=10)
        ax.grid(alpha=0.25, lw=0.45)
        ax.set_xlabel(r"$\Delta(\lambda-\lambda_\odot)$ (deg)")
    axes[0].set_ylabel(r"Ecliptic latitude, $\beta$ (deg)")
    fig.savefig(out, dpi=240)
    plt.close(fig)


def write_sidecar(
    path: Path,
    rows: np.ndarray,
    weights: np.ndarray,
    cos_z: np.ndarray,
    raw_hist: np.ndarray,
    weighted_hist: np.ndarray,
    flux_hist: np.ndarray,
    exposure_hours: np.ndarray,
    xedges: np.ndarray,
    yedges: np.ndarray,
    observation_epoch: np.ndarray,
    observation_hours: np.ndarray,
    exposure_source: str,
    exposure_cadence_minutes: float,
    alpha: float,
    alpha_objective: float,
    apex_north_flux: float,
    apex_south_flux: float,
    min_cos_z: float,
    half_width_deg: float,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h5:
        h5.attrs["selection"] = "intersection of high-quality radiant table with height-velocity fit_sample_idx gate"
        h5.attrs["zenith_weight"] = "1 / max(cos(zenith_angle), min_cos_z)^alpha"
        h5.attrs["velocity_weight"] = "(73 km/s / v_g)^3"
        h5.attrs["exposure_weight"] = "1 / radiant_visibility_hours"
        h5.attrs["exposure_source"] = exposure_source
        h5.attrs["exposure_sampling_minutes"] = np.float32(exposure_cadence_minutes)
        h5.attrs["fitted_zenith_exponent_alpha"] = float(alpha)
        h5.attrs["alpha_fit_objective_log_ratio_squared"] = float(alpha_objective)
        h5.attrs["apex_north_flux_h_inv"] = float(apex_north_flux)
        h5.attrs["apex_south_flux_h_inv"] = float(apex_south_flux)
        h5.attrs["apex_longitude_center_deg"] = APEX_LONGITUDE_DEG
        h5.attrs["apex_longitude_half_width_deg"] = APEX_LONGITUDE_HALF_WIDTH_DEG
        h5.attrs["apex_absolute_beta_range_deg"] = np.asarray([APEX_BETA_MIN_DEG, APEX_BETA_MAX_DEG], dtype=np.float32)
        h5.attrs["min_cos_z"] = float(min_cos_z)
        h5.attrs["snapshot_half_width_solar_longitude_deg"] = float(half_width_deg)
        h5.create_dataset("radiants", data=rows, compression="gzip", compression_opts=3)
        h5.create_dataset("correction_weight", data=weights.astype(np.float32), compression="gzip", compression_opts=3)
        h5.create_dataset("cos_zenith_angle", data=cos_z.astype(np.float32), compression="gzip", compression_opts=3)
        h5.create_dataset("raw_count", data=raw_hist.astype(np.float32), compression="gzip", compression_opts=3)
        h5.create_dataset("weighted_count", data=weighted_hist.astype(np.float32), compression="gzip", compression_opts=3)
        h5.create_dataset("debiased_rate_h_inv", data=flux_hist.astype(np.float32), compression="gzip", compression_opts=3)
        h5.create_dataset("radiant_exposure_hours", data=exposure_hours.astype(np.float32), compression="gzip", compression_opts=3)
        h5.create_dataset("plot_longitude_edges_deg", data=xedges.astype(np.float32))
        h5.create_dataset("ecliptic_latitude_edges_deg", data=yedges.astype(np.float32))
        h5.create_dataset("observation_sample_idx", data=np.rint(observation_epoch * 1e6).astype(np.int64), compression="gzip", compression_opts=3)
        h5.attrs["observation_sample_rate_hz"] = np.int32(1_000_000)
        h5.create_dataset("observation_hours", data=observation_hours.astype(np.float32), compression="gzip", compression_opts=3)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radiants-h5", type=Path, required=True)
    parser.add_argument("--statistics-h5", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-cos-z", type=float, default=0.15)
    parser.add_argument("--mesomode-metadata-dir", type=Path)
    parser.add_argument("--exposure-cadence-minutes", type=float, default=15.0)
    parser.add_argument("--alpha-min", type=float, default=0.0)
    parser.add_argument("--alpha-max", type=float, default=4.0)
    parser.add_argument("--snapshot-half-width-deg", type=float, default=4.0)
    parser.add_argument("--shower-radius-deg", type=float, default=4.0)
    args = parser.parse_args()

    rows = load_filtered_radiants(args.radiants_h5, args.statistics_h5)
    cos_z = zenith_cosine(rows)
    cadence_seconds = float(args.exposure_cadence_minutes) * 60.0
    start_epoch = float(np.nanmin(rows["epoch_unix"]))
    stop_epoch = float(np.nanmax(rows["epoch_unix"]))
    if args.mesomode_metadata_dir is not None:
        intervals = load_mesomode_intervals(args.mesomode_metadata_dir, start_epoch, stop_epoch)
        exposure_source = f"mesomode metadata: {args.mesomode_metadata_dir}"
    else:
        intervals = event_occupied_intervals(load_radiant_epochs(args.radiants_h5), cadence_seconds)
        exposure_source = "15-minute bins containing at least one meteor in the input radiant table"
    observation_epoch, observation_hours = observing_time_samples(intervals, cadence_seconds)
    if len(observation_epoch) == 0:
        raise RuntimeError("No observing-time samples were found for the radiant exposure calculation")
    observation_sun = interpolate_solar_longitude(rows, observation_epoch)
    _, xedges, yedges = radiant_histogram(rows)
    xcenters, ycenters, _, _ = histogram_grid_centers(xedges, yedges)
    _, _, exposure_hours = radiant_exposure_hours_grid(
        observation_epoch,
        observation_sun,
        observation_hours,
        plot_longitude_deg=xcenters,
        beta_deg=ycenters,
    )
    alpha, alpha_objective, apex_north_flux, apex_south_flux = fit_zenith_exponent(
        rows,
        cos_z,
        exposure_hours,
        min_cos_z=args.min_cos_z,
        alpha_bounds=(args.alpha_min, args.alpha_max),
    )
    weights, weighted_hist, flux_hist = corrected_flux_histogram(rows, cos_z, exposure_hours, args.min_cos_z, alpha)
    raw_hist, _, _ = radiant_histogram(rows)
    zero_exposure_count = int(np.sum(raw_hist[exposure_hours <= 0.0]))
    raw_hist = np.where(exposure_hours > 0.0, raw_hist, 0.0)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_all_radiants(
        rows,
        raw_hist,
        flux_hist,
        exposure_hours,
        xedges,
        yedges,
        alpha,
        args.output_dir / "paper_radiant_distribution_corrected.png",
    )
    plot_snapshots(
        rows,
        observation_epoch,
        observation_sun,
        observation_hours,
        xedges,
        yedges,
        args.output_dir / "paper_radiant_snapshots.png",
        half_width_deg=args.snapshot_half_width_deg,
    )
    plot_candidate_showers(rows, args.output_dir / "paper_candidate_showers.png", half_width_deg=args.snapshot_half_width_deg, radius_deg=args.shower_radius_deg)
    write_sidecar(
        args.output_dir / "paper_radiant_results.h5",
        rows,
        weights,
        cos_z,
        raw_hist,
        weighted_hist,
        flux_hist,
        exposure_hours,
        xedges,
        yedges,
        observation_epoch,
        observation_hours,
        exposure_source,
        args.exposure_cadence_minutes,
        alpha,
        alpha_objective,
        apex_north_flux,
        apex_south_flux,
        args.min_cos_z,
        args.snapshot_half_width_deg,
    )
    print(f"paper_radiants {len(rows)}")
    print(f"exposure_samples {len(observation_epoch)} total_mesomode_hours {np.sum(observation_hours):.3f}")
    print(f"detections_in_zero_exposure_bins {zero_exposure_count}")
    print(f"alpha {alpha:.4f} apex_north_flux {apex_north_flux:.6g} apex_south_flux {apex_south_flux:.6g}")
    print(args.output_dir)


if __name__ == "__main__":
    main()
