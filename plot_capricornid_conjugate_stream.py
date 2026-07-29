#!/usr/bin/env python3
"""Plot the Capricornid stream candidate and its opposite-node passage."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, replace
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import FuncFormatter

from plot_omega_eridanids_shower import (
    DEFAULT_EXPOSURE,
    ang2pix_ring,
    interpolate_observation_solar_longitude,
)
from radiant_visibility import (
    centered_plot_longitude_deg,
    radiant_altitude_deg,
    radiant_exposure_hours_points,
)
from shower_selection_windows import WINDOWS, ShowerSelectionWindow, format_window, selection_mask


RADVIEW_ROOT = Path.home() / "src" / "radview"
LIVE_RADVIEW_DATA = Path("/tmp/radiantviz_live_data")
RADVIEW_DATA = LIVE_RADVIEW_DATA if LIVE_RADVIEW_DATA.exists() else RADVIEW_ROOT / "public" / "data"
DEFAULT_CATALOGUE = (
    Path(__file__).resolve().parent
    / "figs"
    / "paper_refresh_20260727_current"
    / "pansy_keplerian_catalogue.h5"
)
DEFAULT_ZENITH_SIDECAR = (
    Path(__file__).resolve().parent
    / "figs"
    / "paper_refresh_20260721_current"
    / "paper_radiant_results.h5"
)
CLUSTER_SOLAR_WINDOW_DEG = 14.0
CLUSTER_VG_RANGE = (14.26, 30.24)
CLUSTER_E_RANGE = (0.664, 0.832)
ORBIT_COLORS = ("C0", "C1", "C2")
COMET_COLOR = "black"
RADIANT_COLOR_SPECS = {
    "vg": (r"Geocentric speed, $v_g$ (km s$^{-1}$)", 10.0, 75.0, "viridis"),
    "inclination": (r"Inclination, $i$ (deg)", 0.0, 20.0, "viridis"),
    "eccentricity": (r"Eccentricity, $e$", 0.65, 0.85, "viridis"),
    "q": (r"Perihelion distance, $q$ (AU)", 0.4, 0.9, "viridis"),
    "a": (r"Semimajor axis, $a$ (AU)", 1.0, 5.0, "viridis"),
    "omega": (r"Argument of perihelion, $\omega$ (deg)", 240.0, 290.0, "viridis"),
}


@dataclass(frozen=True)
class Passage:
    name: str
    solar_lon_deg: float
    sun_centered_lon_deg: float
    beta_deg: float
    vg_km_s: float
    n: int
    mean_kepler: tuple[float, float, float, float, float, float, float] | None = None


@dataclass(frozen=True)
class ActivitySelection:
    short_name: str
    solar_range_deg: tuple[float, float]
    peak_solar_lon_deg: float
    peak_window_deg: float
    sun_centered_lon_deg: float
    beta_deg: float
    vg_range: tuple[float, float] | None
    e_range: tuple[float, float] | None
    healpix_pixels: tuple[int, ...]
    parameter_window: ShowerSelectionWindow | None = None
    healpix_nside: int = 32
    bin_width_deg: float = 1.0
    profile_window_deg: float = 1.0
    minimum_exposure_hours: float = 3.0
    inset_xlim_deg: tuple[float, float] | None = None
    inset_ylim_deg: tuple[float, float] | None = None
    inset_xticks_deg: tuple[float, ...] | None = None
    inset_yticks_deg: tuple[float, ...] | None = None
    inset_solar_lon_deg: float | None = None
    inset_solar_window_deg: float | None = None
    inset_sun_centered_lon_deg: float | None = None
    inset_beta_deg: float | None = None


PASSAGES = (
    Passage(
        name="Daytime Capricornids-Sagittariids (DCS)",
        solar_lon_deg=313.0,
        sun_centered_lon_deg=354.51,
        beta_deg=-8.34,
        vg_km_s=25.02,
        n=19,
    ),
    Passage(
        name=r"$\alpha$ Capricornids (CAP)",
        solar_lon_deg=133.0,
        sun_centered_lon_deg=177.48,
        beta_deg=10.88,
        vg_km_s=21.36,
        n=89,
        mean_kepler=(2.5256, 0.7476, 7.58, 131.47, 264.41, 274.72, 0.6231),
    ),
    Passage(
        name=r"62-Sagittariids (OES)",
        solar_lon_deg=294.2,
        sun_centered_lon_deg=2.59,
        beta_deg=-10.51,
        vg_km_s=21.68,
        n=64,
        mean_kepler=(2.5754, 0.7533, 7.21, 113.06, 277.52, 81.53, 0.6169),
    ),
)

COMET_169P_NEAT = np.asarray([2.604, 0.76796, 11.285, 176.04, 218.13, np.nan, 0.604], dtype=np.float64)
DCS_PROFILE_SOLAR_RANGE_DEG = (300.0, 330.0)
DCS_PROFILE_BIN_WIDTH_DEG = 0.25
DCS_ACTIVITY_SOLAR_RANGE_DEG = (285.0, 335.0)
DCS_ACTIVITY_BIN_WIDTH_DEG = 1.0
DCS_ACTIVITY_HEALPIX_NSIDE = 32
DCS_ACTIVITY_HEALPIX_PIXELS = (6846, 6973, 6974, 7101, 7102, 7103, 7229, 7230)
DCS_ACTIVITY_MEAN_SC_LON_DEG = -5.51
DCS_ACTIVITY_MEAN_BETA_DEG = -8.31
ACTIVITY_SELECTIONS = (
    ActivitySelection(
        short_name="DCS",
        solar_range_deg=DCS_ACTIVITY_SOLAR_RANGE_DEG,
        peak_solar_lon_deg=312.53,
        peak_window_deg=4.0,
        sun_centered_lon_deg=354.07,
        beta_deg=-8.23,
        vg_range=None,
        e_range=None,
        healpix_pixels=DCS_ACTIVITY_HEALPIX_PIXELS,
        parameter_window=WINDOWS["DCS"],
        inset_xlim_deg=(370.0, 345.0),
        inset_ylim_deg=(-14.0, 0.0),
        inset_xticks_deg=(370.0, 365.0, 360.0, 355.0, 350.0, 345.0),
        inset_yticks_deg=(-14.0, -10.0, -5.0, 0.0),
        inset_solar_lon_deg=313.0,
        inset_solar_window_deg=14.0,
        inset_sun_centered_lon_deg=354.51,
        inset_beta_deg=-8.34,
    ),
    ActivitySelection(
        short_name="OES",
        solar_range_deg=DCS_ACTIVITY_SOLAR_RANGE_DEG,
        peak_solar_lon_deg=294.2,
        peak_window_deg=4.0,
        sun_centered_lon_deg=1.63,
        beta_deg=-10.30,
        vg_range=None,
        e_range=None,
        healpix_pixels=(7104, 7105, 7232, 7233, 7234, 7360, 7361),
        parameter_window=WINDOWS["OES"],
    ),
    ActivitySelection(
        short_name="CAP",
        solar_range_deg=(80.0, 150.0),
        peak_solar_lon_deg=116.95,
        peak_window_deg=4.0,
        sun_centered_lon_deg=181.25,
        beta_deg=9.55,
        vg_range=(15.84, 28.74),
        e_range=(0.620, 0.877),
        healpix_pixels=(
            4862,
            4863,
            4864,
            4990,
            4991,
            4992,
            4993,
            5118,
            5119,
            5120,
            5121,
            5247,
            5248,
            5249,
            5375,
            5376,
        ),
        parameter_window=WINDOWS["CAP"],
        inset_xlim_deg=(185.5, 168.0),
        inset_ylim_deg=(5.0, 20.0),
        inset_xticks_deg=(185.0, 180.0, 175.0, 170.0),
        inset_yticks_deg=(5.0, 10.0, 15.0, 20.0),
        inset_solar_lon_deg=133.0,
        inset_solar_window_deg=14.0,
        inset_sun_centered_lon_deg=177.48,
        inset_beta_deg=10.88,
    ),
)


def wrap180(deg):
    return (np.asarray(deg, dtype=np.float64) + 180.0) % 360.0 - 180.0


def wrap360(deg):
    return np.asarray(deg, dtype=np.float64) % 360.0


def angular_separation_deg(lon1, lat1, lon2, lat2):
    lon1 = np.deg2rad(lon1)
    lat1 = np.deg2rad(lat1)
    lon2 = np.deg2rad(lon2)
    lat2 = np.deg2rad(lat2)
    s = np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2)
    return np.rad2deg(np.arccos(np.clip(s, -1.0, 1.0)))


def decode_half_column(buffer: bytes, col: dict, count: int) -> np.ndarray:
    width = int(col.get("width", 1))
    start = int(col.get("offset", 0))
    nbytes = int(col.get("bytes", 0))
    raw = np.frombuffer(buffer[start : start + nbytes], dtype="<f2")
    return raw.astype(np.float64).reshape(count, width)


def load_chunk_rows(data_dir: Path, dataset_id: str, center: float, half_width: float) -> np.ndarray:
    chunk_dir = data_dir / f"{dataset_id}_chunks"
    with (chunk_dir / "chunks.json").open("r") as fh:
        manifest = json.load(fh)
    needed = chunk_indexes(center, half_width)
    rows = []
    dtype = [
        ("dataset", "U8"),
        ("solar_lon", "f8"),
        ("sun_centered_lon", "f8"),
        ("ecliptic_lon", "f8"),
        ("beta", "f8"),
        ("vg", "f8"),
        ("epoch", "f8"),
        ("kepler", "f8", (7,)),
    ]
    for idx in needed:
        chunk = manifest["chunks"][idx]
        if int(chunk.get("count", 0)) <= 0:
            continue
        buffer = (chunk_dir / chunk["url"]).read_bytes()
        count = int(chunk["count"])
        cols = {col["name"]: decode_half_column(buffer, col, count) for col in chunk["columns"]}
        solar = cols["solarLongitude"][:, 0] % 360.0
        keep = np.abs(wrap180(solar - center)) <= half_width
        if not np.any(keep):
            continue
        ecl_lon = cols["eclipticLon"][:, 0] % 360.0
        beta = cols["eclipticLat"][:, 0]
        vg = cols["vg"][:, 0]
        time_days = cols["timeDays"][:, 0]
        kepler = cols["kepler"]
        out = np.empty(int(np.sum(keep)), dtype=dtype)
        out["dataset"] = dataset_id.upper()
        out["solar_lon"] = solar[keep]
        out["sun_centered_lon"] = wrap180(ecl_lon[keep] - solar[keep])
        out["ecliptic_lon"] = ecl_lon[keep]
        out["beta"] = beta[keep]
        out["vg"] = vg[keep]
        out["epoch"] = float(chunk.get("time0Unix", 0.0)) + time_days[keep] * 86400.0
        out["kepler"] = kepler[keep, :]
        rows.append(out)
    return np.concatenate(rows) if rows else np.empty(0, dtype=dtype)


def load_catalogue_rows(
    path: Path,
    center: float | None = None,
    half_width: float | None = None,
    source: str | None = None,
) -> np.ndarray:
    """Load rows from the combined HDF5 catalogue, optionally within a solar slice."""
    dtype = [
        ("dataset", "U8"),
        ("solar_lon", "f8"),
        ("sun_centered_lon", "f8"),
        ("ecliptic_lon", "f8"),
        ("beta", "f8"),
        ("vg", "f8"),
        ("epoch", "f8"),
        ("kepler", "f8", (7,)),
    ]
    with h5py.File(path, "r") as h5:
        solar = np.asarray(h5["solar_longitude_deg"], dtype=np.float64)
        source_id = np.asarray(h5["source_id"], dtype=np.int8)
        if (center is None) != (half_width is None):
            raise ValueError("center and half_width must either both be set or both be omitted")
        keep = np.ones(solar.shape, dtype=bool)
        if center is not None and half_width is not None:
            keep &= np.abs(wrap180(solar - center)) <= half_width
        if source is not None:
            source_lookup = {}
            for value in np.unique(source_id):
                attribute = f"source_id_{value}"
                if attribute in h5.attrs:
                    source_lookup[str(h5.attrs[attribute]).upper()] = int(value)
                elif len(np.unique(source_id)) == 1 and int(value) == 0:
                    source_lookup["PANSY"] = 0
            requested = source.upper()
            if requested not in source_lookup:
                raise ValueError(f"source {source!r} is not present in {path}")
            keep &= source_id == source_lookup[requested]
        index = np.flatnonzero(keep)
        out = np.empty(len(index), dtype=dtype)
        if len(index) == 0:
            return out
        source_names = {}
        for value in np.unique(source_id[index]):
            attribute = f"source_id_{value}"
            source_names[int(value)] = (
                str(h5.attrs[attribute]).upper()
                if attribute in h5.attrs
                else ("PANSY" if int(value) == 0 else f"SOURCE{int(value)}")
            )
        out["dataset"] = np.asarray([source_names[int(value)] for value in source_id[index]])
        out["solar_lon"] = solar[index]
        # Contiguous HDF5 reads followed by a NumPy mask are much faster than
        # million-row h5py fancy-index reads when the full catalogue is used.
        out["sun_centered_lon"] = np.asarray(h5["sun_centered_lon_deg"], dtype=np.float64)[keep]
        out["ecliptic_lon"] = np.asarray(h5["ecliptic_lon_deg"], dtype=np.float64)[keep]
        out["beta"] = np.asarray(h5["ecliptic_lat_deg"], dtype=np.float64)[keep]
        out["vg"] = np.asarray(h5["vg_km_s"], dtype=np.float64)[keep]
        out["epoch"] = np.asarray(h5["epoch_unix"], dtype=np.float64)[keep]
        out["kepler"] = np.asarray(h5["kepler"], dtype=np.float64)[keep, :]
    return out


def load_maarsy_h5_rows(path: Path, center: float, half_width: float) -> np.ndarray:
    import h5py

    dtype = [
        ("dataset", "U8"),
        ("solar_lon", "f8"),
        ("sun_centered_lon", "f8"),
        ("ecliptic_lon", "f8"),
        ("beta", "f8"),
        ("vg", "f8"),
        ("epoch", "f8"),
        ("kepler", "f8", (7,)),
    ]
    if not path.exists():
        return np.empty(0, dtype=dtype)
    with h5py.File(path, "r") as h:
        k = np.asarray(h["kepler"], dtype=np.float64)
        epoch = np.asarray(h["kepler_epoch_unix_second"], dtype=np.float64)
    if k.shape[1] < 6:
        return np.empty(0, dtype=dtype)

    a = k[:, 0]
    e = k[:, 1]
    inc = np.deg2rad(k[:, 2])
    argp = np.deg2rad(k[:, 3])
    raan = np.deg2rad(k[:, 4])
    nu = np.deg2rad(k[:, 5])
    good = np.isfinite(a) & np.isfinite(e) & np.isfinite(epoch)
    good &= np.isfinite(inc) & np.isfinite(argp) & np.isfinite(raan) & np.isfinite(nu)
    good &= (np.abs(a) > 1e-8) & (np.abs(e) < 10.0)
    if not np.any(good):
        return np.empty(0, dtype=dtype)
    a, e, inc, argp, raan, nu, epoch = [x[good] for x in (a, e, inc, argp, raan, nu, epoch)]

    solar = solar_longitude_approx_deg(epoch)
    keep = np.abs(wrap180(solar - center)) <= half_width
    if not np.any(keep):
        return np.empty(0, dtype=dtype)

    a, e, inc, argp, raan, nu, epoch, solar = [x[keep] for x in (a, e, inc, argp, raan, nu, epoch, solar)]
    mu = 0.01720209895**2
    p = a * (1.0 - e * e)
    speed_factor = np.sqrt(mu / np.maximum(np.abs(p), 1e-8))
    vx_p = -speed_factor * np.sin(nu)
    vy_p = speed_factor * (e + np.cos(nu))

    cos_o, sin_o = np.cos(raan), np.sin(raan)
    cos_i, sin_i = np.cos(inc), np.sin(inc)
    cos_w, sin_w = np.cos(argp), np.sin(argp)
    x1 = cos_w * vx_p - sin_w * vy_p
    y1 = sin_w * vx_p + cos_w * vy_p
    vx = cos_o * x1 - sin_o * cos_i * y1
    vy = sin_o * x1 + cos_o * cos_i * y1
    vz = sin_i * y1
    v_met = np.column_stack((vx, vy, vz))

    earth_lon = np.deg2rad((solar + 180.0) % 360.0)
    earth_speed = np.sqrt(mu)
    v_earth = np.column_stack((-earth_speed * np.sin(earth_lon), earth_speed * np.cos(earth_lon), np.zeros_like(earth_lon)))
    v_geo = v_met - v_earth
    vg = np.linalg.norm(v_geo, axis=1) * 149597870.7 / 86400.0
    radiant = -v_geo / np.maximum(np.linalg.norm(v_geo, axis=1)[:, None], 1e-12)
    ecliptic_lon = np.rad2deg(np.arctan2(radiant[:, 1], radiant[:, 0])) % 360.0
    beta = np.rad2deg(np.arcsin(np.clip(radiant[:, 2], -1.0, 1.0)))
    q = a * (1.0 - e)
    kepler = np.column_stack((a, e, np.rad2deg(inc), np.rad2deg(raan), np.rad2deg(argp), np.rad2deg(nu), q))

    finite = np.isfinite(vg) & np.isfinite(ecliptic_lon) & np.isfinite(beta) & (vg > 0.0) & (vg < 120.0)
    out = np.empty(int(np.sum(finite)), dtype=dtype)
    out["dataset"] = "MAARSY"
    out["solar_lon"] = solar[finite]
    out["sun_centered_lon"] = wrap180(ecliptic_lon[finite] - solar[finite])
    out["ecliptic_lon"] = ecliptic_lon[finite]
    out["beta"] = beta[finite]
    out["vg"] = vg[finite]
    out["epoch"] = epoch[finite]
    out["kepler"] = kepler[finite, :]
    return out


def solar_longitude_approx_deg(epoch_unix: np.ndarray) -> np.ndarray:
    jd = np.asarray(epoch_unix, dtype=np.float64) / 86400.0 + 2440587.5
    n = jd - 2451545.0
    mean_lon = (280.46646 + 0.98564736 * n) % 360.0
    mean_anom = np.deg2rad((357.52911 + 0.98560028 * n) % 360.0)
    center = 1.914602 * np.sin(mean_anom) + 0.019993 * np.sin(2 * mean_anom) + 0.000289 * np.sin(3 * mean_anom)
    return (mean_lon + center) % 360.0


def chunk_indexes(center: float, half_width: float) -> list[int]:
    start = int(np.floor(center - half_width))
    stop = int(np.floor(center + half_width))
    return sorted({i % 360 for i in range(start, stop + 1)})


def load_passage_rows(
    data_dir: Path,
    passage: Passage,
    solar_half_width: float,
    catalogue_path: Path | None = None,
    source: str | None = None,
) -> np.ndarray:
    if catalogue_path is not None and catalogue_path.exists():
        return load_catalogue_rows(
            catalogue_path,
            passage.solar_lon_deg,
            solar_half_width,
            source=source,
        )
    datasets = ("pansy", "maarsy") if source is None else (source.lower(),)
    parts = [
        load_chunk_rows(data_dir, dataset, passage.solar_lon_deg, solar_half_width)
        for dataset in datasets
    ]
    return np.concatenate([p for p in parts if len(p)]) if any(len(p) for p in parts) else parts[0]


def validate_passage_sources(
    rows: np.ndarray,
    passage: Passage,
    required_sources: tuple[str, ...] = ("PANSY", "MAARSY"),
) -> None:
    """Reject incomplete exports before they can silently alter a paper figure."""
    present = set(np.unique(rows["dataset"]))
    missing = [source for source in required_sources if source not in present]
    if missing:
        source_counts = {
            source: int(np.sum(rows["dataset"] == source)) for source in required_sources
        }
        raise RuntimeError(
            f"incomplete catalogue coverage for {passage.name} at "
            f"solar longitude {passage.solar_lon_deg:.1f} +/- "
            f"{CLUSTER_SOLAR_WINDOW_DEG / 2.0:.1f} deg: "
            f"missing {', '.join(missing)}; counts={source_counts}"
        )


def select_associated(rows: np.ndarray, passage: Passage, radius_deg: float, velocity_half_width: float) -> np.ndarray:
    sep = angular_separation_deg(rows["sun_centered_lon"], rows["beta"], passage.sun_centered_lon_deg, passage.beta_deg)
    keep = sep <= radius_deg
    keep &= np.abs(rows["vg"] - passage.vg_km_s) <= velocity_half_width
    keep &= np.isfinite(rows["kepler"][:, 0])
    keep &= (rows["kepler"][:, 0] > 0.0) & (rows["kepler"][:, 1] >= 0.0) & (rows["kepler"][:, 1] < 1.0)
    return rows[keep]


def cluster_filter(rows: np.ndarray, vg_range: tuple[float, float], e_range: tuple[float, float]) -> np.ndarray:
    e = rows["kepler"][:, 1]
    keep = np.isfinite(rows["vg"]) & np.isfinite(e)
    keep &= (rows["vg"] >= vg_range[0]) & (rows["vg"] <= vg_range[1])
    keep &= (e >= e_range[0]) & (e <= e_range[1])
    return rows[keep]


def ecliptic_to_equatorial_deg(lon_deg: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Convert J2000 ecliptic radiant coordinates to equatorial coordinates."""
    lon = np.deg2rad(np.asarray(lon_deg, dtype=np.float64))
    lat = np.deg2rad(np.asarray(lat_deg, dtype=np.float64))
    obliquity = np.deg2rad(23.4392911)
    cos_lat = np.cos(lat)
    x = cos_lat * np.cos(lon)
    y = cos_lat * np.sin(lon) * np.cos(obliquity) - np.sin(lat) * np.sin(obliquity)
    z = cos_lat * np.sin(lon) * np.sin(obliquity) + np.sin(lat) * np.cos(obliquity)
    ra_deg = wrap360(np.rad2deg(np.arctan2(y, x)))
    dec_deg = np.rad2deg(np.arcsin(np.clip(z, -1.0, 1.0)))
    return ra_deg, dec_deg


def plot_dcs_solar_longitude_profile(
    rows: np.ndarray,
    out: Path,
    solar_range_deg: tuple[float, float] = DCS_PROFILE_SOLAR_RANGE_DEG,
    bin_width_deg: float = DCS_PROFILE_BIN_WIDTH_DEG,
) -> None:
    """Plot DCS radiant drift and counts in fixed solar-longitude bins."""
    solar_min, solar_max = sorted(map(float, solar_range_deg))
    ra_deg, dec_deg = ecliptic_to_equatorial_deg(rows["ecliptic_lon"], rows["beta"])
    finite = np.isfinite(rows["solar_lon"]) & np.isfinite(ra_deg) & np.isfinite(dec_deg)
    solar_lon = np.asarray(rows["solar_lon"][finite], dtype=np.float64)
    ra_deg = ra_deg[finite]
    dec_deg = dec_deg[finite]

    edges = np.arange(solar_min, solar_max + bin_width_deg, bin_width_deg, dtype=np.float64)
    counts, _ = np.histogram(solar_lon, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.45), constrained_layout=True)
    scatter_style = {
        "s": 8,
        "facecolors": "none",
        "edgecolors": "0.30",
        "linewidths": 0.6,
        "alpha": 0.48,
    }
    axes[0].scatter(solar_lon, ra_deg, **scatter_style)
    axes[1].scatter(solar_lon, dec_deg, **scatter_style)
    axes[2].bar(centers, counts, width=bin_width_deg, color="0.25", edgecolor="0.25", linewidth=0.25)

    for ax in axes:
        ax.set_xlim(solar_min, solar_max)
        ax.set_xlabel(r"Solar longitude, $\lambda_\odot$ (deg)")
        ax.grid(alpha=0.20, linewidth=0.45)
    axes[0].set_ylabel(r"Geocentric right ascension, $\alpha_g$ (deg)")
    axes[1].set_ylabel(r"Geocentric declination, $\delta_g$ (deg)")
    axes[2].set_ylabel(rf"Counts per {bin_width_deg:g}$^\circ$")
    axes[0].set_ylim(295.0, 320.0)
    axes[1].set_ylim(-35.0, -20.0)
    axes[2].set_ylim(bottom=0.0)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=240)
    plt.close(fig)


def select_dcs_activity_pixel(rows: np.ndarray) -> np.ndarray:
    """Select the Radview radiant pixels used for the DCS activity profile."""
    pixel = ang2pix_ring(DCS_ACTIVITY_HEALPIX_NSIDE, rows["sun_centered_lon"], rows["beta"])
    return rows[np.isin(pixel, DCS_ACTIVITY_HEALPIX_PIXELS)]


def select_activity_rows(
    rows: np.ndarray,
    selection: ActivitySelection,
    *,
    require_radiant_pixels: bool = True,
) -> np.ndarray:
    """Apply the CSV-inferred seven-parameter window or the legacy selection."""
    if selection.parameter_window is not None:
        ra_deg, dec_deg = ecliptic_to_equatorial_deg(rows["ecliptic_lon"], rows["beta"])
        return rows[
            selection_mask(
                selection.parameter_window,
                vg_km_s=rows["vg"],
                kepler=rows["kepler"],
                ra_deg=ra_deg,
                dec_deg=dec_deg,
            )
        ]
    keep = np.ones(len(rows), dtype=bool)
    if selection.e_range is not None:
        eccentricity = rows["kepler"][:, 1]
        keep &= np.isfinite(eccentricity)
        keep &= (eccentricity >= selection.e_range[0]) & (eccentricity <= selection.e_range[1])
    if selection.vg_range is not None:
        keep &= np.isfinite(rows["vg"])
        keep &= (rows["vg"] >= selection.vg_range[0]) & (rows["vg"] <= selection.vg_range[1])
    if require_radiant_pixels:
        pixels = ang2pix_ring(selection.healpix_nside, rows["sun_centered_lon"], rows["beta"])
        keep &= np.isin(pixels, selection.healpix_pixels)
    return rows[keep]


def guarded_exposure_mask(
    exposure: np.ndarray,
    minimum_exposure_hours: float,
    minimum_contiguous_bins: int = 3,
) -> np.ndarray:
    """Return conservative exposure support for an activity curve."""
    sufficient = np.asarray(exposure, dtype=np.float64) >= minimum_exposure_hours
    adjacent_to_gap = np.zeros_like(sufficient)
    adjacent_to_gap[1:] |= ~sufficient[:-1]
    adjacent_to_gap[:-1] |= ~sufficient[1:]
    valid = sufficient & ~adjacent_to_gap
    padded = np.concatenate(([False], valid, [False]))
    starts = np.flatnonzero(~padded[:-1] & padded[1:])
    stops = np.flatnonzero(padded[:-1] & ~padded[1:])
    for start, stop in zip(starts, stops, strict=True):
        if stop - start < minimum_contiguous_bins:
            valid[start:stop] = False
    return valid


def zenith_correction_weights(
    rows: np.ndarray,
    *,
    min_cos_z: float,
    alpha: float,
) -> np.ndarray:
    """Return the fitted zenith-angle correction for every detected meteor."""
    ra_deg, dec_deg = ecliptic_to_equatorial_deg(rows["ecliptic_lon"], rows["beta"])
    altitude_deg = radiant_altitude_deg(ra_deg, dec_deg, rows["epoch"])
    cos_z = np.sin(np.deg2rad(altitude_deg))
    weights = np.zeros(len(rows), dtype=np.float64)
    visible = np.isfinite(cos_z) & (cos_z > 0.0)
    weights[visible] = np.maximum(cos_z[visible], float(min_cos_z)) ** (-float(alpha))
    return weights


def load_zenith_correction_parameters(path: Path) -> tuple[float, float]:
    """Read the correction fitted for the paper radiant-distribution analysis."""
    with h5py.File(path, "r") as h5:
        alpha = float(h5.attrs["fitted_zenith_exponent_alpha"])
        min_cos_z = float(h5.attrs["min_cos_z"])
    return alpha, min_cos_z


def activity_profile(
    rows: np.ndarray,
    exposure_path: Path,
    selection: ActivitySelection,
    coverage_rows: np.ndarray | None = None,
    *,
    zenith_alpha: float,
    min_cos_z: float,
) -> dict[str, np.ndarray]:
    """Return raw counts and exposure-plus-zenith-corrected sliding rates."""
    solar_min, solar_max = selection.solar_range_deg
    edges = np.arange(
        solar_min,
        solar_max + selection.bin_width_deg,
        selection.bin_width_deg,
        dtype=np.float64,
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    half_window = 0.5 * selection.profile_window_deg
    zenith_weights = zenith_correction_weights(
        rows,
        min_cos_z=min_cos_z,
        alpha=zenith_alpha,
    )
    window_masks = [
        np.abs(wrap180(rows["solar_lon"] - center)) <= half_window
        for center in centers
    ]
    counts = np.asarray(
        [np.sum(keep) for keep in window_masks],
        dtype=np.float64,
    )
    zenith_weighted_counts = np.asarray(
        [np.sum(zenith_weights[keep]) for keep in window_masks],
        dtype=np.float64,
    )
    zenith_weighted_variance = np.asarray(
        [np.sum(np.square(zenith_weights[keep])) for keep in window_masks],
        dtype=np.float64,
    )
    coverage = np.ones_like(counts, dtype=bool)
    if coverage_rows is not None:
        coverage_counts, _ = np.histogram(coverage_rows["solar_lon"], bins=edges)
        coverage = coverage_counts > 0
        counts[~coverage] = np.nan

    with h5py.File(exposure_path, "r") as h5:
        observation_epoch, observation_solar, observation_hours = interpolate_observation_solar_longitude(h5)
    exposure = np.zeros_like(centers)
    plot_lon = float(centered_plot_longitude_deg(selection.sun_centered_lon_deg))
    for i, center in enumerate(centers):
        keep = np.abs(wrap180(observation_solar - center)) <= half_window
        exposure[i] = float(
            radiant_exposure_hours_points(
                observation_epoch[keep],
                observation_solar[keep],
                observation_hours[keep],
                np.asarray([plot_lon]),
                np.asarray([selection.beta_deg]),
            )[0]
        )

    valid_exposure = guarded_exposure_mask(exposure, selection.minimum_exposure_hours)
    valid = valid_exposure & coverage
    counts[~valid] = np.nan
    zenith_weighted_counts[~valid] = np.nan
    zenith_weighted_variance[~valid] = np.nan
    rate = np.divide(
        zenith_weighted_counts,
        exposure,
        out=np.full_like(exposure, np.nan),
        where=valid,
    )
    uncertainty = np.divide(
        np.sqrt(zenith_weighted_variance),
        exposure,
        out=np.full_like(exposure, np.nan),
        where=valid,
    )
    return {
        "centers": centers,
        "counts": counts,
        "zenith_weighted_counts": zenith_weighted_counts,
        "exposure": exposure,
        "rate": rate,
        "uncertainty": uncertainty,
        "coverage": coverage,
        "valid_exposure": valid_exposure,
    }


def dcs_activity_profile(
    rows: np.ndarray,
    exposure_path: Path,
    solar_range_deg: tuple[float, float] = DCS_ACTIVITY_SOLAR_RANGE_DEG,
    bin_width_deg: float = DCS_ACTIVITY_BIN_WIDTH_DEG,
) -> dict[str, np.ndarray]:
    """Return exposure-corrected PANSY detections in the selected DCS radiant pixel."""
    solar_min, solar_max = sorted(map(float, solar_range_deg))
    edges = np.arange(solar_min, solar_max + bin_width_deg, bin_width_deg, dtype=np.float64)
    centers = 0.5 * (edges[:-1] + edges[1:])
    counts, _ = np.histogram(rows["solar_lon"], bins=edges)
    counts = counts.astype(np.float64)

    with h5py.File(exposure_path, "r") as h5:
        observation_epoch, observation_solar, observation_hours = interpolate_observation_solar_longitude(h5)
    exposure = np.zeros_like(centers)
    plot_lon = float(centered_plot_longitude_deg(DCS_ACTIVITY_MEAN_SC_LON_DEG))
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        keep = (observation_solar >= lo) & (observation_solar < hi)
        exposure[i] = float(
            radiant_exposure_hours_points(
                observation_epoch[keep],
                observation_solar[keep],
                observation_hours[keep],
                np.asarray([plot_lon]),
                np.asarray([DCS_ACTIVITY_MEAN_BETA_DEG]),
            )[0]
        )

    valid_exposure = guarded_exposure_mask(exposure, 3.0)
    counts[~valid_exposure] = np.nan
    rate = np.divide(
        counts,
        exposure,
        out=np.full_like(exposure, np.nan),
        where=valid_exposure,
    )
    uncertainty = np.divide(
        np.sqrt(counts),
        exposure,
        out=np.full_like(exposure, np.nan),
        where=valid_exposure,
    )
    return {
        "centers": centers,
        "counts": counts,
        "exposure": exposure,
        "rate": rate,
        "uncertainty": uncertainty,
    }


def plot_dcs_activity_panel(ax, profile: dict[str, np.ndarray]) -> None:
    """Plot the exposure-corrected DCS activity rate and detected counts."""
    x = profile["centers"]
    y = profile["rate"]
    uncertainty = profile["uncertainty"]
    counts = profile["counts"]
    ax.fill_between(x, y - uncertainty, y + uncertainty, color="C0", alpha=0.20, linewidth=0)
    ax.plot(x, y, color="C0", marker="o", ms=3.2, lw=1.4)
    half_bin = 0.5 * float(np.median(np.diff(x))) if len(x) > 1 else 0.5
    ax.set_xlim(float(x[0] - half_bin), float(x[-1] + half_bin))
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"Solar longitude, $\lambda_\odot$ (deg)")
    ax.set_ylabel(r"Detected rate (h$^{-1}$)")
    ax.grid(alpha=0.22, lw=0.45)
    count_ax = ax.twinx()
    count_ax.step(x, counts, where="mid", color="0.35", lw=1.0, alpha=0.75)
    count_ax.set_ylim(bottom=0.0)
    count_ax.set_ylabel("Raw count")


def plot_activity_inset(
    ax,
    rows: np.ndarray,
    selection: ActivitySelection,
    bounds: tuple[float, float, float, float],
    additional_markers: tuple[tuple[ActivitySelection, str], ...] = (),
) -> None:
    """Show all catalogue meteors within the activity interval and inset limits."""
    inset = ax.inset_axes(bounds)
    inset_lon = (
        selection.sun_centered_lon_deg
        if selection.inset_sun_centered_lon_deg is None
        else selection.inset_sun_centered_lon_deg
    )
    inset_beta = selection.beta_deg if selection.inset_beta_deg is None else selection.inset_beta_deg
    x = inset_lon + wrap180(
        rows["sun_centered_lon"] - inset_lon
    )
    y = rows["beta"]
    xlim = (
        selection.inset_xlim_deg
        if selection.inset_xlim_deg is not None
        else (selection.sun_centered_lon_deg + 2.3, selection.sun_centered_lon_deg - 2.3)
    )
    ylim = (
        selection.inset_ylim_deg
        if selection.inset_ylim_deg is not None
        else (selection.beta_deg - 2.3, selection.beta_deg + 2.3)
    )
    inset.set_xlim(xlim)
    inset.set_ylim(ylim)
    x_lo, x_hi = sorted(xlim)
    y_lo, y_hi = sorted(ylim)
    visible = (x >= x_lo) & (x <= x_hi) & (y >= y_lo) & (y <= y_hi)
    inset.scatter(
        x[visible],
        y[visible],
        s=7.5,
        color="black",
        alpha=0.48,
        linewidths=0,
    )
    inset.add_patch(
        plt.Circle(
            (inset_lon, inset_beta),
            1.0,
            fill=False,
            color="0.25",
            lw=0.8,
            ls="--",
        )
    )
    for marker_selection, marker_color in additional_markers:
        marker_x = inset_lon + float(
            wrap180(marker_selection.sun_centered_lon_deg - inset_lon)
        )
        inset.add_patch(
            plt.Circle(
                (marker_x, marker_selection.beta_deg),
                1.0,
                fill=False,
                color=marker_color,
                lw=0.9,
                ls="--",
            )
        )
    print(
        f"{selection.short_name} inset: {int(np.sum(visible))} contextual "
        f"meteors inside axes from {len(rows)} PANSY meteors in the solar range"
    )
    if selection.inset_xticks_deg is not None:
        inset.set_xticks(selection.inset_xticks_deg)
        inset.set_yticks(selection.inset_yticks_deg)
        inset.set_xticklabels(
            [f"{float(wrap360(value)):g}" for value in selection.inset_xticks_deg]
        )
        inset.tick_params(labelsize=8.0)
    else:
        inset.set_xticks(
            [
                selection.sun_centered_lon_deg + 1.0,
                selection.sun_centered_lon_deg,
                selection.sun_centered_lon_deg - 1.0,
            ]
        )
        inset.set_yticks(
            [
                selection.beta_deg - 1.0,
                selection.beta_deg,
                selection.beta_deg + 1.0,
            ]
        )
        inset.set_xticklabels(
            [
                f"{float(wrap360(selection.sun_centered_lon_deg + 1.0)):.1f}",
                f"{float(wrap360(selection.sun_centered_lon_deg)):.1f}",
                f"{float(wrap360(selection.sun_centered_lon_deg - 1.0)):.1f}",
            ],
            fontsize=8.0,
        )
        inset.set_yticklabels(
            [
                f"{selection.beta_deg - 1.0:.1f}",
                f"{selection.beta_deg:.1f}",
                f"{selection.beta_deg + 1.0:.1f}",
            ],
            fontsize=8.0,
        )
    inset.tick_params(length=2.0, pad=1.0)
    inset.set_xlabel(r"$\lambda_g-\lambda_\odot$ (deg)", fontsize=8.0, labelpad=1.0)
    inset.set_ylabel(r"$\beta_g$ (deg)", fontsize=8.0, labelpad=1.0)
    if selection.short_name == "CAP":
        inset.yaxis.tick_right()
        inset.yaxis.set_label_position("right")
    inset.grid(alpha=0.16, lw=0.35)


def plot_activity_curve(
    ax,
    profile: dict[str, np.ndarray],
    *,
    color: str,
    label: str | None = None,
) -> None:
    x = profile["centers"]
    y = profile["rate"]
    uncertainty = profile["uncertainty"]
    lower = np.maximum(0.0, y - uncertainty)
    upper = y + uncertainty
    ax.fill_between(x, lower, upper, color=color, alpha=0.20, linewidth=0)
    measured = np.isfinite(x) & np.isfinite(y)
    ax.plot(
        x[measured],
        y[measured],
        color=color,
        marker="o",
        ms=3.0,
        lw=1.35,
        label=label,
    )


def plot_activity_panel(
    ax,
    profile: dict[str, np.ndarray],
    selection: ActivitySelection,
    *,
    label_location: str,
    color: str = "C0",
    curve_label: str | None = None,
    panel_label: str | None = None,
) -> None:
    """Plot an exposure-plus-zenith-corrected activity curve."""
    plot_activity_curve(ax, profile, color=color, label=curve_label)
    ax.set_xlim(*selection.solar_range_deg)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"Solar longitude, $\lambda_\odot$ (deg)")
    ax.set_ylabel(r"Zenithal hourly rate (h$^{-1}$)")
    ax.grid(alpha=0.22, lw=0.45)

    label_x = 0.025 if label_location == "left" else 0.975
    ax.text(
        label_x,
        0.97,
        selection.short_name if panel_label is None else panel_label,
        transform=ax.transAxes,
        ha=label_location,
        va="top",
        fontsize=17,
        fontweight="bold",
        color="black",
        zorder=8,
    )
def orbit_xy(kepler: np.ndarray, samples: int = 361) -> tuple[np.ndarray, np.ndarray]:
    a, e, inc, raan, argp = kepler[:5]
    if not np.isfinite(a) or a <= 0.0 or not np.isfinite(e) or e < 0.0 or e >= 1.0:
        return np.asarray([]), np.asarray([])
    ecc_anom = np.linspace(0.0, 2.0 * np.pi, max(samples, int(360 + 720 * e)) + 1)
    x_pf = a * (np.cos(ecc_anom) - e)
    y_pf = a * np.sqrt(max(0.0, 1.0 - e * e)) * np.sin(ecc_anom)
    inc = np.deg2rad(inc)
    raan = np.deg2rad(raan)
    argp = np.deg2rad(argp)
    cos_o, sin_o = np.cos(raan), np.sin(raan)
    cos_w, sin_w = np.cos(argp), np.sin(argp)
    cos_i = np.cos(inc)
    x = (cos_o * cos_w - sin_o * sin_w * cos_i) * x_pf + (-cos_o * sin_w - sin_o * cos_w * cos_i) * y_pf
    y = (sin_o * cos_w + cos_o * sin_w * cos_i) * x_pf + (-sin_o * sin_w + cos_o * cos_w * cos_i) * y_pf
    return x, y


def orbit_xyz(kepler: np.ndarray, samples: int = 361) -> np.ndarray:
    a, e, inc, raan, argp = kepler[:5]
    if not np.isfinite(a) or a <= 0.0 or not np.isfinite(e) or e < 0.0 or e >= 1.0:
        return np.empty((3, 0))
    ecc_anom = np.linspace(0.0, 2.0 * np.pi, max(samples, int(360 + 720 * e)) + 1)
    x_pf = a * (np.cos(ecc_anom) - e)
    y_pf = a * np.sqrt(max(0.0, 1.0 - e * e)) * np.sin(ecc_anom)
    inc = np.deg2rad(inc)
    raan = np.deg2rad(raan)
    argp = np.deg2rad(argp)
    cos_o, sin_o = np.cos(raan), np.sin(raan)
    cos_w, sin_w = np.cos(argp), np.sin(argp)
    cos_i, sin_i = np.cos(inc), np.sin(inc)
    x = (cos_o * cos_w - sin_o * sin_w * cos_i) * x_pf + (-cos_o * sin_w - sin_o * cos_w * cos_i) * y_pf
    y = (sin_o * cos_w + cos_o * sin_w * cos_i) * x_pf + (-sin_o * sin_w + cos_o * cos_w * cos_i) * y_pf
    z = (sin_w * sin_i) * x_pf + (cos_w * sin_i) * y_pf
    return np.vstack((x, y, z))


def node_markers(kepler: np.ndarray) -> list[dict]:
    a, e, _inc, raan, argp = kepler[:5]
    out = []
    for label, nu_deg, node_lon in (
        ("asc.", -argp, raan),
        ("desc.", 180.0 - argp, raan + 180.0),
    ):
        nu = np.deg2rad(nu_deg)
        r = a * (1.0 - e * e) / (1.0 + e * np.cos(nu))
        lon = wrap360(node_lon).item()
        earth_solar = wrap360(lon - 180.0).item()
        out.append(
            {
                "label": label,
                "x": r * np.cos(np.deg2rad(lon)),
                "y": r * np.sin(np.deg2rad(lon)),
                "lon": lon,
                "solar": earth_solar,
            }
        )
    return out


def circular_mean(values: np.ndarray) -> float:
    rad = np.deg2rad(values[np.isfinite(values)])
    return wrap360(np.rad2deg(np.arctan2(np.mean(np.sin(rad)), np.mean(np.cos(rad))))).item() if len(rad) else np.nan


def mean_kepler(rows: np.ndarray) -> np.ndarray:
    kep = rows["kepler"]
    out = np.nanmean(kep, axis=0)
    out[3] = circular_mean(kep[:, 3])
    out[4] = circular_mean(kep[:, 4])
    out[5] = circular_mean(kep[:, 5])
    return out


def radiant_color_values(rows: np.ndarray, color_field: str) -> np.ndarray:
    """Return the requested radiant-panel color quantity."""
    if color_field == "vg":
        return np.asarray(rows["vg"], dtype=np.float64)
    kepler_column = {
        "a": 0,
        "eccentricity": 1,
        "inclination": 2,
        "omega": 4,
        "q": 6,
    }
    if color_field not in kepler_column:
        raise ValueError(f"unsupported radiant color field: {color_field}")
    return np.asarray(rows["kepler"][:, kepler_column[color_field]], dtype=np.float64)


def pansy_display_lon(lon_minus_sun_deg: np.ndarray) -> np.ndarray:
    return wrap180(-lon_minus_sun_deg - 90.0)


def lambert_project(display_lon_deg: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lon = np.deg2rad(wrap180(display_lon_deg))
    lat = np.deg2rad(np.clip(lat_deg, -89.999999, 89.999999))
    cos_lat = np.cos(lat)
    inner = 1.0 + cos_lat * np.cos(lon)
    k = np.sqrt(2.0 / np.maximum(inner, 1e-12))
    return 0.5 * k * cos_lat * np.sin(lon), 0.5 * k * np.sin(lat)


def draw_projection_grid(ax):
    t = np.linspace(0.0, 2.0 * np.pi, 361)
    ax.plot(np.cos(t), np.sin(t), color="0.72", lw=0.55, alpha=0.75)
    for lat in [-60, -30, 0, 30, 60]:
        lon = np.linspace(-180.0, 180.0, 721)
        x, y = lambert_project(lon, np.full_like(lon, lat))
        good = np.isfinite(x) & np.isfinite(y) & (x * x + y * y <= 1.0001)
        ax.plot(x[good], y[good], color="0.78", lw=0.45, alpha=0.55)
    for lon0 in np.arange(-150, 181, 30):
        lat = np.linspace(-89.9, 89.9, 361)
        x, y = lambert_project(np.full_like(lat, lon0), lat)
        good = np.isfinite(x) & np.isfinite(y) & (x * x + y * y <= 1.0001)
        ax.plot(x[good], y[good], color="0.78", lw=0.45, alpha=0.55)


def plot_radiant_panel(
    ax,
    rows: np.ndarray,
    passage: Passage,
    radius_deg: float,
    lon_zoom: float,
    lat_zoom: float,
    solar_half_width: float,
    color_field: str = "vg",
    color_norm: Normalize | None = None,
    color_cmap: str = "viridis",
):
    xoff = wrap180(rows["sun_centered_lon"] - passage.sun_centered_lon_deg)
    x = passage.sun_centered_lon_deg + xoff
    near = (np.abs(xoff) <= lon_zoom) & (np.abs(rows["beta"] - passage.beta_deg) <= lat_zoom)
    if color_field == "black":
        mesh = ax.scatter(
            x[near],
            rows["beta"][near],
            s=8,
            color="black",
            alpha=0.62,
            linewidths=0,
        )
    else:
        if color_norm is None:
            _label, vmin, vmax, _cmap = RADIANT_COLOR_SPECS[color_field]
            color_norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
        mesh = ax.scatter(
            x[near],
            rows["beta"][near],
            s=8,
            c=radiant_color_values(rows, color_field)[near],
            cmap=color_cmap,
            norm=color_norm,
            alpha=0.62,
            linewidths=0,
        )
    ax.add_patch(plt.Circle((passage.sun_centered_lon_deg, passage.beta_deg), radius_deg, fill=False, color="black", lw=0.9, ls="--"))
    ax.set_xlim(passage.sun_centered_lon_deg + lon_zoom, passage.sun_centered_lon_deg - lon_zoom)
    ax.set_ylim(passage.beta_deg - lat_zoom, passage.beta_deg + lat_zoom)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f"{float(wrap360(value)):g}"))
    ax.grid(alpha=0.22, lw=0.45)
    ax.set_title(rf"{passage.name}, $\lambda_\odot={passage.solar_lon_deg:.1f}\pm{solar_half_width:.1f}^\circ$")
    ax.set_xlabel(r"$\lambda'_g=\lambda_g-\lambda_\odot$ (deg)")
    return mesh


def plot_orbits(
    ax,
    selections: list[tuple[Passage, np.ndarray]],
    colors: list[str],
    show_legend: bool = True,
):
    theta = np.linspace(0.0, 2.0 * np.pi, 361)
    ax.plot(np.cos(theta), np.sin(theta), color="0.35", lw=0.9, label="Earth")
    ax.plot(5.204 * np.cos(theta), 5.204 * np.sin(theta), color="0.20", lw=1.0, label="Jupiter")
    ax.scatter([0], [0], marker="o", s=45, color="#f5b342", edgecolor="black", linewidth=0.4, zorder=5)
    cx, cy = orbit_xy(COMET_169P_NEAT)
    comet_good = np.isfinite(cx) & np.isfinite(cy) & (np.hypot(cx, cy) < 6.0)
    for (passage, rows), color in zip(selections, colors, strict=True):
        if len(rows) == 0:
            continue
        orbit_count = 0
        for row in rows:
            a, e = row["kepler"][:2]
            if not np.isfinite(a) or not np.isfinite(e) or a <= 0.0 or e < 0.0 or e >= 1.0:
                continue
            x, y = orbit_xy(row["kepler"])
            if len(x) == 0:
                continue
            ax.plot(x, y, color=color, alpha=0.12, lw=0.7)
            nodes = node_markers(row["kepler"])
            node = min(nodes, key=lambda item: abs(wrap180(item["solar"] - passage.solar_lon_deg)))
            ax.scatter([node["x"]], [node["y"]], marker=".", s=9, color=color, alpha=0.26, linewidths=0, zorder=4)
            orbit_count += 1
            if orbit_count >= 160:
                break
        if "Daytime" in passage.name:
            shower_label = "DCS"
        elif "alpha" in passage.name:
            shower_label = "CAP"
        else:
            shower_label = "OES"
        ax.plot([], [], color=color, alpha=0.75, lw=1.3, label=shower_label)
        earth_lon = np.deg2rad(passage.solar_lon_deg + 180.0)
        ex, ey = np.cos(earth_lon), np.sin(earth_lon)
        ax.scatter([ex], [ey], marker="o", s=34, color=color, edgecolor="black", linewidth=0.35, zorder=6)
    ax.plot(cx[comet_good], cy[comet_good], color=COMET_COLOR, lw=3.0, label="169P/NEAT")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-5.35, 5.35)
    ax.set_ylim(-5.35, 5.35)
    ax.set_xlabel("Ecliptic X (AU)")
    ax.set_ylabel("Ecliptic Y (AU)")
    ax.grid(alpha=0.22, lw=0.45)
    if show_legend:
        legend = ax.legend(
            loc="lower center",
            bbox_to_anchor=(0.5, 1.005),
            ncol=2,
            columnspacing=1.1,
            handletextpad=0.5,
            fontsize=10.5,
            frameon=True,
            framealpha=1.0,
        )
        legend.get_frame().set_facecolor("white")
        legend.get_frame().set_edgecolor("none")


def plot_orbits_side_view(
    ax,
    selections: list[tuple[Passage, np.ndarray]],
    colors: list[str],
    show_legend: bool = True,
):
    reference_rows = [row for _passage, rows in selections for row in rows]
    if reference_rows:
        reference_kepler = mean_kepler(np.asarray(reference_rows, dtype=selections[0][1].dtype))
        reference_xyz = orbit_xyz(reference_kepler, samples=1440)
        radius = np.linalg.norm(reference_xyz, axis=0)
        aphelion = reference_xyz[:, int(np.nanargmax(radius))]
        aphelion[2] = 0.0
        if np.linalg.norm(aphelion[:2]) > 0.0:
            x_unit = aphelion / np.linalg.norm(aphelion)
        else:
            x_unit = np.asarray([1.0, 0.0, 0.0])
    else:
        x_unit = np.asarray([1.0, 0.0, 0.0])
    z_unit = np.asarray([0.0, 0.0, 1.0])

    def project(xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x_unit @ xyz, z_unit @ xyz

    theta = np.linspace(0.0, 2.0 * np.pi, 361)
    for radius, label, color, lw in ((1.0, "Earth", "0.35", 0.9), (5.204, "Jupiter", "0.12", 1.4)):
        xyz = np.vstack((radius * np.cos(theta), radius * np.sin(theta), np.zeros_like(theta)))
        x, z = project(xyz)
        ax.plot(x, z, color=color, lw=lw, label=label)
    for (passage, rows), color in zip(selections, colors, strict=True):
        if len(rows) == 0:
            continue
        orbit_count = 0
        for row in rows:
            a, e = row["kepler"][:2]
            if not np.isfinite(a) or not np.isfinite(e) or a <= 0.0 or e < 0.0 or e >= 1.0:
                continue
            xyz = orbit_xyz(row["kepler"])
            if xyz.shape[1] == 0:
                continue
            x, z = project(xyz)
            ax.plot(x, z, color=color, alpha=0.12, lw=0.7)
            orbit_count += 1
            if orbit_count >= 160:
                break
        if "Daytime" in passage.name:
            shower_label = "DCS"
        elif "alpha" in passage.name:
            shower_label = "CAP"
        else:
            shower_label = "OES"
        ax.plot([], [], color=color, alpha=0.75, lw=1.3, label=shower_label)
        earth_lon = np.deg2rad(passage.solar_lon_deg + 180.0)
        earth_xyz = np.asarray([[np.cos(earth_lon)], [np.sin(earth_lon)], [0.0]])
        ex, ez = project(earth_xyz)
        ax.scatter(ex, ez, marker="o", s=34, color=color, edgecolor="black", linewidth=0.35, zorder=6)
    comet_xyz = orbit_xyz(COMET_169P_NEAT)
    cx, cz = project(comet_xyz)
    comet_good = np.isfinite(cx) & np.isfinite(cz) & (np.hypot(cx, cz) < 6.0)
    ax.plot(cx[comet_good], cz[comet_good], color=COMET_COLOR, lw=3.0, label="169P/NEAT")
    ax.axhline(0.0, color="0.75", lw=0.85)
    ax.scatter([0], [0], marker="o", s=45, color="#f5b342", edgecolor="black", linewidth=0.4, zorder=5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-5.35, 5.35)
    ax.set_ylim(-5.35, 5.35)
    ax.set_xlabel("Mean aphelion direction (AU)")
    ax.set_ylabel("Ecliptic Z (AU)")
    ax.grid(alpha=0.22, lw=0.45)
    if show_legend:
        legend = ax.legend(
            loc="upper right",
            ncol=2,
            columnspacing=0.9,
            handletextpad=0.45,
            fontsize=10.5,
            frameon=True,
            framealpha=1.0,
        )
        legend.get_frame().set_facecolor("white")
        legend.get_frame().set_edgecolor("none")


def plot_orbit_panel_figure(selections: list[tuple[Passage, np.ndarray]], out: Path):
    fig, ax = plt.subplots(figsize=(5.4, 5.4), constrained_layout=True)
    plot_orbits(ax, selections, colors=list(ORBIT_COLORS))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=240)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radview-data", type=Path, default=RADVIEW_DATA)
    parser.add_argument(
        "--catalogue",
        type=Path,
        default=DEFAULT_CATALOGUE,
        help="Complete combined PANSY/MAARSY HDF5 catalogue",
    )
    parser.add_argument("--output", type=Path, default=Path.home() / "src" / "pansy_paper" / "paper_capricornid_conjugate_stream.png")
    parser.add_argument("--orbit-output", type=Path)
    parser.add_argument(
        "--profile-output",
        type=Path,
        default=Path.home() / "src" / "pansy_paper" / "paper_dcs_solar_longitude_profile.png",
    )
    parser.add_argument("--profile-solar-min-deg", type=float, default=DCS_PROFILE_SOLAR_RANGE_DEG[0])
    parser.add_argument("--profile-solar-max-deg", type=float, default=DCS_PROFILE_SOLAR_RANGE_DEG[1])
    parser.add_argument("--profile-bin-width-deg", type=float, default=DCS_PROFILE_BIN_WIDTH_DEG)
    parser.add_argument("--activity-exposure", type=Path, default=DEFAULT_EXPOSURE)
    parser.add_argument(
        "--zenith-correction-sidecar",
        type=Path,
        default=DEFAULT_ZENITH_SIDECAR,
        help="Paper radiant sidecar supplying the fitted zenith exponent and cosine floor.",
    )
    parser.add_argument(
        "--zenith-alpha",
        type=float,
        help="Override the fitted zenith-angle exponent from --zenith-correction-sidecar.",
    )
    parser.add_argument(
        "--min-cos-z",
        type=float,
        help="Override the zenith-correction cosine floor from --zenith-correction-sidecar.",
    )
    parser.add_argument("--activity-solar-min-deg", type=float, default=DCS_ACTIVITY_SOLAR_RANGE_DEG[0])
    parser.add_argument("--activity-solar-max-deg", type=float, default=DCS_ACTIVITY_SOLAR_RANGE_DEG[1])
    parser.add_argument("--activity-bin-width-deg", type=float, default=DCS_ACTIVITY_BIN_WIDTH_DEG)
    parser.add_argument("--solar-half-width-deg", type=float, default=CLUSTER_SOLAR_WINDOW_DEG / 2.0)
    parser.add_argument("--radiant-radius-deg", type=float, default=5.0)
    parser.add_argument("--velocity-half-width-kms", type=float, default=10.0)
    parser.add_argument(
        "--radiant-color",
        choices=("black", *RADIANT_COLOR_SPECS),
        default="black",
        help="Quantity used to color the two radiant panels",
    )
    parser.add_argument("--color-min", type=float)
    parser.add_argument("--color-max", type=float)
    parser.add_argument("--radiant-cmap")
    parser.add_argument("--cluster-filter", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--vg-min", type=float, default=CLUSTER_VG_RANGE[0])
    parser.add_argument("--vg-max", type=float, default=CLUSTER_VG_RANGE[1])
    parser.add_argument("--e-min", type=float, default=CLUSTER_E_RANGE[0])
    parser.add_argument("--e-max", type=float, default=CLUSTER_E_RANGE[1])
    args = parser.parse_args()
    fitted_zenith_alpha, fitted_min_cos_z = load_zenith_correction_parameters(
        args.zenith_correction_sidecar
    )
    zenith_alpha = (
        fitted_zenith_alpha if args.zenith_alpha is None else float(args.zenith_alpha)
    )
    min_cos_z = fitted_min_cos_z if args.min_cos_z is None else float(args.min_cos_z)

    rows_by_passage = [
        (
            p,
            load_passage_rows(
                args.radview_data,
                p,
                args.solar_half_width_deg,
                catalogue_path=args.catalogue,
                source="PANSY",
            ),
        )
        for p in PASSAGES
    ]
    for passage, rows in rows_by_passage:
        validate_passage_sources(rows, passage, required_sources=("PANSY",))
    vg_range = (min(args.vg_min, args.vg_max), max(args.vg_min, args.vg_max))
    e_range = (min(args.e_min, args.e_max), max(args.e_min, args.e_max))
    if args.cluster_filter:
        plot_rows_by_passage = [(p, cluster_filter(rows, vg_range, e_range)) for p, rows in rows_by_passage]
    else:
        plot_rows_by_passage = rows_by_passage
    selections = [(p, select_associated(rows, p, args.radiant_radius_deg, args.velocity_half_width_kms)) for p, rows in plot_rows_by_passage]

    activity_selections = [
        replace(
            ACTIVITY_SELECTIONS[0],
            solar_range_deg=(args.activity_solar_min_deg, args.activity_solar_max_deg),
            bin_width_deg=args.activity_bin_width_deg,
        ),
        ACTIVITY_SELECTIONS[1],
        ACTIVITY_SELECTIONS[2],
    ]
    activity_products = []
    full_catalogue_rows = (
        load_catalogue_rows(args.catalogue, source="PANSY")
        if args.catalogue.exists()
        else None
    )
    activity_row_cache: dict[tuple[float, float], np.ndarray] = {}
    for activity_selection in activity_selections:
        solar_min, solar_max = activity_selection.solar_range_deg
        activity_center = 0.5 * (solar_min + solar_max)
        activity_half_width = 0.5 * (solar_max - solar_min)
        if full_catalogue_rows is not None:
            raw_rows = full_catalogue_rows
        else:
            cache_key = (solar_min, solar_max)
            raw_rows = activity_row_cache.get(cache_key)
            if raw_rows is None:
                raw_rows = load_chunk_rows(
                    args.radview_data,
                    "pansy",
                    activity_center,
                    activity_half_width,
                )
                activity_row_cache[cache_key] = raw_rows
        selected_rows = select_activity_rows(raw_rows, activity_selection)
        profile = activity_profile(
            selected_rows,
            args.activity_exposure,
            activity_selection,
            coverage_rows=raw_rows,
            zenith_alpha=zenith_alpha,
            min_cos_z=min_cos_z,
        )
        activity_products.append((activity_selection, selected_rows, raw_rows, profile))

    plt.rcParams.update(
        {
            "font.size": 13,
            "axes.labelsize": 15,
            "axes.titlesize": 15,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 11,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 9.2), constrained_layout=True)
    dcs_selection, _dcs_rows, _dcs_inset_rows, dcs_profile = activity_products[0]
    oes_selection, _oes_rows, _oes_inset_rows, oes_profile = activity_products[1]
    cap_selection, _cap_rows, _cap_inset_rows, cap_profile = activity_products[2]
    plot_activity_panel(
        axes[0, 0],
        dcs_profile,
        dcs_selection,
        label_location="left",
        curve_label="DCS",
        panel_label="OES / DCS",
    )
    plot_activity_curve(axes[0, 0], oes_profile, color="C2", label="OES")
    combined_upper = max(
        np.nanmax(dcs_profile["rate"] + dcs_profile["uncertainty"]),
        np.nanmax(oes_profile["rate"] + oes_profile["uncertainty"]),
    )
    axes[0, 0].set_ylim(0.0, 1.05 * combined_upper)
    axes[0, 0].legend(loc="lower left", frameon=False, fontsize=11)
    plot_activity_panel(
        axes[0, 1],
        cap_profile,
        cap_selection,
        label_location="right",
        color="C1",
        curve_label="CAP",
    )
    plot_orbits(axes[1, 0], selections, colors=list(ORBIT_COLORS), show_legend=False)
    plot_orbits_side_view(axes[1, 1], selections, colors=list(ORBIT_COLORS), show_legend=True)
    axes[1, 0].set_title("North ecliptic view")
    axes[1, 1].set_title("Ecliptic side view")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=240)
    plt.close(fig)
    if args.orbit_output is not None:
        plot_orbit_panel_figure(selections, args.orbit_output)
    if args.profile_output is not None:
        profile_solar_range = (args.profile_solar_min_deg, args.profile_solar_max_deg)
        profile_center = 0.5 * sum(profile_solar_range)
        profile_half_width = 0.5 * abs(profile_solar_range[1] - profile_solar_range[0])
        if args.catalogue.exists():
            profile_rows = load_catalogue_rows(
                args.catalogue,
                profile_center,
                profile_half_width,
                source="PANSY",
            )
        else:
            profile_rows = load_chunk_rows(
                args.radview_data,
                "pansy",
                profile_center,
                profile_half_width,
            )
        if args.cluster_filter:
            profile_rows = cluster_filter(profile_rows, vg_range, e_range)
        profile_rows = select_associated(
            profile_rows,
            PASSAGES[0],
            args.radiant_radius_deg,
            args.velocity_half_width_kms,
        )
        plot_dcs_solar_longitude_profile(
            profile_rows,
            args.profile_output,
            solar_range_deg=profile_solar_range,
            bin_width_deg=args.profile_bin_width_deg,
        )
        print(f"DCS profile: selected {len(profile_rows)}")
        print(args.profile_output)
    for activity_selection, selected_rows, _inset_rows, profile in activity_products:
        if activity_selection.parameter_window is not None:
            print(format_window(activity_selection.parameter_window))
        finite_rate = np.isfinite(profile["rate"])
        peak_window_count = int(
            np.sum(
                np.abs(
                    wrap180(selected_rows["solar_lon"] - activity_selection.peak_solar_lon_deg)
                )
                <= 0.5 * activity_selection.peak_window_deg
            )
        )
        if np.any(finite_rate):
            peak_index = int(np.nanargmax(profile["rate"]))
            print(
                f"{activity_selection.short_name} activity: selected {len(selected_rows)}; "
                f"{peak_window_count} in peak +/- "
                f"{0.5 * activity_selection.peak_window_deg:.1f} deg; "
                f"peak lambda_sun={profile['centers'][peak_index]:.2f} deg, "
                f"rate={profile['rate'][peak_index]:.3f} h^-1, "
                f"count={profile['counts'][peak_index]}, "
                f"exposure={profile['exposure'][peak_index]:.2f} h"
            )
    print(
        f"activity zenith correction: alpha={zenith_alpha:.6g}, "
        f"min_cos_z={min_cos_z:.6g}, source={args.zenith_correction_sidecar}"
    )
    for passage, rows in selections:
        print(f"{passage.name}: selected {len(rows)}")
        if len(rows):
            print(f"  mean solar {circular_mean(rows['solar_lon']):.2f} mean vg {np.nanmean(rows['vg']):.2f} mean beta {np.nanmean(rows['beta']):.2f}")
    if args.cluster_filter:
        print(f"cluster filter vg={vg_range[0]:.2f}-{vg_range[1]:.2f} km/s e={e_range[0]:.3f}-{e_range[1]:.3f}")
    for passage, rows in rows_by_passage:
        datasets, counts = np.unique(rows["dataset"], return_counts=True)
        detail = ", ".join(f"{name}={count}" for name, count in zip(datasets, counts, strict=True))
        print(f"{passage.name}: radiant rows {len(rows)} ({detail})")
    if args.cluster_filter:
        for passage, rows in plot_rows_by_passage:
            datasets, counts = np.unique(rows["dataset"], return_counts=True) if len(rows) else ([], [])
            detail = ", ".join(f"{name}={count}" for name, count in zip(datasets, counts, strict=True))
            print(f"{passage.name}: filtered radiant rows {len(rows)} ({detail})")
    print(args.output)


if __name__ == "__main__":
    main()
