#!/usr/bin/env python3
"""Plot the Capricornid stream candidate and its opposite-node passage."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import FuncFormatter


RADVIEW_ROOT = Path.home() / "src" / "radview"
LIVE_RADVIEW_DATA = Path("/tmp/radiantviz_live_data")
RADVIEW_DATA = LIVE_RADVIEW_DATA if LIVE_RADVIEW_DATA.exists() else RADVIEW_ROOT / "public" / "data"
CLUSTER_SOLAR_WINDOW_DEG = 14.0
CLUSTER_VG_RANGE = (14.26, 30.24)
CLUSTER_E_RANGE = (0.664, 0.832)
ORBIT_COLORS = ("#0072B2", "#D55E00")
COMET_COLOR = "black"


@dataclass(frozen=True)
class Passage:
    name: str
    solar_lon_deg: float
    sun_centered_lon_deg: float
    beta_deg: float
    vg_km_s: float
    n: int
    mean_kepler: tuple[float, float, float, float, float, float, float] | None = None


PASSAGES = (
    Passage(
        name=r"Daytime $\chi$ Capricornids (DCS)",
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
)

COMET_169P_NEAT = np.asarray([2.604, 0.76796, 11.285, 176.04, 218.13, np.nan, 0.604], dtype=np.float64)
DCS_PROFILE_SOLAR_RANGE_DEG = (300.0, 330.0)
DCS_PROFILE_BIN_WIDTH_DEG = 0.25


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


def load_passage_rows(data_dir: Path, passage: Passage, solar_half_width: float) -> np.ndarray:
    parts = [
        load_chunk_rows(data_dir, "pansy", passage.solar_lon_deg, solar_half_width),
        load_chunk_rows(data_dir, "maarsy", passage.solar_lon_deg, solar_half_width),
    ]
    return np.concatenate([p for p in parts if len(p)]) if any(len(p) for p in parts) else parts[0]


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


def plot_radiant_panel(ax, rows: np.ndarray, passage: Passage, radius_deg: float, lon_zoom: float, lat_zoom: float, solar_half_width: float):
    xoff = wrap180(rows["sun_centered_lon"] - passage.sun_centered_lon_deg)
    x = passage.sun_centered_lon_deg + xoff
    near = (np.abs(xoff) <= lon_zoom) & (np.abs(rows["beta"] - passage.beta_deg) <= lat_zoom)
    norm = Normalize(vmin=10.0, vmax=75.0)
    mesh = ax.scatter(
        x[near],
        rows["beta"][near],
        s=8,
        c=rows["vg"][near],
        cmap="viridis",
        norm=norm,
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


def plot_orbits(ax, selections: list[tuple[Passage, np.ndarray]], colors: list[str]):
    theta = np.linspace(0.0, 2.0 * np.pi, 361)
    ax.plot(np.cos(theta), np.sin(theta), color="0.35", lw=0.9, label="Earth orbit")
    ax.plot(5.204 * np.cos(theta), 5.204 * np.sin(theta), color="0.20", lw=1.0, label="Jupiter orbit")
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
            if a * (1.0 + e) > 5.35:
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
        ax.plot([], [], color=color, alpha=0.75, lw=1.3, label=passage.name)
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
    legend = ax.legend(loc="upper right", fontsize=8, frameon=True, framealpha=1.0)
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
    parser.add_argument("--solar-half-width-deg", type=float, default=CLUSTER_SOLAR_WINDOW_DEG / 2.0)
    parser.add_argument("--radiant-radius-deg", type=float, default=5.0)
    parser.add_argument("--velocity-half-width-kms", type=float, default=10.0)
    parser.add_argument("--cluster-filter", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--vg-min", type=float, default=CLUSTER_VG_RANGE[0])
    parser.add_argument("--vg-max", type=float, default=CLUSTER_VG_RANGE[1])
    parser.add_argument("--e-min", type=float, default=CLUSTER_E_RANGE[0])
    parser.add_argument("--e-max", type=float, default=CLUSTER_E_RANGE[1])
    args = parser.parse_args()

    rows_by_passage = [(p, load_passage_rows(args.radview_data, p, args.solar_half_width_deg)) for p in PASSAGES]
    if args.cluster_filter:
        vg_range = (min(args.vg_min, args.vg_max), max(args.vg_min, args.vg_max))
        e_range = (min(args.e_min, args.e_max), max(args.e_min, args.e_max))
        plot_rows_by_passage = [(p, cluster_filter(rows, vg_range, e_range)) for p, rows in rows_by_passage]
    else:
        plot_rows_by_passage = rows_by_passage
    selections = [(p, select_associated(rows, p, args.radiant_radius_deg, args.velocity_half_width_kms)) for p, rows in plot_rows_by_passage]

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4), constrained_layout=True)
    sc = None
    for ax, (passage, rows), (_p, selected) in zip(axes[:2], plot_rows_by_passage, selections, strict=True):
        sc = plot_radiant_panel(ax, rows, passage, args.radiant_radius_deg, lon_zoom=32.0, lat_zoom=24.0, solar_half_width=args.solar_half_width_deg)
    axes[0].set_ylabel(r"Ecliptic latitude, $\beta$ (deg)")
    if sc is not None:
        cb = fig.colorbar(sc, ax=axes[:2], orientation="horizontal", pad=0.12, fraction=0.06)
        cb.set_label(r"Geocentric speed, $v_g$ (km s$^{-1}$)")
    plot_orbits(axes[2], selections, colors=list(ORBIT_COLORS))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=240)
    plt.close(fig)
    if args.orbit_output is not None:
        plot_orbit_panel_figure(selections, args.orbit_output)
    if args.profile_output is not None:
        profile_solar_range = (args.profile_solar_min_deg, args.profile_solar_max_deg)
        profile_center = 0.5 * sum(profile_solar_range)
        profile_half_width = 0.5 * abs(profile_solar_range[1] - profile_solar_range[0])
        profile_rows = load_chunk_rows(args.radview_data, "pansy", profile_center, profile_half_width)
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
