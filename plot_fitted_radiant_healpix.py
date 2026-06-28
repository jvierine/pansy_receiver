#!/usr/bin/env python3
"""HEALPix histogram of DASST-corrected PANSY meteor radiants."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath
from matplotlib.transforms import Affine2D
import numpy as np

from radiant_visibility import composite_radiant_visibility_points, radiant_visibility_boundary, visibility_samples_from_radiants


PLOT_CENTER_LONGITUDE_DEG = -90.0


def wrap180(deg):
    return (np.asarray(deg, dtype=np.float64) + 180.0) % 360.0 - 180.0


def read_radiants(path: Path):
    with h5py.File(path, "r") as h:
        arr = h["radiants"][()]
    lon = np.asarray(arr["lambda_minus_sun_signed_deg"], dtype=np.float64)
    lat = np.asarray(arr["radiant_beta_ecliptic_deg"], dtype=np.float64)
    speed = np.asarray(arr["speed_km_s"], dtype=np.float64)
    good = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(speed)
    good &= (lat >= -90.0) & (lat <= 90.0)
    return lon[good], lat[good], speed[good]


def read_radiant_rows(path: Path):
    with h5py.File(path, "r") as h:
        return h["radiants"][()]


def healpix_histogram(lon_deg, lat_deg, speed_km_s, nside: int, min_count_for_mean_speed: int):
    if nside <= 0 or nside & (nside - 1):
        raise ValueError("nside must be a positive power of two")
    npix = hp.nside2npix(nside)
    pix = hp.ang2pix(nside, lon_deg, lat_deg, lonlat=True)
    count = np.bincount(pix, minlength=npix).astype(np.float64)
    speed_sum = np.bincount(pix, weights=speed_km_s, minlength=npix).astype(np.float64)
    mean_speed = np.full(npix, np.nan, dtype=np.float64)
    good = count >= min_count_for_mean_speed
    mean_speed[good] = speed_sum[good] / count[good]
    return count, mean_speed


def write_h5(path: Path, input_h5: Path, nside: int, count, mean_speed, n_radiants: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h:
        h.attrs["script"] = Path(__file__).name
        h.attrs["input_h5"] = str(input_h5)
        h.attrs["coordinate_frame"] = "sun-centered ecliptic"
        h.attrs["longitude_dataset"] = "lambda_minus_sun_signed_deg"
        h.attrs["latitude_dataset"] = "radiant_beta_ecliptic_deg"
        h.attrs["longitude_convention"] = "lambda_radiant - lambda_sun in [-180, 180) deg"
        h.attrs["plot_center_longitude_deg"] = PLOT_CENTER_LONGITUDE_DEG
        h.attrs["nside"] = int(nside)
        h.attrs["n_radiants"] = int(n_radiants)
        h.create_dataset("count", data=count)
        h.create_dataset("mean_speed_km_s", data=mean_speed)


def circular_mean_deg(values_deg):
    values = np.asarray(values_deg, dtype=np.float64)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan
    rad = np.deg2rad(values)
    return float(np.rad2deg(np.arctan2(np.mean(np.sin(rad)), np.mean(np.cos(rad)))) % 360.0)


def mean_solar_longitude(path: Path):
    with h5py.File(path, "r") as h:
        arr = h["radiants"][()]
    if "sun_lambda_ecliptic_deg" not in arr.dtype.names:
        return np.nan
    return circular_mean_deg(arr["sun_lambda_ecliptic_deg"])


def add_visibility_boundary_healpix(solar_longitude_deg: float, color: str = "white"):
    if not np.isfinite(float(solar_longitude_deg)):
        return
    plot_lon_deg, beta_deg = radiant_visibility_boundary(solar_longitude_deg)
    lambda_minus_sun = wrap180(PLOT_CENTER_LONGITUDE_DEG - plot_lon_deg)
    hp.projplot(
        lambda_minus_sun,
        beta_deg,
        lonlat=True,
        color=color,
        linewidth=1.35,
        linestyle="--",
        alpha=0.92,
    )


def add_composite_visibility_boundary_healpix(rows, color: str = "white", nside: int = 128):
    sample_epoch, sample_sun = visibility_samples_from_radiants(rows)
    if len(sample_epoch) == 0:
        return False
    lon = np.linspace(-180.0, 180.0, 361)
    lat = np.linspace(-90.0, 90.0, 181)
    lon_grid, lat_grid = np.meshgrid(lon, lat)
    possible = composite_radiant_visibility_points(lon_grid, lat_grid, sample_epoch, sample_sun)
    if not np.any(possible) or np.all(possible):
        return False
    fig_tmp, ax_tmp = plt.subplots()
    contour = ax_tmp.contour(lon, lat, possible.astype(float), levels=[0.5])
    segments = contour.allsegs[0] if contour.allsegs else []
    plt.close(fig_tmp)
    for seg in segments:
        if len(seg) < 2:
            continue
        hp.projplot(seg[:, 0], seg[:, 1], lonlat=True, color=color, linewidth=1.35, linestyle="--", alpha=0.92)
    return True


def mask_healpix_outside_oval(fig):
    """Keep HEALPix image pixels out of the rectangular area around the Mollweide oval."""
    for ax in fig.axes:
        ax.set_facecolor("white")
        if not ax.images:
            continue
        rect = MplPath(
            [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)],
            [MplPath.MOVETO, MplPath.LINETO, MplPath.LINETO, MplPath.LINETO, MplPath.CLOSEPOLY],
        )
        theta = np.linspace(2.0 * np.pi, 0.0, 241)
        ellipse_vertices = np.column_stack((0.5 + 0.5 * np.cos(theta), 0.5 + 0.5 * np.sin(theta)))
        ellipse_codes = np.full(len(ellipse_vertices), MplPath.LINETO, dtype=np.uint8)
        ellipse_codes[0] = MplPath.MOVETO
        ellipse_codes[-1] = MplPath.CLOSEPOLY
        ellipse = MplPath(ellipse_vertices, ellipse_codes)
        outside = MplPath.make_compound_path(rect, ellipse)
        ax.add_patch(PathPatch(outside, transform=ax.transAxes, facecolor="white", edgecolor="none", zorder=20))


def add_healpix_coordinate_labels():
    longitude_labels = [(-90.0, 0.0, "0°"), (0.0, 0.0, "270°"), (90.0, 0.0, "180°")]
    for lon, lat, label in longitude_labels:
        hp.projtext(
            lon,
            lat,
            label,
            lonlat=True,
            color="white",
            fontsize=9,
            ha="center",
            va="center",
            zorder=30,
        )
    for lat in (-30.0, 30.0):
        hp.projtext(
            -170.0,
            lat,
            f"{int(lat)}°",
            lonlat=True,
            color="white",
            fontsize=9,
            ha="left",
            va="center",
            zorder=30,
        )


def plot_healpix(count, output_png: Path, n_radiants: int, solar_longitude_deg: float = np.nan, rows=None):
    output_png.parent.mkdir(parents=True, exist_ok=True)
    plot_count = np.asarray(count, dtype=np.float64).copy()
    plot_count[~np.isfinite(plot_count) | (plot_count <= 0.0)] = 1.0
    cmap = plt.get_cmap("plasma").copy()
    cmap.set_bad(cmap(0.0))
    cmap.set_under(cmap(0.0))
    hp.mollview(
        plot_count,
        fig=1,
        rot=(PLOT_CENTER_LONGITUDE_DEG, 0.0, 0.0),
        flip="astro",
        cmap=cmap,
        badcolor=cmap(0.0),
        bgcolor="white",
        min=1.0,
        norm="log",
        title="",
        unit="Count / pixel",
        cbar=True,
        notext=True,
    )
    hp.graticule(dpar=15, dmer=30, color="0.55", alpha=0.35)
    if rows is None or not add_composite_visibility_boundary_healpix(rows, color="white"):
        add_visibility_boundary_healpix(solar_longitude_deg, color="white")
    add_healpix_coordinate_labels()
    fig = plt.gcf()
    fig.set_size_inches(9.4, 5.6)
    fig.patch.set_facecolor("white")
    mask_healpix_outside_oval(fig)
    fig.text(
        0.5,
        0.035,
        r"Sun-centered ecliptic longitude, $\lambda-\lambda_\odot$ (apex centered at $-90^\circ$)",
        ha="center",
        va="center",
        fontsize=12,
    )
    fig.text(0.035, 0.52, r"Ecliptic latitude, $\beta$", ha="center", va="center", rotation="vertical", fontsize=12)
    fig.text(0.86, 0.92, f"N={n_radiants}", ha="right", va="center", fontsize=10, color="0.2")
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_healpix_with_showers(
    count,
    output_png: Path,
    n_radiants: int,
    showers,
    solar_longitude_deg: float = np.nan,
    rows=None,
):
    plot_count = np.asarray(count, dtype=np.float64).copy()
    plot_count[~np.isfinite(plot_count) | (plot_count <= 0.0)] = 1.0
    cmap = plt.get_cmap("plasma").copy()
    cmap.set_bad(cmap(0.0))
    cmap.set_under(cmap(0.0))
    hp.mollview(
        plot_count,
        fig=1,
        rot=(PLOT_CENTER_LONGITUDE_DEG, 0.0, 0.0),
        flip="astro",
        cmap=cmap,
        badcolor=cmap(0.0),
        bgcolor="white",
        min=1.0,
        norm="log",
        title="",
        unit="Count / pixel",
        cbar=True,
        notext=True,
    )
    hp.graticule(dpar=15, dmer=30, color="0.55", alpha=0.35)
    if rows is None or not add_composite_visibility_boundary_healpix(rows, color="white"):
        add_visibility_boundary_healpix(solar_longitude_deg, color="white")
    add_healpix_coordinate_labels()
    if showers:
        lon = np.asarray([s.radiant_solar_ecliptic_lon_deg for s in showers], dtype=np.float64)
        lon = (lon + 180.0) % 360.0 - 180.0
        lat = np.asarray([s.radiant_ecliptic_lat_deg for s in showers], dtype=np.float64)
        hp.projscatter(
            lon,
            lat,
            lonlat=True,
            marker="o",
            s=150,
            facecolors="none",
            edgecolors="white",
            linewidth=1.5,
            alpha=0.3,
        )
        for shower, lo, la in zip(showers[:18], lon[:18], lat[:18], strict=False):
            hp.projtext(lo, la, shower.code or shower.name[:5], lonlat=True, fontsize=7.5, color="white")
    fig = plt.gcf()
    fig.set_size_inches(9.4, 5.6)
    fig.patch.set_facecolor("white")
    mask_healpix_outside_oval(fig)
    fig.text(
        0.5,
        0.035,
        r"Sun-centered ecliptic longitude, $\lambda-\lambda_\odot$ (apex centered at $-90^\circ$)",
        ha="center",
        va="center",
        fontsize=12,
    )
    fig.text(0.035, 0.52, r"Ecliptic latitude, $\beta$", ha="center", va="center", rotation="vertical", fontsize=12)
    fig.text(0.86, 0.92, f"N={n_radiants}", ha="right", va="center", fontsize=10, color="0.2")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Make a HEALPix histogram of DASST-corrected PANSY radiants.")
    parser.add_argument("--input-h5", type=Path, default=Path("test_plots/fitted_radiant_distribution.h5"))
    parser.add_argument("--output-png", type=Path, default=Path("test_plots/fitted_radiant_distribution_healpix.png"))
    parser.add_argument("--output-h5", type=Path, default=Path("test_plots/fitted_radiant_distribution_healpix.h5"))
    parser.add_argument("--nside", type=int, default=64)
    parser.add_argument("--min-count-for-mean-speed", type=int, default=3)
    parser.add_argument("--show-shower-overlays", action="store_true")
    parser.add_argument("--shower-catalog", type=Path, default=None)
    parser.add_argument("--shower-solar-longitude", type=float, default=None)
    parser.add_argument("--shower-date", default=None)
    parser.add_argument("--shower-peak-tolerance-deg", type=float, default=5.0)
    args = parser.parse_args()

    rows = read_radiant_rows(args.input_h5)
    lon, lat, speed = read_radiants(args.input_h5)
    count, mean_speed = healpix_histogram(lon, lat, speed, args.nside, args.min_count_for_mean_speed)
    write_h5(args.output_h5, args.input_h5, args.nside, count, mean_speed, len(lon))
    if args.show_shower_overlays:
        from shower_radiant_overlay import active_showers
        import h5py

        showers, _query = active_showers(
            rows,
            shower_catalog=args.shower_catalog,
            shower_solar_longitude=args.shower_solar_longitude,
            shower_date=args.shower_date,
            peak_tolerance_deg=args.shower_peak_tolerance_deg,
        )
        plot_healpix_with_showers(count, args.output_png, len(lon), showers, _query, rows=rows)
    else:
        plot_healpix(count, args.output_png, len(lon), mean_solar_longitude(args.input_h5), rows=rows)
    print(f"healpix_radiants {len(lon)} nside {args.nside} occupied_pixels {int(np.sum(count > 0.0))}")
    print(args.output_png)
    print(args.output_h5)


if __name__ == "__main__":
    main()
