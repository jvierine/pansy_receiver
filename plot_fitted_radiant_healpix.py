#!/usr/bin/env python3
"""HEALPix histogram of DASST-corrected PANSY meteor radiants."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


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


def plot_healpix(count, output_png: Path, n_radiants: int):
    output_png.parent.mkdir(parents=True, exist_ok=True)
    masked = hp.ma(count)
    masked.mask = count <= 0.0
    cmap = plt.get_cmap("plasma").copy()
    cmap.set_bad("white")
    cmap.set_under("white")
    hp.mollview(
        masked,
        fig=1,
        rot=(PLOT_CENTER_LONGITUDE_DEG, 0.0, 0.0),
        flip="astro",
        cmap=cmap,
        min=1.0,
        title="",
        unit="Count / pixel",
        cbar=True,
    )
    hp.graticule(dpar=15, dmer=30, color="0.55", alpha=0.35)
    fig = plt.gcf()
    fig.set_size_inches(9.4, 5.6)
    fig.patch.set_facecolor("white")
    for ax in fig.axes:
        ax.set_facecolor("white")
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


def main():
    parser = argparse.ArgumentParser(description="Make a HEALPix histogram of DASST-corrected PANSY radiants.")
    parser.add_argument("--input-h5", type=Path, default=Path("test_plots/fitted_radiant_distribution.h5"))
    parser.add_argument("--output-png", type=Path, default=Path("test_plots/fitted_radiant_distribution_healpix.png"))
    parser.add_argument("--output-h5", type=Path, default=Path("test_plots/fitted_radiant_distribution_healpix.h5"))
    parser.add_argument("--nside", type=int, default=32)
    parser.add_argument("--min-count-for-mean-speed", type=int, default=3)
    args = parser.parse_args()

    lon, lat, speed = read_radiants(args.input_h5)
    count, mean_speed = healpix_histogram(lon, lat, speed, args.nside, args.min_count_for_mean_speed)
    write_h5(args.output_h5, args.input_h5, args.nside, count, mean_speed, len(lon))
    plot_healpix(count, args.output_png, len(lon))
    print(f"healpix_radiants {len(lon)} nside {args.nside} occupied_pixels {int(np.sum(count > 0.0))}")
    print(args.output_png)
    print(args.output_h5)


if __name__ == "__main__":
    main()
