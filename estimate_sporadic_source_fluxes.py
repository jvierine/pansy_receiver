#!/usr/bin/env python3
"""Estimate relative sporadic-source rates from the paper radiant sidecar."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import healpy as hp
import numpy as np


DEFAULT_SIDECAR = Path("figs/paper_radiant_results_current/paper_radiant_results.h5")
DEFAULT_REGIONS = Path("figs/sporadic_source_regions_manual.json")
PLOT_CENTER_LONGITUDE_DEG = -90.0


def wrap180(deg):
    return (np.asarray(deg, dtype=np.float64) + 180.0) % 360.0 - 180.0


def wrap360(deg):
    return np.asarray(deg, dtype=np.float64) % 360.0


def physical_lon_from_plot_x(plot_x_deg):
    return wrap360(PLOT_CENTER_LONGITUDE_DEG - plot_x_deg)


def region_mask(region, lon, beta):
    return (
        (np.abs(wrap180(lon - region["center_lon_deg"])) <= region["half_lon_deg"])
        & (beta >= region["beta_min_deg"])
        & (beta <= region["beta_max_deg"])
    )


def load_regions(path):
    if path.exists():
        with path.open("r", encoding="utf-8") as fh:
            payload = json.load(fh)
        regions = payload["regions"]
        apex_remainder = regions["apex"]
        toroidal_masks = ["northern_toroidal", "southern_toroidal"]
        return regions, toroidal_masks, apex_remainder
    regions = {
        "helion": {"center_lon_deg": 0.0, "half_lon_deg": 30.0, "beta_min_deg": -20.0, "beta_max_deg": 20.0},
        "antihelion": {"center_lon_deg": 180.0, "half_lon_deg": 30.0, "beta_min_deg": -20.0, "beta_max_deg": 20.0},
        "apex": {"center_lon_deg": 270.0, "half_lon_deg": 30.0, "beta_min_deg": -25.0, "beta_max_deg": 25.0},
        "narrow_apex": {"center_lon_deg": 270.0, "half_lon_deg": 15.0, "beta_min_deg": -10.0, "beta_max_deg": 10.0},
        "southern_toroidal": {"center_lon_deg": 270.0, "half_lon_deg": 45.0, "beta_min_deg": -65.0, "beta_max_deg": -35.0},
    }
    return regions, ["southern_toroidal"], regions["apex"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sidecar", type=Path, default=DEFAULT_SIDECAR)
    parser.add_argument("--regions", type=Path, default=DEFAULT_REGIONS)
    args = parser.parse_args()

    with h5py.File(args.sidecar, "r") as h5:
        aperture_rate = np.asarray(h5["healpix_zenith_only_rate_h_inv"], dtype=np.float64)
        velocity_rate = np.asarray(h5["healpix_debiased_rate_h_inv"], dtype=np.float64)
        raw = np.asarray(h5["healpix_raw_count"], dtype=np.float64)
        exposure = np.asarray(h5["healpix_radiant_exposure_hours"], dtype=np.float64)
        nside = int(h5.attrs["healpix_nside"])

    pixel = np.arange(hp.nside2npix(nside), dtype=np.int64)
    lon, beta = hp.pix2ang(nside, pixel, lonlat=True)
    valid = exposure > 0.0

    regions_json, toroidal_region_names, _ = load_regions(args.regions)
    masks = {name: region_mask(region, lon, beta) for name, region in regions_json.items()}
    apex_remainder = masks["apex"] & ~masks["narrow_apex"]
    toroidal = np.zeros_like(valid, dtype=bool)
    for name in toroidal_region_names:
        toroidal |= masks[name]

    regions = [
        ("helion", masks["helion"]),
        ("antihelion", masks["antihelion"]),
        ("apex_remainder", apex_remainder),
        ("narrow_apex", masks["narrow_apex"]),
        ("toroidal", toroidal),
    ]
    values = []
    for name, mask in regions:
        use = mask & valid
        values.append(
            (
                name,
                int(np.nansum(raw[use])),
                float(np.nansum(aperture_rate[use])),
                float(np.nansum(velocity_rate[use])),
            )
        )
    total_raw = sum(raw_count for _, raw_count, _, _ in values)
    total_aperture = sum(aperture for _, _, aperture, _ in values)
    total_velocity = sum(velocity for _, _, _, velocity in values)
    print("region,raw_count,aperture_rate_h_inv,velocity_rate_h_inv,raw_percent,aperture_percent,velocity_percent")
    for name, raw_count, aperture, velocity in values:
        print(
            f"{name},{raw_count},{aperture:.6g},{velocity:.6g},"
            f"{100.0 * raw_count / total_raw:.3f},"
            f"{100.0 * aperture / total_aperture:.3f},"
            f"{100.0 * velocity / total_velocity:.3f}"
        )
    print(f"total,{total_raw},{total_aperture:.6g},{total_velocity:.6g},{100.0:.3f},{100.0:.3f},{100.0:.3f}")
    narrow_aperture_fraction = 100.0 * values[3][2] / (values[2][2] + values[3][2])
    narrow_velocity_fraction = 100.0 * values[3][3] / (values[2][3] + values[3][3])
    print(f"narrow_apex_fraction_of_full_apex_aperture,{narrow_aperture_fraction:.3f}")
    print(f"narrow_apex_fraction_of_full_apex_velocity,{narrow_velocity_fraction:.3f}")


if __name__ == "__main__":
    main()
