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
DEFAULT_FILTER_SPLIT = Path(
    "figs/paper_radiant_results_current/radiant_spherical_harmonic_split.h5"
)
DEFAULT_OUTPUT = Path(
    "figs/paper_radiant_results_current/sporadic_source_flux_estimates.h5"
)
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
    parser.add_argument("--filter-split", type=Path, default=DEFAULT_FILTER_SPLIT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    with h5py.File(args.sidecar, "r") as h5:
        aperture_rate = np.asarray(h5["healpix_zenith_only_rate_h_inv"], dtype=np.float64)
        velocity_rate = np.asarray(h5["healpix_debiased_rate_h_inv"], dtype=np.float64)
        raw = np.asarray(h5["healpix_raw_count"], dtype=np.float64)
        exposure = np.asarray(h5["healpix_radiant_exposure_hours"], dtype=np.float64)
        nside = int(h5.attrs["healpix_nside"])

    with h5py.File(args.filter_split, "r") as h5:
        split_nside = int(h5.attrs["healpix_nside"])
        low_frequency = np.asarray(h5["low_frequency_count"], dtype=np.float64)
        positive_high_frequency = np.asarray(
            h5["positive_high_frequency_count"], dtype=np.float64
        )
        filter_parameters = {
            key: h5.attrs[key]
            for key in (
                "zonal_pass",
                "zonal_taper",
                "meridional_pass",
                "meridional_taper",
                "positive_iterations",
                "positive_fraction",
                "lmax",
            )
        }
    if split_nside != nside or low_frequency.shape != raw.shape:
        raise ValueError("spherical-harmonic split does not match the radiant sidecar")

    pixel = np.arange(hp.nside2npix(nside), dtype=np.int64)
    lon, beta = hp.pix2ang(nside, pixel, lonlat=True)
    valid = exposure > 0.0

    regions_json, _, _ = load_regions(args.regions)
    masks = {name: region_mask(region, lon, beta) for name, region in regions_json.items()}

    component_sum = low_frequency + positive_high_frequency
    low_fraction = np.divide(
        low_frequency,
        component_sum,
        out=np.ones_like(low_frequency),
        where=component_sum > 0.0,
    )
    high_fraction = np.divide(
        positive_high_frequency,
        component_sum,
        out=np.zeros_like(positive_high_frequency),
        where=component_sum > 0.0,
    )

    maps = (raw, aperture_rate, velocity_rate)
    values = []

    def append_components(name, mask):
        use = mask & valid
        values.extend(
            (
                (
                    f"{name}_smooth_lowpass",
                    *(
                        float(np.nansum(data[use] * low_fraction[use]))
                        for data in maps
                    ),
                ),
                (
                    f"{name}_structured_highpass",
                    *(
                        float(np.nansum(data[use] * high_fraction[use]))
                        for data in maps
                    ),
                ),
            )
        )

    append_components("helion", masks["helion"])
    append_components("antihelion", masks["antihelion"])
    append_components("apex", masks["apex"])
    append_components("southern_toroidal", masks["southern_toroidal"])

    source_names = [name for name, *_ in values]
    source_values = np.asarray([row[1:] for row in values], dtype=np.float64)
    totals = np.sum(source_values, axis=0)
    percentages = 100.0 * source_values / totals[np.newaxis, :]
    helion_total = source_values[0] + source_values[1]
    antihelion_total = source_values[2] + source_values[3]
    antihelion_helion_ratio = antihelion_total / helion_total
    apex_component_fraction = (
        100.0
        * source_values[5]
        / (source_values[4] + source_values[5])
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output, "w") as h5:
        h5.attrs["schema"] = "sporadic_source_flux_estimates_v2"
        h5.attrs["source_sidecar"] = str(args.sidecar)
        h5.attrs["source_filter_split"] = str(args.filter_split)
        h5.attrs["source_regions"] = str(args.regions)
        h5.attrs["apex_decomposition"] = (
            "pixelwise normalized low-frequency and positive high-frequency "
            "templates inside each source aperture"
        )
        for key, value in filter_parameters.items():
            h5.attrs[key] = value
        h5.create_dataset(
            "source_name",
            data=np.asarray(source_names, dtype=h5py.string_dtype("utf-8")),
        )
        h5.create_dataset(
            "source_value",
            data=source_values.astype(np.float32),
        )
        h5["source_value"].attrs["columns"] = (
            "raw_count, aperture_rate_h_inv, velocity_rate_h_inv"
        )
        h5.create_dataset("source_percent", data=percentages.astype(np.float32))
        h5.create_dataset(
            "antihelion_helion_ratio",
            data=antihelion_helion_ratio.astype(np.float32),
        )
        h5.create_dataset(
            "narrow_fraction_of_full_apex_percent",
            data=apex_component_fraction.astype(np.float32),
        )

    print(
        f"{'Region':<25} {'Raw':>10} {'Aperture':>12} {'Aperture+vg':>14}"
        f" {'Raw %':>8} {'Aperture %':>12} {'Aperture+vg %':>15}"
    )
    for row_index, (name, raw_count, aperture, velocity) in enumerate(values):
        print(
            f"{name:<25} {raw_count:10.1f} {aperture:12.6g} {velocity:14.6g}"
            f" {percentages[row_index, 0]:8.3f}"
            f" {percentages[row_index, 1]:12.3f}"
            f" {percentages[row_index, 2]:15.3f}"
        )
    print(
        f"{'total':<25} {totals[0]:10.1f} {totals[1]:12.6g} {totals[2]:14.6g}"
        f" {100.0:8.3f} {100.0:12.3f} {100.0:15.3f}"
    )
    print(
        "Narrow-apex fraction of full apex "
        f"(raw/aperture/velocity): {apex_component_fraction[0]:.3f}% / "
        f"{apex_component_fraction[1]:.3f}% / "
        f"{apex_component_fraction[2]:.3f}%"
    )
    print(args.output)


if __name__ == "__main__":
    main()
