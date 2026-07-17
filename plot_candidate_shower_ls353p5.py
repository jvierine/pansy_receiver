#!/usr/bin/env python3
"""Plot and summarize the Radview candidate shower near solar longitude 353.5 deg."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from scipy.stats import norm, poisson

from healpix_hammer import render_healpix_hammer
from plot_paper_radiant_results import PLOT_CENTER_LONGITUDE_DEG, add_source_markers, style_hammer
from radiant_visibility import ecliptic_lonlat_to_radec_deg


DEFAULT_SELECTION = {
    "mode": "radiant-pixels",
    "solarCenterDeg": 353.5,
    "solarWindowDeg": 4.0,
    "solarRangesDeg": [[351.5, 355.5]],
    "healpixOrder": 5,
    "healpixNside": 32,
    "healpixNpix": 12288,
    "selectedRadiantPixelIds": [
        11150,
        11151,
        11243,
        11244,
        11331,
        11332,
        11333,
        11334,
        11418,
        11419,
        11500,
        11577,
    ],
    "distributionFilters": [
        {
            "key": "eccentricity",
            "label": "Eccentricity",
            "unit": "",
            "range": [0.6554599, 0.89862396],
        }
    ],
}
DEFAULT_EXPORTED_COUNT = 38


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--selection-json",
        type=Path,
        default=None,
        help="Optional Radview selection export; defaults to the archived selection in this script.",
    )
    parser.add_argument(
        "--catalog-h5",
        type=Path,
        default=Path("figs/pansy_maarsy_keplerian_catalogue.h5"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("figs/paper_candidate_shower_ls353p5.png"),
    )
    parser.add_argument(
        "--summary-h5",
        type=Path,
        default=Path("figs/paper_candidate_shower_ls353p5.h5"),
    )
    return parser.parse_args()


def circular_interval_mask(values_deg: np.ndarray, ranges_deg: list[list[float]]) -> np.ndarray:
    values = np.asarray(values_deg, dtype=np.float64) % 360.0
    keep = np.zeros(values.shape, dtype=bool)
    for start, stop in ranges_deg:
        start = float(start) % 360.0
        stop = float(stop) % 360.0
        if start <= stop:
            keep |= (values >= start) & (values <= stop)
        else:
            keep |= (values >= start) | (values <= stop)
    return keep


def load_background(catalog_h5: Path, selection: dict) -> dict[str, np.ndarray]:
    filters = {item["key"]: item for item in selection.get("distributionFilters", [])}
    if "eccentricity" not in filters:
        raise ValueError("selection export has no eccentricity filter")
    eccentricity_range = np.asarray(filters["eccentricity"]["range"], dtype=np.float64)

    with h5py.File(catalog_h5, "r") as h5:
        solar_longitude = np.asarray(h5["solar_longitude_deg"], dtype=np.float64)
        eccentricity = np.asarray(h5["e"], dtype=np.float64)
        source_id = np.asarray(h5["source_id"], dtype=np.int8)
        keep = (
            (source_id == 0)
            & circular_interval_mask(solar_longitude, selection["solarRangesDeg"])
            & (eccentricity >= eccentricity_range[0])
            & (eccentricity <= eccentricity_range[1])
        )
        return {
            "epoch_unix": np.asarray(h5["epoch_unix"][keep], dtype=np.float64),
            "solar_longitude_deg": solar_longitude[keep],
            "sun_centered_lon_deg": np.asarray(h5["sun_centered_lon_deg"][keep], dtype=np.float64),
            "ecliptic_lat_deg": np.asarray(h5["ecliptic_lat_deg"][keep], dtype=np.float64),
            "vg_km_s": np.asarray(h5["vg_km_s"][keep], dtype=np.float64),
            "a_AU": np.asarray(h5["a_AU"][keep], dtype=np.float64),
            "e": eccentricity[keep],
            "i_deg": np.asarray(h5["i_deg"][keep], dtype=np.float64),
            "q_AU": np.asarray(h5["q_AU"][keep], dtype=np.float64),
        }


def outer_ring_pixels(nside: int, selected_pixels: np.ndarray) -> np.ndarray:
    selected = set(int(pixel) for pixel in selected_pixels)
    radius = 2.05 * hp.max_pixrad(nside)
    ring: set[int] = set()
    for pixel in selected:
        neighbors = hp.query_disc(
            nside,
            hp.pix2vec(nside, pixel, nest=False),
            radius,
            inclusive=True,
            nest=False,
        )
        ring.update(int(neighbor) for neighbor in neighbors if int(neighbor) not in selected)
    return np.asarray(sorted(ring), dtype=np.int64)


def spherical_mean(lon_deg: np.ndarray, lat_deg: np.ndarray) -> tuple[float, float]:
    lon = np.deg2rad(np.asarray(lon_deg, dtype=np.float64))
    lat = np.deg2rad(np.asarray(lat_deg, dtype=np.float64))
    vector = np.asarray(
        [
            np.mean(np.cos(lat) * np.cos(lon)),
            np.mean(np.cos(lat) * np.sin(lon)),
            np.mean(np.sin(lat)),
        ]
    )
    vector /= np.linalg.norm(vector)
    return float(np.rad2deg(np.arctan2(vector[1], vector[0])) % 360.0), float(np.rad2deg(np.arcsin(vector[2])))


def main() -> None:
    args = parse_args()
    export = (
        json.loads(args.selection_json.read_text())
        if args.selection_json is not None
        else {"selection": DEFAULT_SELECTION, "count": DEFAULT_EXPORTED_COUNT}
    )
    selection = export["selection"]
    nside = int(selection["healpixNside"])
    selected_pixels = np.asarray(selection["selectedRadiantPixelIds"], dtype=np.int64)
    rows = load_background(args.catalog_h5, selection)

    pixel = hp.ang2pix(
        nside,
        rows["sun_centered_lon_deg"],
        rows["ecliptic_lat_deg"],
        lonlat=True,
        nest=False,
    )
    counts = np.bincount(pixel, minlength=hp.nside2npix(nside)).astype(np.float64)
    selected_mask = np.isin(pixel, selected_pixels)
    selected_count = int(np.sum(selected_mask))
    if selected_count != int(export["count"]):
        raise RuntimeError(f"catalogue selection gives {selected_count} meteors, export contains {export['count']}")

    ring_pixels = outer_ring_pixels(nside, selected_pixels)
    background_count = int(np.sum(counts[ring_pixels]))
    expected_background = background_count * len(selected_pixels) / len(ring_pixels)
    p_value = float(poisson.sf(selected_count - 1, expected_background))
    sigma = float(norm.isf(p_value))

    selected = {key: values[selected_mask] for key, values in rows.items()}
    mean_sc_lon, mean_beta = spherical_mean(selected["sun_centered_lon_deg"], selected["ecliptic_lat_deg"])
    mean_solar_longitude = float(np.mean(selected["solar_longitude_deg"]))
    mean_absolute_lon = (mean_sc_lon + mean_solar_longitude) % 360.0
    mean_ra, mean_dec = ecliptic_lonlat_to_radec_deg(mean_absolute_lon, mean_beta)

    fig = plt.figure(figsize=(8.2, 4.7))
    ax = fig.add_subplot(111, projection="hammer")
    vmax = max(float(np.max(counts)), 1.0)
    image = render_healpix_hammer(
        ax,
        counts,
        nside,
        cmap="magma",
        norm=Normalize(vmin=0.0, vmax=vmax),
        center_longitude_deg=PLOT_CENTER_LONGITUDE_DEG,
    )
    style_hammer(ax)
    add_source_markers(ax)
    ax.text(
        0.80,
        0.82,
        (
            rf"$N={selected_count}$" "\n"
            rf"$\overline{{v_g}}={np.mean(selected['vg_km_s']):.1f}\,\mathrm{{km\,s^{{-1}}}}$" "\n"
            rf"Poisson significance: ${sigma:.1f}\sigma$"
        ),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        color="white",
    )
    colorbar = fig.colorbar(image, ax=ax, orientation="horizontal", pad=0.08, fraction=0.055, aspect=35)
    colorbar.set_label(rf"Meteor count per HEALPix pixel ($N_{{\mathrm{{side}}}}={nside}$)")
    fig.subplots_adjust(left=0.035, right=0.985, top=0.985, bottom=0.15)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=240, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    args.summary_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.summary_h5, "w") as h5:
        h5.attrs["description"] = "PANSY candidate shower radiant selection near solar longitude 353.5 deg"
        h5.attrs["selection_json"] = str(args.selection_json) if args.selection_json is not None else "embedded default"
        h5.attrs["catalog_h5"] = str(args.catalog_h5)
        h5.attrs["healpix_ordering"] = "RING"
        h5.create_dataset("healpix_counts", data=counts, compression="gzip")
        h5.create_dataset("selected_pixel_ids", data=selected_pixels)
        h5.create_dataset("outer_ring_pixel_ids", data=ring_pixels)
        summary = h5.create_group("summary")
        for key, value in {
            "selected_count": selected_count,
            "background_count": background_count,
            "background_pixel_count": len(ring_pixels),
            "expected_background": expected_background,
            "poisson_p_value": p_value,
            "poisson_sigma": sigma,
            "mean_solar_longitude_deg": mean_solar_longitude,
            "mean_sun_centered_lon_deg": mean_sc_lon,
            "mean_ecliptic_lat_deg": mean_beta,
            "mean_ra_deg": float(np.asarray(mean_ra)),
            "mean_dec_deg": float(np.asarray(mean_dec)),
            "mean_vg_km_s": float(np.mean(selected["vg_km_s"])),
            "std_vg_km_s": float(np.std(selected["vg_km_s"], ddof=1)),
            "mean_a_AU": float(np.mean(selected["a_AU"])),
            "mean_e": float(np.mean(selected["e"])),
            "mean_i_deg": float(np.mean(selected["i_deg"])),
            "mean_q_AU": float(np.mean(selected["q_AU"])),
            "first_epoch_unix": float(np.min(selected["epoch_unix"])),
            "last_epoch_unix": float(np.max(selected["epoch_unix"])),
        }.items():
            summary.create_dataset(key, data=value)
        selected_group = h5.create_group("selected_meteors")
        for key, values in selected.items():
            selected_group.create_dataset(key, data=values, compression="gzip")

    print(f"wrote {args.output}")
    print(f"wrote {args.summary_h5}")
    print(
        f"N={selected_count}; background={expected_background:.3f}; p={p_value:.3e}; "
        f"significance={sigma:.2f} sigma"
    )
    print(
        f"mean radiant: lambda'={mean_sc_lon:.2f} deg, beta={mean_beta:.2f} deg, "
        f"RA={float(np.asarray(mean_ra)):.2f} deg, Dec={float(np.asarray(mean_dec)):.2f} deg, "
        f"vg={np.mean(selected['vg_km_s']):.2f} km/s"
    )


if __name__ == "__main__":
    main()
