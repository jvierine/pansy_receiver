#!/usr/bin/env python3
"""Build a fixed eta-Aquariid manifest for three-pulse mass fitting."""

from __future__ import annotations

import argparse
import datetime as dt
import os
from pathlib import Path

import h5py
import numpy as np

from plot_omega_eridanids_shower import ang2pix_ring


DEFAULT_CATALOGUE = (
    Path(__file__).resolve().parent
    / "figs"
    / "paper_refresh_20260729_current"
    / "pansy_keplerian_catalogue.h5"
)
ETA_SELECTION_SOLAR_CENTER_DEG = 71.20
ETA_SELECTION_SOLAR_WINDOW_DEG = 53.93
ETA_SELECTION_NSIDE = 64
ETA_SELECTION_PIXELS = np.asarray(
    (
        20045, 20046, 20047, 20301, 20302, 20303, 20304, 20305,
        20556, 20557, 20558, 20559, 20560, 20561, 20562, 20813,
        20814, 20815, 20816, 20817, 20818, 20819, 21068, 21069,
        21070, 21071, 21072, 21073, 21074, 21075, 21325, 21326,
        21327, 21328, 21329, 21330, 21331, 21332, 21581, 21582,
        21583, 21584, 21585, 21586, 21587, 21838, 21839, 21840,
        21841, 21842, 21843, 21844, 22095, 22096, 22097, 22098,
        22099, 22354, 22355, 22356, 22611, 22612,
    ),
    dtype=np.int64,
)


def wrap180(degrees: np.ndarray) -> np.ndarray:
    return (np.asarray(degrees, dtype=float) + 180.0) % 360.0 - 180.0


def diagnostics_path(events_dir: Path, sample_idx: int) -> Path:
    day = dt.datetime.fromtimestamp(
        sample_idx / 1e6, tz=dt.timezone.utc
    ).date().isoformat()
    return events_dir / day / f"pansy_disambiguation_diagnostics_{sample_idx}.h5"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, default=DEFAULT_CATALOGUE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--events-dir",
        type=Path,
        help="If supplied, require a diagnostics file for every selected event.",
    )
    args = parser.parse_args()

    field_names = (
        "epoch_unix",
        "solar_longitude_deg",
        "sun_centered_lon_deg",
        "ecliptic_lat_deg",
        "vg_km_s",
    )
    with h5py.File(args.catalogue, "r") as handle:
        data = {
            name: np.asarray(handle[name], dtype=np.float64) for name in field_names
        }
        kepler = np.asarray(handle["kepler"], dtype=np.float64)

    radiant_pixel = ang2pix_ring(
        ETA_SELECTION_NSIDE,
        data["sun_centered_lon_deg"],
        data["ecliptic_lat_deg"],
    )
    selected = (
        np.isin(radiant_pixel, ETA_SELECTION_PIXELS)
        & (
            np.abs(
                wrap180(
                    data["solar_longitude_deg"] - ETA_SELECTION_SOLAR_CENTER_DEG
                )
            )
            <= ETA_SELECTION_SOLAR_WINDOW_DEG
        )
    )
    sample_idx = np.rint(data["epoch_unix"][selected] * 1e6).astype(np.int64)
    if len(np.unique(sample_idx)) != len(sample_idx):
        raise RuntimeError("ETA selection contains duplicate event identifiers")

    diagnostics = None
    if args.events_dir is not None:
        diagnostics = np.asarray(
            [
                str(diagnostics_path(args.events_dir, int(sample)))
                for sample in sample_idx
            ],
            dtype=object,
        )
        missing = [path for path in diagnostics if not Path(path).is_file()]
        if missing:
            preview = "\n".join(missing[:10])
            raise FileNotFoundError(
                f"missing diagnostics for {len(missing)} of {len(sample_idx)} ETA events:\n"
                f"{preview}"
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(temporary, "w") as output:
        output.attrs["schema"] = "pansy.eta_three_pulse_mass_manifest.v1"
        output.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        output.attrs["catalogue"] = str(args.catalogue)
        output.attrs["selection_source"] = (
            "user-supplied fixed nside=64 HEALPix radiant pixels and solar-longitude window"
        )
        output.attrs["mass_blinding"] = (
            "selection is fixed before mass fitting; no fitted mass or orbital trend "
            "is used for membership or acceptance"
        )
        output.attrs["selection_solar_center_deg"] = ETA_SELECTION_SOLAR_CENTER_DEG
        output.attrs["selection_solar_half_width_deg"] = ETA_SELECTION_SOLAR_WINDOW_DEG
        output.attrs["selection_healpix_ordering"] = "RING"
        output.attrs["selection_healpix_nside"] = ETA_SELECTION_NSIDE
        output.create_dataset("selection_healpix_pixel_ids", data=ETA_SELECTION_PIXELS)
        output.create_dataset("sample_idx", data=sample_idx)
        output.create_dataset(
            "passage_index", data=np.zeros(len(sample_idx), dtype=np.int8)
        )
        output.create_dataset(
            "passage_names",
            data=np.asarray([r"$\eta$-Aquariids (ETA)"], dtype=object),
            dtype=string_dtype,
        )
        output.create_dataset("epoch_unix", data=data["epoch_unix"][selected])
        output.create_dataset(
            "solar_longitude_deg", data=data["solar_longitude_deg"][selected]
        )
        output.create_dataset(
            "sun_centered_lon_deg", data=data["sun_centered_lon_deg"][selected]
        )
        output.create_dataset(
            "ecliptic_latitude_deg", data=data["ecliptic_lat_deg"][selected]
        )
        output.create_dataset(
            "geocentric_speed_km_s", data=data["vg_km_s"][selected]
        )
        output.create_dataset("kepler", data=kepler[selected])
        output["kepler"].attrs["columns"] = (
            "a_AU,e,i_deg,Omega_deg,omega_deg,nu_deg,q_AU"
        )
        if diagnostics is not None:
            output.create_dataset(
                "diagnostics_h5", data=diagnostics, dtype=string_dtype
            )
    os.replace(temporary, args.output)

    years = np.asarray(
        [
            dt.datetime.fromtimestamp(float(epoch), tz=dt.timezone.utc).year
            for epoch in data["epoch_unix"][selected]
        ]
    )
    print(f"manifest={args.output} events={len(sample_idx)}")
    for year in np.unique(years):
        print(f"year={year} events={np.count_nonzero(years == year)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
