#!/usr/bin/env python3
"""Build a fixed eta-Aquariid manifest for three-pulse mass fitting."""

from __future__ import annotations

import argparse
import datetime as dt
import os
from pathlib import Path

import h5py
import numpy as np

from plot_eta_aquariids_snapshot import (
    ETA_BETA_DEG,
    ETA_SC_LON_DEG,
    ETA_SOLAR_LON_DEG,
    ETA_VG_KM_S,
    eta_candidate_mask,
)


DEFAULT_CATALOGUE = (
    Path(__file__).resolve().parent
    / "figs"
    / "paper_refresh_20260729_current"
    / "pansy_keplerian_catalogue.h5"
)


def diagnostics_path(events_dir: Path, sample_idx: int) -> Path:
    day = dt.datetime.fromtimestamp(
        sample_idx / 1e6, tz=dt.timezone.utc
    ).date().isoformat()
    return events_dir / day / f"pansy_disambiguation_diagnostics_{sample_idx}.h5"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, default=DEFAULT_CATALOGUE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--radiant-radius-deg", type=float, default=5.0)
    parser.add_argument("--speed-half-width-kms", type=float, default=10.0)
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

    selected = eta_candidate_mask(
        data, args.radiant_radius_deg, args.speed_half_width_kms
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
            "unchanged eta_candidate_mask selection in plot_eta_aquariids_snapshot.py"
        )
        output.attrs["mass_blinding"] = (
            "selection is fixed before mass fitting; no fitted mass or orbital trend "
            "is used for membership or acceptance"
        )
        output.attrs["selection_solar_center_deg"] = ETA_SOLAR_LON_DEG
        output.attrs["selection_solar_half_width_deg"] = 12.0
        output.attrs["selection_sun_centered_lon_deg"] = ETA_SC_LON_DEG
        output.attrs["selection_ecliptic_latitude_deg"] = ETA_BETA_DEG
        output.attrs["selection_radiant_radius_deg"] = args.radiant_radius_deg
        output.attrs["selection_geocentric_speed_km_s"] = ETA_VG_KM_S
        output.attrs["selection_speed_half_width_kms"] = args.speed_half_width_kms
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
