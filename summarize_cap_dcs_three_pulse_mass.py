#!/usr/bin/env python3
"""Summarize CAP/DCS three-pulse maximum-likelihood initial masses."""

from __future__ import annotations

import argparse
import datetime as dt
import os
from pathlib import Path

import h5py
import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fit-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with h5py.File(args.manifest, "r") as manifest:
        sample_idx = np.asarray(manifest["sample_idx"], dtype=np.int64)
        copied = {
            name: np.asarray(manifest[name])
            for name in (
                "passage_index",
                "epoch_unix",
                "solar_longitude_deg",
                "sun_centered_lon_deg",
                "ecliptic_latitude_deg",
                "geocentric_speed_km_s",
                "kepler",
            )
        }
        passage_names = manifest["passage_names"].asstr()[()]

    n_events = len(sample_idx)
    status = np.zeros(n_events, dtype=np.int8)
    radius_um = np.full(n_events, np.nan)
    mass_kg = np.full(n_events, np.nan)
    profile_grid_radius_um = np.full(n_events, np.nan)
    quality_names = (
        "ew_rms_km",
        "ns_rms_km",
        "up_rms_km",
        "doppler_rms_mps",
        "acceleration_rms_km_s2",
        "variance_weighted_score",
    )
    quality = {name: np.full(n_events, np.nan) for name in quality_names}

    for index, sample in enumerate(sample_idx):
        path = args.fit_dir / f"three_pulse_full_event_{int(sample)}.h5"
        if not path.is_file():
            continue
        try:
            with h5py.File(path, "r") as handle:
                if handle.attrs.get("schema", "") != "pansy.three_pulse_full_event.v1":
                    status[index] = 2
                    continue
                fit = handle["dynamics_refit"]
                event_radius_m = np.asarray(fit["radius_m"], dtype=float)
                event_mass_kg = np.asarray(fit["mass_kg"], dtype=float)
                profile_radius = np.asarray(fit["profile_radius_um"], dtype=float)
                profile_chi2 = np.asarray(fit["profile_chi2"], dtype=float)
                if not len(event_radius_m) or not len(event_mass_kg):
                    status[index] = 2
                    continue
                radius_um[index] = event_radius_m[0] * 1e6
                mass_kg[index] = event_mass_kg[0]
                profile_grid_radius_um[index] = profile_radius[np.nanargmin(profile_chi2)]
                for name in quality_names:
                    quality[name][index] = float(handle["quality"].attrs[name])
                status[index] = 1
        except (OSError, KeyError, ValueError):
            status[index] = 2

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(temporary, "w") as output:
        output.attrs["schema"] = "pansy.cap_dcs_three_pulse_mass_summary.v1"
        output.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        output.attrs["manifest"] = str(args.manifest)
        output.attrs["fit_dir"] = str(args.fit_dir)
        output.attrs["maximum_likelihood_definition"] = (
            "free seven-parameter shrinking-radius trajectory fit; initial mass is "
            "dynamics_refit/mass_kg[0]"
        )
        output.attrs["selection_independence"] = (
            "event membership was fixed before mass fitting and is unchanged in this summary"
        )
        output.attrs["status_codes"] = "0=missing,1=valid,2=invalid"
        output.create_dataset("sample_idx", data=sample_idx)
        output.create_dataset("status", data=status)
        output.create_dataset("maximum_likelihood_initial_radius_um", data=radius_um)
        output.create_dataset("maximum_likelihood_initial_mass_kg", data=mass_kg)
        output.create_dataset("profile_grid_minimum_radius_um", data=profile_grid_radius_um)
        output.create_dataset(
            "passage_names", data=np.asarray(passage_names, dtype=object), dtype=string_dtype
        )
        for name, values in copied.items():
            output.create_dataset(name, data=values)
        for name, values in quality.items():
            output.create_dataset(name, data=values)
    os.replace(temporary, args.output)

    print(
        f"events={n_events} valid={np.count_nonzero(status == 1)} "
        f"missing={np.count_nonzero(status == 0)} invalid={np.count_nonzero(status == 2)}"
    )
    if np.any(status == 1):
        valid_mass = mass_kg[status == 1]
        print(
            "maximum_likelihood_initial_mass_kg "
            f"min={np.min(valid_mass):.6g} median={np.median(valid_mass):.6g} "
            f"max={np.max(valid_mass):.6g}"
        )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
