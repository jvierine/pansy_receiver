#!/usr/bin/env python3
"""Export the flat PANSY Keplerian table used by shower figure scripts."""

import argparse
from pathlib import Path

import h5py
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("radiants_h5", type=Path)
    parser.add_argument("output_h5", type=Path)
    args = parser.parse_args()

    with h5py.File(args.radiants_h5, "r") as source:
        rows = source["radiants"][()]

    args.output_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output_h5, "w") as output:
        output.attrs["description"] = (
            "PANSY Keplerian catalogue derived from the canonical radiant table"
        )
        output.attrs["source_radiant_snapshot"] = str(args.radiants_h5)
        output.attrs["kepler_column_order"] = (
            "a_AU,e,i_deg,Omega_deg,omega_deg,nu_deg,q_AU"
        )
        datasets = {
            "epoch_unix": rows["epoch_unix"],
            "solar_longitude_deg": rows["sun_lambda_ecliptic_deg"],
            "vg_km_s": rows["speed_km_s"],
            "sun_centered_lon_deg": rows["lambda_minus_sun_deg"],
            "ecliptic_lon_deg": rows["radiant_lambda_ecliptic_deg"],
            "ecliptic_lat_deg": rows["radiant_beta_ecliptic_deg"],
            "kepler": rows["kepler"],
            "source_id": np.zeros(len(rows), dtype=np.int8),
        }
        for name, values in datasets.items():
            output.create_dataset(
                name,
                data=values,
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )

    print(f"events {len(rows)}")
    print(args.output_h5)


if __name__ == "__main__":
    main()
