#!/usr/bin/env python3
"""Collect daily radiant rows from compact PANSY orbit metadata tables."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import numpy as np

from plot_fitted_radiant_distribution import (
    centered_plot_longitude_deg,
    wrap180,
    wrap360,
    PLOT_CENTER_LONGITUDE_DEG,
)


RADIANT_DTYPE = np.dtype(
    [
        ("sample_idx", "i8"),
        ("hypothesis", "S8"),
        ("epoch_unix", "f8"),
        ("combined_score", "f8"),
        ("selection_redchi", "f8"),
        ("first_alt_km", "f8"),
        ("source", "S32"),
        ("radiant_ra_gcrs_deg", "f8"),
        ("radiant_dec_gcrs_deg", "f8"),
        ("radiant_lambda_ecliptic_deg", "f8"),
        ("radiant_beta_ecliptic_deg", "f8"),
        ("sun_lambda_ecliptic_deg", "f8"),
        ("sun_beta_ecliptic_deg", "f8"),
        ("lambda_minus_sun_deg", "f8"),
        ("lambda_minus_sun_signed_deg", "f8"),
        ("plot_longitude_deg", "f8"),
        ("speed_km_s", "f8"),
    ]
)

DAILY_COUNT_DTYPE = np.dtype([("utc_day", "S10"), ("count", "i8")])


def rows_from_events(events: np.ndarray) -> np.ndarray:
    good = np.isfinite(events["radiant_ecliptic_lon_deg"])
    good &= np.isfinite(events["radiant_ecliptic_lat_deg"])
    good &= np.isfinite(events["radiant_sun_ecliptic_lon_deg"])
    good &= np.isfinite(events["v_g_km_s"])
    events = events[good]
    rows = np.zeros(len(events), dtype=RADIANT_DTYPE)
    rows["sample_idx"] = events["sample_idx"]
    rows["hypothesis"] = events["selected_hypothesis"]
    rows["epoch_unix"] = events["sample_idx"].astype(np.float64) / 1e6
    rows["combined_score"] = events["combined_score"]
    rows["selection_redchi"] = events["combined_score"]
    rows["first_alt_km"] = events["initial_detection_height_km"]
    rows["source"] = b"orbit_metadata_table"
    rows["radiant_ra_gcrs_deg"] = events["radiant_ra_deg"]
    rows["radiant_dec_gcrs_deg"] = events["radiant_dec_deg"]
    rows["radiant_lambda_ecliptic_deg"] = wrap360(events["radiant_ecliptic_lon_deg"])
    rows["radiant_beta_ecliptic_deg"] = events["radiant_ecliptic_lat_deg"]
    rows["sun_lambda_ecliptic_deg"] = wrap360(events["radiant_sun_ecliptic_lon_deg"])
    rows["sun_beta_ecliptic_deg"] = events["radiant_sun_ecliptic_lat_deg"]
    lambda_minus_sun = wrap360(rows["radiant_lambda_ecliptic_deg"] - rows["sun_lambda_ecliptic_deg"])
    lambda_minus_sun_signed = wrap180(lambda_minus_sun)
    rows["lambda_minus_sun_deg"] = lambda_minus_sun
    rows["lambda_minus_sun_signed_deg"] = lambda_minus_sun_signed
    rows["plot_longitude_deg"] = centered_plot_longitude_deg(lambda_minus_sun_signed)
    rows["speed_km_s"] = events["v_g_km_s"]
    return rows


def collect(orbit_metadata_dir: Path) -> tuple[np.ndarray, int, int]:
    chunks = []
    files_read = 0
    files_skipped = 0
    for path in sorted(orbit_metadata_dir.glob("**/orbit@*.h5")):
        try:
            with h5py.File(path, "r") as h:
                if "events" not in h:
                    files_skipped += 1
                    continue
                chunks.append(rows_from_events(h["events"][()]))
                files_read += 1
        except OSError:
            files_skipped += 1
    if chunks:
        rows = np.concatenate(chunks)
        rows = rows[np.argsort(rows["sample_idx"])]
    else:
        rows = np.zeros(0, dtype=RADIANT_DTYPE)
    return rows, files_read, files_skipped


def daily_counts(rows: np.ndarray) -> np.ndarray:
    if len(rows) == 0:
        return np.zeros(0, dtype=DAILY_COUNT_DTYPE)
    row_days = np.asarray(
        [
            dt.datetime.fromtimestamp(float(epoch_unix), tz=dt.timezone.utc).strftime("%Y-%m-%d").encode()
            for epoch_unix in rows["epoch_unix"]
        ],
        dtype="S10",
    )
    days, counts = np.unique(row_days, return_counts=True)
    out = np.zeros(len(days), dtype=DAILY_COUNT_DTYPE)
    out["utc_day"] = days
    out["count"] = counts
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-metadata-dir", type=Path, default=Path("data/metadata/orbit"))
    parser.add_argument("--output-h5", type=Path, required=True)
    args = parser.parse_args()

    rows, files_read, files_skipped = collect(args.orbit_metadata_dir)
    args.output_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output_h5, "w") as h:
        h.attrs["script"] = Path(__file__).name
        h.attrs["source"] = "compact hourly orbit metadata events tables"
        h.attrs["radiant_definition"] = "DASST zenith-attraction corrected orbit radiant"
        h.attrs["longitude_convention"] = "lambda_radiant - lambda_sun wrapped to [0, 360) deg; signed copy in [-180, 180)"
        h.attrs["plot_center_longitude_deg"] = PLOT_CENTER_LONGITUDE_DEG
        h.attrs["files_read"] = int(files_read)
        h.attrs["files_skipped"] = int(files_skipped)
        h.create_dataset("radiants", data=rows)
        h.create_dataset("daily_counts", data=daily_counts(rows))
    print(f"daily_radiants {len(rows)} files_read {files_read} files_skipped {files_skipped}")
    print(args.output_h5)


if __name__ == "__main__":
    main()
