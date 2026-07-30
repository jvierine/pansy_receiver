#!/usr/bin/env python3
"""Fit a straight trajectory and apparent GCRS radiant from one Level 2 event."""

import argparse
from pathlib import Path

import astropy.units as u
import h5py
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers


SITE = EarthLocation(
    lat=-69.010833 * u.deg,
    lon=39.599722 * u.deg,
    height=100.0 * u.m,
)


def fit_event(level2, requested_event_id=None):
    event_ids = level2["events/event_id"][()]
    starts = level2["events/measurement_start"][()]
    counts = level2["events/measurement_count"][()]
    if requested_event_id is None:
        candidates = range(len(event_ids))
    else:
        candidates = np.flatnonzero(event_ids == requested_event_id)

    for index in candidates:
        start = int(starts[index])
        stop = start + int(counts[index])
        keep = level2["measurements/selection_keep"][start:stop]
        if np.count_nonzero(keep) < 3:
            continue
        time_s = level2["measurements/time_offset_s"][start:stop][keep]
        position_enu_km = np.column_stack(
            [
                level2[f"measurements/{axis}_km"][start:stop][keep]
                for axis in ("east", "north", "up")
            ]
        )
        design = np.column_stack((np.ones(len(time_s)), time_s))
        state = np.linalg.lstsq(design, position_enu_km, rcond=None)[0]
        residual = position_enu_km - design @ state
        rms_km = np.sqrt(np.mean(np.sum(residual**2, axis=1)))
        return int(event_ids[index]), state[0], state[1], rms_km, len(time_s)
    raise RuntimeError("no event with at least three retained measurements")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("level2_h5", type=Path)
    parser.add_argument("--event-id", type=int)
    args = parser.parse_args()
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None

    with h5py.File(args.level2_h5, "r") as level2:
        event_id, position0, velocity, rms_km, count = fit_event(
            level2, args.event_id
        )

    speed = np.linalg.norm(velocity)
    radiant_enu = -velocity / speed
    azimuth = np.arctan2(radiant_enu[0], radiant_enu[1]) * u.rad
    altitude = np.arcsin(radiant_enu[2]) * u.rad
    epoch = Time(event_id / 1e6, format="unix", scale="utc")
    radiant = SkyCoord(
        az=azimuth,
        alt=altitude,
        frame=AltAz(obstime=epoch, location=SITE),
    ).transform_to("gcrs")

    print(f"event_id: {event_id}")
    print(f"epoch_utc: {epoch.utc.isot}")
    print(f"retained_measurements: {count}")
    print(f"straight_line_rms_km: {rms_km:.4f}")
    print(f"speed_km_s: {speed:.4f}")
    print(f"apparent_radiant_ra_deg: {radiant.ra.deg:.4f}")
    print(f"apparent_radiant_dec_deg: {radiant.dec.deg:.4f}")


if __name__ == "__main__":
    main()
