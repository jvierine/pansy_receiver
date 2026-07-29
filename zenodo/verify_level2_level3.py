#!/usr/bin/env python3
"""Derive GCRS states and heliocentric orbits from exported Level 2 data."""

import argparse
from pathlib import Path

import h5py
import numpy as np
from astropy.constants import GM_sun
from astropy.coordinates import (
    AltAz,
    CartesianDifferential,
    CartesianRepresentation,
    EarthLocation,
    GCRS,
    SkyCoord,
    get_body_barycentric_posvel,
)
from astropy.time import Time
import astropy.units as u
from astropy.utils import iers


SITE = EarthLocation(
    lat=-69.010833 * u.deg,
    lon=39.599722 * u.deg,
    height=100.0 * u.m,
)


def fit_state(time_s, position_enu_km):
    """Fit p(t) = p0 + v*t."""
    matrix = np.column_stack((np.ones(len(time_s)), time_s))
    state = np.linalg.lstsq(matrix, position_enu_km, rcond=None)[0]
    residual = position_enu_km - matrix @ state
    rms_km = np.sqrt(np.mean(np.sum(residual**2, axis=1)))
    return state[0], state[1], rms_km


def enu_to_gcrs(position_enu_km, velocity_enu_km_s, epoch):
    # AltAz Cartesian axes are north, east, up.
    position_neu = position_enu_km[[1, 0, 2]]
    velocity_neu = velocity_enu_km_s[[1, 0, 2]]
    cartesian = CartesianRepresentation(
        position_neu * u.km,
        differentials=CartesianDifferential(velocity_neu * u.km / u.s),
    )
    local = SkyCoord(cartesian, frame=AltAz(obstime=epoch, location=SITE))
    gcrs = local.transform_to(GCRS(obstime=epoch)).cartesian
    return (
        gcrs.xyz.to_value(u.km),
        gcrs.differentials["s"].d_xyz.to_value(u.km / u.s),
    )


def heliocentric_ecliptic_state(position_gcrs_km, velocity_gcrs_km_s, epoch):
    earth_p, earth_v = get_body_barycentric_posvel("earth", epoch)
    sun_p, sun_v = get_body_barycentric_posvel("sun", epoch)
    position = (earth_p - sun_p).xyz.to_value(u.km) + position_gcrs_km
    velocity = (earth_v - sun_v).xyz.to_value(u.km / u.s) + velocity_gcrs_km_s

    # ICRS equatorial to mean J2000 ecliptic.
    epsilon = np.deg2rad(23.4392911)
    rotation = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(epsilon), np.sin(epsilon)],
            [0.0, -np.sin(epsilon), np.cos(epsilon)],
        ]
    )
    return rotation @ position, rotation @ velocity


def kepler_elements(position_km, velocity_km_s):
    mu = GM_sun.to_value(u.km**3 / u.s**2)
    radius = np.linalg.norm(position_km)
    speed2 = np.dot(velocity_km_s, velocity_km_s)
    angular_momentum = np.cross(position_km, velocity_km_s)
    h = np.linalg.norm(angular_momentum)
    node = np.cross([0.0, 0.0, 1.0], angular_momentum)
    n = np.linalg.norm(node)
    eccentricity_vector = (
        np.cross(velocity_km_s, angular_momentum) / mu - position_km / radius
    )
    eccentricity = np.linalg.norm(eccentricity_vector)
    semi_major_axis_km = 1.0 / (2.0 / radius - speed2 / mu)
    inclination = np.arccos(np.clip(angular_momentum[2] / h, -1.0, 1.0))
    ascending_node = np.arctan2(node[1], node[0]) % (2.0 * np.pi)
    argument_perihelion = np.arctan2(
        np.dot(np.cross(node, eccentricity_vector), angular_momentum) / (n * eccentricity * h),
        np.dot(node, eccentricity_vector) / (n * eccentricity),
    ) % (2.0 * np.pi)
    true_anomaly = np.arctan2(
        np.dot(np.cross(eccentricity_vector, position_km), angular_momentum)
        / (eccentricity * radius * h),
        np.dot(eccentricity_vector, position_km) / (eccentricity * radius),
    ) % (2.0 * np.pi)
    a_au = (semi_major_axis_km * u.km).to_value(u.AU)
    return np.array(
        [
            a_au,
            eccentricity,
            np.rad2deg(inclination),
            np.rad2deg(ascending_node),
            np.rad2deg(argument_perihelion),
            np.rad2deg(true_anomaly),
            a_au * (1.0 - eccentricity),
        ]
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("release", type=Path)
    parser.add_argument("--month", default="2025-07")
    parser.add_argument("--count", type=int, default=5)
    args = parser.parse_args()
    level2_name = args.release / "level2" / f"pansy_level2_{args.month}.h5"
    level3_name = args.release / "level3" / f"pansy_level3_{args.month}.h5"
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None

    with h5py.File(level2_name, "r") as level2, h5py.File(level3_name, "r") as level3:
        level3_ids = level3["events/event_id"][()]
        level3_row = {int(event_id): i for i, event_id in enumerate(level3_ids)}
        event_ids = level2["events/event_id"][()]
        starts = level2["events/measurement_start"][()]
        counts = level2["events/measurement_count"][()]
        keep_all = level2["measurements/selection_keep"][()]
        time_all = level2["measurements/time_offset_s"][()]
        position_all = np.column_stack(
            [
                level2[f"measurements/{axis}_km"][()]
                for axis in ("east", "north", "up")
            ]
        )
        print("event_id                 RMS  quantity       Level2 fit    Level3")
        shown = 0
        for event_id, start, count in zip(event_ids, starts, counts):
            event_id, start, count = int(event_id), int(start), int(count)
            if count < 20 or event_id not in level3_row:
                continue
            rows = slice(start, start + count)
            keep = keep_all[rows]
            if np.count_nonzero(keep) < 20:
                continue
            time_s = time_all[rows][keep]
            position = position_all[rows][keep]
            position0, velocity0, rms_km = fit_state(time_s, position)
            speed = np.linalg.norm(velocity0)
            if not (45.0 <= speed <= 75.0 and rms_km < 0.25):
                continue

            epoch = Time(event_id / 1e6, format="unix", scale="utc")
            position_gcrs, velocity_gcrs = enu_to_gcrs(position0, velocity0, epoch)
            position_helio, velocity_helio = heliocentric_ecliptic_state(
                position_gcrs, velocity_gcrs, epoch
            )
            elements = kepler_elements(position_helio, velocity_helio)
            # Near-parabolic orbits make a extremely sensitive to tiny velocity
            # differences; use ordinary bound cases for this compact check.
            if not (0.0 < elements[0] < 10.0):
                continue
            j = level3_row[event_id]
            level3_elements = np.array(
                [
                    level3[f"events/{name}"][j]
                    for name in (
                        "semi_major_axis_AU",
                        "eccentricity",
                        "inclination_deg",
                        "longitude_ascending_node_deg",
                        "argument_perihelion_deg",
                        "true_anomaly_deg",
                        "perihelion_distance_AU",
                    )
                ]
            )
            print(f"{event_id:22d} {rms_km:5.3f}  v0 (km/s)   {speed:11.4f} {level3['events/v0_km_s'][j]:9.4f}")
            for name, fit_value, catalogue_value in zip(
                ("a (AU)", "e", "i (deg)", "Omega", "omega", "nu", "q (AU)"),
                elements,
                level3_elements,
            ):
                print(f"{'':29s}  {name:11s} {fit_value:11.4f} {catalogue_value:9.4f}")
            shown += 1
            if shown == args.count:
                break


if __name__ == "__main__":
    main()
