"""Site-independent inertial conversion and incoming Earth-hyperbola asymptote."""

import numpy as np
import csv
import io
import urllib.request
from astropy.table import QTable
from astropy import units as u
from astropy.coordinates import (
    EarthLocation,
    ITRS,
    GCRS,
    CartesianRepresentation,
    CartesianDifferential,
    SkyCoord,
    GeocentricMeanEcliptic,
    get_sun,
)
from astropy.utils import iers
from astropy.time import Time

MU_EARTH = 3.986004418e14
iers.conf.auto_download = False
iers.conf.iers_degraded_accuracy = "warn"


def configure_eop(epoch, path=None):
    """Refresh Earth orientation explicitly, never as an import side effect."""
    payload = urllib.request.urlopen(
        "https://celestrak.org/SpaceData/EOP-All.csv", timeout=30
    ).read()
    rows = list(csv.DictReader(io.StringIO(payload.decode())))
    mjd = Time(epoch, format="unix").mjd
    rows = [r for r in rows if abs(float(r["MJD"]) - mjd) < 5]
    if len(rows) < 3:
        raise ValueError("EOP data does not cover observation")
    data = QTable()
    for name, field, unit in [
        ("MJD", "MJD", u.day),
        ("PM_x", "X", u.arcsec),
        ("PM_y", "Y", u.arcsec),
        ("UT1_UTC", "UT1-UTC", u.s),
    ]:
        data[name] = np.array([float(r[field]) for r in rows]) * unit
    table = iers.IERS(data)
    iers.earth_orientation_table.set(table)
    if path is not None:
        import h5py

        with h5py.File(path, "w") as h:
            for k in data.colnames:
                h[k] = data[k].value
            h.attrs["source"] = "https://celestrak.org/SpaceData/EOP-All.csv"
            h.attrs["data_types"] = " ".join(r["DATA_TYPE"] for r in rows)
    return table


def enu_to_gcrs(station, state, epoch):
    time = Time(epoch, format="unix", scale="utc")
    origin = EarthLocation.from_geodetic(
        station.longitude_deg * u.deg,
        station.latitude_deg * u.deg,
        station.altitude_m * u.m,
    )
    rotation = station.enu_to_ecef()
    xyz = (
        origin.get_itrs(obstime=time).cartesian.xyz.to_value(u.m) + rotation @ state[:3]
    )
    vel = rotation @ state[3:]
    rep = CartesianRepresentation(
        xyz * u.m, differentials=CartesianDifferential(vel * u.m / u.s)
    )
    c = ITRS(rep, obstime=time).transform_to(GCRS(obstime=time))
    return np.r_[
        c.cartesian.xyz.to_value(u.m),
        c.cartesian.differentials["s"].d_xyz.to_value(u.m / u.s),
    ]


def incoming_asymptote(state, mu=MU_EARTH):
    r, v = np.asarray(state[:3]), np.asarray(state[3:])
    radius = np.linalg.norm(r)
    energy = np.dot(v, v) / 2 - mu / radius
    if energy <= 0:
        raise ValueError("Earth-bound state has no hyperbolic geocentric radiant")
    if np.dot(r, v) >= 0:
        raise ValueError("incoming state required")
    speed = np.sqrt(2 * energy)
    h = np.cross(r, v)
    hn = np.linalg.norm(h)
    if hn < 1e-10 * radius * np.linalg.norm(v):
        return speed * v / np.linalg.norm(v)
    evec = np.cross(v, h) / mu - r / radius
    e = np.linalg.norm(evec)
    if e <= 1:
        raise ValueError("nonhyperbolic eccentricity")
    p = evec / e
    q = np.cross(h / hn, p)
    nu = -np.arccos(-1 / e)
    direction = -np.sin(nu) * p + (e + np.cos(nu)) * q
    return speed * direction / np.linalg.norm(direction)


def radiant(station, state, epoch):
    inertial = enu_to_gcrs(station, state, epoch)
    vinf = incoming_asymptote(inertial)
    direction = -vinf / np.linalg.norm(vinf)
    ra = np.degrees(np.arctan2(direction[1], direction[0])) % 360
    dec = np.degrees(np.arcsin(direction[2]))
    time = Time(epoch, format="unix", scale="utc")
    frame = GeocentricMeanEcliptic(equinox=Time("J2000"), obstime=time)
    sky = SkyCoord(
        ra=ra * u.deg, dec=dec * u.deg, frame=GCRS(obstime=time)
    ).transform_to(frame)
    sun = get_sun(time).transform_to(frame)
    return dict(
        geocentric_speed_km_s=np.linalg.norm(vinf) / 1000,
        ra_deg=ra,
        dec_deg=dec,
        sun_centered_ecliptic_longitude_deg=(sky.lon.deg - sun.lon.deg) % 360,
        ecliptic_latitude_deg=sky.lat.deg,
        solar_longitude_deg=sun.lon.deg % 360,
        incoming_asymptotic_velocity_gcrs_mps=vinf,
        state_gcrs_m_mps=inertial,
    )


def hammer(longitude_deg, latitude_deg, centre_deg=270.0):
    lon = -np.radians((np.asarray(longitude_deg) - centre_deg + 180) % 360 - 180)
    lat = np.radians(latitude_deg)
    den = np.sqrt(1 + np.cos(lat) * np.cos(lon / 2))
    return np.stack(
        [
            2 * np.sqrt(2) * np.cos(lat) * np.sin(lon / 2) / den,
            np.sqrt(2) * np.sin(lat) / den,
        ],
        axis=-1,
    )
