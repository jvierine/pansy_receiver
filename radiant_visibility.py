"""Radiant visibility boundary overlays for PANSY radiant maps."""

from __future__ import annotations

import numpy as np


PANSY_SITE_LATITUDE_DEG = -69.010833
PANSY_SITE_LONGITUDE_DEG = 39.599722
PLOT_CENTER_LONGITUDE_DEG = -90.0
MEAN_OBLIQUITY_DEG = 23.4392911


def wrap180(deg):
    return (np.asarray(deg, dtype=np.float64) + 180.0) % 360.0 - 180.0


def wrap360(deg):
    return np.asarray(deg, dtype=np.float64) % 360.0


def centered_plot_longitude_deg(lambda_minus_sun_signed_deg):
    return -wrap180(np.asarray(lambda_minus_sun_signed_deg, dtype=np.float64) - PLOT_CENTER_LONGITUDE_DEG)


def horizon_declination_limit_deg(site_latitude_deg: float = PANSY_SITE_LATITUDE_DEG) -> float:
    """Declination boundary for nonzero above-horizon time at a site latitude."""
    site_latitude_deg = float(site_latitude_deg)
    if site_latitude_deg < 0.0:
        return 90.0 + site_latitude_deg
    return site_latitude_deg - 90.0


def radiant_declination_visible(dec_deg, site_latitude_deg: float = PANSY_SITE_LATITUDE_DEG):
    """Return whether a celestial radiant has nonzero above-horizon time."""
    dec = np.asarray(dec_deg, dtype=np.float64)
    limit = horizon_declination_limit_deg(site_latitude_deg)
    if float(site_latitude_deg) < 0.0:
        return dec <= limit
    return dec >= limit


def radiant_altitude_deg(
    ra_deg,
    dec_deg,
    epoch_unix,
    site_latitude_deg: float = PANSY_SITE_LATITUDE_DEG,
    site_longitude_deg: float = PANSY_SITE_LONGITUDE_DEG,
):
    """Return instantaneous local AltAz altitude for a GCRS radiant."""
    from astropy import units as u
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord
    from astropy.time import Time

    t = Time(np.asarray(epoch_unix, dtype=np.float64), format="unix", scale="utc")
    loc = EarthLocation(lat=float(site_latitude_deg) * u.deg, lon=float(site_longitude_deg) * u.deg, height=0.0 * u.m)
    sky = SkyCoord(
        ra=np.asarray(ra_deg, dtype=np.float64) * u.deg,
        dec=np.asarray(dec_deg, dtype=np.float64) * u.deg,
        frame="gcrs",
        obstime=t,
    )
    alt = sky.transform_to(AltAz(location=loc, obstime=t)).alt.deg
    return np.asarray(alt, dtype=np.float64)


def radiant_above_local_horizon(ra_deg, dec_deg, epoch_unix, min_altitude_deg: float = 0.0):
    """Return whether a GCRS radiant is above the local horizon at event time."""
    alt = radiant_altitude_deg(ra_deg, dec_deg, epoch_unix)
    return np.asarray(alt, dtype=np.float64) > float(min_altitude_deg)


def radiant_visibility_boundary(
    solar_longitude_deg: float,
    site_latitude_deg: float = PANSY_SITE_LATITUDE_DEG,
    n_points: int = 721,
):
    """Return plot longitude and ecliptic latitude for the zero-hours boundary.

    The plotted curve is where the radiant just reaches the local horizon once
    per sidereal day. For PANSY this is declination +20.989 deg; lower
    declinations have nonzero above-horizon time.
    """
    plot_longitude = np.linspace(-180.0, 180.0, int(n_points))
    lambda_minus_sun = wrap180(PLOT_CENTER_LONGITUDE_DEG - plot_longitude)
    absolute_lambda = wrap360(float(solar_longitude_deg) + lambda_minus_sun)
    eps = np.deg2rad(MEAN_OBLIQUITY_DEG)
    decl_limit = np.deg2rad(horizon_declination_limit_deg(site_latitude_deg))
    a = np.cos(eps)
    b = np.sin(eps) * np.sin(np.deg2rad(absolute_lambda))
    phase = np.arctan2(b, a)
    rhs = np.sin(decl_limit) / np.hypot(a, b)
    rhs = np.clip(rhs, -1.0, 1.0)
    beta = np.rad2deg(np.arcsin(rhs) - phase)
    beta = ((beta + 90.0) % 360.0) - 90.0
    return plot_longitude, beta


def add_visibility_boundary_hammer(
    ax,
    solar_longitude_deg: float,
    site_latitude_deg: float = PANSY_SITE_LATITUDE_DEG,
    color: str = "white",
    label: str = "Visibility > 0 h threshold",
):
    x_deg, beta_deg = radiant_visibility_boundary(
        solar_longitude_deg,
        site_latitude_deg=site_latitude_deg,
    )
    ax.plot(
        np.deg2rad(x_deg),
        np.deg2rad(beta_deg),
        color=color,
        linewidth=1.35,
        linestyle="--",
        alpha=0.92,
        zorder=11,
        label=label,
    )
