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


def gmst_deg_from_unix(epoch_unix):
    epoch = np.asarray(epoch_unix, dtype=np.float64)
    jd = epoch / 86400.0 + 2440587.5
    d = jd - 2451545.0
    gmst = 280.46061837 + 360.98564736629 * d
    return wrap360(gmst)


def ecliptic_lonlat_to_radec_deg(lon_deg, lat_deg):
    lon = np.deg2rad(np.asarray(lon_deg, dtype=np.float64))
    lat = np.deg2rad(np.asarray(lat_deg, dtype=np.float64))
    eps = np.deg2rad(MEAN_OBLIQUITY_DEG)
    sin_dec = np.sin(lat) * np.cos(eps) + np.cos(lat) * np.sin(eps) * np.sin(lon)
    dec = np.arcsin(np.clip(sin_dec, -1.0, 1.0))
    ra = np.arctan2(np.sin(lon) * np.cos(eps) - np.tan(lat) * np.sin(eps), np.cos(lon))
    return wrap360(np.rad2deg(ra)), np.rad2deg(dec)


def altitude_from_radec_deg(
    ra_deg,
    dec_deg,
    epoch_unix,
    site_latitude_deg: float = PANSY_SITE_LATITUDE_DEG,
    site_longitude_deg: float = PANSY_SITE_LONGITUDE_DEG,
):
    lat = np.deg2rad(float(site_latitude_deg))
    lst = np.deg2rad(wrap360(gmst_deg_from_unix(epoch_unix) + float(site_longitude_deg)))
    ra = np.deg2rad(np.asarray(ra_deg, dtype=np.float64))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=np.float64))
    hour_angle = lst - ra
    sin_alt = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(hour_angle)
    return np.rad2deg(np.arcsin(np.clip(sin_alt, -1.0, 1.0)))


def visibility_samples_from_radiants(rows, step_seconds: float = 3600.0, max_samples: int = 1500):
    """Return epoch and solar-longitude samples spanning the plotted radiants."""
    if rows is None or len(rows) == 0:
        return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64)
    epoch = np.asarray(rows["epoch_unix"], dtype=np.float64) if "epoch_unix" in rows.dtype.names else np.asarray([], dtype=np.float64)
    sun = (
        np.asarray(rows["sun_lambda_ecliptic_deg"], dtype=np.float64)
        if "sun_lambda_ecliptic_deg" in rows.dtype.names
        else np.asarray([], dtype=np.float64)
    )
    good = np.isfinite(epoch) & np.isfinite(sun)
    epoch = epoch[good]
    sun = sun[good]
    if len(epoch) == 0:
        return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64)
    order = np.argsort(epoch)
    epoch = epoch[order]
    sun = sun[order]
    start = float(epoch[0])
    stop = float(epoch[-1])
    if stop <= start:
        return epoch[:1], sun[:1]
    sample_count = int(np.floor((stop - start) / float(step_seconds))) + 1
    sample_count = max(2, min(int(max_samples), sample_count))
    sample_epoch = np.linspace(start, stop, sample_count)
    unwrapped_sun = np.rad2deg(np.unwrap(np.deg2rad(sun)))
    sample_sun = wrap360(np.interp(sample_epoch, epoch, unwrapped_sun))
    return sample_epoch, sample_sun


def composite_radiant_visibility_grid(
    epoch_unix,
    solar_longitude_deg,
    plot_longitude_deg=None,
    beta_deg=None,
    min_altitude_deg: float = 0.0,
):
    """Return a mask of sun-centered radiant positions visible at least once."""
    epoch = np.asarray(epoch_unix, dtype=np.float64)
    sun = np.asarray(solar_longitude_deg, dtype=np.float64)
    good = np.isfinite(epoch) & np.isfinite(sun)
    epoch = epoch[good]
    sun = sun[good]
    if plot_longitude_deg is None:
        plot_longitude_deg = np.linspace(-180.0, 180.0, 361)
    if beta_deg is None:
        beta_deg = np.linspace(-90.0, 90.0, 181)
    plot_longitude_deg = np.asarray(plot_longitude_deg, dtype=np.float64)
    beta_deg = np.asarray(beta_deg, dtype=np.float64)
    if len(epoch) == 0:
        return plot_longitude_deg, beta_deg, np.zeros((len(beta_deg), len(plot_longitude_deg)), dtype=bool)
    lon_grid, beta_grid = np.meshgrid(plot_longitude_deg, beta_deg)
    lambda_minus_sun = wrap180(PLOT_CENTER_LONGITUDE_DEG - lon_grid)
    possible = composite_radiant_visibility_points(
        lambda_minus_sun,
        beta_grid,
        epoch,
        sun,
        min_altitude_deg=min_altitude_deg,
    )
    return plot_longitude_deg, beta_deg, possible


def radiant_exposure_hours_grid(
    epoch_unix,
    solar_longitude_deg,
    observing_hours,
    plot_longitude_deg=None,
    beta_deg=None,
    min_altitude_deg: float = 0.0,
    chunk_size: int = 512,
):
    """Accumulate observing hours for each Sun-centered radiant direction.

    Each input sample represents an observing interval centered at
    ``epoch_unix`` with duration ``observing_hours``. The local zenith is
    rotated into Sun-centered ecliptic coordinates, allowing the full grid to
    be evaluated with matrix products rather than repeated coordinate-object
    construction.
    """
    epoch = np.asarray(epoch_unix, dtype=np.float64)
    sun = np.asarray(solar_longitude_deg, dtype=np.float64)
    hours = np.asarray(observing_hours, dtype=np.float64)
    epoch, sun, hours = np.broadcast_arrays(epoch, sun, hours)
    good = np.isfinite(epoch) & np.isfinite(sun) & np.isfinite(hours) & (hours > 0.0)
    epoch = epoch[good]
    sun = sun[good]
    hours = hours[good]

    if plot_longitude_deg is None:
        plot_longitude_deg = np.linspace(-180.0, 180.0, 361)
    if beta_deg is None:
        beta_deg = np.linspace(-90.0, 90.0, 181)
    plot_longitude_deg = np.asarray(plot_longitude_deg, dtype=np.float64)
    beta_deg = np.asarray(beta_deg, dtype=np.float64)
    exposure = np.zeros((len(beta_deg), len(plot_longitude_deg)), dtype=np.float64)
    if len(epoch) == 0:
        return plot_longitude_deg, beta_deg, exposure

    plot_lon_grid, beta_grid = np.meshgrid(plot_longitude_deg, beta_deg)
    relative_lon = np.deg2rad(wrap180(PLOT_CENTER_LONGITUDE_DEG - plot_lon_grid.ravel()))
    beta_rad = np.deg2rad(beta_grid.ravel())
    cos_beta = np.cos(beta_rad)
    radiant_vectors = np.column_stack(
        (
            cos_beta * np.cos(relative_lon),
            cos_beta * np.sin(relative_lon),
            np.sin(beta_rad),
        )
    )

    site_lat = np.deg2rad(PANSY_SITE_LATITUDE_DEG)
    lst = np.deg2rad(wrap360(gmst_deg_from_unix(epoch) + PANSY_SITE_LONGITUDE_DEG))
    zenith_equatorial = np.column_stack(
        (
            np.cos(site_lat) * np.cos(lst),
            np.cos(site_lat) * np.sin(lst),
            np.full(len(epoch), np.sin(site_lat)),
        )
    )
    eps = np.deg2rad(MEAN_OBLIQUITY_DEG)
    zenith_ecliptic = np.column_stack(
        (
            zenith_equatorial[:, 0],
            np.cos(eps) * zenith_equatorial[:, 1] + np.sin(eps) * zenith_equatorial[:, 2],
            -np.sin(eps) * zenith_equatorial[:, 1] + np.cos(eps) * zenith_equatorial[:, 2],
        )
    )
    sun_rad = np.deg2rad(sun)
    zenith_sun_centered = np.column_stack(
        (
            np.cos(sun_rad) * zenith_ecliptic[:, 0] + np.sin(sun_rad) * zenith_ecliptic[:, 1],
            -np.sin(sun_rad) * zenith_ecliptic[:, 0] + np.cos(sun_rad) * zenith_ecliptic[:, 1],
            zenith_ecliptic[:, 2],
        )
    )

    threshold = np.sin(np.deg2rad(float(min_altitude_deg)))
    flat_exposure = exposure.ravel()
    chunk_size = max(1, int(chunk_size))
    for start in range(0, len(epoch), chunk_size):
        stop = min(start + chunk_size, len(epoch))
        sin_altitude = zenith_sun_centered[start:stop] @ radiant_vectors.T
        flat_exposure += (sin_altitude > threshold).T @ hours[start:stop]
    return plot_longitude_deg, beta_deg, exposure


def composite_radiant_visibility_points(
    lambda_minus_sun_deg,
    beta_deg,
    epoch_unix,
    solar_longitude_deg,
    min_altitude_deg: float = 0.0,
):
    """Return whether each sun-centered ecliptic radiant was visible at least once."""
    epoch = np.asarray(epoch_unix, dtype=np.float64)
    sun = np.asarray(solar_longitude_deg, dtype=np.float64)
    good = np.isfinite(epoch) & np.isfinite(sun)
    epoch = epoch[good]
    sun = sun[good]
    lambda_minus_sun = np.asarray(lambda_minus_sun_deg, dtype=np.float64)
    beta = np.asarray(beta_deg, dtype=np.float64)
    lambda_minus_sun, beta = np.broadcast_arrays(lambda_minus_sun, beta)
    possible = np.zeros(lambda_minus_sun.shape, dtype=bool)
    for ep, sun_lon in zip(epoch, sun, strict=False):
        absolute_lambda = wrap360(float(sun_lon) + lambda_minus_sun)
        ra, dec = ecliptic_lonlat_to_radec_deg(absolute_lambda, beta)
        possible |= altitude_from_radec_deg(ra, dec, float(ep)) > float(min_altitude_deg)
    return possible


def add_composite_visibility_boundary_hammer(
    ax,
    epoch_unix,
    solar_longitude_deg,
    color: str = "white",
    label: str | None = "Composite visibility > 0 h",
    min_altitude_deg: float = 0.0,
):
    x_deg, beta_deg, possible = composite_radiant_visibility_grid(
        epoch_unix,
        solar_longitude_deg,
        min_altitude_deg=min_altitude_deg,
    )
    if not np.any(possible) or np.all(possible):
        return
    ax.contour(
        np.deg2rad(x_deg),
        np.deg2rad(beta_deg),
        possible.astype(float),
        levels=[0.5],
        colors=color,
        linewidths=1.35,
        linestyles="--",
        alpha=0.92,
        zorder=11,
    )


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
