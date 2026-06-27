"""Optional IAU meteor shower overlays for sun-centered radiant plots."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from plot_fitted_radiant_distribution import centered_plot_longitude_deg, wrap180


def circular_mean_deg(values_deg):
    values = np.asarray(values_deg, dtype=np.float64)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return math.nan
    rad = np.deg2rad(values)
    return float(np.rad2deg(np.arctan2(np.mean(np.sin(rad)), np.mean(np.cos(rad)))) % 360.0)


def resolve_shower_query(rows, shower_solar_longitude=None, shower_date=None):
    if shower_solar_longitude is not None:
        return float(shower_solar_longitude)
    if shower_date is not None:
        return shower_date
    arr = rows
    if isinstance(rows, list):
        values = [row.get("sun_lambda_ecliptic_deg", np.nan) for row in rows]
    else:
        values = arr["sun_lambda_ecliptic_deg"] if "sun_lambda_ecliptic_deg" in arr.dtype.names else []
    return circular_mean_deg(values)


def active_showers(rows, shower_catalog=None, shower_solar_longitude=None, shower_date=None, peak_tolerance_deg=5.0):
    try:
        import iau_meteor_showers as ims
    except ImportError as exc:
        raise RuntimeError(
            "IAU shower overlays require the iau_meteor_showers package. "
            "Install it or add its src directory to PYTHONPATH."
        ) from exc
    query = resolve_shower_query(rows, shower_solar_longitude=shower_solar_longitude, shower_date=shower_date)
    query_lon = ims.solar_longitude(query)
    if not np.isfinite(float(query_lon)):
        return [], query
    showers = ims.get_showers(
        query,
        path=Path(shower_catalog) if shower_catalog is not None else None,
        peak_tolerance_deg=peak_tolerance_deg,
    )
    return showers, query_lon


def shower_plot_points(showers):
    if not showers:
        return np.asarray([]), np.asarray([]), []
    lon = np.asarray([shower.radiant_solar_ecliptic_lon_deg for shower in showers], dtype=np.float64)
    lat = np.asarray([shower.radiant_ecliptic_lat_deg for shower in showers], dtype=np.float64)
    labels = [shower.code or shower.name[:5] for shower in showers]
    x = centered_plot_longitude_deg(wrap180(lon))
    return x, lat, labels


def add_shower_overlay_hammer(ax, showers, max_labels=18, label_color="black"):
    x_deg, lat_deg, labels = shower_plot_points(showers)
    if len(x_deg) == 0:
        return
    ax.scatter(
        np.deg2rad(x_deg),
        np.deg2rad(lat_deg),
        marker="D",
        s=34,
        facecolor="#fff2a8",
        edgecolor="black",
        linewidth=0.55,
        zorder=12,
        label="IAU showers",
    )
    for x, y, label in zip(x_deg[:max_labels], lat_deg[:max_labels], labels[:max_labels], strict=False):
        ax.annotate(
            label,
            (np.deg2rad(x), np.deg2rad(y)),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=7.5,
            color=label_color,
            zorder=13,
        )
