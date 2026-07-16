"""Render equal-area HEALPix maps on Matplotlib Hammer axes."""

from __future__ import annotations

import healpy as hp
import numpy as np
from matplotlib.collections import PolyCollection


def wrap180(deg):
    return (np.asarray(deg, dtype=np.float64) + 180.0) % 360.0 - 180.0


def centered_plot_longitude_deg(lambda_minus_sun_deg, center_longitude_deg: float = -90.0):
    signed = wrap180(lambda_minus_sun_deg)
    return -wrap180(signed - float(center_longitude_deg))


def healpix_centers(nside: int, center_longitude_deg: float = -90.0):
    pixel = np.arange(hp.nside2npix(int(nside)), dtype=np.int64)
    lon_deg, beta_deg = hp.pix2ang(int(nside), pixel, lonlat=True)
    plot_lon_deg = centered_plot_longitude_deg(lon_deg, center_longitude_deg=center_longitude_deg)
    return pixel, np.asarray(lon_deg), np.asarray(beta_deg), np.asarray(plot_lon_deg)


def _pixel_polygons(nside: int, pixels: np.ndarray, center_longitude_deg: float):
    boundaries = hp.boundaries(int(nside), np.asarray(pixels, dtype=np.int64), step=1)
    lon_deg = np.rad2deg(np.arctan2(boundaries[:, 1, :], boundaries[:, 0, :])) % 360.0
    beta_deg = np.rad2deg(np.arcsin(np.clip(boundaries[:, 2, :], -1.0, 1.0)))
    plot_lon_deg = centered_plot_longitude_deg(lon_deg, center_longitude_deg=center_longitude_deg)
    center_lon_deg, _ = hp.pix2ang(int(nside), np.asarray(pixels, dtype=np.int64), lonlat=True)
    center_plot_deg = centered_plot_longitude_deg(center_lon_deg, center_longitude_deg=center_longitude_deg)

    polygons = []
    source_indices = []
    for source_idx, (x_deg, y_deg, center_deg) in enumerate(zip(plot_lon_deg, beta_deg, center_plot_deg, strict=True)):
        x_deg = x_deg + 360.0 * np.round((center_deg - x_deg) / 360.0)
        polygon = np.column_stack((np.deg2rad(x_deg), np.deg2rad(y_deg)))
        polygons.append(polygon)
        source_indices.append(source_idx)
        if np.max(x_deg) > 180.0:
            polygons.append(np.column_stack((np.deg2rad(x_deg - 360.0), np.deg2rad(y_deg))))
            source_indices.append(source_idx)
        if np.min(x_deg) < -180.0:
            polygons.append(np.column_stack((np.deg2rad(x_deg + 360.0), np.deg2rad(y_deg))))
            source_indices.append(source_idx)
    return polygons, np.asarray(source_indices, dtype=np.int64)


def render_healpix_hammer(
    ax,
    values,
    nside: int,
    *,
    cmap="magma",
    norm=None,
    center_longitude_deg: float = -90.0,
):
    """Draw finite positive HEALPix cells as their true spherical boundaries."""
    values = np.asarray(values, dtype=np.float64)
    expected = hp.nside2npix(int(nside))
    if values.shape != (expected,):
        raise ValueError(f"expected {expected} HEALPix values for nside={nside}, got {values.shape}")
    keep = np.isfinite(values) & (values > 0.0)
    pixels = np.flatnonzero(keep)
    polygons, source_indices = _pixel_polygons(int(nside), pixels, float(center_longitude_deg))
    collection = PolyCollection(
        polygons,
        array=values[pixels][source_indices],
        cmap=cmap,
        norm=norm,
        edgecolors="none",
        linewidths=0.0,
        antialiased=False,
        rasterized=True,
        zorder=1,
    )
    collection.set_clip_path(ax.patch)
    ax.add_collection(collection)
    ax.set_facecolor("black")
    return collection
