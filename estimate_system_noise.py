#!/usr/bin/env python3
"""Estimate PANSY sky/system noise with antenna-pattern weighting."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u

import pansy_config as pc
import pansy_gain as pgain


BEAM_NAMES = ["Zenith", "North", "East", "South", "West"]


def parse_utc(value: str) -> float:
    """Parse a unix timestamp or ISO UTC string."""
    try:
        return float(value)
    except ValueError:
        text = value.strip().replace("Z", "+00:00")
        dt = datetime.fromisoformat(text)
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return dt.timestamp()


def local_sky_grid(n_az: int, n_el: int, min_elevation_deg: float) -> dict[str, np.ndarray]:
    """Return a midpoint alt/az grid, direction cosines, and solid-angle weights."""
    az_edges = np.linspace(0.0, 360.0, n_az + 1)
    el_edges = np.linspace(min_elevation_deg, 90.0, n_el + 1)
    az_deg = 0.5 * (az_edges[:-1] + az_edges[1:])
    el_deg = 0.5 * (el_edges[:-1] + el_edges[1:])
    az_grid, el_grid = np.meshgrid(az_deg, el_deg)

    az_rad = np.deg2rad(az_grid)
    el_rad = np.deg2rad(el_grid)
    east = np.cos(el_rad) * np.sin(az_rad)
    north = np.cos(el_rad) * np.cos(az_rad)
    up = np.sin(el_rad)
    uvw = np.column_stack([east.ravel(), north.ravel(), -up.ravel()])

    az_step = np.deg2rad(az_edges[1] - az_edges[0])
    el_step = np.deg2rad(el_edges[1] - el_edges[0])
    solid_angle = np.cos(el_rad).ravel() * az_step * el_step
    return {
        "az_deg": az_grid.ravel(),
        "el_deg": el_grid.ravel(),
        "uvw": uvw,
        "solid_angle_sr": solid_angle,
    }


def gsm_temperature(az_deg: np.ndarray, el_deg: np.ndarray, unix_time: float, freq_mhz: float) -> np.ndarray:
    """Evaluate the Global Sky Model on a local alt/az grid."""
    from pygdsm import GlobalSkyModel

    location = EarthLocation(lat=pc.lat * u.deg, lon=pc.lon * u.deg, height=100.0 * u.m)
    frame = AltAz(obstime=Time(unix_time, format="unix"), location=location)
    coords = SkyCoord(az=az_deg * u.deg, alt=el_deg * u.deg, frame=frame)
    gsm = GlobalSkyModel()
    gsm.generate(float(freq_mhz))
    try:
        return np.asarray(gsm.get_sky_temperature(coords), dtype=np.float64).reshape(-1)
    except Exception:
        return np.asarray([float(gsm.get_sky_temperature(coord)) for coord in coords], dtype=np.float64)


def gain_weights(
    uvw: np.ndarray,
    beam_id: int,
    model: str,
    rx_channel: int | str | None = None,
    include_element_pattern: bool = True,
) -> np.ndarray:
    """Return linear gain weights for one beam."""
    beam_vecs = pgain.tx_beam_unit_vectors()
    if model == "rx":
        return pgain.rx_power_gain(uvw, channel=rx_channel, steer=beam_vecs[beam_id], include_element_pattern=include_element_pattern)
    if model == "tx":
        return pgain.tx_power_gain(uvw, beam_id, beam_vecs=beam_vecs)
    if model == "two_way":
        return pgain.two_way_power_gain(uvw, beam_id, beam_vecs=beam_vecs, rx_channel=rx_channel)
    raise ValueError(f"unknown gain model {model!r}")


def weighted_temperature(
    sky_temp_k: np.ndarray,
    gain: np.ndarray,
    solid_angle_sr: np.ndarray,
) -> tuple[float, float]:
    """Return gain-weighted sky temperature and beam solid angle."""
    weight = np.asarray(gain, dtype=np.float64) * np.asarray(solid_angle_sr, dtype=np.float64)
    good = np.isfinite(sky_temp_k) & np.isfinite(weight) & (weight > 0.0)
    if not np.any(good):
        return np.nan, np.nan
    beam_solid_angle = float(np.sum(weight[good]))
    temp = float(np.sum(sky_temp_k[good] * weight[good]) / beam_solid_angle)
    return temp, beam_solid_angle


def estimate(
    unix_time: float,
    model: str,
    receiver_temp_k: float,
    n_az: int,
    n_el: int,
    min_elevation_deg: float,
    rx_channel: int | str | None,
    freq_mhz: float,
) -> list[dict[str, float | int | str]]:
    grid = local_sky_grid(n_az=n_az, n_el=n_el, min_elevation_deg=min_elevation_deg)
    sky_temp = gsm_temperature(grid["az_deg"], grid["el_deg"], unix_time=unix_time, freq_mhz=freq_mhz)
    rows = []
    for beam_id, beam_name in enumerate(BEAM_NAMES):
        gain = gain_weights(grid["uvw"], beam_id=beam_id, model=model, rx_channel=rx_channel)
        sky_k, beam_sr = weighted_temperature(sky_temp, gain, grid["solid_angle_sr"])
        rows.append(
            {
                "unix_time": float(unix_time),
                "utc": datetime.fromtimestamp(unix_time, tz=timezone.utc).isoformat().replace("+00:00", "Z"),
                "beam_id": beam_id,
                "beam_name": beam_name,
                "gain_model": model,
                "freq_mhz": float(freq_mhz),
                "sky_temp_k": sky_k,
                "receiver_temp_k": float(receiver_temp_k),
                "system_temp_k": sky_k + float(receiver_temp_k),
                "beam_solid_angle_sr": beam_sr,
            }
        )
    return rows


def write_csv(rows: list[dict[str, float | int | str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--time", default=datetime.now(timezone.utc).isoformat(), help="Unix seconds or ISO UTC time.")
    parser.add_argument("--gain-model", choices=["rx", "tx", "two_way"], default="rx")
    parser.add_argument("--receiver-temp-k", type=float, default=0.0, help="Additive receiver temperature.")
    parser.add_argument("--freq-mhz", type=float, default=pc.freq / 1e6)
    parser.add_argument("--n-az", type=int, default=360)
    parser.add_argument("--n-el", type=int, default=90)
    parser.add_argument("--min-elevation-deg", type=float, default=0.5)
    parser.add_argument("--rx-channel", default=None, help="Optional receive channel index or module name.")
    parser.add_argument("--output", type=Path, default=Path("figs/system_noise_estimate.csv"))
    args = parser.parse_args()

    rx_channel: int | str | None
    if args.rx_channel is None:
        rx_channel = None
    else:
        try:
            rx_channel = int(args.rx_channel)
        except ValueError:
            rx_channel = args.rx_channel

    rows = estimate(
        unix_time=parse_utc(args.time),
        model=args.gain_model,
        receiver_temp_k=args.receiver_temp_k,
        n_az=args.n_az,
        n_el=args.n_el,
        min_elevation_deg=args.min_elevation_deg,
        rx_channel=rx_channel,
        freq_mhz=args.freq_mhz,
    )
    write_csv(rows, args.output)
    for row in rows:
        print(
            f"{row['utc']} beam {row['beam_id']} {row['beam_name']:>6s} "
            f"{row['gain_model']} sky={row['sky_temp_k']:.1f} K "
            f"system={row['system_temp_k']:.1f} K"
        )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
