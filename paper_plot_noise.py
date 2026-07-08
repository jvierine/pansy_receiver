#!/usr/bin/env python3
"""Plot the PANSY paper system-noise figure with analytic sky-noise modeling."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import multiprocessing as mp
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u

import estimate_system_noise as esn
import pansy_config as pc


BEAM_NAMES = ["Zenith", "North", "East", "South", "West"]
DEFAULT_METADATA_PATHS = (
    Path("/mnt/data/juha/pansy/metadata/cut"),
    Path("/media/analysis/metadata/cut"),
)


def utc_day_bounds(day: str) -> tuple[int, int, datetime]:
    date = datetime.strptime(day, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    start_us = int(date.timestamp() * 1_000_000)
    return start_us, start_us + int(86_400 * 1_000_000), date


def default_metadata_path() -> Path:
    for path in DEFAULT_METADATA_PATHS:
        if path.exists():
            return path
    return DEFAULT_METADATA_PATHS[0]


def read_xc2_day(metadata_path: Path, day: str) -> dict[str, np.ndarray | int | str]:
    import digital_rf as drf

    start_us, end_us, date = utc_day_bounds(day)
    dm = drf.DigitalMetadataReader(str(metadata_path))
    records = dm.read(start_us, end_us)
    if not records:
        raise RuntimeError(f"no xc2 metadata records found in {metadata_path} for {day}")

    times_s = []
    noise_floor = []
    range_noise = []
    r0 = r1 = rdec = None
    for key in sorted(records):
        record = records[key]
        xc_arr = np.asarray(record["xc_arr"])
        if xc_arr.ndim != 4 or xc_arr.shape[1] < len(BEAM_NAMES):
            continue
        if r0 is None:
            r0 = int(record["r0"])
            r1 = int(record["r1"])
            rdec = int(record["rdec"])
        times_s.append(float(key) / 1e6)
        noise_floor.append([float(np.nanmedian(np.real(xc_arr[0, beam_i, :, :]))) for beam_i in range(len(BEAM_NAMES))])
        range_noise.append(np.nanmax(np.real(xc_arr[0, 0, :, :]), axis=0))

    if len(times_s) == 0 or r0 is None or r1 is None or rdec is None:
        raise RuntimeError(f"xc2 records in {metadata_path} for {day} did not contain usable xc_arr data")

    return {
        "day": day,
        "date": date,
        "times_s": np.asarray(times_s, dtype=np.float64),
        "noise_floor": np.asarray(noise_floor, dtype=np.float64),
        "range_noise": np.asarray(range_noise, dtype=np.float64),
        "r0": r0,
        "r1": r1,
        "rdec": rdec,
    }


def _read_cut_minute(task: tuple[str, int, int, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read one minute of cut metadata and return pulse noise samples."""
    import digital_rf as drf

    metadata_path, start_us, end_us, guard_samples = task
    dm = drf.DigitalMetadataReader(metadata_path)
    try:
        records = dm.read(start_us, end_us)
    except Exception:
        return (
            np.asarray([], dtype=np.float64),
            np.asarray([], dtype=np.int16),
            np.asarray([], dtype=np.float64),
        )

    times = []
    beams = []
    powers = []
    for _key, record in records.items():
        try:
            zrx_re = np.asarray(record["zrx_echoes_re"], dtype=np.float32)
            zrx_im = np.asarray(record["zrx_echoes_im"], dtype=np.float32)
            beam_id = np.asarray(record["beam_id"], dtype=np.int16).reshape(-1)
            tx_idx = np.asarray(record["tx_idx"], dtype=np.float64).reshape(-1)
            pad = int(np.asarray(record["pad"]).reshape(()))
            txlen = int(np.asarray(record["txlen"]).reshape(()))
        except Exception:
            continue
        if zrx_re.ndim != 3 or zrx_re.shape != zrx_im.shape:
            continue
        n_pulse = min(zrx_re.shape[0], len(beam_id), len(tx_idx))
        if n_pulse == 0:
            continue
        n_sample = zrx_re.shape[2]
        pre_end = max(0, min(pad - guard_samples, n_sample))
        post_start = max(0, min(pad + txlen + guard_samples, n_sample))
        pre = np.arange(0, pre_end, dtype=np.int64)
        post = np.arange(post_start, n_sample, dtype=np.int64)
        noise_idx = np.concatenate([pre, post])
        if len(noise_idx) == 0:
            continue
        power = zrx_re[:n_pulse, :, :][:, :, noise_idx] ** 2 + zrx_im[:n_pulse, :, :][:, :, noise_idx] ** 2
        pulse_power = np.nanmean(power, axis=(1, 2))
        good = np.isfinite(pulse_power) & np.isfinite(tx_idx[:n_pulse]) & (beam_id[:n_pulse] >= 0) & (beam_id[:n_pulse] < len(BEAM_NAMES))
        if np.any(good):
            times.append(tx_idx[:n_pulse][good] / 1e6)
            beams.append(beam_id[:n_pulse][good])
            powers.append(pulse_power[good])

    if not times:
        return (
            np.asarray([], dtype=np.float64),
            np.asarray([], dtype=np.int16),
            np.asarray([], dtype=np.float64),
        )
    return (
        np.concatenate(times).astype(np.float64),
        np.concatenate(beams).astype(np.int16),
        np.concatenate(powers).astype(np.float64),
    )


def read_cut_day(
    metadata_path: Path,
    day: str,
    workers: int,
    chunk_seconds: int = 60,
    guard_samples: int = 25,
) -> dict[str, np.ndarray | str | int]:
    start_us, end_us, _date = utc_day_bounds(day)
    step_us = int(chunk_seconds * 1_000_000)
    tasks = [
        (str(metadata_path), start, min(start + step_us, end_us), int(guard_samples))
        for start in range(start_us, end_us, step_us)
    ]
    if workers <= 1:
        parts = [_read_cut_minute(task) for task in tasks]
    else:
        with mp.Pool(processes=workers) as pool:
            parts = list(pool.imap_unordered(_read_cut_minute, tasks, chunksize=1))
    times = [part[0] for part in parts if len(part[0])]
    beams = [part[1] for part in parts if len(part[1])]
    powers = [part[2] for part in parts if len(part[2])]
    if not times:
        raise RuntimeError(f"no cut metadata noise samples found in {metadata_path} for {day}")
    return {
        "day": day,
        "times_s": np.concatenate(times),
        "beam_id": np.concatenate(beams),
        "noise_power": np.concatenate(powers),
        "guard_samples": int(guard_samples),
    }


def sky_temperature_matrix(
    times_s: np.ndarray,
    gain_model: str,
    n_az: int,
    n_el: int,
    min_elevation_deg: float,
    rx_channel: int | str | None,
    freq_mhz: float,
) -> np.ndarray:
    """Return beam-weighted T_sky with shape (time, beam)."""
    from pygdsm import GlobalSkyModel

    grid = esn.local_sky_grid(n_az=n_az, n_el=n_el, min_elevation_deg=min_elevation_deg)
    gains = [
        esn.gain_weights(grid["uvw"], beam_id=beam_i, model=gain_model, rx_channel=rx_channel)
        for beam_i in range(len(BEAM_NAMES))
    ]

    location = EarthLocation(lat=pc.lat * u.deg, lon=pc.lon * u.deg, height=100.0 * u.m)
    gsm = GlobalSkyModel()
    gsm.generate(float(freq_mhz))
    out = np.full((len(times_s), len(BEAM_NAMES)), np.nan, dtype=np.float64)
    for ti, unix_time in enumerate(times_s):
        frame = AltAz(obstime=Time(float(unix_time), format="unix"), location=location)
        coords = SkyCoord(az=grid["az_deg"] * u.deg, alt=grid["el_deg"] * u.deg, frame=frame)
        try:
            sky_temp = np.asarray(gsm.get_sky_temperature(coords), dtype=np.float64).reshape(-1)
        except Exception:
            sky_temp = np.asarray([float(gsm.get_sky_temperature(coord)) for coord in coords], dtype=np.float64)
        for beam_i, gain in enumerate(gains):
            out[ti, beam_i], _beam_sr = esn.weighted_temperature(sky_temp, gain, grid["solid_angle_sr"])
    return out


def smooth_trace(values: np.ndarray, kernel: int) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if kernel <= 1 or len(values) < 3:
        return values
    kernel = min(int(kernel), len(values) if len(values) % 2 == 1 else len(values) - 1)
    if kernel < 3:
        return values
    if kernel % 2 == 0:
        kernel -= 1
    pad = kernel // 2
    padded = np.pad(values, pad, mode="edge")
    return np.asarray([np.nanmedian(padded[i : i + kernel]) for i in range(len(values))], dtype=np.float64)


def fit_receiver_temperature(measured_power: np.ndarray, sky_temp_k: np.ndarray) -> dict[str, float | np.ndarray]:
    good = np.isfinite(measured_power) & np.isfinite(sky_temp_k)
    if np.count_nonzero(good) < 3:
        raise RuntimeError("not enough finite samples to fit receiver temperature")
    slope, intercept = np.linalg.lstsq(
        np.column_stack([sky_temp_k[good], np.ones(np.count_nonzero(good))]),
        measured_power[good],
        rcond=None,
    )[0]
    if slope <= 0.0:
        raise RuntimeError(f"non-positive fitted power-to-temperature slope: {slope}")
    return {
        "power_per_kelvin": float(slope),
        "receiver_temp_k": float(intercept / slope),
        "measured_tsys_k": measured_power / slope,
    }


def plot_noise_figure(
    day_data: dict[str, np.ndarray | int | str],
    sky_temp_k: np.ndarray,
    output: Path,
    gain_model: str,
    fit_beam: int,
    smooth_kernel: int,
) -> dict[str, float]:
    times_s = np.asarray(day_data["times_s"], dtype=np.float64)
    measured_power = np.asarray(day_data["noise_floor"], dtype=np.float64)
    range_noise = np.asarray(day_data["range_noise"], dtype=np.float64)
    fit_power = smooth_trace(measured_power[:, fit_beam], smooth_kernel)
    fit = fit_receiver_temperature(fit_power, sky_temp_k[:, fit_beam])
    receiver_temp_k = float(fit["receiver_temp_k"])
    scale = 1.0 / float(fit["power_per_kelvin"])
    measured_tsys = measured_power * scale
    measured_fit_tsys = np.asarray(fit["measured_tsys_k"], dtype=np.float64)
    model_tsys = sky_temp_k + receiver_temp_k

    tms = np.asarray(times_s, dtype="datetime64[s]")
    r0 = int(day_data["r0"])
    r1 = int(day_data["r1"])
    rdec = int(day_data["rdec"])
    rvec = np.arange(1600) * 0.15
    range_axis = rvec[r0:r1:rdec]

    output.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(6.9, 5.3))
    gs = gridspec.GridSpec(nrows=4, ncols=1, height_ratios=[1.0, 1.0, 0.20, 0.45], hspace=0.05)
    ax_range = fig.add_subplot(gs[0])
    ax_temp = fig.add_subplot(gs[1], sharex=ax_range)
    ax_pad = fig.add_subplot(gs[2])
    ax_cb = fig.add_subplot(gs[3])
    ax_pad.axis("off")
    ax_cb.axis("off")

    image = np.log10(np.maximum(scale * range_noise.T, 1.0))
    mesh = ax_range.pcolormesh(tms, range_axis, image, vmin=3.0, vmax=5.0, shading="auto")
    ax_range.set_ylabel("Range (km)")
    ax_range.tick_params(labelbottom=False)
    ax_range.set_title(str(day_data["day"]))
    axr = ax_range.twinx()
    axr.set_ylabel(BEAM_NAMES[fit_beam])
    axr.set_yticks([])

    colors = ["black", "tab:blue", "tab:orange", "tab:green", "tab:red"]
    ax_temp.plot(tms, measured_fit_tsys, color="0.55", lw=1.7, label=f"{BEAM_NAMES[fit_beam]} measured")
    for beam_i, beam_name in enumerate(BEAM_NAMES):
        lw = 2.1 if beam_i == fit_beam else 1.0
        alpha = 1.0 if beam_i == fit_beam else 0.7
        ax_temp.plot(tms, model_tsys[:, beam_i], color=colors[beam_i], lw=lw, alpha=alpha, label=f"{beam_name} model")
    ax_temp.set_ylabel(r"$T_{\mathrm{sys}}$ (K)")
    ax_temp.set_xlabel("Time (UTC)")
    ax_temp.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ylim_values = np.concatenate([measured_fit_tsys.reshape(-1), model_tsys.reshape(-1)])
    ax_temp.set_ylim(bottom=max(0.0, 0.85 * np.nanmin(ylim_values)), top=1.12 * np.nanmax(ylim_values))
    ax_temp.text(
        0.02,
        0.93,
        rf"$T_{{\mathrm{{rec}}}}={receiver_temp_k:.0f}$ K, {gain_model} beam model",
        transform=ax_temp.transAxes,
        ha="left",
        va="top",
    )
    ax_temp.legend(loc="upper right", fontsize=7, ncol=2, frameon=False)

    cb = fig.colorbar(mesh, ax=ax_cb, orientation="horizontal")
    cb.set_label(r"Noise temperature ($\log_{10}$ K)")
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {
        "receiver_temp_k": receiver_temp_k,
        "power_per_kelvin": float(fit["power_per_kelvin"]),
        "fit_beam_median_tsys_k": float(np.nanmedian(measured_fit_tsys)),
        "fit_beam_median_tsky_k": float(np.nanmedian(sky_temp_k[:, fit_beam])),
    }


def parse_rx_channel(value: str | None) -> int | str | None:
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--day", default="2025-05-10", help="UTC day to plot, YYYY-MM-DD.")
    parser.add_argument("--metadata", type=Path, default=None, help="Digital RF xc2 metadata path.")
    parser.add_argument("--output", type=Path, default=Path("figs/noise_temp.png"))
    parser.add_argument("--gain-model", choices=["rx", "tx", "two_way"], default="rx")
    parser.add_argument("--fit-beam", type=int, default=0)
    parser.add_argument("--freq-mhz", type=float, default=pc.freq / 1e6)
    parser.add_argument("--n-az", type=int, default=180)
    parser.add_argument("--n-el", type=int, default=60)
    parser.add_argument("--min-elevation-deg", type=float, default=0.5)
    parser.add_argument("--rx-channel", default=None)
    parser.add_argument("--smooth-kernel", type=int, default=5)
    args = parser.parse_args()

    metadata_path = args.metadata or default_metadata_path()
    day_data = read_xc2_day(metadata_path, args.day)
    sky_temp = sky_temperature_matrix(
        np.asarray(day_data["times_s"], dtype=np.float64),
        gain_model=args.gain_model,
        n_az=args.n_az,
        n_el=args.n_el,
        min_elevation_deg=args.min_elevation_deg,
        rx_channel=parse_rx_channel(args.rx_channel),
        freq_mhz=args.freq_mhz,
    )
    stats = plot_noise_figure(
        day_data,
        sky_temp,
        output=args.output,
        gain_model=args.gain_model,
        fit_beam=args.fit_beam,
        smooth_kernel=args.smooth_kernel,
    )
    print(f"metadata {metadata_path}")
    print(f"output {args.output}")
    for key, value in stats.items():
        print(f"{key} {value:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
