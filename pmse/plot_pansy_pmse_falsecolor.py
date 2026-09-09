#!/usr/bin/env python3
"""Create false-color PANSY PMSE RTDI plots from xc2 autocorrelation spectra."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
import numpy as np


SAMPLE_RATE_HZ = 1_000_000.0
RANGE_GATE_KM = 0.15
AUTOCORRELATION_COUNT = 7
RADAR_FREQUENCY_HZ = 47_000_000.0
# f_D = 2 f_radar v_radial / c; widths scale by the same positive factor.
MPS_PER_HZ = 299_792_458.0 / (2.0 * RADAR_FREQUENCY_HZ)


def utc_timestamp(day: str) -> float:
    return dt.datetime.fromisoformat(day).replace(tzinfo=dt.timezone.utc).timestamp()


def read_blocks(root: Path, day: str):
    blocks = []
    for hour in range(24):
        directory = root / f"{day}T{hour:02d}-00-00"
        for path in sorted(directory.glob("xc@*.h5")):
            with h5py.File(path, "r") as handle:
                for key in handle:
                    group = handle[key]
                    xc = group["xc_arr"]
                    power = np.sum(
                        np.abs(xc[:AUTOCORRELATION_COUNT, :, :, :]), axis=0
                    )
                    n_fft = int(group["n_fft"][()])
                    ipp = int(group["ipp"][()])
                    f0 = int(group["f0"][()])
                    f1 = int(group["f1"][()])
                    frequency = np.fft.fftshift(
                        np.fft.fftfreq(n_fft, d=5 * ipp / SAMPLE_RATE_HZ)
                    )[f0:f1]
                    outer = np.abs(frequency) >= 0.75 * np.max(np.abs(frequency))
                    noise = np.nanmedian(power[:, outer, :], axis=(1, 2))
                    ratio = power / noise[:, None, None]
                    snr_db = 10.0 * np.log10(
                        np.maximum(np.nanmax(ratio, axis=1), 1e-12)
                    )
                    weights = np.maximum(ratio - 1.0, 0.0)
                    weight_sum = np.sum(weights, axis=1)
                    mean = np.sum(weights * frequency[None, :, None], axis=1)
                    mean = np.divide(
                        mean, weight_sum, out=np.zeros_like(mean), where=weight_sum > 0
                    )
                    variance = np.sum(
                        weights * (frequency[None, :, None] - mean[:, None, :]) ** 2,
                        axis=1,
                    )
                    variance = np.divide(
                        variance,
                        weight_sum,
                        out=np.zeros_like(variance),
                        where=weight_sum > 0,
                    )
                    width = np.sqrt(np.maximum(variance, 0.0))
                    r0 = int(group["r0"][()])
                    r1 = int(group["r1"][()])
                    rdec = int(group["rdec"][()])
                    ranges = np.arange(r0, r1, rdec) * RANGE_GATE_KM
                    blocks.append(
                        (
                            int(group["i0"][()]) / SAMPLE_RATE_HZ,
                            int(group["i1"][()]) / SAMPLE_RATE_HZ,
                            ranges,
                            snr_db,
                            mean * MPS_PER_HZ,
                            width * MPS_PER_HZ,
                        )
                    )
    return sorted(blocks, key=lambda item: item[0])


def rgb_image(snr_db, mean, width, valid, snr_max, doppler_max, width_max):
    hue = np.clip((doppler_max - mean) / (2.0 * doppler_max), 0.0, 1.0) * (2.0 / 3.0)
    saturation = np.clip(width / width_max, 0.0, 1.0)
    value = np.clip(snr_db / snr_max, 0.0, 1.0)
    hsv = np.nan_to_num(
        np.stack((hue, saturation, value), axis=-1), nan=0.0, posinf=1.0, neginf=0.0
    )
    rgb = hsv_to_rgb(hsv)
    rgb[~valid, :] = 1.0
    return rgb


def add_mapping_bar(ax, kind, maximum):
    x = np.linspace(0.0, 1.0, 512)
    hsv = np.zeros((1, len(x), 3))
    if kind == "doppler":
        hsv[0, :, 0] = (1.0 - x) * (2.0 / 3.0)
        hsv[0, :, 1:] = 1.0
        extent = [-maximum, maximum, 0, 1]
        label = "Hue: radial velocity centroid (m/s)"
    elif kind == "snr":
        hsv[0, :, 1] = 0.0
        hsv[0, :, 2] = x
        extent = [0, maximum, 0, 1]
        label = "Intensity/value: peak SNR (dB)"
    else:
        hsv[0, :, 0] = 0.78
        hsv[0, :, 1] = x
        hsv[0, :, 2] = 1.0
        extent = [0, maximum, 0, 1]
        label = "Saturation: velocity RMS width (m/s)"
    ax.imshow(hsv_to_rgb(hsv), aspect="auto", extent=extent, origin="lower")
    ax.set_yticks([])
    ax.set_xlabel(label, fontsize=11)
    ax.tick_params(axis="x", labelsize=11, pad=1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--day", default="2025-01-30")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--minute-resolution", type=float, default=1.0)
    parser.add_argument("--snr-max", type=float, default=35.0)
    parser.add_argument("--doppler-max", type=float, default=3.0 * MPS_PER_HZ,
                        help="Symmetric radial velocity hue limit in m/s")
    parser.add_argument("--width-max", type=float, default=2.0 * MPS_PER_HZ,
                        help="RMS velocity width saturation limit in m/s")
    parser.add_argument("--hide-provenance", action="store_true")
    args = parser.parse_args()

    blocks = read_blocks(args.root, args.day)
    if not blocks:
        raise SystemExit(f"no xc2 blocks found for {args.day} under {args.root}")
    ranges = blocks[0][2]
    beam_count = blocks[0][3].shape[0]
    day_start = utc_timestamp(args.day)
    day_end = day_start + 86400.0
    cadence_s = args.minute_resolution * 60.0
    time_edges = np.arange(day_start, day_end + cadence_s, cadence_s)
    time_numbers = mdates.date2num(
        [dt.datetime.fromtimestamp(value, dt.timezone.utc) for value in time_edges]
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    for beam in range(beam_count):
        shape = (len(ranges), len(time_edges) - 1)
        snr = np.full(shape, np.nan)
        mean = np.full(shape, np.nan)
        width = np.full(shape, np.nan)
        for start, end, _ranges, block_snr, block_mean, block_width in blocks:
            lo = max(0, int(np.floor((start - day_start) / cadence_s)))
            hi = min(shape[1], int(np.ceil((end - day_start) / cadence_s)))
            for index in range(lo, hi):
                replace = np.isnan(snr[:, index]) | (block_snr[beam] > snr[:, index])
                snr[replace, index] = block_snr[beam, replace]
                mean[replace, index] = block_mean[beam, replace]
                width[replace, index] = block_width[beam, replace]
        available = np.flatnonzero(np.any(np.isfinite(snr), axis=0))
        if len(available) == 0:
            raise RuntimeError(f"beam {beam} has no valid columns")
        missing = np.flatnonzero(~np.any(np.isfinite(snr), axis=0))
        for index in missing:
            nearest = available[np.argmin(np.abs(available - index))]
            snr[:, index] = snr[:, nearest]
            mean[:, index] = mean[:, nearest]
            width[:, index] = width[:, nearest]
        valid = np.isfinite(snr)
        rgb = rgb_image(
            snr, mean, width, valid,
            args.snr_max, args.doppler_max, args.width_max,
        )

        plt.rcParams.update({'font.size': 11, 'axes.titlesize': 11,
                             'axes.labelsize': 11, 'xtick.labelsize': 11,
                             'ytick.labelsize': 11})
        fig = plt.figure(figsize=(166/25.4, 6.7), constrained_layout=True)
        grid = fig.add_gridspec(4, 1, height_ratios=[12, 1, 1, 1])
        ax = fig.add_subplot(grid[0, 0])
        ax.imshow(
            rgb,
            origin="lower",
            aspect="auto",
            interpolation="nearest",
            extent=[time_numbers[0], time_numbers[-1], ranges[0] - 0.3, ranges[-1] + 0.3],
        )
        ax.set_xlim(time_numbers[0], time_numbers[-1])
        ax.set_ylim(75, 100)
        ax.xaxis_date()
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=4, tz=dt.timezone.utc))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=dt.timezone.utc))
        ax.set_xlabel("Time (UTC)")
        ax.set_ylabel("Range (km)")
        label = "zenith beam 0" if beam == 0 else f"TX beam {beam}"
        ax.set_title(f"{args.day} | {label}")
        if args.day == "2025-01-30":
            imaging_time = dt.datetime(2025, 1, 30, 4, 15, 28, tzinfo=dt.timezone.utc)
            ax.axvline(imaging_time, color="white", linewidth=1.0, linestyle="--")

        add_mapping_bar(fig.add_subplot(grid[1, 0]), "doppler", args.doppler_max)
        add_mapping_bar(fig.add_subplot(grid[2, 0]), "snr", args.snr_max)
        add_mapping_bar(fig.add_subplot(grid[3, 0]), "width", args.width_max)
        if not args.hide_provenance:
            ax.text(
                0.995, 0.008,
                "plot_pansy_pmse_falsecolor.py; xc2 autocorrelations 0:7; nearest-time display fill",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5,
                color="white", alpha=0.8,
            )
        output = args.output.with_name(f"{args.output.stem}_beam{beam}{args.output.suffix}")
        fig.savefig(output, dpi=300)
        plt.close(fig)
        print(f"beam={beam} output={output}")
    print(f"blocks={len(blocks)}")


if __name__ == "__main__":
    main()
