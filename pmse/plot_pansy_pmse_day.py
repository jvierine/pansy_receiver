#!/usr/bin/env python3
"""Plot one UTC day of PANSY PMSE peak spectral intensity from xc2 metadata."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np


SAMPLE_RATE_HZ = 1_000_000.0
RANGE_GATE_KM = 0.15
AUTOCORRELATION_COUNT = 7


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
                    peak_db = 10.0 * np.log10(
                        np.maximum(np.nanmax(power, axis=1) / noise[:, None], 1e-12)
                    )
                    r0 = int(group["r0"][()])
                    r1 = int(group["r1"][()])
                    rdec = int(group["rdec"][()])
                    ranges = np.arange(r0, r1, rdec) * RANGE_GATE_KM
                    blocks.append(
                        (
                            int(group["i0"][()]) / SAMPLE_RATE_HZ,
                            int(group["i1"][()]) / SAMPLE_RATE_HZ,
                            ranges,
                            peak_db,
                            str(path),
                            str(key),
                        )
                    )
    return sorted(blocks, key=lambda item: item[0])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--day", default="2025-01-30")
    parser.add_argument("--beam", default="all", help="beam number or 'all'")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--minute-resolution", type=float, default=1.0)
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float, default=35.0)
    parser.add_argument("--hide-provenance", action="store_true")
    args = parser.parse_args()

    blocks = read_blocks(args.root, args.day)
    if not blocks:
        raise SystemExit(f"no xc2 blocks found for {args.day} under {args.root}")

    reference_ranges = blocks[0][2]
    for block in blocks:
        if not np.array_equal(block[2], reference_ranges):
            raise RuntimeError("range grids differ within the requested day")

    day_start = utc_timestamp(args.day)
    day_end = day_start + 86400.0
    cadence_s = args.minute_resolution * 60.0
    time_edges = np.arange(day_start, day_end + cadence_s, cadence_s)
    range_step = float(np.median(np.diff(reference_ranges)))
    range_edges = np.r_[
        reference_ranges - range_step / 2.0, reference_ranges[-1] + range_step / 2.0
    ]
    time_numbers = mdates.date2num(
        [dt.datetime.fromtimestamp(value, dt.timezone.utc) for value in time_edges]
    )

    beam_count = blocks[0][3].shape[0]
    beams = list(range(beam_count)) if args.beam == "all" else [int(args.beam)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for beam in beams:
        if not 0 <= beam < beam_count:
            raise ValueError(f"beam {beam} is unavailable; valid beams are 0..{beam_count - 1}")
        image = np.full((len(reference_ranges), len(time_edges) - 1), np.nan)
        hits = np.zeros(len(time_edges) - 1, dtype=np.int32)
        for start, end, _ranges, peak_db, _path, _key in blocks:
            lo = max(0, int(np.floor((start - day_start) / cadence_s)))
            hi = min(len(hits), int(np.ceil((end - day_start) / cadence_s)))
            for index in range(lo, hi):
                if hits[index] == 0:
                    image[:, index] = peak_db[beam]
                else:
                    image[:, index] = np.maximum(image[:, index], peak_db[beam])
                hits[index] += 1

        # Fill duty-cycle gaps for display, matching the false-color RTDI.
        available = np.flatnonzero(hits > 0)
        for index in np.flatnonzero(hits == 0):
            nearest = available[np.argmin(np.abs(available - index))]
            image[:, index] = image[:, nearest]

        plt.rcParams.update({'font.size': 11, 'axes.titlesize': 11,
                             'axes.labelsize': 11, 'xtick.labelsize': 11,
                             'ytick.labelsize': 11})
        fig, ax = plt.subplots(figsize=(166/25.4, 3.8), constrained_layout=True)
        mesh = ax.pcolormesh(
            time_numbers, range_edges, image, shading="flat", cmap="plasma",
            vmin=args.vmin, vmax=args.vmax,
        )
        colorbar = fig.colorbar(mesh, ax=ax, pad=0.012)
        colorbar.set_label("Peak spectral power / noise (dB)")
        ax.set_ylabel("Range (km)")
        ax.set_xlabel("Time (UTC)")
        ax.set_ylim(75, 100)
        ax.set_xlim(
            mdates.date2num([
                dt.datetime.fromtimestamp(day_start, dt.timezone.utc),
                dt.datetime.fromtimestamp(day_end, dt.timezone.utc),
            ])
        )
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=4, tz=dt.timezone.utc))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=dt.timezone.utc))
        imaging_time = dt.datetime(2025, 1, 30, 4, 15, 28, tzinfo=dt.timezone.utc)
        if args.day == "2025-01-30":
            ax.axvline(imaging_time, color="cyan", linewidth=1.1, linestyle="--")
            ax.text(mdates.date2num(imaging_time), 99.5, " imaging block 04:15:28",
                    color="cyan", fontsize=11, va="top", ha="left")
        label = "zenith beam 0" if beam == 0 else f"TX beam {beam}"
        ax.set_title(f"{args.day} | {label}")
        if not args.hide_provenance:
            fig.text(
                0.995, 0.003,
                "plot_pansy_pmse_day.py; xc2 autocorrelations 0:7; max over Doppler; nearest-time display fill",
                ha="right", va="bottom", fontsize=6.5, color="0.35",
            )
        output = args.output
        if len(beams) > 1:
            output = output.with_name(f"{output.stem}_beam{beam}{output.suffix}")
        fig.savefig(output, dpi=300)
        plt.close(fig)
        print(f"beam={beam} covered_minutes={np.count_nonzero(hits)} output={output}")
    print(f"blocks={len(blocks)}")
    print(f"range_km={reference_ranges[0]:.2f}..{reference_ranges[-1]:.2f}")


if __name__ == "__main__":
    main()
