#!/usr/bin/env python3
"""Plot representative PANSY PMSE range-Doppler spectra for all TX beams."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


SAMPLE_RATE_HZ = 1_000_000.0
RANGE_GATE_KM = 0.15
AUTOCORRELATION_COUNT = 7
RADAR_FREQUENCY_HZ = 47_000_000.0
# f_D = 2 f_radar v_radial / c; retain the stored Doppler sign.
MPS_PER_HZ = 299_792_458.0 / (2.0 * RADAR_FREQUENCY_HZ)


def decode(group):
    xc = group["xc_arr"]
    power = np.sum(np.abs(xc[:AUTOCORRELATION_COUNT, :, :, :]), axis=0)
    n_fft = int(group["n_fft"][()])
    ipp = int(group["ipp"][()])
    f0 = int(group["f0"][()])
    f1 = int(group["f1"][()])
    frequency = np.fft.fftshift(
        np.fft.fftfreq(n_fft, d=5 * ipp / SAMPLE_RATE_HZ)
    )[f0:f1]
    outer = np.abs(frequency) >= 0.75 * np.max(np.abs(frequency))
    noise = np.nanmedian(power[:, outer, :], axis=(1, 2))
    db = 10.0 * np.log10(np.maximum(power / noise[:, None, None], 1e-12))
    r0 = int(group["r0"][()])
    r1 = int(group["r1"][()])
    rdec = int(group["rdec"][()])
    ranges = np.arange(r0, r1, rdec) * RANGE_GATE_KM
    timestamp = int(group["i0"][()]) / SAMPLE_RATE_HZ
    return timestamp, frequency, ranges, db


def inventory(root: Path, day: str):
    records = []
    for hour in range(24):
        for path in sorted((root / f"{day}T{hour:02d}-00-00").glob("xc@*.h5")):
            with h5py.File(path, "r") as handle:
                for key in handle:
                    timestamp, frequency, ranges, db = decode(handle[key])
                    band = (ranges >= 75.0) & (ranges <= 100.0)
                    score = float(np.nanmax(db[:, :, band]))
                    records.append((score, timestamp, str(path), str(key)))
    return records


def choose(records, count: int):
    imaging = min(records, key=lambda row: abs(row[1] - 1738210528.0))
    chosen = [imaging]
    for record in sorted(records, reverse=True):
        if all(abs(record[1] - item[1]) >= 2 * 3600 for item in chosen):
            chosen.append(record)
            if len(chosen) == count:
                break
    return sorted(chosen, key=lambda row: row[1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--day", default="2025-01-30")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--examples", type=int, default=3)
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float, default=35.0)
    parser.add_argument("--hide-provenance", action="store_true")
    args = parser.parse_args()

    selected = choose(inventory(args.root, args.day), args.examples)
    decoded = []
    for _score, _timestamp, path, key in selected:
        with h5py.File(path, "r") as handle:
            decoded.append(decode(handle[key]))

    beam_count = decoded[0][3].shape[0]
    plt.rcParams.update({'font.size':11,'axes.titlesize':11,'axes.labelsize':11,
                         'xtick.labelsize':11,'ytick.labelsize':11})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for beam in range(beam_count):
        fig, axes = plt.subplots(len(decoded),1,figsize=(166/25.4,8.0),
                                 sharex=True,sharey=True,layout='constrained')
        axes = np.atleast_1d(axes)
        for row, (timestamp, frequency, ranges, db) in enumerate(decoded):
            stamp = dt.datetime.fromtimestamp(timestamp, dt.timezone.utc)
            ax = axes[row]
            mesh = ax.pcolormesh(
                frequency * MPS_PER_HZ, ranges, db[beam].T, shading="auto", cmap="plasma",
                vmin=args.vmin, vmax=args.vmax,
            )
            ax.set_xlim(-10 * MPS_PER_HZ, 10 * MPS_PER_HZ)
            ax.set_ylim(75, 100)
            ax.set_title(f"{stamp:%H:%M:%S} UTC")
            ax.set_ylabel("Range (km)")
            if row == len(decoded) - 1:
                ax.set_xlabel("Radial velocity (m/s)")
            fig.colorbar(mesh,ax=ax,label='Spectral power / noise (dB)')
        fig.suptitle(f'{args.day} | TX beam {beam}',fontsize=11)
        output=args.output.with_name(f'{args.output.stem}_beam{beam}{args.output.suffix}')
        fig.savefig(output,dpi=300)
        plt.close(fig)
    for score, timestamp, path, key in selected:
        stamp = dt.datetime.fromtimestamp(timestamp, dt.timezone.utc)
        print(f"{stamp.isoformat()} score_db={score:.2f} {path} group={key}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
