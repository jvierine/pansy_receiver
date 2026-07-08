#!/usr/bin/env python3
"""Debug plots for PANSY cut-metadata noise estimates."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

import paper_plot_noise as ppn


def read_example_cut(metadata_path: Path, day: str):
    import digital_rf as drf

    start_us, end_us, _date = ppn.utc_day_bounds(day)
    dm = drf.DigitalMetadataReader(str(metadata_path))
    for start in range(start_us, end_us, int(60e6)):
        records = dm.read(start, min(start + int(60e6), end_us))
        for key in sorted(records):
            record = records[key]
            if "zrx_echoes_re" in record and "zrx_echoes_im" in record:
                return key, record
    return None, None


def example_power_and_mask(record, guard_samples: int) -> tuple[np.ndarray, np.ndarray, int, int]:
    zrx_re = np.asarray(record["zrx_echoes_re"], dtype=np.float32)
    zrx_im = np.asarray(record["zrx_echoes_im"], dtype=np.float32)
    pad = int(np.asarray(record["pad"]).reshape(()))
    txlen = int(np.asarray(record["txlen"]).reshape(()))
    sample_power = np.nanmean(zrx_re**2 + zrx_im**2, axis=(0, 1))
    n_sample = len(sample_power)
    used = np.zeros(n_sample, dtype=bool)
    used[: max(0, min(pad - guard_samples, n_sample))] = True
    post_start = max(0, min(pad + txlen + guard_samples, n_sample))
    used[post_start:] = True
    return sample_power, used, pad, txlen


def minute_bins(day: str) -> tuple[np.ndarray, np.ndarray]:
    start_us, _end_us, date = ppn.utc_day_bounds(day)
    edges_s = start_us / 1e6 + np.arange(1441, dtype=np.float64) * 60.0
    centers = edges_s[:-1] + 30.0
    return edges_s, centers.astype("datetime64[s]")


def binned_counts_and_medians(times_s: np.ndarray, beam_id: np.ndarray, power: np.ndarray, day: str):
    edges_s, centers = minute_bins(day)
    minute = np.searchsorted(edges_s, times_s, side="right") - 1
    valid = (minute >= 0) & (minute < 1440) & np.isfinite(power) & (beam_id >= 0) & (beam_id < len(ppn.BEAM_NAMES))
    minute = minute[valid]
    beam_id = beam_id[valid]
    power = power[valid]
    counts = np.zeros((len(ppn.BEAM_NAMES), 1440), dtype=np.int64)
    medians = np.full((len(ppn.BEAM_NAMES), 1440), np.nan, dtype=np.float64)
    for beam_i in range(len(ppn.BEAM_NAMES)):
        for minute_i in np.unique(minute[beam_id == beam_i]):
            idx = (beam_id == beam_i) & (minute == minute_i)
            counts[beam_i, minute_i] = int(np.count_nonzero(idx))
            medians[beam_i, minute_i] = float(np.nanmedian(power[idx]))
    return centers, counts, medians, valid


def power_db(power: np.ndarray) -> np.ndarray:
    return 10.0 * np.log10(np.maximum(np.asarray(power, dtype=np.float64), 1.0))


def plot_debug(day_data: dict[str, np.ndarray | str], output: Path, example_record=None) -> dict[str, float]:
    times_s = np.asarray(day_data["times_s"], dtype=np.float64)
    beam_id = np.asarray(day_data["beam_id"], dtype=np.int16)
    power = np.asarray(day_data["noise_power"], dtype=np.float64)
    day = str(day_data["day"])
    centers, counts, medians, valid = binned_counts_and_medians(times_s, beam_id, power, day)
    tms = centers
    guard_samples = int(day_data.get("guard_samples", 25))

    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 1, figsize=(10.0, 8.0), constrained_layout=True)
    colors = ["black", "tab:blue", "tab:orange", "tab:green", "tab:red"]
    used_fraction = np.nan
    if example_record is not None:
        sample_power, used, pad, txlen = example_power_and_mask(example_record, guard_samples)
        x = np.arange(len(sample_power))
        sample_db = power_db(sample_power)
        axes[0].plot(x, sample_db, color="0.15", lw=1.2)
        axes[0].fill_between(x, sample_db, where=used, color="tab:green", alpha=0.30, label="used padding")
        axes[0].fill_between(x, sample_db, where=~used, color="tab:red", alpha=0.18, label="excluded echo window")
        axes[0].axvline(pad, color="tab:red", lw=0.8)
        axes[0].axvline(pad + txlen, color="tab:red", lw=0.8)
        axes[0].set_ylabel("Raw power (dB)")
        used_fraction = float(np.count_nonzero(used) / len(used))
        axes[0].set_title(
            f"Used samples in one cut: {np.count_nonzero(used)}/{len(used)} "
            f"({100.0 * used_fraction:.1f}%), guard {guard_samples} samples"
        )
        axes[0].legend(loc="upper right", fontsize=8)
    else:
        axes[0].text(0.5, 0.5, "No example cut found", ha="center", va="center", transform=axes[0].transAxes)
        axes[0].set_axis_off()

    for beam_i, name in enumerate(ppn.BEAM_NAMES):
        axes[1].plot(tms, counts[beam_i], color=colors[beam_i], lw=1.0, label=name)
        axes[2].plot(tms, power_db(medians[beam_i]), color=colors[beam_i], lw=1.0, label=name)

    axes[1].set_title(f"Cut padding noise debug, {day}; guard {guard_samples} samples; per-pulse mean power")
    axes[1].set_ylabel("Pulses / minute")
    axes[2].set_ylabel("Binned power (dB)")
    axes[1].legend(loc="upper right", ncol=5, fontsize=8)
    for ax in axes[1:3]:
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax.set_xlabel("Time (UTC)")
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return {
        "samples": float(len(power)),
        "valid_samples": float(np.count_nonzero(valid)),
        "used_fraction_per_cut_sample": float(used_fraction),
        "median_power": float(np.nanmedian(power)),
        "output": str(output),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--day", default="2025-05-10")
    parser.add_argument("--metadata", type=Path, default=Path("/mnt/data/juha/pansy/metadata/cut"))
    parser.add_argument("--output", type=Path, default=Path("figs/cut_noise_debug_2025-05-10.png"))
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--chunk-seconds", type=int, default=60)
    parser.add_argument("--guard-samples", type=int, default=25)
    args = parser.parse_args()

    example_key, example_record = read_example_cut(args.metadata, args.day)
    data = ppn.read_cut_day(
        args.metadata,
        args.day,
        workers=args.workers,
        chunk_seconds=args.chunk_seconds,
        guard_samples=args.guard_samples,
    )
    stats = plot_debug(data, args.output, example_record=example_record)
    print(f"metadata {args.metadata}")
    print(f"example_key {example_key}")
    for key, value in stats.items():
        print(f"{key} {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
