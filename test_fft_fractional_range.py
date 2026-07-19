#!/usr/bin/env python3
"""Compare repeated-sample and FFT fractional-delay range estimates for one event."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from interferometer_alias_diagnostics import amp_scale, load_cut
import plot_interferometric_disambiguation as pid
from range_interpolation_study import (
    DRG_KM,
    FFTFractionalRangeDopplerSearch,
    PowerOnlyRangeDopplerSearch,
    quadratic_peak_2d,
)


def estimate(search, z_tx, z_rx, delays, pulses, range_scale):
    ranges = np.full(len(pulses), np.nan)
    dopplers = np.full(len(pulses), np.nan)
    for output_index, pulse in enumerate(pulses):
        if hasattr(search, "peak"):
            range_index, doppler = search.peak(z_tx[pulse], z_rx[pulse])
            ranges[output_index] = (delays[pulse] + range_index / range_scale) * DRG_KM
            dopplers[output_index] = doppler
            continue
        power = search.mf(z_tx[pulse], z_rx[pulse])
        noise = float(np.nanmedian(power))
        range_profile = np.nanmax(power, axis=1)
        ri = int(np.nanargmax(range_profile))
        di = int(np.nanargmax(power[ri]))
        dr, dd = quadratic_peak_2d(
            np.log(np.maximum(power, max(noise, 1e-12))), ri, di
        )
        ranges[output_index] = (delays[pulse] + (ri + dr) / range_scale) * DRG_KM
        doppler_step = float(search.dopv[1] - search.dopv[0])
        dopplers[output_index] = search.dopv[di] + dd * doppler_step
    return ranges, dopplers


def rms(values, keep):
    selected = np.asarray(values)[np.asarray(keep, dtype=bool)]
    selected = selected[np.isfinite(selected)]
    return float(np.sqrt(np.mean(selected**2)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-idx", type=int, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fft-pad", type=int, default=16)
    args = parser.parse_args()

    day = dt.datetime.fromtimestamp(
        args.sample_idx / 1e6, tz=dt.timezone.utc
    ).strftime("%Y-%m-%d")
    diagnostics = (
        args.base
        / "events"
        / day
        / f"pansy_disambiguation_diagnostics_{args.sample_idx}.h5"
    )
    with h5py.File(diagnostics, "r") as handle:
        label = handle.attrs["selected_hypothesis"]
        if isinstance(label, bytes):
            label = label.decode()
        group = handle["hypotheses"][label]
        time_s = np.asarray(group["t_rel_s"], dtype=float)
        model_range = np.linalg.norm(
            np.asarray(group["physics_ceplecha_model"], dtype=float), axis=1
        )
        keep = np.asarray(
            group.get("physics_ceplecha_keep", np.ones(len(time_s))), dtype=bool
        )
        selected_range_time_component = int(
            handle.attrs.get("selected_range_time_component", -1)
        )

    cut = load_cut(args.base / "metadata/cut", args.sample_idx)
    cut_observations = pid.recompute_cut_observables(cut, interp=1)
    snr_gate = np.asarray(cut_observations["snr"], dtype=float) > 7.0
    observations_for_clock = {
        key: value[snr_gate]
        if isinstance(value, np.ndarray) and len(value) == len(snr_gate)
        else value
        for key, value in cut_observations.items()
    }
    if selected_range_time_component >= 0:
        segments = pid.split_observations_by_range_time(
            observations_for_clock, min_points=3
        )
        if selected_range_time_component < len(segments):
            observations_for_clock = pid.subset_pulse_observations(
                observations_for_clock, segments[selected_range_time_component]
            )
    diagnostic_zero_s = float(
        np.asarray(observations_for_clock["tx_idx"], dtype=float)[0]
    ) / 1e6
    measurement_absolute_s = diagnostic_zero_s + time_s
    raw_absolute_s = np.asarray(cut["tx_idx"], dtype=float) / 1e6
    pulses = np.asarray(
        [int(np.argmin(np.abs(raw_absolute_s - value))) for value in measurement_absolute_s],
        dtype=int,
    )
    z_rx = np.asarray(cut["zrx_echoes_re"], dtype=np.complex64) + 1j * np.asarray(
        cut["zrx_echoes_im"], dtype=np.complex64
    )
    z_tx = np.asarray(cut["ztx_pulses_re"], dtype=np.complex64) + 1j * np.asarray(
        cut["ztx_pulses_im"], dtype=np.complex64
    )
    z_rx *= amp_scale()[None, :, None]
    delays = np.asarray(cut["delays"], dtype=int)

    methods = []
    repeated = PowerOnlyRangeDopplerSearch(
        z_tx.shape[1], z_rx.shape[2], z_rx.shape[1], interp=2, fft_pad=args.fft_pad
    )
    repeated_range, _ = estimate(repeated, z_tx, z_rx, delays, pulses, 2.0)
    methods.append(("Repeated samples x2", repeated_range))
    for oversample in (2, 4, 8):
        search = FFTFractionalRangeDopplerSearch(
            z_tx.shape[1],
            z_rx.shape[2],
            z_rx.shape[1],
            range_oversample=oversample,
            fft_pad=args.fft_pad,
        )
        fft_range, _ = estimate(
            search, z_tx, z_rx, delays, pulses, float(oversample)
        )
        methods.append((f"FFT fractional delay x{oversample}", fft_range))

    fig, axes = plt.subplots(2, 1, figsize=(8.5, 7.0), sharex=True, constrained_layout=True)
    for label, measured_range in methods:
        residual = measured_range - model_range
        value = rms(residual, keep)
        axes[0].plot(time_s[keep], measured_range[keep], ".", ms=3, label=label)
        axes[1].plot(
            time_s[keep], residual[keep] * 1e3, ".", ms=3, label=f"{label}: {value * 1e3:.1f} m"
        )
        print(label, f"range_rms_m={value * 1e3:.6g}", flush=True)
    axes[0].plot(time_s[keep], model_range[keep], color="black", lw=1.0, label="best-fit range")
    axes[0].set_ylabel("Range (km)")
    axes[1].axhline(0.0, color="black", lw=0.7)
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Range residual (m)")
    for axis in axes:
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, fontsize=8)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=190)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
