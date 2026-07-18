#!/usr/bin/env python3
"""Append decoded beat-phase acceleration to an existing event diagnostic."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

import pansy_config as pc


def load_nonoverlapping_phase_acceleration(path: Path) -> dict:
    with h5py.File(path, "r") as handle:
        group = handle["baud_averaged_beat_pairs"]
        pair = {name: np.asarray(values) for name, values in group.items()}
        measurement_tx_s = np.asarray(handle["measurement/tx_idx"], dtype=float) / 1e6
    previous = pair["previous_index"].astype(int)
    current = pair["current_index"].astype(int)
    beam = pair["beam_id"].astype(int)
    beat = pair["beat_real"] + 1j * pair["beat_imag"]
    rows = []
    for beam_id in np.unique(beam):
        indices = np.flatnonzero(beam == beam_id)
        indices = indices[np.argsort(pair["time_s"][indices])]
        candidates = []
        for first, second in zip(indices[:-1], indices[1:]):
            if previous[second] != current[first]:
                continue
            if abs(pair["delta_t_s"][second] - pair["delta_t_s"][first]) > 5e-6:
                continue
            candidates.append(
                (
                    pair["time_s"][first],
                    pair["time_s"][second],
                    pair["delta_t_s"][first],
                    pair["delta_t_s"][second],
                    np.angle(beat[second] * np.conj(beat[first])),
                    np.hypot(pair["beat_phase_std_rad"][first], pair["beat_phase_std_rad"][second]),
                )
            )
        rows.extend(candidates[::3])
    dtype = [
        ("first_time_s", "f8"),
        ("second_time_s", "f8"),
        ("first_dt_s", "f8"),
        ("second_dt_s", "f8"),
        ("observed_delta_phase_rad", "f8"),
        ("formal_phase_std_rad", "f8"),
    ]
    return {"samples": np.asarray(rows, dtype=dtype), "measurement_tx_s": measurement_tx_s}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-plot", type=Path, required=True)
    parser.add_argument("--beat-h5", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    phase_data = load_nonoverlapping_phase_acceleration(args.beat_h5)
    samples = phase_data["samples"]
    mean_dt = 0.5 * (samples["first_dt_s"] + samples["second_dt_s"])
    acceleration = pc.wavelength * samples["observed_delta_phase_rad"] / (4.0 * np.pi * mean_dt**2)
    acceleration_std = pc.wavelength * samples["formal_phase_std_rad"] / (4.0 * np.pi * mean_dt**2)
    time_s = 0.5 * (samples["first_time_s"] + samples["second_time_s"])
    time_s -= phase_data["measurement_tx_s"][0]

    with h5py.File(args.beat_h5, "r") as handle:
        model_velocity = np.asarray(handle["measurement/trajectory_model_doppler_mps"], dtype=float)
        model_time = np.asarray(handle["measurement/tx_idx"], dtype=float) / 1e6
    absolute_zero = phase_data["measurement_tx_s"][0]
    first_velocity = np.interp(samples["first_time_s"], model_time, model_velocity)
    second_velocity = np.interp(samples["second_time_s"], model_time, model_velocity)
    model_acceleration = (second_velocity - first_velocity) / (
        samples["second_time_s"] - samples["first_time_s"]
    )
    ambiguity_half_span = pc.wavelength / (4.0 * np.nanmedian(mean_dt) ** 2)
    finite_fit = np.isfinite(time_s) & np.isfinite(acceleration) & np.isfinite(acceleration_std)
    fit_keep = finite_fit.copy()
    fit_center = float(np.nanmedian(time_s[finite_fit]))
    fit_x = time_s - fit_center
    fit_coefficients = np.asarray([0.0, np.nanmedian(acceleration[finite_fit])])
    for _ in range(5):
        weights = 1.0 / np.maximum(acceleration_std[fit_keep], 1.0)
        fit_coefficients = np.polyfit(fit_x[fit_keep], acceleration[fit_keep], 1, w=weights)
        fit_residual = acceleration - np.polyval(fit_coefficients, fit_x)
        residual_median = np.nanmedian(fit_residual[fit_keep])
        residual_sigma = 1.4826 * np.nanmedian(np.abs(fit_residual[fit_keep] - residual_median))
        if not np.isfinite(residual_sigma) or residual_sigma <= 0.0:
            break
        fit_keep = finite_fit & (np.abs(fit_residual - residual_median) < 4.0 * residual_sigma)

    event_image = mpimg.imread(args.event_plot)
    fig = plt.figure(figsize=(16.0, 8.0), dpi=140)
    grid = fig.add_gridspec(1, 2, width_ratios=(3.15, 1.0), wspace=0.04)
    image_ax = fig.add_subplot(grid[0, 0])
    image_ax.imshow(event_image)
    image_ax.axis("off")

    ax = fig.add_subplot(grid[0, 1])
    ax.errorbar(
        time_s,
        acceleration / 1e3,
        yerr=acceleration_std / 1e3,
        fmt=".",
        ms=4,
        color="black",
        ecolor="0.78",
        elinewidth=0.6,
        capsize=0,
        label="decoded beat phase",
    )
    order = np.argsort(time_s)
    ax.plot(
        time_s[order],
        model_acceleration[order] / 1e3,
        color="C0",
        lw=1.2,
        label="shrinking-radius model",
    )
    ax.plot(
        time_s[order],
        np.polyval(fit_coefficients, fit_x[order]) / 1e3,
        color="C3",
        lw=1.3,
        label="weighted fit to phase measurements",
    )
    ax.axhline(ambiguity_half_span / 1e3, color="0.45", ls="--", lw=0.8)
    ax.axhline(-ambiguity_half_span / 1e3, color="0.45", ls="--", lw=0.8, label="phase ambiguity limits")
    finite = np.isfinite(acceleration)
    median = float(np.nanmedian(acceleration[finite]))
    weighted_mean = float(
        np.average(
            acceleration[finite],
            weights=1.0 / np.maximum(acceleration_std[finite], 1.0) ** 2,
        )
    )
    model_median = float(np.nanmedian(model_acceleration))
    fitted_median = float(np.polyval(fit_coefficients, 0.0))
    ax.text(
        0.04,
        0.97,
        f"N {np.count_nonzero(finite)}\n"
        + f"median {median / 1e3:.2f} km s$^{{-2}}$\n"
        + f"weighted {weighted_mean / 1e3:.2f} km s$^{{-2}}$\n"
        + f"phase fit {fitted_median / 1e3:.2f} km s$^{{-2}}$\n"
        + f"model {model_median / 1e3:.2f} km s$^{{-2}}$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(r"Radial acceleration (km s$^{-2}$)")
    ax.set_title("Decoded beat-phase deceleration", fontsize=11)
    ax.grid(alpha=0.2, lw=0.5)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
