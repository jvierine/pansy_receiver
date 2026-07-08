#!/usr/bin/env python3
"""Fit cut-metadata noise power to PyGDSM sky noise plus receiver temperature."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

import paper_plot_noise as ppn
import pansy_config as pc


def load_measurements(path: Path) -> dict[str, np.ndarray | str | int]:
    data = np.load(path)
    return {
        "day": str(data["day"]),
        "times_s": np.asarray(data["times_s"], dtype=np.float64),
        "beam_id": np.asarray(data["beam_id"], dtype=np.int16),
        "noise_power": np.asarray(data["noise_power"], dtype=np.float64),
        "guard_samples": int(data["guard_samples"]),
    }


def minute_edges(day: str) -> tuple[np.ndarray, np.ndarray]:
    start_us, _end_us, _date = ppn.utc_day_bounds(day)
    edges = start_us / 1e6 + np.arange(1441, dtype=np.float64) * 60.0
    centers = edges[:-1] + 30.0
    return edges, centers


def bin_cut_power(day_data: dict[str, np.ndarray | str | int]) -> dict[str, np.ndarray]:
    times_s = np.asarray(day_data["times_s"], dtype=np.float64)
    beam_id = np.asarray(day_data["beam_id"], dtype=np.int16)
    power = np.asarray(day_data["noise_power"], dtype=np.float64)
    edges, centers = minute_edges(str(day_data["day"]))
    minute = np.searchsorted(edges, times_s, side="right") - 1
    good = (
        (minute >= 0)
        & (minute < 1440)
        & (beam_id >= 0)
        & (beam_id < len(ppn.BEAM_NAMES))
        & np.isfinite(power)
        & (power > 0.0)
    )
    flat_index = beam_id[good].astype(np.int64) * 1440 + minute[good].astype(np.int64)
    counts = np.bincount(flat_index, minlength=len(ppn.BEAM_NAMES) * 1440).reshape(len(ppn.BEAM_NAMES), 1440)
    sums = np.bincount(flat_index, weights=power[good], minlength=len(ppn.BEAM_NAMES) * 1440).reshape(len(ppn.BEAM_NAMES), 1440)
    sums2 = np.bincount(flat_index, weights=power[good] ** 2, minlength=len(ppn.BEAM_NAMES) * 1440).reshape(len(ppn.BEAM_NAMES), 1440)
    mean = np.full_like(sums, np.nan, dtype=np.float64)
    var = np.full_like(sums, np.nan, dtype=np.float64)
    nonzero = counts > 0
    mean[nonzero] = sums[nonzero] / counts[nonzero]
    var[nonzero] = np.maximum(0.0, sums2[nonzero] / counts[nonzero] - mean[nonzero] ** 2)
    sem = np.full_like(mean, np.nan)
    sem[nonzero] = np.sqrt(var[nonzero] / counts[nonzero])
    return {"centers_s": centers, "counts": counts, "mean_power": mean, "sem_power": sem}


def interpolate_sky_model(model_times_s: np.ndarray, model_tsky: np.ndarray, target_times_s: np.ndarray) -> np.ndarray:
    out = np.full((len(ppn.BEAM_NAMES), len(target_times_s)), np.nan, dtype=np.float64)
    for beam_i in range(len(ppn.BEAM_NAMES)):
        out[beam_i] = np.interp(target_times_s, model_times_s, model_tsky[:, beam_i])
    return out


def fit_linear_trec(tsky: np.ndarray, power: np.ndarray, weights: np.ndarray) -> dict[str, float | np.ndarray]:
    good = np.isfinite(tsky) & np.isfinite(power) & np.isfinite(weights) & (weights > 0.0)
    if np.count_nonzero(good) < 3:
        raise RuntimeError("not enough finite samples for fit")
    x = np.asarray(tsky[good], dtype=np.float64)
    y = np.asarray(power[good], dtype=np.float64)
    w = np.asarray(weights[good], dtype=np.float64)
    xmat = np.column_stack([x, np.ones_like(x)])
    sw = np.sqrt(w / np.nanmedian(w))
    beta = np.linalg.lstsq(xmat * sw[:, None], y * sw, rcond=None)[0]
    slope, intercept = float(beta[0]), float(beta[1])
    if slope <= 0.0:
        raise RuntimeError(f"non-positive fitted slope {slope}")
    model = slope * x + intercept
    resid = y - model
    dof = max(1, len(y) - 2)
    sigma2 = float(np.sum((sw * resid) ** 2) / dof)
    cov = sigma2 * np.linalg.inv((xmat * sw[:, None]).T @ (xmat * sw[:, None]))
    grad = np.asarray([-intercept / slope**2, 1.0 / slope], dtype=np.float64)
    trec_var = float(grad @ cov @ grad)
    return {
        "slope": slope,
        "intercept": intercept,
        "receiver_temp_k": intercept / slope,
        "receiver_temp_sigma_k": np.sqrt(max(0.0, trec_var)),
        "used": good,
    }


def rolling_nanmedian(values: np.ndarray, half_width: int) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    out = np.full_like(values, np.nan)
    for i in range(len(values)):
        vals = values[max(0, i - half_width) : min(len(values), i + half_width + 1)]
        vals = vals[np.isfinite(vals)]
        if len(vals):
            out[i] = float(np.nanmedian(vals))
    return out


def robust_jump_outlier_mask(
    power_db: np.ndarray,
    counts: np.ndarray,
    min_count: int,
    half_width: int = 20,
    sigma_threshold: float = 6.0,
    db_threshold: float = 1.2,
    jump_db: float = 2.0,
) -> tuple[np.ndarray, np.ndarray]:
    mask = np.zeros_like(power_db, dtype=bool)
    score = np.full_like(power_db, np.nan, dtype=np.float64)
    for beam_i in range(power_db.shape[0]):
        series = np.asarray(power_db[beam_i], dtype=np.float64)
        baseline = rolling_nanmedian(series, half_width)
        residual = series - baseline
        local_mad = rolling_nanmedian(np.abs(residual), half_width)
        sigma = 1.4826 * np.maximum(local_mad, 0.08)
        score[beam_i] = residual / sigma
        jump = np.full_like(series, np.nan)
        jump[1:-1] = np.maximum(np.abs(series[1:-1] - series[:-2]), np.abs(series[1:-1] - series[2:]))
        valid = np.isfinite(series) & (counts[beam_i] >= min_count)
        mask[beam_i] = valid & (
            (np.abs(residual) > np.maximum(db_threshold, sigma_threshold * sigma))
            | (jump > jump_db)
        )
    return mask, score


def bootstrap_trec(
    minute_index: np.ndarray,
    tsky: np.ndarray,
    power: np.ndarray,
    weights: np.ndarray,
    n_bootstrap: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    finite = np.isfinite(tsky) & np.isfinite(power) & np.isfinite(weights) & (weights > 0.0)
    minute_index = minute_index[finite]
    tsky = tsky[finite]
    power = power[finite]
    weights = weights[finite]
    unique_minutes = np.unique(minute_index)
    trec = np.full(n_bootstrap, np.nan, dtype=np.float64)
    by_minute = {minute: np.where(minute_index == minute)[0] for minute in unique_minutes}
    for bi in range(n_bootstrap):
        sampled = rng.choice(unique_minutes, size=len(unique_minutes), replace=True)
        idx = np.concatenate([by_minute[minute] for minute in sampled])
        try:
            trec[bi] = float(fit_linear_trec(tsky[idx], power[idx], weights[idx])["receiver_temp_k"])
        except Exception:
            pass
    return trec[np.isfinite(trec)]


def plot_fit(
    output: Path,
    day: str,
    centers_s: np.ndarray,
    counts: np.ndarray,
    mean_power: np.ndarray,
    tsky: np.ndarray,
    fit: dict[str, float | np.ndarray],
    trec_bootstrap: np.ndarray,
    min_count: int,
    outlier_mask: np.ndarray,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    colors = ["black", "tab:blue", "tab:orange", "tab:green", "tab:red"]
    tms = centers_s.astype(np.int64).astype("datetime64[s]")
    slope = float(fit["slope"])
    intercept = float(fit["intercept"])
    trec = float(fit["receiver_temp_k"])
    if len(trec_bootstrap):
        lo, hi = np.nanpercentile(trec_bootstrap, [16.0, 84.0])
        err_lo, err_hi = trec - lo, hi - trec
    else:
        err_lo = err_hi = float(fit["receiver_temp_sigma_k"])

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.4), constrained_layout=True, sharex=True)
    for beam_i, name in enumerate(ppn.BEAM_NAMES):
        good = (counts[beam_i] >= min_count) & np.isfinite(mean_power[beam_i]) & ~outlier_mask[beam_i]
        bad = outlier_mask[beam_i]
        axes[0].plot(tms[good], 10.0 * np.log10(mean_power[beam_i, good]), ".", ms=2.5, color=colors[beam_i], alpha=0.45)
        if np.any(bad):
            axes[0].plot(tms[bad], 10.0 * np.log10(mean_power[beam_i, bad]), "x", ms=4.0, color=colors[beam_i], alpha=0.9)
        model_power = slope * tsky[beam_i] + intercept
        axes[0].plot(tms, 10.0 * np.log10(np.maximum(model_power, 1.0)), "-", color=colors[beam_i], lw=1.6, label=name)
    axes[0].set_title(
        f"PANSY cut-noise power fit, {day}: "
        rf"$T_{{rec}}={trec:.0f}^{{+{err_hi:.0f}}}_{{-{err_lo:.0f}}}$ K"
    )
    axes[0].set_ylabel("Raw noise power (dB)")
    def db_to_kelvin(power_db):
        return np.power(10.0, np.asarray(power_db) / 10.0) / slope

    def kelvin_to_db(temp_k):
        return 10.0 * np.log10(np.maximum(np.asarray(temp_k) * slope, 1.0))

    temp_axis = axes[0].secondary_yaxis("right", functions=(db_to_kelvin, kelvin_to_db))
    temp_axis.set_ylabel(r"Equivalent $T_{\mathrm{sys}}$ (K)")
    axes[0].legend(loc="upper right", ncol=5, fontsize=8)

    for beam_i, name in enumerate(ppn.BEAM_NAMES):
        axes[1].plot(tms, tsky[beam_i] + trec, "-", color=colors[beam_i], lw=1.5, label=f"{name} $T_{{sky}}+T_{{rec}}$")
    axes[1].axhline(trec, color="black", lw=1.6, label=rf"$T_{{rec}}={trec:.0f}$ K")
    axes[1].axhspan(trec - err_lo, trec + err_hi, color="0.3", alpha=0.16, label="68% bootstrap")
    axes[1].set_ylabel("Temperature (K)")
    axes[1].set_xlabel("Time (UTC)")
    axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axes[1].legend(loc="upper right", fontsize=7, ncol=2)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--day", default="2025-05-10")
    parser.add_argument("--metadata", type=Path, default=Path("/mnt/data/juha/pansy/metadata/cut"))
    parser.add_argument("--measurements", type=Path, default=None, help="Existing cut noise measurement NPZ.")
    parser.add_argument("--output", type=Path, default=Path("figs/cut_noise_pygdsm_fit_2025-05-10.png"))
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--chunk-seconds", type=int, default=60)
    parser.add_argument("--guard-samples", type=int, default=25)
    parser.add_argument("--min-count", type=int, default=20)
    parser.add_argument("--gain-model", choices=["rx", "tx", "two_way"], default="rx")
    parser.add_argument("--sky-step-minutes", type=float, default=10.0)
    parser.add_argument("--n-az", type=int, default=180)
    parser.add_argument("--n-el", type=int, default=60)
    parser.add_argument("--min-elevation-deg", type=float, default=0.5)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20250708)
    args = parser.parse_args()

    if args.measurements is None:
        day_data = ppn.read_cut_day(
            args.metadata,
            args.day,
            workers=args.workers,
            chunk_seconds=args.chunk_seconds,
            guard_samples=args.guard_samples,
        )
    else:
        day_data = load_measurements(args.measurements)
        args.day = str(day_data["day"])
    binned = bin_cut_power(day_data)
    centers_s = binned["centers_s"]
    sky_step_s = float(args.sky_step_minutes) * 60.0
    model_times = np.arange(centers_s[0], centers_s[-1] + 0.5 * sky_step_s, sky_step_s)
    model_tsky = ppn.sky_temperature_matrix(
        model_times,
        gain_model=args.gain_model,
        n_az=args.n_az,
        n_el=args.n_el,
        min_elevation_deg=args.min_elevation_deg,
        rx_channel=None,
        freq_mhz=pc.freq / 1e6,
    )
    tsky = interpolate_sky_model(model_times, model_tsky, centers_s)
    counts = np.asarray(binned["counts"], dtype=np.float64)
    mean_power = np.asarray(binned["mean_power"], dtype=np.float64)
    power_db = 10.0 * np.log10(np.maximum(mean_power, 1.0))
    outlier_mask, robust_score = robust_jump_outlier_mask(power_db, counts, args.min_count)
    keep = (counts >= args.min_count) & np.isfinite(mean_power) & ~outlier_mask
    minute_grid = np.broadcast_to(np.arange(1440, dtype=np.int64), mean_power.shape)
    fit = fit_linear_trec(tsky[keep], mean_power[keep], counts[keep])
    trec_bootstrap = bootstrap_trec(
        minute_grid[keep],
        tsky[keep],
        mean_power[keep],
        counts[keep],
        n_bootstrap=args.bootstrap,
        seed=args.seed,
    )
    plot_fit(args.output, args.day, centers_s, counts, mean_power, tsky, fit, trec_bootstrap, args.min_count, outlier_mask)
    print(f"metadata {args.metadata}")
    print(f"output {args.output}")
    print(f"samples {len(day_data['noise_power'])}")
    print(f"fit_bins {np.count_nonzero(keep)}")
    print(f"outlier_bins {np.count_nonzero(outlier_mask)}")
    print(f"slope_power_per_k {float(fit['slope']):.9g}")
    print(f"intercept_power {float(fit['intercept']):.9g}")
    print(f"receiver_temp_k {float(fit['receiver_temp_k']):.6g}")
    print(f"receiver_temp_sigma_k {float(fit['receiver_temp_sigma_k']):.6g}")
    if len(trec_bootstrap):
        lo, med, hi = np.nanpercentile(trec_bootstrap, [16.0, 50.0, 84.0])
        print(f"receiver_temp_bootstrap_p16_p50_p84_k {lo:.6g} {med:.6g} {hi:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
