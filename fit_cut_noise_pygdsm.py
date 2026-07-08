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


def receive_sky_temperature_matrix(
    times_s: np.ndarray,
    n_az: int,
    n_el: int,
    min_elevation_deg: float,
) -> np.ndarray:
    """Beam-weighted sky temperature using the one-way receive pattern.

    Thermal sky noise enters through the receiving antenna pattern only. This is
    deliberately not configurable in the system-noise fit.
    """
    return ppn.sky_temperature_matrix(
        times_s,
        gain_model="rx",
        n_az=n_az,
        n_el=n_el,
        min_elevation_deg=min_elevation_deg,
        rx_channel=None,
        freq_mhz=pc.freq / 1e6,
    )


def load_measurements(path: Path) -> dict[str, np.ndarray | str | int]:
    import h5py

    if path.suffix.lower() == ".npz":
        raise ValueError("use HDF5 measurement files, not NPZ")
    with h5py.File(path, "r") as data:
        day = data.attrs.get("day", "")
        guard_samples = int(data.attrs.get("guard_samples", 25))
        times_s = np.asarray(data["times_s"][()], dtype=np.float64)
        beam_id = np.asarray(data["beam_id"][()], dtype=np.int16)
        noise_power = np.asarray(data["noise_power"][()], dtype=np.float64)
        per_module = int(data.attrs.get("per_module", 1 if noise_power.ndim == 2 else 0))
    return {
        "day": str(day),
        "times_s": times_s,
        "beam_id": beam_id,
        "noise_power": noise_power,
        "guard_samples": guard_samples,
        "per_module": per_module,
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
    if power.ndim == 1:
        power = power[:, None]
    if power.ndim != 2:
        raise ValueError("noise_power must have shape (pulse,) or (pulse, module)")
    n_module = power.shape[1]
    edges, centers = minute_edges(str(day_data["day"]))
    minute = np.searchsorted(edges, times_s, side="right") - 1
    base_good = (
        (minute >= 0)
        & (minute < 1440)
        & (beam_id >= 0)
        & (beam_id < len(ppn.BEAM_NAMES))
    )
    counts = np.zeros((n_module, len(ppn.BEAM_NAMES), 1440), dtype=np.int64)
    sums = np.zeros_like(counts, dtype=np.float64)
    sums2 = np.zeros_like(counts, dtype=np.float64)
    for module_i in range(n_module):
        good = base_good & np.isfinite(power[:, module_i]) & (power[:, module_i] > 0.0)
        flat_index = beam_id[good].astype(np.int64) * 1440 + minute[good].astype(np.int64)
        counts[module_i] = np.bincount(flat_index, minlength=len(ppn.BEAM_NAMES) * 1440).reshape(len(ppn.BEAM_NAMES), 1440)
        sums[module_i] = np.bincount(flat_index, weights=power[good, module_i], minlength=len(ppn.BEAM_NAMES) * 1440).reshape(
            len(ppn.BEAM_NAMES), 1440
        )
        sums2[module_i] = np.bincount(
            flat_index,
            weights=power[good, module_i] ** 2,
            minlength=len(ppn.BEAM_NAMES) * 1440,
        ).reshape(len(ppn.BEAM_NAMES), 1440)
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


def fit_receiver_temperature(tsky: np.ndarray, power: np.ndarray, weights: np.ndarray, max_trec_k: float = 20000.0) -> dict[str, float | np.ndarray]:
    good = np.isfinite(tsky) & np.isfinite(power) & np.isfinite(weights) & (weights > 0.0)
    if np.count_nonzero(good) < 3:
        raise RuntimeError("not enough finite samples for fit")
    x = np.asarray(tsky[good], dtype=np.float64)
    y = np.asarray(power[good], dtype=np.float64)
    w = np.asarray(weights[good], dtype=np.float64)

    def slope_for(trec: float) -> float:
        xt = x + trec
        return float(np.sum(w * xt * y) / np.sum(w * xt * xt))

    def objective(trec: float) -> float:
        slope = slope_for(trec)
        resid = y - slope * (x + trec)
        return float(np.sum(w * resid * resid))

    # Bounded one-dimensional search. The coarse pass prevents the golden
    # search from wandering when the optimum is pinned to the physical boundary.
    coarse = np.linspace(0.0, float(max_trec_k), 401)
    coarse_obj = np.asarray([objective(float(trec)) for trec in coarse])
    best_i = int(np.nanargmin(coarse_obj))
    if best_i == 0:
        trec = 0.0
    elif best_i == len(coarse) - 1:
        trec = float(max_trec_k)
    else:
        lo = float(coarse[best_i - 1])
        hi = float(coarse[best_i + 1])
        gr = (np.sqrt(5.0) - 1.0) / 2.0
        c = hi - gr * (hi - lo)
        d = lo + gr * (hi - lo)
        fc = objective(c)
        fd = objective(d)
        for _ in range(80):
            if abs(hi - lo) < 1e-3:
                break
            if fc < fd:
                hi = d
                d = c
                fd = fc
                c = hi - gr * (hi - lo)
                fc = objective(c)
            else:
                lo = c
                c = d
                fc = fd
                d = lo + gr * (hi - lo)
                fd = objective(d)
        trec = 0.5 * (lo + hi)
    slope = slope_for(trec)
    intercept = slope * trec
    return {
        "slope": slope,
        "intercept": intercept,
        "receiver_temp_k": trec,
        "receiver_temp_sigma_k": np.nan,
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
    flat_power = power_db.reshape((-1, power_db.shape[-1]))
    flat_counts = counts.reshape((-1, counts.shape[-1]))
    flat_mask = mask.reshape((-1, mask.shape[-1]))
    flat_score = score.reshape((-1, score.shape[-1]))
    for series_i in range(flat_power.shape[0]):
        series = np.asarray(flat_power[series_i], dtype=np.float64)
        baseline = rolling_nanmedian(series, half_width)
        residual = series - baseline
        local_mad = rolling_nanmedian(np.abs(residual), half_width)
        sigma = 1.4826 * np.maximum(local_mad, 0.08)
        flat_score[series_i] = residual / sigma
        jump = np.full_like(series, np.nan)
        jump[1:-1] = np.maximum(np.abs(series[1:-1] - series[:-2]), np.abs(series[1:-1] - series[2:]))
        valid = np.isfinite(series) & (flat_counts[series_i] >= min_count)
        flat_mask[series_i] = valid & (
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
            trec[bi] = float(fit_receiver_temperature(tsky[idx], power[idx], weights[idx])["receiver_temp_k"])
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
    if mean_power.ndim == 3:
        plot_module_fits(output, day, centers_s, counts, mean_power, tsky, fit, trec_bootstrap, min_count, outlier_mask)
        return
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

    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    for beam_i, name in enumerate(ppn.BEAM_NAMES):
        good = (counts[beam_i] >= min_count) & np.isfinite(mean_power[beam_i]) & ~outlier_mask[beam_i]
        bad = outlier_mask[beam_i]
        ax.plot(tms[good], 10.0 * np.log10(mean_power[beam_i, good]), ".", ms=3.2, color=colors[beam_i], alpha=0.40)
        if np.any(bad):
            ax.plot(tms[bad], 10.0 * np.log10(mean_power[beam_i, bad]), "x", ms=5.5, mew=1.2, color=colors[beam_i], alpha=0.95)
        model_power = slope * tsky[beam_i] + intercept
        ax.plot(tms, 10.0 * np.log10(np.maximum(model_power, 1.0)), "-", color=colors[beam_i], lw=2.0, label=name)
    ax.set_title(f"PANSY cut-noise power fit, {day}")
    ax.set_ylabel("Raw noise power (dB)")
    ax.set_xlabel("Time (UTC)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    def db_to_kelvin(power_db):
        return np.power(10.0, np.asarray(power_db) / 10.0) / slope

    def kelvin_to_db(temp_k):
        return 10.0 * np.log10(np.maximum(np.asarray(temp_k) * slope, 1.0))

    temp_axis = ax.secondary_yaxis("right", functions=(db_to_kelvin, kelvin_to_db))
    temp_axis.set_ylabel(r"Equivalent $T_{\mathrm{sys}}$ (K)")
    if trec > 0.0:
        ax.axhline(kelvin_to_db(trec), color="0.2", lw=1.0, alpha=0.6)
    ax.text(
        0.015,
        0.96,
        rf"$T_{{rec}}={trec:.0f}^{{+{err_hi:.0f}}}_{{-{err_lo:.0f}}}$ K",
        transform=ax.transAxes,
        ha="left",
        va="top",
    )
    ax.legend(loc="upper right", ncol=5, fontsize=8, frameon=False, handlelength=1.8, markerscale=1.8)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_module_fits(
    output: Path,
    day: str,
    centers_s: np.ndarray,
    counts: np.ndarray,
    mean_power: np.ndarray,
    tsky: np.ndarray,
    fits: list[dict[str, float | np.ndarray]],
    trec_bootstrap: list[np.ndarray],
    min_count: int,
    outlier_mask: np.ndarray,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    colors = ["black", "tab:blue", "tab:orange", "tab:green", "tab:red"]
    tms = centers_s.astype(np.int64).astype("datetime64[s]")
    n_module = mean_power.shape[0]
    ncol = 2
    nrow = int(np.ceil(n_module / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.6, 2.0 * nrow), sharex=True, constrained_layout=True)
    axes = np.ravel(axes)
    for module_i in range(n_module):
        ax = axes[module_i]
        fit = fits[module_i]
        slope = float(fit["slope"])
        trec = float(fit["receiver_temp_k"])
        if len(trec_bootstrap[module_i]):
            lo, hi = np.nanpercentile(trec_bootstrap[module_i], [16.0, 84.0])
            err_lo, err_hi = trec - lo, hi - trec
        else:
            err_lo = err_hi = np.nan
        for beam_i, name in enumerate(ppn.BEAM_NAMES):
            good = (
                (counts[module_i, beam_i] >= min_count)
                & np.isfinite(mean_power[module_i, beam_i])
                & ~outlier_mask[module_i, beam_i]
            )
            if np.any(good):
                ax.plot(
                    tms[good],
                    10.0 * np.log10(mean_power[module_i, beam_i, good]),
                    ".",
                    ms=2.0,
                    color=colors[beam_i],
                    alpha=0.22,
                )
            model_power = slope * (tsky[beam_i] + trec)
            ax.plot(tms, 10.0 * np.log10(np.maximum(model_power, 1.0)), "-", color=colors[beam_i], lw=1.15, label=name)
        ax.text(
            0.02,
            0.94,
            rf"M{module_i}: $T_{{rec}}={trec:.0f}^{{+{err_hi:.0f}}}_{{-{err_lo:.0f}}}$ K",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
        )
        ax.set_ylabel("dB")
        ax.tick_params(labelsize=8)
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    handles, labels = axes[0].get_legend_handles_labels()
    for ax in axes[n_module:]:
        ax.axis("off")
    if len(axes) > n_module:
        axes[n_module].legend(handles, labels, loc="center", ncol=1, frameon=False, markerscale=1.8, fontsize=10)
    else:
        fig.legend(handles, labels, loc="upper center", ncol=5, frameon=False, markerscale=1.8, fontsize=8)
    for ax in axes[-ncol:]:
        if ax.get_visible():
            ax.set_xlabel("Time (UTC)")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_single_module_fit(
    output: Path,
    centers_s: np.ndarray,
    counts: np.ndarray,
    mean_power: np.ndarray,
    tsky: np.ndarray,
    fit: dict[str, float | np.ndarray],
    trec_bootstrap: np.ndarray,
    min_count: int,
    outlier_mask: np.ndarray,
    module_i: int,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    colors = ["black", "tab:blue", "tab:orange", "tab:green", "tab:red"]
    tms = centers_s.astype(np.int64).astype("datetime64[s]")
    slope = float(fit["slope"])
    trec = float(fit["receiver_temp_k"])
    if len(trec_bootstrap):
        lo, hi = np.nanpercentile(trec_bootstrap, [16.0, 84.0])
        err_lo, err_hi = trec - lo, hi - trec
    else:
        err_lo = err_hi = np.nan

    fig, ax = plt.subplots(figsize=(7.2, 4.1), constrained_layout=True)
    for beam_i, name in enumerate(ppn.BEAM_NAMES):
        good = (
            (counts[module_i, beam_i] >= min_count)
            & np.isfinite(mean_power[module_i, beam_i])
            & ~outlier_mask[module_i, beam_i]
        )
        ax.plot(
            tms[good],
            10.0 * np.log10(mean_power[module_i, beam_i, good]),
            ".",
            ms=2.8,
            color=colors[beam_i],
            alpha=0.32,
        )
        model_power = slope * (tsky[beam_i] + trec)
        ax.plot(tms, 10.0 * np.log10(np.maximum(model_power, 1.0)), "-", color=colors[beam_i], lw=1.9, label=name)
    ax.set_ylabel("Raw noise power (dB)")
    ax.set_xlabel("Time (UTC)")
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

    def db_to_kelvin(power_db):
        return np.power(10.0, np.asarray(power_db) / 10.0) / slope

    def kelvin_to_db(temp_k):
        return 10.0 * np.log10(np.maximum(np.asarray(temp_k) * slope, 1.0))

    temp_axis = ax.secondary_yaxis("right", functions=(db_to_kelvin, kelvin_to_db))
    temp_axis.set_ylabel(r"Equivalent $T_{\mathrm{sys}}$ (K)")
    ax.axhline(kelvin_to_db(trec), color="0.25", lw=1.0, alpha=0.55)
    ax.text(
        0.015,
        0.96,
        rf"Module {module_i}: $T_{{rec}}={trec:.0f}^{{+{err_hi:.0f}}}_{{-{err_lo:.0f}}}$ K",
        transform=ax.transAxes,
        ha="left",
        va="top",
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=5, fontsize=8, frameon=False, markerscale=2.0, handlelength=1.8)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def fit_all_modules(
    mean_power: np.ndarray,
    counts: np.ndarray,
    tsky: np.ndarray,
    min_count: int,
    outlier_mask: np.ndarray,
    n_bootstrap: int,
    seed: int,
) -> tuple[list[dict[str, float | np.ndarray]], list[np.ndarray]]:
    fits = []
    bootstraps = []
    minute_grid = np.broadcast_to(np.arange(1440, dtype=np.int64), mean_power.shape[1:])
    for module_i in range(mean_power.shape[0]):
        keep = (counts[module_i] >= min_count) & np.isfinite(mean_power[module_i]) & ~outlier_mask[module_i]
        fit = fit_receiver_temperature(tsky[keep], mean_power[module_i][keep], counts[module_i][keep])
        boot = bootstrap_trec(
            minute_grid[keep],
            tsky[keep],
            mean_power[module_i][keep],
            counts[module_i][keep],
            n_bootstrap=n_bootstrap,
            seed=seed + module_i,
        )
        fits.append(fit)
        bootstraps.append(boot)
    return fits, bootstraps


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--day", default="2025-05-10")
    parser.add_argument("--metadata", type=Path, default=Path("/mnt/data/juha/pansy/metadata/cut"))
    parser.add_argument("--measurements", type=Path, default=None, help="Existing cut noise measurement HDF5 sidecar.")
    parser.add_argument("--output", type=Path, default=Path("figs/cut_noise_pygdsm_fit_2025-05-10.png"))
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--chunk-seconds", type=int, default=60)
    parser.add_argument("--guard-samples", type=int, default=25)
    parser.add_argument("--min-count", type=int, default=20)
    parser.add_argument("--sky-step-minutes", type=float, default=10.0)
    parser.add_argument("--n-az", type=int, default=180)
    parser.add_argument("--n-el", type=int, default=60)
    parser.add_argument("--min-elevation-deg", type=float, default=0.5)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20250708)
    parser.add_argument("--paper-module", type=int, default=None, help="Write a one-panel paper plot for this receiver module index.")
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
    model_tsky = receive_sky_temperature_matrix(
        model_times,
        n_az=args.n_az,
        n_el=args.n_el,
        min_elevation_deg=args.min_elevation_deg,
    )
    tsky = interpolate_sky_model(model_times, model_tsky, centers_s)
    counts = np.asarray(binned["counts"], dtype=np.float64)
    mean_power = np.asarray(binned["mean_power"], dtype=np.float64)
    power_db = 10.0 * np.log10(np.maximum(mean_power, 1.0))
    outlier_mask, robust_score = robust_jump_outlier_mask(power_db, counts, args.min_count)
    fits, trec_bootstrap = fit_all_modules(
        mean_power,
        counts,
        tsky,
        args.min_count,
        outlier_mask,
        n_bootstrap=args.bootstrap,
        seed=args.seed,
    )
    if args.paper_module is None:
        plot_fit(args.output, args.day, centers_s, counts, mean_power, tsky, fits, trec_bootstrap, args.min_count, outlier_mask)
    else:
        plot_single_module_fit(
            args.output,
            centers_s,
            counts,
            mean_power,
            tsky,
            fits[args.paper_module],
            trec_bootstrap[args.paper_module],
            args.min_count,
            outlier_mask,
            args.paper_module,
        )
    print(f"metadata {args.metadata}")
    print(f"output {args.output}")
    print(f"samples {len(day_data['noise_power'])}")
    print(f"modules {mean_power.shape[0]}")
    print(f"outlier_bins {np.count_nonzero(outlier_mask)}")
    slopes = np.asarray([float(fit["slope"]) for fit in fits], dtype=np.float64)
    ref_slope = slopes[0]
    power_equalization = ref_slope / slopes
    voltage_equalization = np.sqrt(power_equalization)
    for module_i, fit in enumerate(fits):
        keep_i = (counts[module_i] >= args.min_count) & np.isfinite(mean_power[module_i]) & ~outlier_mask[module_i]
        print(f"module {module_i}")
        print(f"  fit_bins {np.count_nonzero(keep_i)}")
        print(f"  slope_power_per_k {float(fit['slope']):.9g}")
        print(f"  receiver_temp_k {float(fit['receiver_temp_k']):.6g}")
        print(f"  power_equalization_to_module0 {power_equalization[module_i]:.9g}")
        print(f"  voltage_equalization_to_module0 {voltage_equalization[module_i]:.9g}")
        if len(trec_bootstrap[module_i]):
            lo, med, hi = np.nanpercentile(trec_bootstrap[module_i], [16.0, 50.0, 84.0])
            print(f"  receiver_temp_bootstrap_p16_p50_p84_k {lo:.6g} {med:.6g} {hi:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
