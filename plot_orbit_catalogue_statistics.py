#!/usr/bin/env python3
"""Make paper-ready catalogue statistics plots from compact orbit metadata."""

from __future__ import annotations

import argparse
import datetime as dt
import warnings
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from erfa import ErfaWarning
from astropy.coordinates import get_sun
from astropy.time import Time
from matplotlib.colors import LogNorm

import orbit_metadata_table as omt


DEFAULT_MIN_SAMPLE_IDX = int(dt.datetime(2025, 1, 1, tzinfo=dt.timezone.utc).timestamp() * 1_000_000)
DEFAULT_MAX_SAMPLE_IDX = int(dt.datetime(2027, 1, 1, tzinfo=dt.timezone.utc).timestamp() * 1_000_000)
HEIGHT_VELOCITY_HEIGHT_MIN_KM = 20.0
HEIGHT_VELOCITY_HEIGHT_MAX_KM = 180.0
HEIGHT_VELOCITY_SPEED_MIN_KM_S = 0.0
HEIGHT_VELOCITY_SPEED_MAX_KM_S = 80.0
LOW_HEIGHT_AUDIT_DTYPE = np.dtype(
    [
        ("sample_idx", "<i8"),
        ("n_selected_measurements", "<i4"),
        ("n_below_threshold", "<i4"),
        ("frac_below_threshold", "<f4"),
        ("min_height_km", "<f4"),
        ("p05_height_km", "<f4"),
        ("median_height_km", "<f4"),
        ("max_height_km", "<f4"),
        ("v_g_km_s", "<f4"),
        ("combined_score", "<f4"),
    ]
)


def wrap360(deg: np.ndarray) -> np.ndarray:
    return np.asarray(deg, dtype=np.float64) % 360.0


def sun_ecliptic_longitude_deg(unix_time: np.ndarray) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        time = Time(np.asarray(unix_time, dtype=np.float64), format="unix", scale="utc")
        sun = get_sun(time).transform_to("geocentricmeanecliptic")
    return wrap360(sun.lon.deg)


def date_labels_for_solar_longitude(year: int, ticks_deg: np.ndarray) -> list[str]:
    start = dt.datetime(year, 1, 1, tzinfo=dt.timezone.utc)
    stop = dt.datetime(year + 1, 1, 1, tzinfo=dt.timezone.utc)
    step_seconds = 6 * 3600
    times = np.arange(int(start.timestamp()), int(stop.timestamp()) + step_seconds, step_seconds, dtype=np.float64)
    lon = np.rad2deg(np.unwrap(np.deg2rad(sun_ecliptic_longitude_deg(times))))
    labels = []
    for tick in np.asarray(ticks_deg, dtype=np.float64):
        target = float(tick)
        while target < lon[0]:
            target += 360.0
        while target > lon[-1]:
            target -= 360.0
        if target < lon[0] or target > lon[-1]:
            labels.append("")
            continue
        unix = float(np.interp(target, lon, times))
        date = dt.datetime.fromtimestamp(unix, tz=dt.timezone.utc)
        labels.append(date.strftime("%b %d"))
    return labels


def add_calendar_axes(ax, ticks_deg: np.ndarray) -> None:
    ax.set_xticks(ticks_deg)
    for year, offset, pad in ((2025, 0, 2), (2026, 34, 2)):
        top = ax.twiny()
        top.set_xlim(ax.get_xlim())
        top.set_xticks(ticks_deg)
        top.set_xticklabels(date_labels_for_solar_longitude(year, ticks_deg), fontsize=8)
        top.set_xlabel(f"{year} UTC date", labelpad=pad)
        top.tick_params(axis="x", direction="out", pad=pad)
        top.spines["top"].set_position(("outward", offset))


def month_tick_positions_deg(year: int) -> tuple[np.ndarray, list[str]]:
    dates = [dt.datetime(year, month, 1, 12, tzinfo=dt.timezone.utc) for month in range(1, 13)]
    unix = np.asarray([date.timestamp() for date in dates], dtype=np.float64)
    labels = [date.strftime("%b") for date in dates]
    return sun_ecliptic_longitude_deg(unix), labels


def add_month_axis(ax, year: int) -> None:
    positions, labels = month_tick_positions_deg(year)
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    top.set_xticks(positions)
    top.set_xticklabels(labels, fontsize=8)
    top.set_xlabel("")
    top.tick_params(axis="x", direction="out", pad=2)


def coerce_structured(arr: np.ndarray, dtype: np.dtype) -> np.ndarray:
    arr = np.asarray(arr)
    out = np.zeros(arr.shape, dtype=dtype)
    if arr.dtype.names is None:
        return out
    for name in dtype.names:
        if name in arr.dtype.names:
            out[name] = arr[name]
        elif dtype[name].kind == "f":
            out[name] = np.nan
    return out


def iter_orbit_tables(orbit_metadata_dir: Path, include_paths: bool):
    for path in sorted(orbit_metadata_dir.glob("**/orbit@*.h5")):
        try:
            with h5py.File(path, "r") as h:
                if "events" not in h:
                    continue
                events = coerce_structured(h["events"][()], omt.EVENT_DTYPE)
                paths = coerce_structured(h["paths"][()], omt.PATH_DTYPE) if include_paths and "paths" in h else None
                yield path, events, paths
        except OSError:
            continue


def collect_events_and_paths(orbit_metadata_dir: Path, include_paths: bool) -> tuple[np.ndarray, np.ndarray, int]:
    event_chunks = []
    path_chunks = []
    files_read = 0
    for _path, events, paths in iter_orbit_tables(orbit_metadata_dir, include_paths=include_paths):
        event_chunks.append(events)
        if paths is not None and len(paths):
            path_chunks.append(paths)
        files_read += 1
    events = np.concatenate(event_chunks) if event_chunks else np.zeros(0, dtype=omt.EVENT_DTYPE)
    paths = np.concatenate(path_chunks) if path_chunks else np.zeros(0, dtype=omt.PATH_DTYPE)
    if len(events) and "sample_idx" in events.dtype.names:
        events = events[np.argsort(events["sample_idx"])]
        _, unique_i = np.unique(events["sample_idx"], return_index=True)
        if len(unique_i) != len(events):
            events = events[np.sort(unique_i)]
    if len(paths) and "sample_idx" in paths.dtype.names:
        paths = paths[np.argsort(paths["sample_idx"])]
    return events, paths, files_read


def finite_field(events: np.ndarray, name: str) -> np.ndarray:
    if events.dtype.names is None or name not in events.dtype.names:
        return np.full(len(events), np.nan, dtype=np.float32)
    values = np.asarray(events[name], dtype=np.float64)
    out = np.full(values.shape, np.nan, dtype=np.float32)
    finite = np.isfinite(values)
    finite &= np.abs(values) <= np.finfo(np.float32).max
    out[finite] = values[finite].astype(np.float32)
    return out


def quality_mask(
    events: np.ndarray,
    max_radiant_sigma_deg: float | None,
    min_sample_idx: int,
    max_sample_idx: int,
) -> np.ndarray:
    good = np.ones(len(events), dtype=bool)
    if events.dtype.names is None or "sample_idx" not in events.dtype.names:
        return np.zeros(len(events), dtype=bool)
    sample_idx = np.asarray(events["sample_idx"], dtype=np.int64)
    good &= (sample_idx >= int(min_sample_idx)) & (sample_idx < int(max_sample_idx))
    for name in ("initial_detection_height_km", "v_g_km_s", "radiant_sun_ecliptic_lon_deg"):
        good &= np.isfinite(finite_field(events, name))
    good &= finite_field(events, "v_g_km_s") > 0.0
    good &= finite_field(events, "initial_detection_height_km") > 0.0
    if max_radiant_sigma_deg is not None and "fit_parameter_covariance" in events.dtype.names:
        cov = np.asarray(events["fit_parameter_covariance"], dtype=np.float64)
        vel_var = np.trace(cov[:, 3:6, 3:6], axis1=1, axis2=2)
        vel_sigma = np.sqrt(np.where(vel_var >= 0.0, vel_var, np.nan))
        radiant_sigma = np.rad2deg(np.arctan2(vel_sigma, np.maximum(finite_field(events, "v_g_km_s") * 1e3, 1.0)))
        good &= np.isfinite(radiant_sigma) & (radiant_sigma <= float(max_radiant_sigma_deg))
    return good


def fit_quality_mask(
    events: np.ndarray,
    max_radiant_sigma_deg: float | None,
    min_sample_idx: int,
    max_sample_idx: int,
    max_combined_score: float | None,
    max_frac_e_gt_1: float | None,
    min_uncertainty_samples: int,
    max_initial_state_position_sigma_m: float | None,
    max_initial_state_radiant_angle_sigma_deg: float | None,
) -> np.ndarray:
    good = quality_mask(events, max_radiant_sigma_deg, min_sample_idx, max_sample_idx)
    if events.dtype.names is None:
        return good
    if "orbit_solution_type" in events.dtype.names:
        good &= np.asarray(events["orbit_solution_type"]) == b"dasst_winning_alias"
    if max_combined_score is not None and "combined_score" in events.dtype.names:
        combined_score = finite_field(events, "combined_score")
        good &= np.isfinite(combined_score) & (combined_score <= float(max_combined_score))
    if max_frac_e_gt_1 is not None and "frac_e_gt_1" in events.dtype.names:
        frac_e_gt_1 = finite_field(events, "frac_e_gt_1")
        good &= np.isfinite(frac_e_gt_1) & (frac_e_gt_1 <= float(max_frac_e_gt_1))
    if "n_uncertainty_samples" in events.dtype.names:
        good &= np.asarray(events["n_uncertainty_samples"], dtype=np.int64) >= int(min_uncertainty_samples)
    if (
        max_initial_state_position_sigma_m is not None
        or max_initial_state_radiant_angle_sigma_deg is not None
    ) and "fit_parameter_covariance" in events.dtype.names:
        cov = np.asarray(events["fit_parameter_covariance"], dtype=np.float64)
        pos_var = np.trace(cov[:, :3, :3], axis1=1, axis2=2)
        vel_var = np.trace(cov[:, 3:6, 3:6], axis1=1, axis2=2)
        pos_sigma = np.sqrt(np.where(pos_var >= 0.0, pos_var, np.nan))
        vel_sigma = np.sqrt(np.where(vel_var >= 0.0, vel_var, np.nan))
        angle_sigma = np.rad2deg(
            np.arctan2(vel_sigma, np.maximum(finite_field(events, "v_g_km_s") * 1e3, 1.0))
        )
        if max_initial_state_position_sigma_m is not None:
            good &= np.isfinite(pos_sigma) & (pos_sigma <= float(max_initial_state_position_sigma_m))
        if max_initial_state_radiant_angle_sigma_deg is not None:
            good &= np.isfinite(angle_sigma) & (angle_sigma <= float(max_initial_state_radiant_angle_sigma_deg))
    return good


def catalogue_arrays(
    events: np.ndarray,
    max_radiant_sigma_deg: float | None,
    min_sample_idx: int,
    max_sample_idx: int,
) -> dict[str, np.ndarray]:
    keep = quality_mask(events, max_radiant_sigma_deg, min_sample_idx, max_sample_idx)
    kept = events[keep]
    out = {
        "sample_idx": np.asarray(kept["sample_idx"], dtype=np.int64),
        "solar_longitude_deg": wrap360(finite_field(kept, "radiant_sun_ecliptic_lon_deg")).astype(np.float32),
        "initial_detection_height_km": finite_field(kept, "initial_detection_height_km").astype(np.float32),
        "v_g_km_s": finite_field(kept, "v_g_km_s").astype(np.float32),
    }
    if "combined_score" in kept.dtype.names:
        out["combined_score"] = finite_field(kept, "combined_score").astype(np.float32)
    return out


def fit_catalogue_arrays(
    events: np.ndarray,
    max_radiant_sigma_deg: float | None,
    min_sample_idx: int,
    max_sample_idx: int,
    max_combined_score: float | None,
    max_frac_e_gt_1: float | None,
    min_uncertainty_samples: int,
    max_initial_state_position_sigma_m: float | None,
    max_initial_state_radiant_angle_sigma_deg: float | None,
) -> dict[str, np.ndarray]:
    keep = fit_quality_mask(
        events,
        max_radiant_sigma_deg,
        min_sample_idx,
        max_sample_idx,
        max_combined_score,
        max_frac_e_gt_1,
        min_uncertainty_samples,
        max_initial_state_position_sigma_m,
        max_initial_state_radiant_angle_sigma_deg,
    )
    kept = events[keep]
    return {
        "fit_sample_idx": np.asarray(kept["sample_idx"], dtype=np.int64),
        "fit_initial_detection_height_km": finite_field(kept, "initial_detection_height_km").astype(np.float32),
        "fit_v_g_km_s": finite_field(kept, "v_g_km_s").astype(np.float32),
    }


def empty_measurement_arrays() -> dict[str, np.ndarray]:
    return {
        "measurement_event_sample_idx": np.zeros(0, dtype=np.int64),
        "measurement_height_km": np.zeros(0, dtype=np.float32),
        "measurement_event_v_g_km_s": np.zeros(0, dtype=np.float32),
    }


def measurement_height_velocity_arrays(
    events: np.ndarray,
    paths: np.ndarray,
    max_radiant_sigma_deg: float | None,
    min_sample_idx: int,
    max_sample_idx: int,
    max_combined_score: float | None,
    max_frac_e_gt_1: float | None,
    min_uncertainty_samples: int,
    max_initial_state_position_sigma_m: float | None,
    max_initial_state_radiant_angle_sigma_deg: float | None,
) -> dict[str, np.ndarray]:
    if len(events) == 0 or len(paths) == 0 or paths.dtype.names is None:
        return empty_measurement_arrays()
    keep = fit_quality_mask(
        events,
        max_radiant_sigma_deg,
        min_sample_idx,
        max_sample_idx,
        max_combined_score,
        max_frac_e_gt_1,
        min_uncertainty_samples,
        max_initial_state_position_sigma_m,
        max_initial_state_radiant_angle_sigma_deg,
    )
    kept = events[keep]
    sample_idx = np.asarray(kept["sample_idx"], dtype=np.int64)
    speeds = finite_field(kept, "v_g_km_s")
    order = np.argsort(sample_idx)
    sample_idx = sample_idx[order]
    speeds = speeds[order]

    path_sample_idx = np.asarray(paths["sample_idx"], dtype=np.int64)
    match = np.searchsorted(sample_idx, path_sample_idx)
    in_range = match < len(sample_idx)
    matched = np.zeros(len(paths), dtype=bool)
    matched[in_range] = sample_idx[match[in_range]] == path_sample_idx[in_range]
    if "position_enu_km" not in paths.dtype.names:
        matched &= False
    height = np.asarray(paths["position_enu_km"], dtype=np.float32)[:, 2] if "position_enu_km" in paths.dtype.names else np.zeros(len(paths), dtype=np.float32)
    good = matched & np.isfinite(height)
    if "selection_keep" in paths.dtype.names:
        good &= np.asarray(paths["selection_keep"], dtype=bool)
    else:
        good &= False
    good &= height > 0.0
    out_speed = np.full(len(paths), np.nan, dtype=np.float32)
    out_speed[matched] = speeds[match[matched]]
    good &= np.isfinite(out_speed)
    good &= height >= HEIGHT_VELOCITY_HEIGHT_MIN_KM
    good &= height < HEIGHT_VELOCITY_HEIGHT_MAX_KM
    good &= out_speed >= HEIGHT_VELOCITY_SPEED_MIN_KM_S
    good &= out_speed < HEIGHT_VELOCITY_SPEED_MAX_KM_S
    return {
        "measurement_event_sample_idx": path_sample_idx[good].astype(np.int64),
        "measurement_height_km": height[good].astype(np.float32),
        "measurement_event_v_g_km_s": out_speed[good].astype(np.float32),
    }


def measurement_low_height_diagnostics(
    events: np.ndarray,
    paths: np.ndarray,
    max_radiant_sigma_deg: float | None,
    min_sample_idx: int,
    max_sample_idx: int,
    max_combined_score: float | None,
    max_frac_e_gt_1: float | None,
    min_uncertainty_samples: int,
    max_initial_state_position_sigma_m: float | None,
    max_initial_state_radiant_angle_sigma_deg: float | None,
    low_height_threshold_km: float,
) -> np.ndarray:
    if len(events) == 0 or len(paths) == 0 or paths.dtype.names is None:
        return np.zeros(0, dtype=LOW_HEIGHT_AUDIT_DTYPE)
    keep = fit_quality_mask(
        events,
        max_radiant_sigma_deg,
        min_sample_idx,
        max_sample_idx,
        max_combined_score,
        max_frac_e_gt_1,
        min_uncertainty_samples,
        max_initial_state_position_sigma_m,
        max_initial_state_radiant_angle_sigma_deg,
    )
    kept = events[keep]
    if len(kept) == 0 or "position_enu_km" not in paths.dtype.names:
        return np.zeros(0, dtype=LOW_HEIGHT_AUDIT_DTYPE)

    sample_idx = np.asarray(kept["sample_idx"], dtype=np.int64)
    speeds = finite_field(kept, "v_g_km_s")
    combined_score = finite_field(kept, "combined_score") if "combined_score" in kept.dtype.names else np.full(len(kept), np.nan, dtype=np.float32)
    order = np.argsort(sample_idx)
    sample_idx = sample_idx[order]
    speeds = speeds[order]
    combined_score = combined_score[order]

    path_sample_idx = np.asarray(paths["sample_idx"], dtype=np.int64)
    match = np.searchsorted(sample_idx, path_sample_idx)
    in_range = match < len(sample_idx)
    matched = np.zeros(len(paths), dtype=bool)
    matched[in_range] = sample_idx[match[in_range]] == path_sample_idx[in_range]
    height = np.asarray(paths["position_enu_km"], dtype=np.float32)[:, 2]
    path_speed = np.full(len(paths), np.nan, dtype=np.float32)
    path_speed[matched] = speeds[match[matched]]

    good = matched & np.isfinite(height) & np.isfinite(path_speed)
    if "selection_keep" in paths.dtype.names:
        good &= np.asarray(paths["selection_keep"], dtype=bool)
    else:
        good &= False
    good &= height > 0.0
    good &= height < HEIGHT_VELOCITY_HEIGHT_MAX_KM
    good &= path_speed >= HEIGHT_VELOCITY_SPEED_MIN_KM_S
    good &= path_speed < HEIGHT_VELOCITY_SPEED_MAX_KM_S
    if not np.any(good):
        return np.zeros(0, dtype=LOW_HEIGHT_AUDIT_DTYPE)

    good_sample = path_sample_idx[good]
    good_height = height[good]
    sort_i = np.argsort(good_sample)
    good_sample = good_sample[sort_i]
    good_height = good_height[sort_i]
    unique_sample, start, n_selected = np.unique(good_sample, return_index=True, return_counts=True)
    n_below = np.add.reduceat((good_height < float(low_height_threshold_km)).astype(np.int32), start)
    suspicious = n_below > 0
    if not np.any(suspicious):
        return np.zeros(0, dtype=LOW_HEIGHT_AUDIT_DTYPE)

    unique_sample = unique_sample[suspicious]
    start = start[suspicious]
    n_selected = n_selected[suspicious].astype(np.int32)
    n_below = n_below[suspicious].astype(np.int32)
    event_match = np.searchsorted(sample_idx, unique_sample)

    out = np.zeros(len(unique_sample), dtype=LOW_HEIGHT_AUDIT_DTYPE)
    out["sample_idx"] = unique_sample.astype(np.int64)
    out["n_selected_measurements"] = n_selected
    out["n_below_threshold"] = n_below
    out["frac_below_threshold"] = (n_below / np.maximum(n_selected, 1)).astype(np.float32)
    out["v_g_km_s"] = speeds[event_match].astype(np.float32)
    out["combined_score"] = combined_score[event_match].astype(np.float32)
    for i, (s0, n) in enumerate(zip(start, n_selected)):
        h = good_height[s0 : s0 + n]
        out["min_height_km"][i] = np.min(h)
        out["p05_height_km"][i] = np.percentile(h, 5)
        out["median_height_km"][i] = np.median(h)
        out["max_height_km"][i] = np.max(h)
    return out


def histogram_solar_longitude(solar_longitude_deg: np.ndarray, bin_width_deg: float) -> tuple[np.ndarray, np.ndarray]:
    edges = np.arange(0.0, 360.0 + float(bin_width_deg), float(bin_width_deg), dtype=np.float32)
    counts, edges = np.histogram(solar_longitude_deg, bins=edges)
    return counts.astype(np.int32), edges.astype(np.float32)


def histogram_solar_longitude_by_year(
    sample_idx: np.ndarray,
    solar_longitude_deg: np.ndarray,
    bin_width_deg: float,
    years: tuple[int, ...] = (2025, 2026),
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    edges = np.arange(0.0, 360.0 + float(bin_width_deg), float(bin_width_deg), dtype=np.float32)
    epoch_s = np.asarray(sample_idx, dtype=np.float64) / 1_000_000.0
    event_year = np.asarray([dt.datetime.fromtimestamp(float(t), tz=dt.timezone.utc).year for t in epoch_s], dtype=np.int16)
    counts = []
    solar_longitude = np.asarray(solar_longitude_deg, dtype=np.float64)
    for year in years:
        hist, _ = np.histogram(solar_longitude[event_year == int(year)], bins=edges)
        counts.append(hist.astype(np.int32))
    return np.asarray(years, dtype=np.int16), np.asarray(counts, dtype=np.int32), edges.astype(np.float32)


def count_density_from_counts(
    counts_by_year: np.ndarray,
    edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    counts = np.asarray(counts_by_year, dtype=np.float64)
    bin_width_deg = np.diff(np.asarray(edges, dtype=np.float64))
    density_by_year = (counts / bin_width_deg[np.newaxis, :]).astype(np.float32)
    all_density = (np.sum(counts, axis=0) / bin_width_deg).astype(np.float32)
    return density_by_year, all_density


def histogram_height_velocity(height_km: np.ndarray, speed_km_s: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    height_edges = np.arange(HEIGHT_VELOCITY_HEIGHT_MIN_KM, HEIGHT_VELOCITY_HEIGHT_MAX_KM + 1.0, 1.0, dtype=np.float32)
    speed_edges = np.arange(HEIGHT_VELOCITY_SPEED_MIN_KM_S, HEIGHT_VELOCITY_SPEED_MAX_KM_S + 1.0, 1.0, dtype=np.float32)
    hist, h_edges, v_edges = np.histogram2d(height_km, speed_km_s, bins=[height_edges, speed_edges])
    return hist.astype(np.int32), h_edges.astype(np.float32), v_edges.astype(np.float32)


def write_statistics_h5(
    path: Path,
    arrays: dict[str, np.ndarray],
    fit_arrays: dict[str, np.ndarray],
    measurement_count: int,
    solar_counts,
    solar_years,
    solar_year_counts,
    solar_year_count_density,
    solar_all_count_density,
    solar_edges,
    hv_counts,
    measurement_hv_counts,
    low_height_diagnostics,
    height_edges,
    speed_edges,
    files_read: int,
    low_height_threshold_km: float,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h:
        h.attrs["script"] = Path(__file__).name
        h.attrs["source"] = "compact orbit metadata events tables"
        h.attrs["height_velocity_source"] = "orbit metadata paths.position_enu_km[:,2] joined to events.v_g_km_s by sample_idx"
        h.attrs["event_height_velocity_source"] = "orbit metadata events.initial_detection_height_km and events.v_g_km_s"
        h.attrs["no_detection_metadata_used"] = True
        h.attrs["files_read"] = int(files_read)
        h.attrs["event_count"] = int(len(arrays["sample_idx"]))
        h.attrs["fit_event_count"] = int(len(fit_arrays["fit_sample_idx"]))
        h.attrs["measurement_count"] = int(measurement_count)
        h.attrs["low_height_audit_threshold_km"] = float(low_height_threshold_km)
        h.attrs["low_height_audit_event_count"] = int(len(low_height_diagnostics))
        h.attrs["solar_longitude_count_density_unit"] = "count per degree"
        h.attrs["height_velocity_quality_filter"] = (
            "static radiant monitor high-quality gate: valid DASST winning alias, "
            "combined_score, n_uncertainty_samples, initial_state_position_sigma_m, "
            "and initial_state_radiant_angle_sigma_deg"
        )
        for name, values in arrays.items():
            h.create_dataset(name, data=values, compression="gzip", shuffle=True)
        for name, values in fit_arrays.items():
            h.create_dataset(name, data=values, compression="gzip", shuffle=True)
        h.create_dataset("solar_longitude_count", data=solar_counts, compression="gzip", shuffle=True)
        h.create_dataset("solar_longitude_year", data=solar_years)
        h.create_dataset("solar_longitude_count_by_year", data=solar_year_counts, compression="gzip", shuffle=True)
        h.create_dataset("solar_longitude_count_per_degree_by_year", data=solar_year_count_density, compression="gzip", shuffle=True)
        h.create_dataset("solar_longitude_count_per_degree", data=solar_all_count_density, compression="gzip", shuffle=True)
        h.create_dataset("solar_longitude_edges_deg", data=solar_edges)
        h.create_dataset("height_velocity_count", data=hv_counts, compression="gzip", shuffle=True)
        h.create_dataset("measurement_height_velocity_count", data=measurement_hv_counts, compression="gzip", shuffle=True)
        h.create_dataset("measurement_low_height_diagnostics", data=low_height_diagnostics, compression="gzip", shuffle=True)
        h.create_dataset("height_edges_km", data=height_edges)
        h.create_dataset("speed_edges_km_s", data=speed_edges)


def plot_solar_counts(
    path: Path,
    years: np.ndarray,
    count_density_by_year: np.ndarray,
    all_count_density: np.ndarray,
    edges: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    fig, ax = plt.subplots(figsize=(8.0, 3.6), constrained_layout=True)
    ax.step(centers, all_count_density, where="mid", color="black", linewidth=1.55, label="All")
    colors = ("#2f5f8f", "#c04b37", "#4f7f3f")
    for year, count_density, color in zip(years, count_density_by_year, colors):
        ax.step(centers, count_density, where="mid", color=color, linewidth=1.2, label=str(int(year)))
    ax.set_xlim(0, 360)
    finite_density = np.concatenate((np.ravel(all_count_density), np.ravel(count_density_by_year)))
    finite_density = finite_density[np.isfinite(finite_density)]
    ax.set_ylim(0, max(1.0, float(np.max(finite_density)) if len(finite_density) else 1.0) * 1.08)
    ax.set_xlabel(r"Solar longitude, $\lambda_\odot$ (deg)")
    ax.set_ylabel(r"Catalogue count (deg$^{-1}$)")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, loc="upper right", ncol=3)
    add_month_axis(ax, 2025)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_height_velocity(
    path: Path,
    hv_counts: np.ndarray,
    height_edges: np.ndarray,
    speed_edges: np.ndarray,
    ylabel: str,
    colorbar_label: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plot_count = np.asarray(hv_counts, dtype=np.float64)
    plot_count[plot_count <= 0.0] = np.nan
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("black")
    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    mesh = ax.pcolormesh(
        speed_edges,
        height_edges,
        plot_count,
        shading="auto",
        cmap=cmap,
        norm=LogNorm(vmin=1.0, vmax=max(1.0, float(np.nanmax(plot_count)) if np.any(np.isfinite(plot_count)) else 1.0)),
    )
    ax.set_xlabel(r"Geocentric velocity, $v_g$ (km s$^{-1}$)")
    ax.set_ylabel(ylabel)
    ax.set_xlim(float(speed_edges[0]), float(speed_edges[-1]))
    ax.set_ylim(float(height_edges[0]), float(height_edges[-1]))
    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label(colorbar_label)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-metadata-dir", type=Path, default=Path("data/metadata/orbit"))
    parser.add_argument("--output-h5", type=Path, default=Path("paper_orbit_catalogue_statistics.h5"))
    parser.add_argument("--solar-output", type=Path, default=Path("paper_counts_solar_longitude.png"))
    parser.add_argument("--height-velocity-output", type=Path, default=Path("paper_height_velocity.png"))
    parser.add_argument("--measurement-height-velocity-output", type=Path, default=Path("paper_measurement_height_velocity.png"))
    parser.add_argument("--solar-bin-width-deg", type=float, default=1.0)
    parser.add_argument("--max-radiant-sigma-deg", type=float, default=None)
    parser.add_argument("--min-sample-idx", type=int, default=DEFAULT_MIN_SAMPLE_IDX)
    parser.add_argument("--max-sample-idx", type=int, default=DEFAULT_MAX_SAMPLE_IDX)
    parser.add_argument("--include-measurements", action="store_true")
    parser.add_argument("--height-max-combined-score", type=float, default=1.5)
    parser.add_argument("--height-max-frac-e-gt-1", type=float, default=0.5)
    parser.add_argument("--height-min-uncertainty-samples", type=int, default=3)
    parser.add_argument("--height-max-initial-state-position-sigma-m", type=float, default=1000.0)
    parser.add_argument("--height-max-initial-state-radiant-angle-sigma-deg", type=float, default=3.0)
    parser.add_argument("--low-height-audit-threshold-km", type=float, default=60.0)
    args = parser.parse_args()

    events, paths, files_read = collect_events_and_paths(args.orbit_metadata_dir, include_paths=args.include_measurements)
    arrays = catalogue_arrays(events, args.max_radiant_sigma_deg, args.min_sample_idx, args.max_sample_idx)
    fit_arrays = fit_catalogue_arrays(
        events,
        args.max_radiant_sigma_deg,
        args.min_sample_idx,
        args.max_sample_idx,
        args.height_max_combined_score,
        args.height_max_frac_e_gt_1,
        args.height_min_uncertainty_samples,
        args.height_max_initial_state_position_sigma_m,
        args.height_max_initial_state_radiant_angle_sigma_deg,
    )
    measurement_arrays = (
        measurement_height_velocity_arrays(
            events,
            paths,
            args.max_radiant_sigma_deg,
            args.min_sample_idx,
            args.max_sample_idx,
            args.height_max_combined_score,
            args.height_max_frac_e_gt_1,
            args.height_min_uncertainty_samples,
            args.height_max_initial_state_position_sigma_m,
            args.height_max_initial_state_radiant_angle_sigma_deg,
        )
        if args.include_measurements
        else empty_measurement_arrays()
    )
    low_height_diagnostics = (
        measurement_low_height_diagnostics(
            events,
            paths,
            args.max_radiant_sigma_deg,
            args.min_sample_idx,
            args.max_sample_idx,
            args.height_max_combined_score,
            args.height_max_frac_e_gt_1,
            args.height_min_uncertainty_samples,
            args.height_max_initial_state_position_sigma_m,
            args.height_max_initial_state_radiant_angle_sigma_deg,
            args.low_height_audit_threshold_km,
        )
        if args.include_measurements
        else np.zeros(0, dtype=LOW_HEIGHT_AUDIT_DTYPE)
    )
    solar_counts, solar_edges = histogram_solar_longitude(arrays["solar_longitude_deg"], args.solar_bin_width_deg)
    solar_years, solar_year_counts, solar_edges = histogram_solar_longitude_by_year(
        arrays["sample_idx"],
        arrays["solar_longitude_deg"],
        args.solar_bin_width_deg,
    )
    solar_year_count_density, solar_all_count_density = count_density_from_counts(solar_year_counts, solar_edges)
    hv_counts, height_edges, speed_edges = histogram_height_velocity(
        fit_arrays["fit_initial_detection_height_km"],
        fit_arrays["fit_v_g_km_s"],
    )
    measurement_hv_counts, _, _ = histogram_height_velocity(
        measurement_arrays["measurement_height_km"],
        measurement_arrays["measurement_event_v_g_km_s"],
    )
    write_statistics_h5(
        args.output_h5,
        arrays,
        fit_arrays,
        len(measurement_arrays["measurement_height_km"]),
        solar_counts,
        solar_years,
        solar_year_counts,
        solar_year_count_density,
        solar_all_count_density,
        solar_edges,
        hv_counts,
        measurement_hv_counts,
        low_height_diagnostics,
        height_edges,
        speed_edges,
        files_read,
        args.low_height_audit_threshold_km,
    )
    plot_solar_counts(args.solar_output, solar_years, solar_year_count_density, solar_all_count_density, solar_edges)
    plot_height_velocity(
        args.height_velocity_output,
        hv_counts,
        height_edges,
        speed_edges,
        "Initial detection height (km)",
        "Catalogue meteor count",
    )
    if args.include_measurements:
        plot_height_velocity(
            args.measurement_height_velocity_output,
            measurement_hv_counts,
            height_edges,
            speed_edges,
            "Measurement height (km)",
            "Measurement count",
        )
    print(
        f"orbit_catalogue_statistics events={len(arrays['sample_idx'])} "
        f"fit_events={len(fit_arrays['fit_sample_idx'])} "
        f"measurements={len(measurement_arrays['measurement_height_km'])} "
        f"low_height_events={len(low_height_diagnostics)} files_read={files_read}"
    )
    print(args.output_h5)
    print(args.solar_output)
    print(args.height_velocity_output)
    if args.include_measurements:
        print(args.measurement_height_velocity_output)


if __name__ == "__main__":
    main()
