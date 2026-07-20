#!/usr/bin/env python3
"""Build and process a reproducible random subset of PANSY mass profiles."""

from __future__ import annotations

import argparse
import datetime as dt
import os
import signal
import sys
import traceback
from pathlib import Path

import h5py
import numpy as np
import scipy.optimize as opt
from scipy.stats import chi2 as chi2_distribution

import fit_best_alias_physics_models as physics


METEOROID_DENSITY_KG_M3 = 3000.0
DEFAULT_SEED = 20260717


class PointTimeout(Exception):
    pass


def timeout_handler(_signum, _frame):
    raise PointTimeout()


def radius_um_to_mass_kg(radius_um):
    radius_m = np.asarray(radius_um, dtype=np.float64) * 1e-6
    return (4.0 / 3.0) * np.pi * METEOROID_DENSITY_KG_M3 * radius_m**3


def weighted_quantile(values, weights, quantiles):
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    quantiles = np.asarray(quantiles, dtype=np.float64)
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(good):
        return np.full(quantiles.shape, np.nan)
    values = values[good]
    weights = weights[good]
    order = np.argsort(values)
    values = values[order]
    cumulative = np.cumsum(weights[order])
    cumulative /= cumulative[-1]
    return np.interp(quantiles, cumulative, values, left=values[0], right=values[-1])


def log_grid_probability(radius_um, delta_chi2):
    radius_um = np.asarray(radius_um, dtype=np.float64)
    probability_density = np.exp(-0.5 * np.clip(delta_chi2, 0.0, 700.0))
    probability_density[~np.isfinite(delta_chi2)] = 0.0
    log_radius = np.log(radius_um)
    edges = np.empty(len(log_radius) + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (log_radius[:-1] + log_radius[1:])
    edges[0] = log_radius[0] - 0.5 * (log_radius[1] - log_radius[0])
    edges[-1] = log_radius[-1] + 0.5 * (log_radius[-1] - log_radius[-2])
    weights = probability_density * np.diff(edges)
    total = np.sum(weights)
    if total > 0.0:
        weights /= total
        probability_density /= total
    return probability_density, weights


def threshold_crossing(log_radius, delta_chi2, inside_index, outside_index, threshold):
    x1 = float(log_radius[inside_index])
    x2 = float(log_radius[outside_index])
    y1 = float(delta_chi2[inside_index])
    y2 = float(delta_chi2[outside_index])
    if not np.isfinite(y1) or not np.isfinite(y2) or y1 == y2:
        return float(np.exp(x1))
    fraction = np.clip((threshold - y1) / (y2 - y1), 0.0, 1.0)
    return float(np.exp(x1 + fraction * (x2 - x1)))


def profile_interval(radius_um, delta_chi2, threshold):
    radius_um = np.asarray(radius_um, dtype=np.float64)
    delta_chi2 = np.asarray(delta_chi2, dtype=np.float64)
    finite = np.isfinite(delta_chi2)
    if not np.any(finite):
        return np.nan, np.nan, False, False, "no_finite_profile", "no_finite_profile"
    best_index = int(np.nanargmin(delta_chi2))
    accepted = finite & (delta_chi2 <= float(threshold))
    if not accepted[best_index]:
        return np.nan, np.nan, False, False, "no_accepted_profile", "no_accepted_profile"

    left = best_index
    while left > 0 and accepted[left - 1]:
        left -= 1
    right = best_index
    while right + 1 < len(radius_um) and accepted[right + 1]:
        right += 1

    log_radius = np.log(radius_um)
    lower_bounded = left > 0 and np.isfinite(delta_chi2[left - 1])
    upper_bounded = right + 1 < len(radius_um) and np.isfinite(delta_chi2[right + 1])
    lower_status = "bounded" if lower_bounded else ("open_grid" if left == 0 else "fit_failure")
    upper_status = "bounded" if upper_bounded else ("open_grid" if right + 1 == len(radius_um) else "fit_failure")
    lower = (
        threshold_crossing(log_radius, delta_chi2, left, left - 1, threshold)
        if lower_bounded
        else 0.0
    )
    upper = (
        threshold_crossing(log_radius, delta_chi2, right, right + 1, threshold)
        if upper_bounded
        else np.inf
    )
    return lower, upper, lower_bounded, upper_bounded, lower_status, upper_status


def iter_diagnostics(events_dir: Path):
    for day_entry in sorted(os.scandir(events_dir), key=lambda item: item.name):
        if not day_entry.is_dir(follow_symlinks=False):
            continue
        for event_entry in sorted(os.scandir(day_entry.path), key=lambda item: item.name):
            if (
                event_entry.is_file(follow_symlinks=False)
                and event_entry.name.startswith("pansy_disambiguation_diagnostics_")
                and event_entry.name.endswith(".h5")
            ):
                yield Path(event_entry.path)


def build_manifest(events_dir: Path, output_h5: Path, sample_size: int, seed: int):
    rng = np.random.default_rng(int(seed))
    reservoir: list[str] = []
    total = 0
    for path in iter_diagnostics(events_dir):
        total += 1
        path_string = str(path)
        if len(reservoir) < sample_size:
            reservoir.append(path_string)
        else:
            replacement = int(rng.integers(0, total))
            if replacement < sample_size:
                reservoir[replacement] = path_string
        if total % 100000 == 0:
            print(f"manifest_scan diagnostics={total}", flush=True)
    rng.shuffle(reservoir)
    sample_idx = np.asarray(
        [int(Path(path).stem.rsplit("_", 1)[1]) for path in reservoir],
        dtype=np.int64,
    )
    string_dtype = h5py.string_dtype(encoding="utf-8")
    output_h5.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_h5.with_suffix(output_h5.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.create_dataset("sample_idx", data=sample_idx)
        handle.create_dataset("diagnostics_h5", data=np.asarray(reservoir, dtype=object), dtype=string_dtype)
        handle.attrs["schema"] = "pansy.catalogue_mass_profile_manifest.v1"
        handle.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        handle.attrs["events_dir"] = str(events_dir)
        handle.attrs["random_seed"] = int(seed)
        handle.attrs["available_diagnostics"] = int(total)
        handle.attrs["sample_size"] = int(len(reservoir))
        handle.attrs["sampling"] = "uniform reservoir sample without replacement over diagnostics files"
    os.replace(temporary, output_h5)
    print(f"manifest_complete available={total} selected={len(reservoir)} output={output_h5}", flush=True)


def load_selected_fit(diagnostics_h5: Path):
    observations = physics.load_best_alias(diagnostics_h5)
    with h5py.File(diagnostics_h5, "r") as handle:
        label = handle.attrs["selected_hypothesis"]
        if isinstance(label, bytes):
            label = label.decode()
        group = handle["hypotheses"][label]
        required = (
            "physics_ceplecha_params",
            "physics_ceplecha_model",
            "physics_ceplecha_velocity_km_s",
        )
        if any(name not in group for name in required):
            raise RuntimeError("selected hypothesis has no shrinking-radius physics fit")
        fit = {
            "label": str(label),
            "params": np.asarray(group["physics_ceplecha_params"], dtype=np.float64),
            "model_km": np.asarray(group["physics_ceplecha_model"], dtype=np.float64),
            "velocity_km_s": np.asarray(group["physics_ceplecha_velocity_km_s"], dtype=np.float64),
            "keep": np.asarray(group.get("physics_ceplecha_keep", np.ones(len(observations["t_s"]))), dtype=bool),
            "physics_reduced_chi2": float(group.attrs.get("physics_ceplecha_reduced_chi2", np.nan)),
            "selection_reduced_chi2": float(group.attrs.get("selection_reduced_chi2", np.nan)),
            "combined_good_fit": bool(group.attrs.get("combined_good_fit", False)),
            "trajectory_quality_reject": bool(group.attrs.get("trajectory_quality_reject", False)),
        }
    return observations, fit


def process_profile(
    diagnostics_h5: Path,
    output_h5: Path,
    grid_n: int,
    radius_min_um: float,
    radius_max_um: float,
    point_timeout_s: int,
    max_nfev: int,
):
    observations, stored_fit = load_selected_fit(diagnostics_h5)
    t_s = observations["t_s"]
    points_km = observations["points_km"]
    doppler_km_s = observations["doppler_km_s"]
    finite = np.isfinite(t_s) & np.all(np.isfinite(points_km), axis=1) & np.isfinite(doppler_km_s)
    keep = stored_fit["keep"] & finite
    if np.count_nonzero(keep) < 10:
        keep = finite
    if np.count_nonzero(keep) < 10:
        raise RuntimeError("fewer than 10 finite position-Doppler measurements")

    stored_prediction = physics.predicted_doppler(stored_fit["model_km"], stored_fit["velocity_km_s"])
    position_residual = points_km[keep] - stored_fit["model_km"][keep]
    doppler_residual = doppler_km_s[keep] - stored_prediction[keep]
    sigma_position_km = float(np.sqrt(np.mean(position_residual**2)))
    sigma_doppler_km_s = float(np.sqrt(np.mean(doppler_residual**2)))
    if not np.isfinite(sigma_position_km) or sigma_position_km <= 0.0:
        raise RuntimeError("invalid fitted position residual RMS")
    if not np.isfinite(sigma_doppler_km_s) or sigma_doppler_km_s <= 0.0:
        raise RuntimeError("invalid fitted Doppler residual RMS")

    density, density_metadata = physics.pbal.density_interpolator(observations["sample_epoch_unix"])
    radius_um = np.geomspace(float(radius_min_um), float(radius_max_um), int(grid_n))
    log_radius_m = np.log10(radius_um * 1e-6)
    lower6 = np.asarray([-np.inf, -np.inf, 20e3, -90e3, -90e3, -90e3])
    upper6 = np.asarray([np.inf, np.inf, 220e3, 90e3, 90e3, 90e3])
    scale6 = np.asarray([1e5, 1e5, 1e5, 7e4, 7e4, 7e4])

    def residual6(parameters6, fixed_log_radius):
        parameters7 = np.r_[parameters6, fixed_log_radius]
        position, velocity, _radius, _mass, _success, _message = physics.propagate_shrinking_radius_model(
            parameters7, t_s, density
        )
        prediction = physics.predicted_doppler(position, velocity)
        return np.r_[
            ((points_km[keep] - position[keep]) / sigma_position_km).ravel(),
            (doppler_km_s[keep] - prediction[keep]) / sigma_doppler_km_s,
        ]

    def residual7(parameters7):
        return residual6(parameters7[:6], parameters7[6])

    initial7 = stored_fit["params"].copy()
    initial7[6] = np.clip(initial7[6], log_radius_m[0], log_radius_m[-1])
    free_result = opt.least_squares(
        residual7,
        initial7,
        bounds=(np.r_[lower6, log_radius_m[0]], np.r_[upper6, log_radius_m[-1]]),
        x_scale=np.r_[scale6, 1.0],
        loss="linear",
        max_nfev=max(2 * int(max_nfev), 80),
    )
    free_parameters = np.asarray(free_result.x, dtype=np.float64)
    free_chi2 = float(np.sum(residual7(free_parameters) ** 2))

    profile_chi2 = np.full(len(radius_um), np.nan)
    profile_parameters6 = np.full((len(radius_um), 6), np.nan)
    profile_success = np.zeros(len(radius_um), dtype=bool)
    evaluation_order = np.argsort(np.abs(log_radius_m - free_parameters[6]))
    previous = free_parameters[:6].copy()
    old_handler = signal.signal(signal.SIGALRM, timeout_handler)
    try:
        for index in evaluation_order:
            try:
                signal.alarm(int(point_timeout_s))
                candidates = []
                for start in (previous, free_parameters[:6]):
                    result = opt.least_squares(
                        lambda parameters: residual6(parameters, log_radius_m[index]),
                        start,
                        bounds=(lower6, upper6),
                        x_scale=scale6,
                        loss="linear",
                        max_nfev=int(max_nfev),
                    )
                    candidates.append(result)
                best = min(candidates, key=lambda result: result.cost)
                residual = residual6(best.x, log_radius_m[index])
                profile_chi2[index] = float(np.sum(residual**2))
                profile_parameters6[index] = best.x
                profile_success[index] = bool(best.success)
                previous = best.x.copy()
            except Exception:
                profile_success[index] = False
            finally:
                signal.alarm(0)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)

    finite_profile = np.isfinite(profile_chi2)
    if np.count_nonzero(finite_profile) < max(5, grid_n // 2):
        raise RuntimeError(f"only {np.count_nonzero(finite_profile)} of {grid_n} profile points completed")
    profile_best_index = int(np.nanargmin(profile_chi2))
    profile_best_radius_um = float(radius_um[profile_best_index])
    grid_seed = np.r_[profile_parameters6[profile_best_index], log_radius_m[profile_best_index]]
    refined_free_result = opt.least_squares(
        residual7,
        grid_seed,
        bounds=(np.r_[lower6, log_radius_m[0]], np.r_[upper6, log_radius_m[-1]]),
        x_scale=np.r_[scale6, 1.0],
        loss="linear",
        max_nfev=max(2 * int(max_nfev), 80),
    )
    refined_free_parameters = np.asarray(refined_free_result.x, dtype=np.float64)
    refined_free_chi2 = float(np.sum(residual7(refined_free_parameters) ** 2))
    if refined_free_chi2 < free_chi2:
        free_parameters = refined_free_parameters
        free_chi2 = refined_free_chi2
    minimum_chi2 = float(min(free_chi2, np.nanmin(profile_chi2)))
    delta_chi2 = profile_chi2 - minimum_chi2
    probability_density, probability_weights = log_grid_probability(radius_um, delta_chi2)
    marginal_radius_quantiles_um = weighted_quantile(radius_um, probability_weights, [0.025, 0.5, 0.975])
    free_best_radius_um = float(10.0 ** free_parameters[6] * 1e6)
    confidence_threshold = float(chi2_distribution.ppf(0.95, 1))
    profile_lower_um, profile_upper_um, lower_bounded, upper_bounded, lower_status, upper_status = profile_interval(
        radius_um, delta_chi2, confidence_threshold
    )

    free_position, free_velocity, _free_radius, _free_mass, _success, _message = (
        physics.propagate_shrinking_radius_model(free_parameters, t_s, density)
    )
    free_prediction = physics.predicted_doppler(free_position, free_velocity)
    path_length_km = float(np.sum(np.linalg.norm(np.diff(free_position[keep], axis=0), axis=1)))
    initial_speed_km_s = float(np.linalg.norm(free_parameters[3:6]) / 1e3)
    fitted_speed_mean_km_s = float(np.mean(np.linalg.norm(free_velocity[keep], axis=1)))
    free_position_rms_km = float(np.sqrt(np.mean(np.sum((points_km[keep] - free_position[keep]) ** 2, axis=1))))
    free_doppler_rms_km_s = float(np.sqrt(np.mean((doppler_km_s[keep] - free_prediction[keep]) ** 2)))

    output_h5.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_h5.with_suffix(output_h5.suffix + f".{os.getpid()}.tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["schema"] = "pansy.catalogue_mass_profile.v1"
        handle.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        handle.attrs["diagnostics_h5"] = str(diagnostics_h5)
        handle.attrs["sample_idx"] = int(observations["sample_idx"])
        handle.attrs["sample_epoch_unix"] = float(observations["sample_epoch_unix"])
        handle.attrs["selected_hypothesis"] = stored_fit["label"]
        handle.attrs["meteoroid_density_kg_m3"] = METEOROID_DENSITY_KG_M3
        handle.attrs["sigma_model"] = "RMS of stored free-fit residuals; one pooled position-component sigma and one Doppler sigma fixed across all models"
        handle.attrs["confidence_interval"] = "connected profile-likelihood interval containing the best fit; Delta chi-square threshold for one parameter"
        handle.attrs["marginal_interval"] = "2.5, 50, and 97.5 percent quantiles with density proportional to exp(-Delta chi-square/2) and a flat prior in log radius on the stored grid"
        handle.attrs["upper_bound_convention"] = "+inf when the profile-likelihood interval remains open at the maximum radius"
        profile = handle.create_group("profile")
        profile.create_dataset("radius_um", data=radius_um)
        profile.create_dataset("mass_kg", data=radius_um_to_mass_kg(radius_um))
        profile.create_dataset("chi2", data=profile_chi2)
        profile.create_dataset("delta_chi2", data=delta_chi2)
        profile.create_dataset("relative_probability_density_log_radius", data=probability_density)
        profile.create_dataset("probability_weight", data=probability_weights)
        profile.create_dataset("success", data=profile_success)
        profile.create_dataset("parameters6", data=profile_parameters6)
        result = handle.create_group("result")
        result.create_dataset("free_best_parameters7", data=free_parameters)
        result.create_dataset("free_best_radius_um", data=free_best_radius_um)
        result.create_dataset("free_best_mass_kg", data=radius_um_to_mass_kg(free_best_radius_um))
        result.create_dataset("profile_grid_best_radius_um", data=profile_best_radius_um)
        result.create_dataset("profile_grid_best_mass_kg", data=radius_um_to_mass_kg(profile_best_radius_um))
        result.create_dataset("profile_ci95_lower_radius_um", data=profile_lower_um)
        result.create_dataset("profile_ci95_upper_radius_um", data=profile_upper_um)
        result.create_dataset("profile_ci95_lower_mass_kg", data=radius_um_to_mass_kg(profile_lower_um))
        result.create_dataset("profile_ci95_upper_mass_kg", data=radius_um_to_mass_kg(profile_upper_um))
        result.create_dataset("profile_ci95_lower_bounded", data=lower_bounded)
        result.create_dataset("profile_ci95_upper_bounded", data=upper_bounded)
        result.attrs["profile_ci95_lower_status"] = lower_status
        result.attrs["profile_ci95_upper_status"] = upper_status
        result.create_dataset("profile_ci95_delta_chi2_threshold", data=confidence_threshold)
        result.create_dataset("marginal_radius_quantiles_um", data=marginal_radius_quantiles_um)
        result.create_dataset("marginal_mass_quantiles_kg", data=radius_um_to_mass_kg(marginal_radius_quantiles_um))
        result.create_dataset("free_chi2", data=free_chi2)
        result.create_dataset("minimum_chi2", data=minimum_chi2)
        quality = handle.create_group("quality")
        quality.create_dataset("n_measurements", data=int(np.count_nonzero(keep)))
        quality.create_dataset("sigma_position_component_km", data=sigma_position_km)
        quality.create_dataset("sigma_doppler_km_s", data=sigma_doppler_km_s)
        quality.create_dataset("free_position_3d_rms_km", data=free_position_rms_km)
        quality.create_dataset("free_doppler_rms_km_s", data=free_doppler_rms_km_s)
        quality.create_dataset("path_length_km", data=path_length_km)
        quality.create_dataset("initial_speed_km_s", data=initial_speed_km_s)
        quality.create_dataset("fitted_speed_mean_km_s", data=fitted_speed_mean_km_s)
        quality.create_dataset("stored_physics_reduced_chi2", data=stored_fit["physics_reduced_chi2"])
        quality.create_dataset("selection_reduced_chi2", data=stored_fit["selection_reduced_chi2"])
        quality.create_dataset("combined_good_fit", data=stored_fit["combined_good_fit"])
        quality.create_dataset("trajectory_quality_reject", data=stored_fit["trajectory_quality_reject"])
        atmosphere = handle.create_group("atmosphere")
        atmosphere.create_dataset("altitude_km", data=density_metadata["msis_alt_grid_km"])
        atmosphere.create_dataset("density_kg_m3", data=density_metadata["msis_density_kg_m3"])
    os.replace(temporary, output_h5)


def run_worker(args):
    with h5py.File(args.manifest_h5, "r") as handle:
        sample_idx = np.asarray(handle["sample_idx"], dtype=np.int64)
        diagnostics = handle["diagnostics_h5"].asstr()[()]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_dir = args.output_dir / "worker_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    summary_path = log_dir / f"worker_{args.worker_index:03d}.tsv"
    assigned = np.arange(args.worker_index, len(sample_idx), args.worker_count, dtype=np.int64)
    with summary_path.open("a", buffering=1) as summary:
        for manifest_index in assigned:
            event_sample = int(sample_idx[manifest_index])
            output_h5 = args.output_dir / "profiles" / f"mass_profile_{event_sample}.h5"
            if output_h5.exists():
                continue
            try:
                process_profile(
                    Path(diagnostics[manifest_index]),
                    output_h5,
                    args.grid_n,
                    args.radius_min_um,
                    args.radius_max_um,
                    args.point_timeout_s,
                    args.max_nfev,
                )
                status = "ok"
            except Exception as error:
                status = f"ERR {type(error).__name__}: {error}"
                traceback.print_exc()
            summary.write(f"{manifest_index}\t{event_sample}\t{status}\n")
            print(f"worker={args.worker_index:03d} index={manifest_index} sample={event_sample} {status}", flush=True)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    manifest = subparsers.add_parser("build-manifest")
    manifest.add_argument("--events-dir", type=Path, required=True)
    manifest.add_argument("--output-h5", type=Path, required=True)
    manifest.add_argument("--sample-size", type=int, default=50000)
    manifest.add_argument("--seed", type=int, default=DEFAULT_SEED)
    worker = subparsers.add_parser("worker")
    worker.add_argument("--manifest-h5", type=Path, required=True)
    worker.add_argument("--output-dir", type=Path, required=True)
    worker.add_argument("--worker-index", type=int, required=True)
    worker.add_argument("--worker-count", type=int, required=True)
    worker.add_argument("--grid-n", type=int, default=41)
    worker.add_argument("--radius-min-um", type=float, default=1.0)
    worker.add_argument("--radius-max-um", type=float, default=10000.0)
    worker.add_argument("--point-timeout-s", type=int, default=8)
    worker.add_argument("--max-nfev", type=int, default=45)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.command == "build-manifest":
        build_manifest(args.events_dir, args.output_h5, args.sample_size, args.seed)
    elif args.command == "worker":
        run_worker(args)
    else:
        raise RuntimeError(f"unknown command {args.command}")


if __name__ == "__main__":
    main()
