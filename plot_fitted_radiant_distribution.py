#!/usr/bin/env python3
"""Plot DASST zenith-attraction corrected radiant distribution for PANSY meteors."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from radiant_visibility import radiant_above_local_horizon
from dasst_orbits_from_candidate_states import (
    DEFAULT_KERNEL,
    dasst_orbit_from_gcrs,
    first_radiant_column,
    radiant_payload_from_dasst,
)
from radiant_visibility import add_composite_visibility_boundary_hammer, add_visibility_boundary_hammer, visibility_samples_from_radiants


PLOT_CENTER_LONGITUDE_DEG = -90.0
SCATTER_DIAMETER_CONSTANT = 7_500.0
SCATTER_DIAMETER_MIN_PT = 2.2
SCATTER_DIAMETER_MAX_PT = 7.0
HEAD_ECHO_MIN_SPEED_KM_S = float(os.environ.get("PANSY_HEAD_ECHO_MIN_SPEED_KM_S", "6.0"))


def wrap180(deg):
    return (np.asarray(deg, dtype=np.float64) + 180.0) % 360.0 - 180.0


def wrap360(deg):
    return np.asarray(deg, dtype=np.float64) % 360.0


def centered_plot_longitude_deg(lambda_minus_sun_signed_deg):
    return -wrap180(np.asarray(lambda_minus_sun_signed_deg, dtype=np.float64) - PLOT_CENTER_LONGITUDE_DEG)


def centered_tick_labels():
    tick_positions_deg = np.arange(-150.0, 180.0, 30.0)
    labels = []
    for tick in tick_positions_deg:
        label_value = int(wrap180(PLOT_CENTER_LONGITUDE_DEG - tick))
        labels.append("" if label_value == -75 else f"{label_value}°")
    return tick_positions_deg, labels


def scaled_scatter_area(n_points: int) -> float:
    """Return matplotlib scatter marker area with roughly constant total coverage."""
    n = max(1, int(n_points))
    diameter_pt = np.sqrt(SCATTER_DIAMETER_CONSTANT / n)
    diameter_pt = float(np.clip(diameter_pt, SCATTER_DIAMETER_MIN_PT, SCATTER_DIAMETER_MAX_PT))
    return diameter_pt**2


def circular_mean_deg(values_deg):
    values = np.asarray(values_deg, dtype=np.float64)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan
    rad = np.deg2rad(values)
    return float(np.rad2deg(np.arctan2(np.mean(np.sin(rad)), np.mean(np.cos(rad)))) % 360.0)


def accepted_hypothesis_from_diagnostics(input_dir: Path, sample_idx: int):
    path = input_dir / f"pansy_disambiguation_diagnostics_{sample_idx}.h5"
    if not path.exists():
        return None
    try:
        with h5py.File(path, "r") as h:
            selected = str(h.attrs.get("selected_hypothesis", ""))
            if not selected or "hypotheses" not in h or selected not in h["hypotheses"]:
                return None
            grp = h["hypotheses"][selected]
            if bool(grp.attrs.get("combined_reject", True)):
                return None
            if not bool(grp.attrs.get("combined_good_fit", False)):
                return None
            return selected
    except OSError:
        return None


def first_existing_hypothesis(h5: h5py.File, accepted: str | None, include_rejected_best: bool):
    if accepted is not None and accepted in h5:
        return accepted
    if not include_rejected_best:
        return None
    winners = [name for name, grp in h5.items() if int(grp.attrs.get("combined_rank", -1)) == 0]
    return sorted(winners)[0] if winners else None


def radiant_row_from_payload(sample_idx: int, label: str, grp, payload: dict[str, np.ndarray], source: str):
    radiant_gcrs = first_radiant_column(payload, "radiant_orbit_GCRS", 2)
    radiant_gme = first_radiant_column(payload, "radiant_orbit_GeocentricMeanEcliptic", 2)
    radiant_sun_gme = first_radiant_column(payload, "radiant_sun_GeocentricMeanEcliptic", 2)
    radiant_state_gme = first_radiant_column(payload, "radiant_orbit_states_GeocentricMeanEcliptic", 6)
    speed_km_s = float(np.linalg.norm(radiant_state_gme[3:6]) / 1e3)
    lambda_minus_sun = wrap360(radiant_gme[0] - radiant_sun_gme[0])
    lambda_minus_sun_signed = wrap180(lambda_minus_sun)
    return {
        "sample_idx": sample_idx,
        "hypothesis": label,
        "epoch_unix": float(grp.attrs.get("epoch_unix", sample_idx / 1e6)),
        "combined_score": float(grp.attrs.get("combined_score", np.nan)),
        "selection_redchi": float(grp.attrs.get("fixed_velocity_reduced_chi2", np.nan)),
        "first_alt_km": float(grp.attrs.get("first_alt_km", np.nan)),
        "source": source,
        "radiant_ra_gcrs_deg": float(radiant_gcrs[0]),
        "radiant_dec_gcrs_deg": float(radiant_gcrs[1]),
        "radiant_lambda_ecliptic_deg": float(wrap360(radiant_gme[0])),
        "radiant_beta_ecliptic_deg": float(radiant_gme[1]),
        "sun_lambda_ecliptic_deg": float(wrap360(radiant_sun_gme[0])),
        "sun_beta_ecliptic_deg": float(radiant_sun_gme[1]),
        "lambda_minus_sun_deg": float(lambda_minus_sun),
        "lambda_minus_sun_signed_deg": float(lambda_minus_sun_signed),
        "plot_longitude_deg": float(centered_plot_longitude_deg(lambda_minus_sun_signed)),
        "speed_km_s": speed_km_s,
    }


def radiant_row_is_physically_possible(row: dict) -> bool:
    ra = float(row.get("radiant_ra_gcrs_deg", np.nan))
    dec = float(row.get("radiant_dec_gcrs_deg", np.nan))
    epoch = float(row.get("epoch_unix", np.nan))
    speed = float(row.get("speed_km_s", np.nan))
    return bool(
        np.isfinite(ra)
        and np.isfinite(dec)
        and np.isfinite(epoch)
        and np.isfinite(speed)
        and radiant_above_local_horizon(ra, dec, epoch)
        and speed >= HEAD_ECHO_MIN_SPEED_KM_S
    )


def stored_radiant_payload(grp):
    payload = {}
    for frame in ("GCRS", "GeocentricMeanEcliptic"):
        for prefix in ("radiant_obs", "radiant_orbit", "radiant_sun", "radiant_obs_states", "radiant_orbit_states"):
            key = f"{prefix}_{frame}"
            if key in grp:
                payload[key] = np.asarray(grp[key], dtype=np.float64)
    required = ("radiant_orbit_GeocentricMeanEcliptic", "radiant_sun_GeocentricMeanEcliptic")
    return payload if all(key in payload for key in required) else None


def collect_one(task):
    input_dir, sample_idx, accepted, include_rejected_best, kernel, recompute_missing = task
    orbit_path = input_dir / f"pansy_candidate_orbits_dasst_{sample_idx}.h5"
    state_path = input_dir / f"pansy_candidate_orbit_states_{sample_idx}.h5"
    if not orbit_path.exists() or not state_path.exists():
        return None
    label = accepted
    try:
        with h5py.File(orbit_path, "r") as ho:
            label = first_existing_hypothesis(ho, accepted, include_rejected_best)
            if label is None:
                label = accepted
            if label is not None and label in ho:
                grp = ho[label]
                payload = stored_radiant_payload(grp)
                if payload is not None:
                    row = radiant_row_from_payload(sample_idx, label, grp, payload, "stored_dasst_h5")
                    return row if radiant_row_is_physically_possible(row) else None
    except OSError:
        return None
    if not recompute_missing:
        return None
    try:
        with h5py.File(state_path, "r") as hs, h5py.File(orbit_path, "r") as ho:
            if label is None:
                label = first_existing_hypothesis(hs, accepted, include_rejected_best)
            if label is None or label not in hs:
                return None
            state_grp = hs[label]
            attr_grp = ho[label] if label in ho else state_grp
            state = np.asarray(state_grp["state_gcrs_m_mps"][()], dtype=np.float64)
            if state.shape != (6,) or not np.all(np.isfinite(state)):
                return None
            epoch = float(state_grp.attrs.get("epoch_unix", attr_grp.attrs.get("epoch_unix", sample_idx / 1e6)))
            od = dasst_orbit_from_gcrs(state[None, :], epoch, Path(kernel))
            payload = radiant_payload_from_dasst(od)
            row = radiant_row_from_payload(sample_idx, label, attr_grp, payload, "recomputed_nominal_dasst")
            return row if radiant_row_is_physically_possible(row) else None
    except Exception as exc:
        return {"sample_idx": sample_idx, "error": repr(exc)}


def collect_fitted_radiants(input_dir: Path, include_rejected_best: bool, recompute_missing: bool, jobs: int, kernel: Path):
    tasks = []
    for path in sorted(input_dir.glob("pansy_candidate_orbits_dasst_*.h5")):
        sample_idx = int(path.stem.split("_")[-1])
        accepted = accepted_hypothesis_from_diagnostics(input_dir, sample_idx)
        if accepted is None and not include_rejected_best:
            continue
        tasks.append((input_dir, sample_idx, accepted, include_rejected_best, str(kernel), recompute_missing))
    if jobs <= 1:
        results = [collect_one(task) for task in tasks]
    else:
        with mp.get_context("spawn").Pool(processes=jobs) as pool:
            results = list(pool.imap_unordered(collect_one, tasks, chunksize=4))
    rows = [r for r in results if isinstance(r, dict) and "error" not in r]
    errors = [r for r in results if isinstance(r, dict) and "error" in r]
    return rows, errors, len(tasks)


def rows_to_arrays(rows):
    dtype = [
        ("sample_idx", "i8"),
        ("hypothesis", "S8"),
        ("epoch_unix", "f8"),
        ("combined_score", "f8"),
        ("selection_redchi", "f8"),
        ("first_alt_km", "f8"),
        ("source", "S32"),
        ("radiant_ra_gcrs_deg", "f8"),
        ("radiant_dec_gcrs_deg", "f8"),
        ("radiant_lambda_ecliptic_deg", "f8"),
        ("radiant_beta_ecliptic_deg", "f8"),
        ("sun_lambda_ecliptic_deg", "f8"),
        ("sun_beta_ecliptic_deg", "f8"),
        ("lambda_minus_sun_deg", "f8"),
        ("lambda_minus_sun_signed_deg", "f8"),
        ("plot_longitude_deg", "f8"),
        ("speed_km_s", "f8"),
    ]
    arr = np.zeros(len(rows), dtype=dtype)
    for i, row in enumerate(rows):
        for name in arr.dtype.names:
            if name in ("hypothesis", "source"):
                arr[name][i] = row[name].encode("ascii", "ignore")
            else:
                arr[name][i] = row[name]
    return arr


def arrays_to_rows(arr):
    rows = []
    for rec in arr:
        row = {}
        for name in arr.dtype.names:
            val = rec[name]
            if name in ("hypothesis", "source"):
                val = val.decode("ascii", "ignore") if isinstance(val, bytes) else str(val)
            elif np.issubdtype(arr.dtype[name], np.integer):
                val = int(val)
            else:
                val = float(val)
            row[name] = val
        rows.append(row)
    return rows


def write_h5(path: Path, rows, errors, task_count: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    arr = rows_to_arrays(rows)
    with h5py.File(path, "w") as h:
        h.attrs["script"] = Path(__file__).name
        h.attrs["source"] = "DASST rebound_od radiant_orbit_GeocentricMeanEcliptic"
        h.attrs["selection"] = "accepted selected hypotheses unless include_rejected_best was requested"
        h.attrs["physical_filter"] = (
            "radiant GCRS RA/Dec must be above the PANSY local horizon at event time "
            f"and speed_km_s >= {HEAD_ECHO_MIN_SPEED_KM_S:g}"
        )
        h.attrs["radiant_definition"] = "DASST zenith-attraction corrected orbit radiant"
        h.attrs["longitude_convention"] = "lambda_radiant - lambda_sun wrapped to [0, 360) deg; signed copy in [-180, 180)"
        h.attrs["plot_center_longitude_deg"] = PLOT_CENTER_LONGITUDE_DEG
        h.attrs["task_count"] = int(task_count)
        h.attrs["error_count"] = int(len(errors))
        h.create_dataset("radiants", data=arr)
        if errors:
            err_dtype = [("sample_idx", "i8"), ("error", h5py.string_dtype("utf-8"))]
            err = np.empty(len(errors), dtype=err_dtype)
            for i, row in enumerate(errors):
                err["sample_idx"][i] = int(row["sample_idx"])
                err["error"][i] = row["error"]
            h.create_dataset("errors", data=err)


def add_source_markers(ax):
    for label, lon_deg, lat_deg, marker, color, size in [
        ("Helion", 0.0, 0.0, "o", "#ffd21f", 58),
        ("Apex", -90.0, 0.0, r"$\otimes$", "black", 120),
        ("Antihelion", 180.0, 0.0, "o", "black", 58),
    ]:
        x = np.deg2rad(centered_plot_longitude_deg(wrap180(lon_deg)))
        y = np.deg2rad(lat_deg)
        ax.scatter(
            x,
            y,
            marker=marker,
            s=size,
            color=color,
            edgecolor="black" if marker == "o" else None,
            linewidth=0.35 if marker == "o" else 0.0,
            zorder=10,
        )


def plot_radiants(rows, output_png: Path):
    output_png.parent.mkdir(parents=True, exist_ok=True)
    arr = rows_to_arrays(rows)
    good = (
        np.isfinite(arr["plot_longitude_deg"])
        & np.isfinite(arr["radiant_beta_ecliptic_deg"])
        & np.isfinite(arr["speed_km_s"])
    )
    arr = arr[good]
    fig = plt.figure(figsize=(9.4, 6.2), constrained_layout=True)
    ax = fig.add_subplot(111, projection="hammer")
    sc = ax.scatter(
        np.deg2rad(arr["plot_longitude_deg"]),
        np.deg2rad(arr["radiant_beta_ecliptic_deg"]),
        c=arr["speed_km_s"],
        s=scaled_scatter_area(len(arr)),
        cmap="turbo",
        vmin=10.0,
        vmax=80.0,
        alpha=0.72,
        linewidths=0.15,
        edgecolors="white",
        zorder=3,
    )
    add_source_markers(ax)
    sample_epoch, sample_sun = visibility_samples_from_radiants(arr)
    if len(sample_epoch):
        add_composite_visibility_boundary_hammer(ax, sample_epoch, sample_sun, color="black")
    query = circular_mean_deg(arr["sun_lambda_ecliptic_deg"])
    tick_pos, tick_labels = centered_tick_labels()
    ax.set_xticks(np.deg2rad(tick_pos))
    ax.set_xticklabels(tick_labels)
    ax.grid(True, alpha=0.42)
    ax.set_xlabel(r"Sun-centered ecliptic longitude, $\lambda-\lambda_\odot$ (apex centered at $-90^\circ$)", labelpad=18)
    ax.set_ylabel(r"Ecliptic latitude, $\beta$")
    cb = fig.colorbar(sc, ax=ax, orientation="horizontal", pad=0.14, fraction=0.046)
    cb.set_label("Geocentric velocity (km/s)")

    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def plot_radiants_with_options(
    rows,
    output_png: Path,
    show_shower_overlays: bool = False,
    shower_catalog: Path | None = None,
    shower_solar_longitude: float | None = None,
    shower_date: str | None = None,
    shower_peak_tolerance_deg: float = 5.0,
):
    if not show_shower_overlays:
        plot_radiants(rows, output_png)
        return
    output_png.parent.mkdir(parents=True, exist_ok=True)
    arr = rows_to_arrays(rows)
    good = (
        np.isfinite(arr["plot_longitude_deg"])
        & np.isfinite(arr["radiant_beta_ecliptic_deg"])
        & np.isfinite(arr["speed_km_s"])
    )
    arr = arr[good]
    fig = plt.figure(figsize=(9.4, 6.2), constrained_layout=True)
    ax = fig.add_subplot(111, projection="hammer")
    sc = ax.scatter(
        np.deg2rad(arr["plot_longitude_deg"]),
        np.deg2rad(arr["radiant_beta_ecliptic_deg"]),
        c=arr["speed_km_s"],
        s=scaled_scatter_area(len(arr)),
        cmap="turbo",
        vmin=10.0,
        vmax=80.0,
        alpha=0.72,
        linewidths=0.15,
        edgecolors="white",
        zorder=3,
    )
    add_source_markers(ax)
    from shower_radiant_overlay import active_showers, add_shower_overlay_hammer

    showers, query = active_showers(
        arr,
        shower_catalog=shower_catalog,
        shower_solar_longitude=shower_solar_longitude,
        shower_date=shower_date,
        peak_tolerance_deg=shower_peak_tolerance_deg,
    )
    sample_epoch, sample_sun = visibility_samples_from_radiants(arr)
    if len(sample_epoch):
        add_composite_visibility_boundary_hammer(ax, sample_epoch, sample_sun, color="black")
    elif np.isfinite(query):
        add_visibility_boundary_hammer(ax, query, color="black")
    add_shower_overlay_hammer(ax, showers)
    tick_pos, tick_labels = centered_tick_labels()
    ax.set_xticks(np.deg2rad(tick_pos))
    ax.set_xticklabels(tick_labels)
    ax.grid(True, alpha=0.42)
    ax.set_xlabel(r"Sun-centered ecliptic longitude, $\lambda-\lambda_\odot$ (apex centered at $-90^\circ$)", labelpad=18)
    ax.set_ylabel(r"Ecliptic latitude, $\beta$")
    cb = fig.colorbar(sc, ax=ax, orientation="horizontal", pad=0.14, fraction=0.046)
    cb.set_label("Geocentric velocity (km/s)")
    ax.legend(loc="lower left", fontsize=8, frameon=True)
    ax.set_title(f"IAU established showers at solar longitude {float(query):.1f} deg", fontsize=10)
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot PANSY radiant distribution from DASST orbit-determination products.")
    parser.add_argument("--input-dir", type=Path, default=Path("test_plots"))
    parser.add_argument("--output-png", type=Path, default=Path("test_plots/fitted_radiant_distribution.png"))
    parser.add_argument("--output-h5", type=Path, default=Path("test_plots/fitted_radiant_distribution.h5"))
    parser.add_argument("--kernel", type=Path, default=DEFAULT_KERNEL)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--include-rejected-best", action="store_true", help="Include rank-0 best states even when the event was rejected.")
    parser.add_argument("--no-recompute-missing", action="store_true", help="Use only radiant datasets already present in DASST HDF5 files.")
    parser.add_argument("--reuse-h5", action="store_true", help="Only redraw the PNG from an existing output HDF5.")
    parser.add_argument("--show-shower-overlays", action="store_true", help="Overlay active IAU established meteor showers.")
    parser.add_argument("--shower-catalog", type=Path, default=None, help="IAU established shower catalogue path.")
    parser.add_argument("--shower-solar-longitude", type=float, default=None, help="Solar longitude for shower overlay selection.")
    parser.add_argument("--shower-date", default=None, help="UTC date/datetime for shower overlay selection.")
    parser.add_argument("--shower-peak-tolerance-deg", type=float, default=5.0)
    args = parser.parse_args()

    if args.reuse_h5:
        with h5py.File(args.output_h5, "r") as h:
            rows = arrays_to_rows(h["radiants"][()])
            errors = []
            task_count = int(h.attrs.get("task_count", len(rows)))
    else:
        rows, errors, task_count = collect_fitted_radiants(
            args.input_dir,
            include_rejected_best=args.include_rejected_best,
            recompute_missing=not args.no_recompute_missing,
            jobs=max(1, args.jobs),
            kernel=args.kernel,
        )
        write_h5(args.output_h5, rows, errors, task_count)
    plot_radiants_with_options(
        rows,
        args.output_png,
        show_shower_overlays=args.show_shower_overlays,
        shower_catalog=args.shower_catalog,
        shower_solar_longitude=args.shower_solar_longitude,
        shower_date=args.shower_date,
        shower_peak_tolerance_deg=args.shower_peak_tolerance_deg,
    )
    print(f"dasst_radiants {len(rows)} tasks {task_count} errors {len(errors)}")
    print(args.output_png)
    print(args.output_h5)


if __name__ == "__main__":
    main()
