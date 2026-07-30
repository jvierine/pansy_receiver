#!/usr/bin/env python3
"""Export a versioned, monthly PANSY catalogue snapshot for Zenodo."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np

from orbit_metadata_table import EVENT_DTYPE, PATH_DTYPE, _coerce_structured


SCHEMA_VERSION = "pansy_catalogue_v1"
SAMPLE_RATE_HZ = 1_000_000
DEFAULT_EXAMPLE_EVENT_ID = 1757458600402300
KEPLER_NAMES = (
    "semi_major_axis_AU",
    "eccentricity",
    "inclination_deg",
    "longitude_ascending_node_deg",
    "argument_perihelion_deg",
    "true_anomaly_deg",
    "perihelion_distance_AU",
)
FIT_PARAMETER_NAMES = (
    "east_m",
    "north_m",
    "up_m",
    "east_velocity_mps",
    "north_velocity_mps",
    "up_velocity_mps",
    "log10_beta_kg_m2",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-metadata-dir", required=True, type=Path)
    parser.add_argument("--cut-metadata-dir", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--release-version", default="1")
    parser.add_argument("--start-date", type=dt.date.fromisoformat)
    parser.add_argument("--stop-date", type=dt.date.fromisoformat)
    parser.add_argument("--example-event-id", type=int, default=DEFAULT_EXAMPLE_EVENT_ID)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parent,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def _set_common_attrs(h: h5py.File, release_version: str, product_level: str, commit: str) -> None:
    h.attrs["title"] = f"PANSY meteor head-echo catalogue Level {product_level}"
    h.attrs["release_version"] = str(release_version)
    h.attrs["schema_version"] = SCHEMA_VERSION
    h.attrs["created_utc"] = utc_now()
    h.attrs["processing_git_commit"] = commit
    h.attrs["event_id_definition"] = "Unix epoch sample index at 1 MHz"
    h.attrs["sample_rate_hz"] = SAMPLE_RATE_HZ
    h.attrs["site_name"] = "PANSY radar, Syowa Station, Antarctica"
    h.attrs["site_latitude_deg"] = -69.010833
    h.attrs["site_longitude_deg"] = 39.599722
    h.attrs["site_altitude_m"] = 100.0
    h.attrs["missing_float_value"] = "NaN"


def _create_dataset(
    group: h5py.Group,
    name: str,
    shape: tuple[int, ...],
    dtype,
    units: str = "",
    description: str = "",
) -> h5py.Dataset:
    trailing = shape[1:]
    chunk_rows = max(1, min(8192, 1_048_576 // max(np.dtype(dtype).itemsize * int(np.prod(trailing or (1,))), 1)))
    ds = group.create_dataset(
        name,
        shape=shape,
        maxshape=(None, *trailing),
        dtype=dtype,
        chunks=(chunk_rows, *trailing),
        compression="gzip",
        compression_opts=4,
        shuffle=True,
    )
    if units:
        ds.attrs["units"] = units
    if description:
        ds.attrs["description"] = description
    return ds


def _append(ds: h5py.Dataset, values) -> None:
    values = np.asarray(values, dtype=ds.dtype)
    if not len(values):
        return
    start = len(ds)
    ds.resize(start + len(values), axis=0)
    ds[start:] = values


def _init_level2(path: Path, release_version: str, commit: str) -> h5py.File:
    h = h5py.File(path, "w")
    _set_common_attrs(h, release_version, "2", commit)
    h.attrs["coordinate_frame"] = "local east-north-up (ENU)"
    events = h.create_group("events")
    _create_dataset(events, "event_id", (0,), "<i8", "sample", "Stable event join key")
    _create_dataset(events, "measurement_start", (0,), "<i8", "row", "First row in /measurements for this event")
    _create_dataset(events, "measurement_count", (0,), "<i4", "count", "Number of path measurements")

    measurements = h.create_group("measurements")
    _create_dataset(measurements, "event_id", (0,), "<i8", "sample", "Parent event identifier")
    _create_dataset(
        measurements,
        "time_sample_idx",
        (0,),
        "<i8",
        "sample",
        "Measurement epoch as Unix sample index at 1 MHz",
    )
    _create_dataset(measurements, "time_offset_s", (0,), "<f4", "s", "Time relative to event epoch")
    _create_dataset(measurements, "east_km", (0,), "<f4", "km")
    _create_dataset(measurements, "north_km", (0,), "<f4", "km")
    _create_dataset(measurements, "up_km", (0,), "<f4", "km")
    _create_dataset(measurements, "range_km", (0,), "<f4", "km", "Slant range, norm of ENU position")
    _create_dataset(measurements, "snr_db", (0,), "<f4", "dB", "10 log10 of stored linear SNR")
    _create_dataset(measurements, "doppler_velocity_mps", (0,), "<f4", "m s-1")
    _create_dataset(measurements, "beam_id", (0,), "<i1", "index", "Transmit beam identifier")
    _create_dataset(
        measurements,
        "selection_keep",
        (0,),
        "?",
        description="True where the measurement was retained by the final trajectory fit",
    )
    return h


def _init_level3(path: Path, release_version: str, commit: str) -> h5py.File:
    h = h5py.File(path, "w")
    _set_common_attrs(h, release_version, "3", commit)
    h.attrs["radiant_frame"] = "GCRS for RA/Dec; geocentric true ecliptic for ecliptic coordinates"
    h.attrs["initial_state_frame"] = "GCRS"
    h.attrs["initial_state_component_order"] = "x_m,y_m,z_m,vx_mps,vy_mps,vz_mps"
    h.attrs["fit_parameter_covariance_frame"] = "local ENU"
    h.attrs["fit_parameter_covariance_order"] = ",".join(FIT_PARAMETER_NAMES)
    h.attrs["kepler_order"] = ",".join(KEPLER_NAMES)
    h.attrs["orbital_uncertainty_method"] = (
        "The Kepler covariance is estimated by Monte Carlo sampling of the fitted trajectory "
        "covariance and propagating the sampled states through the geocentric and heliocentric "
        "transformations."
    )
    h.attrs["orbital_uncertainty_interpretation"] = (
        "First-order covariance summary, not a sampled posterior distribution or a precise "
        "estimate of tail probabilities."
    )
    h.attrs["covariance_caveat"] = (
        "fit_parameter_covariance is the covariance of the local ENU trajectory fit; "
        "it has not been transformed into a covariance of initial_state_gcrs_m_mps"
    )
    g = h.create_group("events")
    scalar_fields = [
        ("event_id", "<i8", "sample", "Stable event join key"),
        ("epoch_unix_s", "<f8", "s", "UTC Unix epoch"),
        ("solar_longitude_deg", "<f4", "deg", "Geocentric ecliptic longitude of Earth"),
        ("sun_ecliptic_longitude_deg", "<f4", "deg", "Apparent geocentric ecliptic longitude of the Sun"),
        ("v0_km_s", "<f4", "km s-1", "Initial atmospheric speed from the local trajectory fit"),
        ("geocentric_speed_km_s", "<f4", "km s-1", "Geocentric speed after attraction correction"),
        ("radiant_ra_deg", "<f4", "deg", "Geocentric radiant right ascension"),
        ("radiant_dec_deg", "<f4", "deg", "Geocentric radiant declination"),
        ("radiant_ecliptic_longitude_deg", "<f4", "deg", ""),
        ("radiant_ecliptic_latitude_deg", "<f4", "deg", ""),
        ("sun_centered_ecliptic_longitude_deg", "<f4", "deg", "Radiant longitude minus Sun longitude"),
        ("sun_centered_ecliptic_latitude_deg", "<f4", "deg", ""),
        ("radiant_angular_std_deg", "<f4", "deg", "Scalar angular uncertainty derived from velocity covariance"),
        ("initial_position_std_m", "<f4", "m", "sqrt(trace(C_position)) in the local fit frame"),
        ("initial_velocity_std_mps", "<f4", "m s-1", "sqrt(trace(C_velocity)) in the local fit frame"),
        ("combined_fit_score", "<f4", "", ""),
        (
            "n_uncertainty_samples",
            "<i4",
            "count",
            "Number of jointly propagated states used to estimate the Kepler uncertainty summary",
        ),
        ("fraction_eccentricity_gt_1", "<f4", "", ""),
        ("candidate_number", "<i4", "index", ""),
        ("combined_rank", "<i4", "rank", ""),
        ("n_aliases_orbit_tested", "<i4", "count", ""),
    ]
    for name, dtype, units, description in scalar_fields:
        _create_dataset(g, name, (0,), dtype, units, description)
    for name, units in zip(KEPLER_NAMES, ("AU", "", "deg", "deg", "deg", "deg", "AU")):
        _create_dataset(g, name, (0,), "<f4", units)
        _create_dataset(g, f"{name}_std", (0,), "<f4", units)
    _create_dataset(g, "initial_state_gcrs_m_mps", (0, 6), "<f4", "m and m s-1")
    _create_dataset(
        g,
        "geocentric_state_gcrs_m_mps",
        (0, 6),
        "<f4",
        "m and m s-1",
        "Initial position with the geocentric radiant velocity after zenith-attraction correction",
    )
    _create_dataset(g, "fit_parameters_enu", (0, 7), "<f4", "m, m s-1, and log10(kg m-2)")
    _create_dataset(g, "fit_parameter_covariance_enu", (0, 7, 7), "<f4", "mixed; see parameter order")
    _create_dataset(
        g,
        "kepler_covariance",
        (0, 7, 7),
        "<f4",
        "mixed; see Kepler order",
        "Sample covariance of the jointly propagated Kepler elements",
    )
    string_dtype = h5py.string_dtype("utf-8")
    _create_dataset(g, "selected_hypothesis", (0,), string_dtype)
    _create_dataset(g, "selection_model_type", (0,), string_dtype)
    _create_dataset(g, "orbit_solution_type", (0,), string_dtype)
    return h


def _decode_strings(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [v.decode("utf-8", "replace").rstrip("\x00") if isinstance(v, bytes) else str(v) for v in values],
        dtype=object,
    )


def append_level2(h: h5py.File, events: np.ndarray, paths: np.ndarray) -> None:
    if not len(events):
        return
    order = np.argsort(paths["sample_idx"], kind="stable") if len(paths) else np.asarray([], dtype=int)
    paths = paths[order]
    event_ids = events["sample_idx"].astype(np.int64)
    path_ids = paths["sample_idx"].astype(np.int64)
    starts = np.searchsorted(path_ids, event_ids, side="left")
    stops = np.searchsorted(path_ids, event_ids, side="right")
    base = len(h["measurements/event_id"])
    _append(h["events/event_id"], event_ids)
    _append(h["events/measurement_start"], base + starts)
    _append(h["events/measurement_count"], stops - starts)
    if not len(paths):
        return
    t_rel = paths["t_rel_s"].astype(np.float64)
    pos = paths["position_enu_km"].astype(np.float64)
    snr = paths["snr"].astype(np.float64)
    snr_db = np.full(len(snr), np.nan, dtype=np.float32)
    good_snr = np.isfinite(snr) & (snr > 0.0)
    snr_db[good_snr] = 10.0 * np.log10(snr[good_snr])
    _append(h["measurements/event_id"], path_ids)
    _append(h["measurements/time_sample_idx"], path_ids + np.rint(t_rel * SAMPLE_RATE_HZ).astype(np.int64))
    _append(h["measurements/time_offset_s"], t_rel)
    _append(h["measurements/east_km"], pos[:, 0])
    _append(h["measurements/north_km"], pos[:, 1])
    _append(h["measurements/up_km"], pos[:, 2])
    _append(h["measurements/range_km"], np.linalg.norm(pos, axis=1))
    _append(h["measurements/snr_db"], snr_db)
    _append(h["measurements/doppler_velocity_mps"], paths["doppler_mps"])
    _append(h["measurements/beam_id"], paths["beam_id"])
    _append(h["measurements/selection_keep"], paths["selection_keep"])


def append_level3(h: h5py.File, events: np.ndarray) -> None:
    if not len(events):
        return
    g = h["events"]
    event_id = events["sample_idx"].astype(np.int64)
    state = events["initial_state_gcrs_m_mps"].astype(np.float64)
    fit_cov = events["fit_parameter_covariance"].astype(np.float64)
    pos_var = np.trace(fit_cov[:, :3, :3], axis1=1, axis2=2)
    vel_var = np.trace(fit_cov[:, 3:6, 3:6], axis1=1, axis2=2)
    pos_std = np.sqrt(np.where(pos_var >= 0.0, pos_var, np.nan))
    vel_std = np.sqrt(np.where(vel_var >= 0.0, vel_var, np.nan))
    v0 = np.linalg.norm(events["fit_parameters"][:, 3:6].astype(np.float64), axis=1) / 1e3
    radiant_std = np.rad2deg(np.arctan2(vel_std, np.maximum(v0 * 1e3, 1.0)))
    sun_lon = np.mod(events["radiant_sun_ecliptic_lon_deg"].astype(np.float64), 360.0)
    ra = np.deg2rad(events["radiant_ra_deg"].astype(np.float64))
    dec = np.deg2rad(events["radiant_dec_deg"].astype(np.float64))
    vg_mps = events["v_g_km_s"].astype(np.float64) * 1e3
    geocentric_state = state.copy()
    geocentric_state[:, 3] = -vg_mps * np.cos(dec) * np.cos(ra)
    geocentric_state[:, 4] = -vg_mps * np.cos(dec) * np.sin(ra)
    geocentric_state[:, 5] = -vg_mps * np.sin(dec)
    values = {
        "event_id": event_id,
        "epoch_unix_s": event_id.astype(np.float64) / SAMPLE_RATE_HZ,
        "solar_longitude_deg": np.mod(sun_lon + 180.0, 360.0),
        "sun_ecliptic_longitude_deg": sun_lon,
        "v0_km_s": v0,
        "geocentric_speed_km_s": events["v_g_km_s"],
        "radiant_ra_deg": events["radiant_ra_deg"],
        "radiant_dec_deg": events["radiant_dec_deg"],
        "radiant_ecliptic_longitude_deg": events["radiant_ecliptic_lon_deg"],
        "radiant_ecliptic_latitude_deg": events["radiant_ecliptic_lat_deg"],
        "sun_centered_ecliptic_longitude_deg": np.mod(
            events["radiant_ecliptic_lon_deg"] - sun_lon, 360.0
        ),
        "sun_centered_ecliptic_latitude_deg": events["radiant_sun_ecliptic_lat_deg"],
        "radiant_angular_std_deg": radiant_std,
        "initial_position_std_m": pos_std,
        "initial_velocity_std_mps": vel_std,
        "combined_fit_score": events["combined_score"],
        "n_uncertainty_samples": events["n_uncertainty_samples"],
        "fraction_eccentricity_gt_1": events["frac_e_gt_1"],
        "candidate_number": events["candidate_number"],
        "combined_rank": events["combined_rank"],
        "n_aliases_orbit_tested": events["n_aliases_orbit_tested"],
        "initial_state_gcrs_m_mps": state,
        "geocentric_state_gcrs_m_mps": geocentric_state,
        "fit_parameters_enu": events["fit_parameters"],
        "fit_parameter_covariance_enu": fit_cov,
        "kepler_covariance": events["kepler_covariance"],
        "selected_hypothesis": _decode_strings(events["selected_hypothesis"]),
        "selection_model_type": _decode_strings(events["selection_model_type"]),
        "orbit_solution_type": _decode_strings(events["orbit_solution_type"]),
    }
    for i, name in enumerate(KEPLER_NAMES):
        values[name] = events["kepler"][:, i]
        values[f"{name}_std"] = events["kepler_std"][:, i]
    for name, value in values.items():
        _append(g[name], value)


def orbit_days(root: Path, start: dt.date | None, stop: dt.date | None) -> list[dt.date]:
    days = sorted(
        {
            dt.date.fromisoformat(path.parent.name[:10])
            for path in root.glob("*/orbit@*.h5")
            if len(path.parent.name) >= 10
        }
    )
    if start is not None:
        days = [day for day in days if day >= start]
    if stop is not None:
        days = [day for day in days if day <= stop]
    return days


def read_day(root: Path, day: dt.date) -> tuple[np.ndarray, np.ndarray, int]:
    event_chunks = []
    path_chunks = []
    skipped = 0
    for path in sorted(root.glob(f"{day.isoformat()}T*/orbit@*.h5")):
        try:
            with h5py.File(path, "r") as h:
                if "events" not in h:
                    skipped += 1
                    continue
                event_chunks.append(_coerce_structured(h["events"][()], EVENT_DTYPE))
                if "paths" in h:
                    path_chunks.append(_coerce_structured(h["paths"][()], PATH_DTYPE))
        except OSError:
            skipped += 1
    events = np.concatenate(event_chunks) if event_chunks else np.zeros(0, dtype=EVENT_DTYPE)
    paths = np.concatenate(path_chunks) if path_chunks else np.zeros(0, dtype=PATH_DTYPE)
    if len(events):
        order = np.argsort(events["sample_idx"], kind="stable")
        events = events[order]
        keep = np.r_[True, np.diff(events["sample_idx"]) != 0]
        events = events[keep]
        paths = paths[np.isin(paths["sample_idx"], events["sample_idx"])]
    return events, paths, skipped


def write_example_cut(
    cut_metadata_dir: Path,
    output_path: Path,
    event_id: int,
    release_version: str,
    commit: str,
) -> None:
    import digital_rf as drf

    reader = drf.DigitalMetadataReader(str(cut_metadata_dir))
    records = reader.read(event_id - 1, event_id + 1)
    if event_id not in records:
        raise KeyError(f"example event {event_id} is not in {cut_metadata_dir}")
    cut = records[event_id]
    with h5py.File(output_path, "w") as h:
        _set_common_attrs(h, release_version, "1", commit)
        h.attrs["event_id"] = event_id
        h.attrs["description"] = (
            "Finite raw-voltage cut for the example event in the PANSY catalogue article; "
            "this is not a continuous voltage recording"
        )
        h.attrs["voltage_units"] = "uncalibrated signed ADC counts"
        h.attrs["complex_storage"] = "real and imaginary components are stored separately"
        h.attrs["same_beam_revisit_s"] = 0.008
        h.attrs["inter_pulse_period_s"] = 0.0016
        h.attrs["receive_module_count"] = 7
        for name, value in cut.items():
            arr = np.asarray(value)
            dataset_kwargs = {}
            if arr.dtype.kind in "UO":
                string_dtype = h5py.string_dtype("utf-8")
                data = str(arr.item()) if arr.ndim == 0 else arr.astype(object)
                dataset_kwargs["dtype"] = string_dtype
            else:
                data = arr
            if arr.ndim:
                dataset_kwargs.update(compression="gzip", compression_opts=4, shuffle=True)
            ds = h.create_dataset(
                name,
                data=data,
                **dataset_kwargs,
            )
            if name in {"tx_idx", "c_tx_idx"}:
                ds.attrs["units"] = "Unix sample index at 1 MHz"
            elif name in {"delays"}:
                ds.attrs["units"] = "sample"
            elif name in {"c_range_km"}:
                ds.attrs["units"] = "km"
            elif name in {"c_doppler"}:
                ds.attrs["units"] = "m s-1"
            elif name in {"c_snr"}:
                ds.attrs["units"] = "linear power ratio"
            elif name.endswith("_re") or name.endswith("_im"):
                ds.attrs["units"] = "uncalibrated signed ADC count"


def validate_month(level2_path: Path, level3_path: Path) -> dict[str, int]:
    with h5py.File(level2_path, "r") as l2, h5py.File(level3_path, "r") as l3:
        ids2 = l2["events/event_id"][()]
        ids3 = l3["events/event_id"][()]
        measurement_ids = l2["measurements/event_id"][()]
        starts = l2["events/measurement_start"][()]
        counts = l2["events/measurement_count"][()]
        if len(ids2) != len(np.unique(ids2)) or len(ids3) != len(np.unique(ids3)):
            raise ValueError(f"duplicate event ID in {level2_path.name} or {level3_path.name}")
        if not np.array_equal(ids2, ids3):
            raise ValueError(f"Level 2/3 event IDs differ for {level2_path.stem}")
        if len(ids2) and not np.array_equal(starts, np.r_[0, np.cumsum(counts[:-1])]):
            raise ValueError(f"invalid Level 2 event offsets in {level2_path.name}")
        if int(np.sum(counts)) != len(measurement_ids):
            raise ValueError(f"invalid Level 2 measurement counts in {level2_path.name}")
        if len(measurement_ids) and not np.all(np.isin(measurement_ids, ids2)):
            raise ValueError(f"orphan Level 2 measurement in {level2_path.name}")
        return {"events": len(ids2), "measurements": len(measurement_ids)}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_release(args: argparse.Namespace) -> Path:
    output = args.output_dir.resolve()
    staging = output.with_name(f".{output.name}.building")
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"{output} exists; use --overwrite for an atomic rebuild")
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    (staging / "level2").mkdir()
    (staging / "level3").mkdir()
    commit = git_commit()
    days = orbit_days(args.orbit_metadata_dir, args.start_date, args.stop_date)
    if not days:
        raise RuntimeError("no orbit metadata days found")

    handles: dict[str, tuple[h5py.File, h5py.File]] = {}
    skipped_files = 0
    try:
        for index, day in enumerate(days, start=1):
            events, paths, skipped = read_day(args.orbit_metadata_dir, day)
            skipped_files += skipped
            month = day.strftime("%Y-%m")
            if month not in handles:
                level2_path = staging / "level2" / f"pansy_level2_{month}.h5"
                level3_path = staging / "level3" / f"pansy_level3_{month}.h5"
                handles[month] = (
                    _init_level2(level2_path, args.release_version, commit),
                    _init_level3(level3_path, args.release_version, commit),
                )
            append_level2(handles[month][0], events, paths)
            append_level3(handles[month][1], events)
            print(
                f"[{index:04d}/{len(days):04d}] {day}: "
                f"{len(events)} events, {len(paths)} measurements",
                flush=True,
            )
    finally:
        for level2, level3 in handles.values():
            level2.close()
            level3.close()

    if args.cut_metadata_dir is not None:
        write_example_cut(
            args.cut_metadata_dir,
            staging / f"pansy_level1_example_{args.example_event_id}.h5",
            args.example_event_id,
            args.release_version,
            commit,
        )

    summary = {
        "release_version": str(args.release_version),
        "schema_version": SCHEMA_VERSION,
        "created_utc": utc_now(),
        "processing_git_commit": commit,
        "coverage_start_utc_day": days[0].isoformat(),
        "coverage_stop_utc_day": days[-1].isoformat(),
        "orbit_files_skipped": skipped_files,
        "months": {},
    }
    for month in sorted(handles):
        stats = validate_month(
            staging / "level2" / f"pansy_level2_{month}.h5",
            staging / "level3" / f"pansy_level3_{month}.h5",
        )
        summary["months"][month] = stats
    summary["event_count"] = sum(item["events"] for item in summary["months"].values())
    summary["measurement_count"] = sum(item["measurements"] for item in summary["months"].values())
    (staging / "release_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    readme_source = Path(__file__).resolve().parent / "zenodo" / "catalogue_v1" / "README.md"
    shutil.copy2(readme_source, staging / "README.md")
    verifier_source = Path(__file__).resolve().parent / "zenodo" / "verify_level2_level3.py"
    shutil.copy2(verifier_source, staging / "verify_level2_level3.py")
    example_source = Path(__file__).resolve().parent / "zenodo" / "example_level2_radiant.py"
    shutil.copy2(example_source, staging / "example_level2_radiant.py")
    metadata_source = Path(__file__).resolve().parent / "zenodo" / "catalogue_v1" / "zenodo_metadata.json"
    shutil.copy2(metadata_source, staging / "zenodo_metadata.json")
    release_files = sorted(path for path in staging.rglob("*") if path.is_file())
    checksum_lines = [f"{sha256(path)}  {path.relative_to(staging)}" for path in release_files]
    (staging / "SHA256SUMS").write_text("\n".join(checksum_lines) + "\n")

    backup = output.with_name(f".{output.name}.previous")
    if backup.exists():
        shutil.rmtree(backup)
    if output.exists():
        os.replace(output, backup)
    try:
        os.replace(staging, output)
    except Exception:
        if backup.exists() and not output.exists():
            os.replace(backup, output)
        raise
    if backup.exists():
        shutil.rmtree(backup)
    print(f"validated release written to {output}")
    return output


def main() -> None:
    args = parse_args()
    build_release(args)


if __name__ == "__main__":
    main()
