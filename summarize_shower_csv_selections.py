#!/usr/bin/env python3
"""Summarize radview meteor-shower CSV selections from orbit metadata."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

import h5py
import numpy as np


EVENT_FIELDS = (
    "sample_idx",
    "radiant_ra_deg",
    "radiant_dec_deg",
    "v_g_km_s",
    "radiant_ecliptic_lon_deg",
    "radiant_ecliptic_lat_deg",
    "radiant_sun_ecliptic_lon_deg",
    "kepler",
)


def circular_mean_deg(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if not len(values):
        return math.nan
    radians = np.deg2rad(values)
    return float(np.rad2deg(np.arctan2(np.mean(np.sin(radians)), np.mean(np.cos(radians)))) % 360.0)


def read_selection(path: Path) -> tuple[list[int], dict[int, dict[str, float]]]:
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    header_index = next(i for i, line in enumerate(lines) if line.startswith("CurNum,"))
    sample_ids: list[int] = []
    csv_values: dict[int, dict[str, float]] = {}
    for row in csv.DictReader(lines[header_index:]):
        try:
            sample_id = int(row["MetCod"])
            values = {
                "SolLon": float(row["SolLon"]),
                "SCELoG": float(row["SCELoG"]),
                "ELaG": float(row["ELaG"]),
                "VG": float(row["VG"]),
            }
        except (KeyError, TypeError, ValueError):
            continue
        sample_ids.append(sample_id)
        csv_values[sample_id] = values
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError(f"{path}: duplicate MetCod values")
    return sample_ids, csv_values


def candidate_hour_files(metadata: Path, sample_id: int) -> list[Path]:
    event_time = datetime.fromtimestamp(sample_id / 1e6, tz=timezone.utc)
    paths: list[Path] = []
    for delta_hours in (-1, 0, 1):
        hour = (event_time + timedelta(hours=delta_hours)).replace(minute=0, second=0, microsecond=0)
        directory = metadata / hour.strftime("%Y-%m-%dT%H-00-00")
        paths.extend(sorted(directory.glob("*.h5")))
    return paths


def load_events(metadata: Path, wanted: set[int]) -> dict[int, np.void]:
    ids_by_file: dict[Path, set[int]] = defaultdict(set)
    for sample_id in wanted:
        candidates = candidate_hour_files(metadata, sample_id)
        if not candidates:
            continue
        for path in candidates:
            ids_by_file[path].add(sample_id)

    found: dict[int, np.void] = {}
    for path, possible_ids in sorted(ids_by_file.items()):
        with h5py.File(path, "r") as handle:
            if "events" not in handle:
                continue
            events = handle["events"][:]
        if not len(events):
            continue
        event_ids = np.asarray(events["sample_idx"], dtype=np.int64)
        keep = np.isin(event_ids, np.fromiter(possible_ids, dtype=np.int64))
        for event in events[keep]:
            sample_id = int(event["sample_idx"])
            if sample_id in found:
                raise RuntimeError(f"sample {sample_id} occurs in more than one orbit file")
            found[sample_id] = event
    return found


def jopek_dh(reference: np.ndarray, candidate: np.ndarray) -> float:
    """Return the Jopek (1993) hybrid orbital dissimilarity D_H."""
    q1, e1, i1, w1, node1 = reference
    q2, e2, i2, w2, node2 = candidate
    i1, w1, node1, i2, w2, node2 = np.deg2rad([i1, w1, node1, i2, w2, node2])

    def normal(inc: float, node: float) -> np.ndarray:
        return np.array([np.sin(inc) * np.sin(node), -np.sin(inc) * np.cos(node), np.cos(inc)])

    def perihelion(inc: float, argp: float, node: float) -> np.ndarray:
        return np.array(
            [
                np.cos(node) * np.cos(argp) - np.sin(node) * np.sin(argp) * np.cos(inc),
                np.sin(node) * np.cos(argp) + np.cos(node) * np.sin(argp) * np.cos(inc),
                np.sin(argp) * np.sin(inc),
            ]
        )

    cos_i = float(np.clip(np.dot(normal(i1, node1), normal(i2, node2)), -1.0, 1.0))
    cos_pi = float(np.clip(np.dot(perihelion(i1, w1, node1), perihelion(i2, w2, node2)), -1.0, 1.0))
    return float(
        np.sqrt(
            (e1 - e2) ** 2
            + ((q1 - q2) / (q1 + q2)) ** 2
            + (2.0 * np.sin(0.5 * np.arccos(cos_i))) ** 2
            + (((e1 + e2) / 2.0) * 2.0 * np.sin(0.5 * np.arccos(cos_pi))) ** 2
        )
    )


def summarize(
    path: Path,
    sample_ids: list[int],
    csv_values: dict[int, dict[str, float]],
    events: dict[int, np.void],
) -> dict[str, object]:
    missing = sorted(set(sample_ids) - events.keys())
    if missing:
        raise RuntimeError(f"{path}: missing {len(missing)} orbit events; first missing IDs: {missing[:10]}")
    rows = [events[sample_id] for sample_id in sample_ids]
    values = {field: np.asarray([row[field] for row in rows]) for field in EVENT_FIELDS[1:]}
    kepler = np.asarray(values["kepler"], dtype=np.float64)
    sun_centered_lon = (
        values["radiant_ecliptic_lon_deg"] - values["radiant_sun_ecliptic_lon_deg"]
    ) % 360.0
    csv_columns = {
        key: np.asarray([csv_values[sample_id][key] for sample_id in sample_ids], dtype=np.float64)
        for key in ("SolLon", "SCELoG", "ELaG", "VG")
    }
    metadata_columns = {
        "SolLon": np.asarray(values["radiant_sun_ecliptic_lon_deg"], dtype=np.float64),
        "SCELoG": sun_centered_lon,
        "ELaG": np.asarray(values["radiant_ecliptic_lat_deg"], dtype=np.float64),
        "VG": np.asarray(values["v_g_km_s"], dtype=np.float64),
    }
    csv_columns["EclipticLon"] = (csv_columns["SolLon"] + csv_columns["SCELoG"]) % 360.0
    metadata_columns["EclipticLon"] = np.asarray(values["radiant_ecliptic_lon_deg"], dtype=np.float64)
    validation: dict[str, dict[str, float]] = {}
    for key in csv_columns:
        residual = metadata_columns[key] - csv_columns[key]
        if key in {"SolLon", "SCELoG", "EclipticLon"}:
            residual = (residual + 180.0) % 360.0 - 180.0
        absolute = np.abs(residual)
        validation[key] = {
            "max_abs": float(np.nanmax(absolute)),
            "median_abs": float(np.nanmedian(absolute)),
            "rms": float(np.sqrt(np.nanmean(residual**2))),
            "p95_abs": float(np.nanpercentile(absolute, 95.0)),
        }
    selection_parameters = {
        "vg_km_s": np.asarray(values["v_g_km_s"], dtype=np.float64),
        "e": kepler[:, 1],
        "a_au": kepler[:, 0],
        "i_deg": kepler[:, 2],
        "q_au": kepler[:, 6],
        "ra_deg": np.asarray(values["radiant_ra_deg"], dtype=np.float64),
        "dec_deg": np.asarray(values["radiant_dec_deg"], dtype=np.float64),
    }
    parameter_distributions: dict[str, dict[str, float]] = {}
    for key, parameter_values in selection_parameters.items():
        finite = parameter_values[np.isfinite(parameter_values)]
        quantiles = np.nanpercentile(finite, [0.0, 1.0, 2.5, 5.0, 50.0, 95.0, 97.5, 99.0, 100.0])
        parameter_distributions[key] = {
            label: float(value)
            for label, value in zip(
                ("min", "p01", "p025", "p05", "median", "p95", "p975", "p99", "max"),
                quantiles,
                strict=True,
            )
        }
    mean_kepler = np.nanmean(kepler, axis=0)
    for column in (3, 4, 5):
        mean_kepler[column] = circular_mean_deg(kepler[:, column])

    # 169P/NEAT: q, e, i, omega, Omega.
    neat = np.array([0.604, 0.76796, 11.285, 218.13, 176.04])
    candidate = np.array([mean_kepler[6], mean_kepler[1], mean_kepler[2], mean_kepler[4], mean_kepler[3]])
    return {
        "selection_file": str(path),
        "count": len(rows),
        "solar_longitude_deg": circular_mean_deg(values["radiant_sun_ecliptic_lon_deg"]),
        "radiant_ra_deg": circular_mean_deg(values["radiant_ra_deg"]),
        "radiant_dec_deg": float(np.nanmean(values["radiant_dec_deg"])),
        "sun_centered_ecliptic_lon_deg": circular_mean_deg(sun_centered_lon),
        "sun_centered_ecliptic_lat_deg": float(np.nanmean(values["radiant_ecliptic_lat_deg"])),
        "v_g_km_s": float(np.nanmean(values["v_g_km_s"])),
        "a_au": float(mean_kepler[0]),
        "e": float(mean_kepler[1]),
        "i_deg": float(mean_kepler[2]),
        "node_deg": float(mean_kepler[3]),
        "arg_peri_deg": float(mean_kepler[4]),
        "true_anomaly_deg": float(mean_kepler[5]),
        "q_au": float(mean_kepler[6]),
        "d_h_169p_neat": jopek_dh(neat, candidate),
        "csv_orbit_validation": validation,
        "selection_parameter_distributions": parameter_distributions,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("selections", type=Path, nargs="+")
    args = parser.parse_args()

    selections: dict[Path, tuple[list[int], dict[int, dict[str, float]]]] = {
        path: read_selection(path) for path in args.selections
    }
    wanted = {sample_id for sample_ids, _ in selections.values() for sample_id in sample_ids}
    events = load_events(args.metadata, wanted)
    result = [
        summarize(path, sample_ids, csv_values, events)
        for path, (sample_ids, csv_values) in selections.items()
    ]
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
