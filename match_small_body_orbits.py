#!/usr/bin/env python3
"""Rank MPC NEA and comet orbits by the Jopek D_H criterion."""

from __future__ import annotations

import argparse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np


NEA_URL = "https://www.minorplanetcenter.net/iau/MPCORB/NEA.txt"
COMET_URL = "https://www.minorplanetcenter.net/iau/MPCORB/CometEls.txt"
DEFAULT_CACHE = Path.home() / ".cache" / "pansy_receiver" / "mpc"
DEFAULT_OUTPUT = Path("figs/ser_parent_body_jopek_matches.h5")


@dataclass(frozen=True)
class Orbit:
    q_au: float
    e: float
    i_deg: float
    omega_deg: float
    node_deg: float


def perihelion_direction(orbit: Orbit) -> np.ndarray:
    """Return the unit vector from the Sun toward perihelion."""
    inclination = np.deg2rad(orbit.i_deg)
    omega = np.deg2rad(orbit.omega_deg)
    node = np.deg2rad(orbit.node_deg)
    return np.asarray(
        (
            np.cos(node) * np.cos(omega)
            - np.sin(node) * np.sin(omega) * np.cos(inclination),
            np.sin(node) * np.cos(omega)
            + np.cos(node) * np.sin(omega) * np.cos(inclination),
            np.sin(omega) * np.sin(inclination),
        ),
        dtype=np.float64,
    )


def orbit_normal(orbit: Orbit) -> np.ndarray:
    inclination = np.deg2rad(orbit.i_deg)
    node = np.deg2rad(orbit.node_deg)
    return np.asarray(
        (
            np.sin(inclination) * np.sin(node),
            -np.sin(inclination) * np.cos(node),
            np.cos(inclination),
        ),
        dtype=np.float64,
    )


def vector_angle(first: np.ndarray, second: np.ndarray) -> float:
    return float(np.arctan2(np.linalg.norm(np.cross(first, second)), np.dot(first, second)))


def jopek_dh(reference: Orbit, candidate: Orbit) -> float:
    """Return the Jopek (1993) hybrid orbital dissimilarity criterion."""
    plane_angle = vector_angle(orbit_normal(reference), orbit_normal(candidate))
    perihelion_angle = vector_angle(
        perihelion_direction(reference),
        perihelion_direction(candidate),
    )

    q_scale = reference.q_au + candidate.q_au
    if q_scale <= 0.0:
        return np.inf
    distance_squared = (
        (reference.e - candidate.e) ** 2
        + ((reference.q_au - candidate.q_au) / q_scale) ** 2
        + (2.0 * np.sin(0.5 * plane_angle)) ** 2
        + (
            0.5
            * (reference.e + candidate.e)
            * 2.0
            * np.sin(0.5 * perihelion_angle)
        )
        ** 2
    )
    return float(np.sqrt(distance_squared))


def download_catalog(path: Path, url: str, refresh: bool) -> None:
    if path.exists() and not refresh:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".download")
    urllib.request.urlretrieve(url, temporary)
    temporary.replace(path)


def parse_nea_catalog(path: Path) -> list[dict]:
    records = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 11:
            continue
        try:
            omega = float(parts[5])
            node = float(parts[6])
            inclination = float(parts[7])
            eccentricity = float(parts[8])
            semi_major_axis = float(parts[10])
        except ValueError:
            continue
        packed_designation = parts[0]
        name = line[166:194].strip() if len(line) >= 194 else packed_designation
        records.append(
            {
                "designation": packed_designation,
                "name": name or packed_designation,
                "epoch": parts[3],
                "a_au": semi_major_axis,
                "q_au": semi_major_axis * (1.0 - eccentricity),
                "e": eccentricity,
                "i_deg": inclination,
                "omega_deg": omega,
                "node_deg": node,
            }
        )
    return records


def parse_comet_catalog(path: Path) -> list[dict]:
    records = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 12:
            continue
        try:
            q_au = float(parts[4])
            eccentricity = float(parts[5])
            omega = float(parts[6])
            node = float(parts[7])
            inclination = float(parts[8])
        except ValueError:
            continue
        designation = parts[0]
        name = line[101:157].strip() if len(line) >= 157 else designation
        semi_major_axis = (
            q_au / (1.0 - eccentricity) if eccentricity < 1.0 else np.nan
        )
        records.append(
            {
                "designation": designation,
                "name": name or designation,
                "epoch": parts[9],
                "a_au": semi_major_axis,
                "q_au": q_au,
                "e": eccentricity,
                "i_deg": inclination,
                "omega_deg": omega,
                "node_deg": node,
            }
        )
    return records


def rank_catalog(reference: Orbit, records: list[dict]) -> list[dict]:
    for record in records:
        candidate = Orbit(
            q_au=record["q_au"],
            e=record["e"],
            i_deg=record["i_deg"],
            omega_deg=record["omega_deg"],
            node_deg=record["node_deg"],
        )
        record["D_H"] = jopek_dh(reference, candidate)
    return sorted(records, key=lambda record: record["D_H"])


def write_ranked_group(group: h5py.Group, records: list[dict], top: int) -> None:
    records = records[:top]
    text_type = h5py.string_dtype("utf-8")
    group.create_dataset(
        "rank",
        data=np.arange(1, len(records) + 1, dtype=np.int32),
    )
    for field in ("designation", "name", "epoch"):
        group.create_dataset(
            field,
            data=np.asarray([record[field] for record in records], dtype=text_type),
        )
    for field in (
        "D_H",
        "a_au",
        "q_au",
        "e",
        "i_deg",
        "omega_deg",
        "node_deg",
    ):
        group.create_dataset(
            field,
            data=np.asarray([record[field] for record in records], dtype=np.float32),
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q-au", type=float, default=0.899)
    parser.add_argument("--eccentricity", type=float, default=0.712)
    parser.add_argument("--inclination-deg", type=float, default=91.22)
    parser.add_argument("--omega-deg", type=float, default=317.12)
    parser.add_argument("--node-deg", type=float, default=290.12)
    parser.add_argument("--nea-catalog", type=Path, default=DEFAULT_CACHE / "NEA.txt")
    parser.add_argument(
        "--comet-catalog",
        type=Path,
        default=DEFAULT_CACHE / "CometEls.txt",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--top", type=int, default=25)
    parser.add_argument("--refresh", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    download_catalog(args.nea_catalog, NEA_URL, args.refresh)
    download_catalog(args.comet_catalog, COMET_URL, args.refresh)
    reference = Orbit(
        q_au=args.q_au,
        e=args.eccentricity,
        i_deg=args.inclination_deg,
        omega_deg=args.omega_deg,
        node_deg=args.node_deg,
    )
    nea = rank_catalog(reference, parse_nea_catalog(args.nea_catalog))
    comets = rank_catalog(reference, parse_comet_catalog(args.comet_catalog))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output, "w") as h5:
        h5.attrs["schema"] = "small_body_jopek_matches_v1"
        h5.attrs["created_utc"] = datetime.now(timezone.utc).isoformat()
        h5.attrs["criterion"] = "Jopek (1993) D_H"
        h5.attrs["nea_catalog_url"] = NEA_URL
        h5.attrs["comet_catalog_url"] = COMET_URL
        for field in ("q_au", "e", "i_deg", "omega_deg", "node_deg"):
            h5.attrs[f"reference_{field}"] = np.float32(getattr(reference, field))
        write_ranked_group(h5.create_group("nea"), nea, args.top)
        write_ranked_group(h5.create_group("comet"), comets, args.top)

    print(
        f"NEA: {nea[0]['name']} (D_H={nea[0]['D_H']:.6f}); "
        f"comet: {comets[0]['name']} (D_H={comets[0]['D_H']:.6f})"
    )
    print(args.output)


if __name__ == "__main__":
    main()
