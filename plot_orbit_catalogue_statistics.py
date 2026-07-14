#!/usr/bin/env python3
"""Make paper-ready catalogue statistics plots from compact orbit metadata."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


def wrap360(deg: np.ndarray) -> np.ndarray:
    return np.asarray(deg, dtype=np.float64) % 360.0


def iter_orbit_event_tables(orbit_metadata_dir: Path):
    for path in sorted(orbit_metadata_dir.glob("**/orbit@*.h5")):
        try:
            with h5py.File(path, "r") as h:
                if "events" not in h:
                    continue
                yield path, h["events"][()]
        except OSError:
            continue


def collect_events(orbit_metadata_dir: Path) -> tuple[np.ndarray, int]:
    chunks = []
    files_read = 0
    for _path, events in iter_orbit_event_tables(orbit_metadata_dir):
        chunks.append(events)
        files_read += 1
    if not chunks:
        return np.zeros(0, dtype=[]), files_read
    events = np.concatenate(chunks)
    if "sample_idx" in events.dtype.names:
        events = events[np.argsort(events["sample_idx"])]
    return events, files_read


def finite_field(events: np.ndarray, name: str) -> np.ndarray:
    if events.dtype.names is None or name not in events.dtype.names:
        return np.full(len(events), np.nan, dtype=np.float32)
    values = np.asarray(events[name], dtype=np.float64)
    out = np.full(values.shape, np.nan, dtype=np.float32)
    finite = np.isfinite(values)
    finite &= np.abs(values) <= np.finfo(np.float32).max
    out[finite] = values[finite].astype(np.float32)
    return out


def quality_mask(events: np.ndarray, max_radiant_sigma_deg: float | None) -> np.ndarray:
    good = np.ones(len(events), dtype=bool)
    if events.dtype.names is None or "sample_idx" not in events.dtype.names:
        return np.zeros(len(events), dtype=bool)
    good &= np.isfinite(np.asarray(events["sample_idx"], dtype=np.float64))
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


def catalogue_arrays(events: np.ndarray, max_radiant_sigma_deg: float | None) -> dict[str, np.ndarray]:
    keep = quality_mask(events, max_radiant_sigma_deg)
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


def histogram_solar_longitude(solar_longitude_deg: np.ndarray, bin_width_deg: float) -> tuple[np.ndarray, np.ndarray]:
    edges = np.arange(0.0, 360.0 + float(bin_width_deg), float(bin_width_deg), dtype=np.float32)
    counts, edges = np.histogram(solar_longitude_deg, bins=edges)
    return counts.astype(np.int32), edges.astype(np.float32)


def histogram_height_velocity(height_km: np.ndarray, speed_km_s: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    height_edges = np.arange(70.0, 140.0 + 1.0, 1.0, dtype=np.float32)
    speed_edges = np.arange(0.0, 80.0 + 1.0, 1.0, dtype=np.float32)
    hist, h_edges, v_edges = np.histogram2d(height_km, speed_km_s, bins=[height_edges, speed_edges])
    return hist.astype(np.int32), h_edges.astype(np.float32), v_edges.astype(np.float32)


def write_statistics_h5(path: Path, arrays: dict[str, np.ndarray], solar_counts, solar_edges, hv_counts, height_edges, speed_edges, files_read: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h:
        h.attrs["script"] = Path(__file__).name
        h.attrs["source"] = "compact orbit metadata events tables"
        h.attrs["files_read"] = int(files_read)
        h.attrs["event_count"] = int(len(arrays["sample_idx"]))
        for name, values in arrays.items():
            h.create_dataset(name, data=values, compression="gzip", shuffle=True)
        h.create_dataset("solar_longitude_count", data=solar_counts, compression="gzip", shuffle=True)
        h.create_dataset("solar_longitude_edges_deg", data=solar_edges)
        h.create_dataset("height_velocity_count", data=hv_counts, compression="gzip", shuffle=True)
        h.create_dataset("height_edges_km", data=height_edges)
        h.create_dataset("speed_edges_km_s", data=speed_edges)


def plot_solar_counts(path: Path, counts: np.ndarray, edges: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    fig, ax = plt.subplots(figsize=(8.0, 3.4), constrained_layout=True)
    ax.step(centers, counts, where="mid", color="#2f5f8f", linewidth=1.6)
    ax.fill_between(centers, counts, step="mid", color="#2f5f8f", alpha=0.24)
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"Solar longitude, $\lambda_\odot$ (deg)")
    ax.set_ylabel("Catalogue meteor count")
    ax.grid(True, alpha=0.25)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_height_velocity(path: Path, hv_counts: np.ndarray, height_edges: np.ndarray, speed_edges: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plot_count = np.asarray(hv_counts, dtype=np.float64).T
    plot_count[plot_count <= 0.0] = np.nan
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("white")
    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    mesh = ax.pcolormesh(
        height_edges,
        speed_edges,
        plot_count,
        shading="auto",
        cmap=cmap,
        norm=LogNorm(vmin=1.0, vmax=max(1.0, float(np.nanmax(plot_count)) if np.any(np.isfinite(plot_count)) else 1.0)),
    )
    ax.set_xlabel("Initial detection height (km)")
    ax.set_ylabel(r"Geocentric velocity, $v_g$ (km s$^{-1}$)")
    ax.set_xlim(float(height_edges[0]), float(height_edges[-1]))
    ax.set_ylim(float(speed_edges[0]), float(speed_edges[-1]))
    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label("Catalogue meteor count")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-metadata-dir", type=Path, default=Path("data/metadata/orbit"))
    parser.add_argument("--output-h5", type=Path, default=Path("paper_orbit_catalogue_statistics.h5"))
    parser.add_argument("--solar-output", type=Path, default=Path("paper_counts_solar_longitude.png"))
    parser.add_argument("--height-velocity-output", type=Path, default=Path("paper_height_velocity.png"))
    parser.add_argument("--solar-bin-width-deg", type=float, default=1.0)
    parser.add_argument("--max-radiant-sigma-deg", type=float, default=None)
    args = parser.parse_args()

    events, files_read = collect_events(args.orbit_metadata_dir)
    arrays = catalogue_arrays(events, args.max_radiant_sigma_deg)
    solar_counts, solar_edges = histogram_solar_longitude(arrays["solar_longitude_deg"], args.solar_bin_width_deg)
    hv_counts, height_edges, speed_edges = histogram_height_velocity(arrays["initial_detection_height_km"], arrays["v_g_km_s"])
    write_statistics_h5(args.output_h5, arrays, solar_counts, solar_edges, hv_counts, height_edges, speed_edges, files_read)
    plot_solar_counts(args.solar_output, solar_counts, solar_edges)
    plot_height_velocity(args.height_velocity_output, hv_counts, height_edges, speed_edges)
    print(f"orbit_catalogue_statistics events={len(arrays['sample_idx'])} files_read={files_read}")
    print(args.output_h5)
    print(args.solar_output)
    print(args.height_velocity_output)


if __name__ == "__main__":
    main()
