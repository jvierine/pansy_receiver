#!/usr/bin/env python3
"""Plot meteor angular-position histograms with the PANSY zenith TX gain pattern."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

import plot_interferometric_disambiguation as disamb


def direction_cosines_to_offsets_deg(uvw: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Convert PANSY direction cosines to east/north angular offsets from zenith."""
    uvw = np.asarray(uvw, dtype=np.float64)
    down = -uvw[:, 2]
    east_deg = np.rad2deg(np.arctan2(uvw[:, 0], down))
    north_deg = np.rad2deg(np.arctan2(uvw[:, 1], down))
    return east_deg, north_deg


def interferometry_grid_offsets_deg(grid_n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    u, v, w, valid = disamb.horizon_grid(grid_n)
    uvw = np.column_stack([u.ravel(), v.ravel(), w.ravel()])
    east, north = direction_cosines_to_offsets_deg(uvw)
    return u, v, east.reshape(u.shape), north.reshape(u.shape), valid


def read_selected_grid_cells(path: Path, use_selection_keep: bool = True) -> dict | None:
    try:
        with h5py.File(path, "r") as handle:
            label = str(handle.attrs.get("selected_hypothesis", ""))
            if not label or "hypotheses" not in handle or label not in handle["hypotheses"]:
                return None
            group = handle["hypotheses"][label]
            if bool(group.attrs.get("combined_reject", False)):
                return None
            if "candidate_indices" not in group or "beam_id" not in group or "candidates" not in handle:
                return None
            cand_group = handle["candidates"]
            if "grid_row" not in cand_group or "grid_col" not in cand_group:
                return None
            candidate_indices = np.asarray(group["candidate_indices"], dtype=np.int64)
            beam_id = np.asarray(group["beam_id"], dtype=np.int64)
            n = min(len(candidate_indices), len(beam_id))
            candidate_indices = candidate_indices[:n]
            beam_id = beam_id[:n]
            keep = np.ones(n, dtype=bool)
            if use_selection_keep and "selection_keep" in group:
                keep &= np.asarray(group["selection_keep"], dtype=bool)[:n]
            keep &= np.isfinite(beam_id)
            grid_row_all = np.asarray(cand_group["grid_row"], dtype=np.int64)
            grid_col_all = np.asarray(cand_group["grid_col"], dtype=np.int64)
            valid_idx = (candidate_indices >= 0) & (candidate_indices < len(grid_row_all))
            keep &= valid_idx
            selected_idx = candidate_indices[keep]
            return {
                "sample_idx": int(handle.attrs.get("sample_idx", -1)),
                "hypothesis": label,
                "grid_row": grid_row_all[selected_idx],
                "grid_col": grid_col_all[selected_idx],
                "beam_id": beam_id[keep],
            }
    except OSError:
        return None


def collect_grid_counts(diagnostics_dir: Path, grid_n: int, use_selection_keep: bool = True) -> dict:
    files = sorted(diagnostics_dir.glob("pansy_disambiguation_diagnostics_*.h5"))
    zenith_counts = np.zeros((grid_n, grid_n), dtype=np.int64)
    all_counts = np.zeros((grid_n, grid_n), dtype=np.int64)
    event_count = 0
    skipped = 0
    total_points = 0
    for path in files:
        row = read_selected_grid_cells(path, use_selection_keep=use_selection_keep)
        if row is None or len(row["beam_id"]) == 0:
            skipped += 1
            continue
        grid_row = np.asarray(row["grid_row"], dtype=np.int64)
        grid_col = np.asarray(row["grid_col"], dtype=np.int64)
        beam_id = np.asarray(row["beam_id"], dtype=np.int64)
        valid = (grid_row >= 0) & (grid_row < grid_n) & (grid_col >= 0) & (grid_col < grid_n)
        if not np.any(valid):
            skipped += 1
            continue
        event_count += 1
        grid_row = grid_row[valid]
        grid_col = grid_col[valid]
        beam_id = beam_id[valid]
        np.add.at(all_counts, (grid_row, grid_col), 1)
        zenith = beam_id == 0
        np.add.at(zenith_counts, (grid_row[zenith], grid_col[zenith]), 1)
        total_points += len(grid_row)
    return {
        "zenith_counts": zenith_counts,
        "all_counts": all_counts,
        "event_count": event_count,
        "skipped_count": skipped,
        "file_count": len(files),
        "total_points": total_points,
    }


def zenith_tx_gain_grid(u: np.ndarray, v: np.ndarray, valid: np.ndarray) -> np.ndarray:
    w = np.full_like(u, np.nan, dtype=np.float64)
    w[valid] = -np.sqrt(np.maximum(0.0, 1.0 - u[valid] ** 2 - v[valid] ** 2))
    uvw = np.column_stack([u.ravel(), v.ravel(), w.ravel()])
    beam_id = np.zeros(len(uvw), dtype=np.int64)
    gain = disamb.tx_array_gain_db(uvw, beam_id, model="module_center_coherent").reshape(u.shape)
    gain -= np.nanmax(gain)
    return gain


def add_histogram_panel(ax, east_deg, north_deg, counts, title, cmap="magma"):
    image = np.asarray(counts, dtype=np.float64)
    image = np.where(np.isfinite(image), np.maximum(image, 1.0), np.nan)
    mesh = ax.pcolormesh(
        east_deg,
        north_deg,
        image,
        shading="auto",
        cmap=cmap,
        norm=LogNorm(vmin=1, vmax=max(1, float(np.nanmax(image)))),
    )
    ax.set_title(title)
    return mesh, int(np.sum(image))


def format_direction_axis(ax, extent_deg):
    ax.set_xlim(-extent_deg, extent_deg)
    ax.set_ylim(-extent_deg, extent_deg)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("East of zenith (deg)")
    ax.set_ylabel("North of zenith (deg)")
    ax.grid(True, color="white", alpha=0.18, linewidth=0.5)


def plot_position_histograms(
    grid_counts: dict,
    output: Path,
    extent_deg: float,
    grid_n: int,
):
    output.parent.mkdir(parents=True, exist_ok=True)
    u, v, east_grid, north_grid, valid_grid = interferometry_grid_offsets_deg(grid_n)
    in_view = (
        np.isfinite(east_grid)
        & np.isfinite(north_grid)
        & (east_grid >= -extent_deg)
        & (east_grid <= extent_deg)
        & (north_grid >= -extent_deg)
        & (north_grid <= extent_deg)
    )
    zenith_counts = np.where(in_view, grid_counts["zenith_counts"], 0)
    all_counts = np.where(in_view, grid_counts["all_counts"], 0)
    n_zenith = int(np.sum(zenith_counts))
    n_all = int(np.sum(all_counts))
    row_use = np.flatnonzero(np.any(in_view, axis=1))
    col_use = np.flatnonzero(np.any(in_view, axis=0))
    if len(row_use) == 0 or len(col_use) == 0:
        raise RuntimeError("no interferometry grid cells inside requested angular extent")
    rs = slice(int(row_use[0]), int(row_use[-1]) + 1)
    cs = slice(int(col_use[0]), int(col_use[-1]) + 1)
    east_plot = east_grid[rs, cs]
    north_plot = north_grid[rs, cs]
    zenith_plot = zenith_counts[rs, cs]
    all_plot = all_counts[rs, cs]

    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.9), constrained_layout=True)
    im0, _ = add_histogram_panel(
        axes[0],
        east_plot,
        north_plot,
        zenith_plot,
        f"Zenith TX meteor positions\nbeam 0, N={n_zenith} points",
        cmap="magma",
    )
    cb0 = fig.colorbar(im0, ax=axes[0], shrink=0.82)
    cb0.set_label("Meteor detections per bin")
    format_direction_axis(axes[0], extent_deg)

    gain = zenith_tx_gain_grid(u, v, valid_grid)
    im1 = axes[1].pcolormesh(east_plot, north_plot, gain[rs, cs], shading="auto", cmap="viridis", vmin=-30.0, vmax=0.0)
    axes[1].set_title("Zenith TX gain pattern\nrelative array factor (dB)")
    cb1 = fig.colorbar(im1, ax=axes[1], shrink=0.82)
    cb1.set_label("Relative gain (dB)")
    format_direction_axis(axes[1], extent_deg)

    im2, _ = add_histogram_panel(
        axes[2],
        east_plot,
        north_plot,
        all_plot,
        f"All TX beams meteor positions\nN={n_all} points",
        cmap="magma",
    )
    cb2 = fig.colorbar(im2, ax=axes[2], shrink=0.82)
    cb2.set_label("Meteor detections per bin")
    format_direction_axis(axes[2], extent_deg)

    fig.savefig(output, dpi=220)
    plt.close(fig)
    return {
        "zenith_points_in_view": n_zenith,
        "all_points_in_view": n_all,
        "total_points": int(grid_counts["total_points"]),
    }


def plot_stacked_position_histogram(input_h5: Path, output: Path, extent_deg: float):
    """Plot a current all-days position histogram from a daily histogram stack."""
    with h5py.File(input_h5, "r") as h:
        counts = np.asarray(h["all_counts"], dtype=np.int64)
        east_grid = np.asarray(h["east_grid_deg"], dtype=np.float64)
        north_grid = np.asarray(h["north_grid_deg"], dtype=np.float64)
        metrics = h["metrics"][()] if "metrics" in h else np.zeros(0, dtype=[])
    if counts.ndim != 3:
        raise ValueError(f"{input_h5} all_counts must have shape (day, row, col)")
    total_counts = np.sum(counts, axis=0)
    in_view = (
        np.isfinite(east_grid)
        & np.isfinite(north_grid)
        & (east_grid >= -extent_deg)
        & (east_grid <= extent_deg)
        & (north_grid >= -extent_deg)
        & (north_grid <= extent_deg)
    )
    total_counts = np.where(in_view, total_counts, 0)
    row_use = np.flatnonzero(np.any(in_view, axis=1))
    col_use = np.flatnonzero(np.any(in_view, axis=0))
    if len(row_use) == 0 or len(col_use) == 0:
        raise RuntimeError("no interferometry grid cells inside requested angular extent")
    rs = slice(int(row_use[0]), int(row_use[-1]) + 1)
    cs = slice(int(col_use[0]), int(col_use[-1]) + 1)
    n_points = int(np.sum(total_counts))
    n_days = int(len(metrics)) if len(metrics) else int(counts.shape[0])

    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    image = np.where(np.isfinite(total_counts[rs, cs]), np.maximum(total_counts[rs, cs], 1), np.nan)
    mesh = ax.pcolormesh(
        east_grid[rs, cs],
        north_grid[rs, cs],
        image,
        shading="auto",
        cmap="magma",
        norm=LogNorm(vmin=1, vmax=max(1, float(np.nanmax(image)))),
    )
    ax.set_title(f"All TX beams meteor positions\n{n_days} days, N={n_points} points")
    format_direction_axis(ax, extent_deg)
    cb = fig.colorbar(mesh, ax=ax, shrink=0.86)
    cb.set_label("Meteor detections per bin")
    fig.savefig(output, dpi=220)
    plt.close(fig)
    return {"days": n_days, "points_in_view": n_points}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics-dir", type=Path, default=Path("test_plots"))
    parser.add_argument("--input-h5", type=Path, default=None, help="Daily histogram stack HDF5 to plot instead of scanning diagnostics.")
    parser.add_argument("--output", type=Path, default=Path("test_plots/pansy_meteor_position_histograms.png"))
    parser.add_argument("--extent-deg", type=float, default=15.0)
    parser.add_argument("--grid-n", type=int, default=501, help="Interferometry search grid size used to form grid_row/grid_col.")
    parser.add_argument(
        "--include-clipped",
        action="store_true",
        help="Use all selected-hypothesis points instead of only trajectory-fit kept points.",
    )
    args = parser.parse_args()

    if args.input_h5 is not None:
        stats = plot_stacked_position_histogram(args.input_h5, args.output, args.extent_deg)
        print(args.output)
        print(f"daily_histogram_days {stats['days']}")
        print(f"points_in_view {stats['points_in_view']}")
        return

    grid_counts = collect_grid_counts(args.diagnostics_dir, args.grid_n, use_selection_keep=not args.include_clipped)
    stats = plot_position_histograms(grid_counts, args.output, args.extent_deg, args.grid_n)
    print(args.output)
    print(f"diagnostic_files {grid_counts['file_count']}")
    print(f"fitted_events {grid_counts['event_count']}")
    print(f"skipped_events {grid_counts['skipped_count']}")
    print(f"total_points {stats['total_points']}")
    print(f"zenith_points_in_view {stats['zenith_points_in_view']}")
    print(f"all_points_in_view {stats['all_points_in_view']}")


if __name__ == "__main__":
    main()
