#!/usr/bin/env python3
"""Generate PANSY two-way gain pattern figures for the paper memo."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import pansy_gain as pgain


BEAM_NAMES = ["Zenith", "North", "East", "South", "West"]


def direction_grid(extent_deg: float, grid_n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    east_deg = np.linspace(-extent_deg, extent_deg, grid_n)
    north_deg = np.linspace(-extent_deg, extent_deg, grid_n)
    east_grid, north_grid = np.meshgrid(east_deg, north_deg)
    u = np.sin(np.deg2rad(east_grid))
    v = np.sin(np.deg2rad(north_grid))
    rho2 = u**2 + v**2
    valid = rho2 <= 1.0
    w = np.full_like(u, np.nan, dtype=np.float64)
    w[valid] = -np.sqrt(np.maximum(0.0, 1.0 - rho2[valid]))
    return east_grid, north_grid, u, v, w


def relative_db(power: np.ndarray) -> np.ndarray:
    db = pgain.power_to_db(power)
    finite = np.isfinite(db)
    if np.any(finite):
        db = db - np.nanmax(db)
    return db


def format_axis(ax, extent_deg: float) -> None:
    ax.set_xlim(-extent_deg, extent_deg)
    ax.set_ylim(-extent_deg, extent_deg)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("East of zenith (deg)")
    ax.set_ylabel("North of zenith (deg)")
    ax.grid(True, color="white", alpha=0.18, linewidth=0.5)


def add_gain_panel(fig, ax, east, north, gain_db, title, extent_deg, vmin=-35.0):
    im = ax.pcolormesh(east, north, gain_db, shading="auto", cmap="viridis", vmin=vmin, vmax=0.0)
    ax.set_title(title)
    format_axis(ax, extent_deg)
    return im


def plot_zenith_components(output: Path, extent_deg: float, grid_n: int) -> None:
    east, north, u, v, w = direction_grid(extent_deg, grid_n)
    valid = np.isfinite(w)
    uvw = np.column_stack([u[valid], v[valid], w[valid]])
    beam_vecs = pgain.tx_beam_unit_vectors()
    steer = beam_vecs[0]
    modules = pgain.tx_module_positions()
    tx_rel_module = modules[0] - np.mean(modules[0], axis=0, keepdims=True)
    module_centers = pgain.tx_module_center_positions(modules)

    panels: list[tuple[str, np.ndarray]] = []
    tx_module = np.full(u.shape, np.nan, dtype=np.float64)
    tx_sparse = np.full(u.shape, np.nan, dtype=np.float64)
    rx_module = np.full(u.shape, np.nan, dtype=np.float64)
    two_way = np.full(u.shape, np.nan, dtype=np.float64)

    tx_module[valid] = pgain.module_power_gain(uvw, steer=steer, module_positions=tx_rel_module)
    tx_sparse[valid] = pgain.sparse_tx_array_power_gain(uvw, 0, beam_vecs=beam_vecs, module_centers=module_centers)
    rx_module[valid] = pgain.rx_power_gain(uvw, steer=steer)
    two_way[valid] = pgain.two_way_power_gain(uvw, 0, beam_vecs=beam_vecs)

    panels.append(("TX module", relative_db(tx_module)))
    panels.append(("TX sparse array", relative_db(tx_sparse)))
    panels.append(("RX module", relative_db(rx_module)))
    panels.append(("Two-way gain", relative_db(two_way)))

    output.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({"font.size": 13, "axes.titlesize": 15, "axes.labelsize": 13}):
        fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.0), constrained_layout=True)
        im = None
        for ax, (title, gain_db) in zip(axes.ravel(), panels):
            im = add_gain_panel(fig, ax, east, north, gain_db, title, extent_deg)
        fig.colorbar(im, ax=axes.ravel().tolist(), label="Relative power gain (dB)", shrink=0.88)
        fig.suptitle("PANSY zenith beam gain components")
        fig.savefig(output, dpi=220)
        plt.close(fig)


def plot_five_beam_two_way(output: Path, extent_deg: float, grid_n: int) -> None:
    east, north, u, v, w = direction_grid(extent_deg, grid_n)
    valid = np.isfinite(w)
    gain_maps = pgain.precompute_two_way_gain_maps(u, v, w, valid)
    for beam_i in range(gain_maps.shape[0]):
        finite = np.isfinite(gain_maps[beam_i])
        if np.any(finite):
            gain_maps[beam_i] -= np.nanmax(gain_maps[beam_i])
    linear_sum = np.nansum(10.0 ** (gain_maps / 10.0), axis=0)
    composite = np.full_like(linear_sum, np.nan, dtype=np.float64)
    good = np.isfinite(linear_sum) & (linear_sum > 0)
    composite[good] = 10.0 * np.log10(linear_sum[good])
    composite -= np.nanmax(composite)

    output.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({"font.size": 12, "axes.titlesize": 14, "axes.labelsize": 12}):
        fig, axes = plt.subplots(2, 3, figsize=(13.6, 8.4), constrained_layout=True)
        im = None
        for beam_i, ax in enumerate(axes.ravel()[:5]):
            im = add_gain_panel(fig, ax, east, north, gain_maps[beam_i], BEAM_NAMES[beam_i], extent_deg)
        im = add_gain_panel(fig, axes.ravel()[5], east, north, composite, "Five-beam composite", extent_deg)
        fig.colorbar(im, ax=axes.ravel().tolist(), label="Relative two-way power gain (dB)", shrink=0.88)
        fig.suptitle("PANSY two-way antenna gain patterns")
        fig.savefig(output, dpi=220)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/Users/jvi019/src/pansy_paper/memos/figures/pansy_gain_patterns"),
    )
    parser.add_argument("--extent-deg", type=float, default=20.0)
    parser.add_argument("--grid-n", type=int, default=501)
    args = parser.parse_args()

    plot_zenith_components(args.output_dir / "pansy_zenith_gain_components.png", args.extent_deg, args.grid_n)
    plot_five_beam_two_way(args.output_dir / "pansy_two_way_gain_five_beams.png", args.extent_deg, args.grid_n)


if __name__ == "__main__":
    main()
