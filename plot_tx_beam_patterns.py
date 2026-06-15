#!/usr/bin/env python3
"""Plot PANSY transmit beam patterns from antenna module centers."""

from __future__ import annotations

from pathlib import Path

import scipy.constants as c
import matplotlib.pyplot as plt
import numpy as np

import pansy_config as pc


def module_centers() -> np.ndarray:
    """Use one coherent array element per PANSY antenna module center."""
    centers = []
    for name, pos in sorted(pc.module_center.items()):
        if name == "RFTX":
            continue
        centers.append(np.asarray(pos, dtype=np.float64))
    centers = np.asarray(centers, dtype=np.float64)
    centers -= np.mean(centers, axis=0, keepdims=True)
    return centers


def direction_from_zenith_angles(ew_deg, ns_deg):
    """Local direction vector using EW/NS zenith-angle coordinates."""
    u = np.sin(np.deg2rad(ew_deg))
    v = np.sin(np.deg2rad(ns_deg))
    rho2 = u**2 + v**2
    valid = rho2 <= 1.0
    w = np.full_like(u, np.nan, dtype=np.float64)
    w[valid] = -np.sqrt(1.0 - rho2[valid])
    return u, v, w, valid


def module_array_power(module_pos, u, v, w, valid):
    """Coherent module-center array factor for zenith, north, east, south, west."""
    wavelength = c.c / pc.freq
    k0 = 2.0 * np.pi / wavelength
    beam_ew = np.array([0.0, 0.0, 10.0, 0.0, -10.0])
    beam_ns = np.array([0.0, 10.0, 0.0, -10.0, 0.0])
    bu, bv, bw, _ = direction_from_zenith_angles(beam_ew, beam_ns)
    beam_dirs = np.column_stack([bu, bv, bw])

    uvw = np.column_stack([u[valid], v[valid], w[valid]])
    power = np.zeros((5,) + u.shape, dtype=np.float64)
    for beam_i, steer in enumerate(beam_dirs):
        phase = k0 * (module_pos @ (uvw - steer).T)
        af = np.abs(np.mean(np.exp(1j * phase), axis=0)) ** 2
        panel = np.full(u.shape, np.nan, dtype=np.float64)
        panel[valid] = af / np.nanmax(af)
        power[beam_i] = panel
    return power, beam_ew, beam_ns


def main() -> int:
    out = Path("/Users/jvi019/src/pansy_paper/memos/figures/pansy_tx_beam_patterns_zenith_angle.png")
    out_panels = Path("/Users/jvi019/src/pansy_paper/memos/figures/pansy_tx_beam_patterns_zenith_angle_panels.png")
    out_single_dir = Path("/Users/jvi019/src/pansy_paper/memos/figures/pansy_tx_beam_patterns")
    out.parent.mkdir(parents=True, exist_ok=True)
    out_single_dir.mkdir(parents=True, exist_ok=True)

    lim = 20.0
    n_grid = 501
    ew_deg = np.linspace(-lim, lim, n_grid)
    ns_deg = np.linspace(-lim, lim, n_grid)
    ew, ns = np.meshgrid(ew_deg, ns_deg)
    u, v, w, valid = direction_from_zenith_angles(ew, ns)
    modules = module_centers()
    power, beam_ew, beam_ns = module_array_power(modules, u, v, w, valid)

    summed = np.sum(power, axis=0)
    summed /= np.nanmax(summed)

    fig, ax = plt.subplots(figsize=(8.0, 6.4), constrained_layout=True)
    im = ax.pcolormesh(ew, ns, summed, shading="auto", cmap="turbo", vmin=0.0, vmax=1.0)
    ax.contour(ew, ns, summed, levels=[0.1, 0.25, 0.5, 0.75], colors="white", linewidths=0.6, alpha=0.8)
    ax.plot(beam_ew, beam_ns, "rx", ms=8, mew=1.5)
    for i, (x, y) in enumerate(zip(beam_ew, beam_ns)):
        ax.text(x + 0.5, y + 0.5, str(i), color="white", fontsize=9, weight="bold")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("East-West zenith angle (deg)")
    ax.set_ylabel("North-South zenith angle (deg)")
    ax.set_title("PANSY transmit beam patterns from module centers, five-beam sum")
    fig.colorbar(im, ax=ax, label="Normalized linear transmit power")
    fig.savefig(out, dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(3, 3, figsize=(10.5, 10.0), constrained_layout=True)
    beam_names = ["Zenith", "North", "East", "South", "West"]
    panel_locs = {
        1: (0, 1),
        4: (1, 0),
        0: (1, 1),
        2: (1, 2),
        3: (2, 1),
    }
    for ax in axes.ravel():
        ax.axis("off")
    im = None
    for i, (row, col) in panel_locs.items():
        ax = axes[row, col]
        ax.axis("on")
        panel = power[i] / np.nanmax(power[i])
        im = ax.pcolormesh(ew, ns, panel, shading="auto", cmap="turbo", vmin=0.0, vmax=1.0)
        ax.contour(ew, ns, panel, levels=[0.1, 0.25, 0.5, 0.75], colors="white", linewidths=0.5, alpha=0.75)
        ax.plot(beam_ew[i], beam_ns[i], "rx", ms=7, mew=1.4)
        ax.set_title(f"{i}: {beam_names[i]}")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel("EW ZA (deg)")
        if col == 0:
            ax.set_ylabel("NS ZA (deg)")
    fig.colorbar(im, ax=axes.ravel().tolist(), label="Beam-normalized linear transmit power")
    fig.suptitle("PANSY transmit beam patterns from module centers, individual beams")
    fig.savefig(out_panels, dpi=220)
    plt.close(fig)

    # Full visible hemisphere in direction-cosine coordinates for each pointing.
    dc_n = 701
    axis = np.linspace(-1.0, 1.0, dc_n)
    dc_u, dc_v = np.meshgrid(axis, axis)
    dc_rho2 = dc_u**2 + dc_v**2
    dc_valid = dc_rho2 <= 1.0
    dc_w = np.full_like(dc_u, np.nan)
    dc_w[dc_valid] = -np.sqrt(1.0 - dc_rho2[dc_valid])
    dc_power, beam_ew_dc, beam_ns_dc = module_array_power(modules, dc_u, dc_v, dc_w, dc_valid)
    beam_u = np.sin(np.deg2rad(beam_ew_dc))
    beam_v = np.sin(np.deg2rad(beam_ns_dc))
    for i, name in enumerate(beam_names):
        panel = dc_power[i] / np.nanmax(dc_power[i])
        fig, ax = plt.subplots(figsize=(7.0, 6.4), constrained_layout=True)
        im = ax.pcolormesh(dc_u, dc_v, panel, shading="auto", cmap="turbo", vmin=0.0, vmax=1.0)
        ax.contour(dc_u, dc_v, panel, levels=[0.1, 0.25, 0.5, 0.75], colors="white", linewidths=0.6, alpha=0.75)
        ax.add_patch(plt.Circle((0.0, 0.0), 1.0, color="black", fill=False, linewidth=1.0))
        ax.plot(beam_u[i], beam_v[i], "rx", ms=8, mew=1.5)
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(-1.0, 1.0)
        ax.set_aspect("equal")
        ax.set_xlabel("u direction cosine")
        ax.set_ylabel("v direction cosine")
        ax.set_title(f"PANSY TX beam {i}: {name}, full horizon")
        fig.colorbar(im, ax=ax, label="Beam-normalized linear transmit power")
        single_out = out_single_dir / f"pansy_tx_beam_{i}_{name.lower()}_zenith_angle.png"
        fig.savefig(single_out, dpi=220)
        plt.close(fig)
        print(single_out)

        fig, ax = plt.subplots(figsize=(7.0, 6.4), constrained_layout=True)
        im = ax.pcolormesh(dc_u, dc_v, panel, shading="auto", cmap="turbo", vmin=0.0, vmax=1.0)
        ax.contour(dc_u, dc_v, panel, levels=[0.1, 0.25, 0.5, 0.75], colors="white", linewidths=0.6, alpha=0.75)
        ax.add_patch(plt.Circle((0.0, 0.0), 1.0, color="black", fill=False, linewidth=1.0))
        ax.plot(beam_u[i], beam_v[i], "rx", ms=8, mew=1.5)
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(-1.0, 1.0)
        ax.set_aspect("equal")
        ax.set_xlabel("u direction cosine")
        ax.set_ylabel("v direction cosine")
        ax.set_title(f"PANSY TX beam {i}: {name}, full horizon")
        fig.colorbar(im, ax=ax, label="Beam-normalized linear transmit power")
        dc_out = out_single_dir / f"pansy_tx_beam_{i}_{name.lower()}_direction_cosine.png"
        fig.savefig(dc_out, dpi=220)
        plt.close(fig)
        print(dc_out)
    print(out)
    print(out_panels)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
