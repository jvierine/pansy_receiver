#!/usr/bin/env python3
"""Plot DASST-corrected PANSY orbit samples in heliocentric ecliptic x-y."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


PLANETS = [
    ("Mercury", 0.387, "#9a9a9a"),
    ("Venus", 0.723, "#b08a5b"),
    ("Earth", 1.000, "#4f8fd6"),
    ("Mars", 1.524, "#c46852"),
    ("Jupiter", 5.203, "#8d7658"),
]


def rotation_matrix_313(node_deg: float, inc_deg: float, argp_deg: float) -> np.ndarray:
    node = np.deg2rad(node_deg)
    inc = np.deg2rad(inc_deg)
    argp = np.deg2rad(argp_deg)
    cn, sn = np.cos(node), np.sin(node)
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(argp), np.sin(argp)
    rz_node = np.array([[cn, -sn, 0.0], [sn, cn, 0.0], [0.0, 0.0, 1.0]])
    rx_inc = np.array([[1.0, 0.0, 0.0], [0.0, ci, -si], [0.0, si, ci]])
    rz_argp = np.array([[cw, -sw, 0.0], [sw, cw, 0.0], [0.0, 0.0, 1.0]])
    return rz_node @ rx_inc @ rz_argp


def orbit_xy_from_elements(elements: np.ndarray, n_points: int = 900) -> tuple[np.ndarray, np.ndarray]:
    """Return heliocentric ecliptic x/y coordinates in AU for a conic section."""
    a_au, e, inc_deg, node_deg, argp_deg, _nu_deg, q_au = [float(x) for x in elements]
    if not np.all(np.isfinite(elements[:6])) or not np.isfinite(q_au) or e < 0:
        return np.array([]), np.array([])
    if e < 1.0:
        f = np.linspace(0.0, 2.0 * np.pi, n_points)
        p = max(abs(a_au) * max(1.0 - e * e, 1e-12), 1e-9)
    else:
        fmax = np.arccos(np.clip(-1.0 / max(e, 1.0 + 1e-9), -1.0, 1.0)) - 1e-3
        f = np.linspace(-fmax, fmax, n_points)
        p = max(abs(q_au) * (1.0 + e), 1e-9)
    denom = 1.0 + e * np.cos(f)
    good = denom > 1e-6
    r = p / denom[good]
    f = f[good]
    perifocal = np.vstack([r * np.cos(f), r * np.sin(f), np.zeros_like(r)])
    xyz = rotation_matrix_313(node_deg, inc_deg, argp_deg) @ perifocal
    return xyz[0], xyz[1]


def format_element_box(kep: np.ndarray, std: np.ndarray, frac_e_gt_1: float) -> str:
    a, e, inc, node, argp, nu, q = kep
    sa, se, sinc, snode, sargp, snu, sq = std
    q_aphelion = a * (1.0 + e) if e < 1.0 else np.nan
    return "\n".join(
        [
            "DASST-corrected H03",
            f"a = {a:.1f} +/- {sa:.1f} AU",
            f"q = {q:.3f} +/- {sq:.3f} AU",
            f"Q = {q_aphelion:.1f} AU" if np.isfinite(q_aphelion) else "Q = unbound",
            f"e = {e:.4f} +/- {se:.4f}",
            f"i = {inc:.2f} +/- {sinc:.2f} deg",
            f"omega = {argp:.2f} +/- {sargp:.2f} deg",
            f"Omega = {node:.3f} +/- {snode:.3f} deg",
            f"nu = {nu:.2f} +/- {snu:.2f} deg",
            f"P(e > 1) = {frac_e_gt_1:.2f}",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot PANSY DASST-corrected orbit samples.")
    parser.add_argument("dasst_orbit_h5", type=Path)
    parser.add_argument("--hypothesis", default="H03")
    parser.add_argument("--output", type=Path, default=Path("../pansy_paper/memos/figures/pansy_orbit_xy_h03.png"))
    parser.add_argument("--max-samples", type=int, default=300)
    args = parser.parse_args()

    with h5py.File(args.dasst_orbit_h5, "r") as h:
        grp = h[args.hypothesis]
        kep = grp["kepler"][()]
        std = grp["kepler_std"][()]
        samples = grp["kepler_samples"][()]
        frac_e_gt_1 = float(grp.attrs["frac_e_gt_1"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.0, 9.0), constrained_layout=True)
    theta = np.linspace(0.0, 2.0 * np.pi, 720)
    for name, radius, color in PLANETS:
        ax.plot(radius * np.cos(theta), radius * np.sin(theta), color=color, lw=1.1, alpha=0.75, label=name)
        ax.scatter([radius], [0.0], s=48, color=color, zorder=4)

    finite_samples = samples[np.all(np.isfinite(samples), axis=1)]
    finite_samples = finite_samples[: args.max_samples]
    for sample in finite_samples:
        x, y = orbit_xy_from_elements(sample)
        if len(x):
            ax.plot(x, y, color="0.25", alpha=0.055, lw=0.8)

    x, y = orbit_xy_from_elements(kep)
    ax.plot(x, y, color="black", alpha=0.82, lw=2.8, label=f"{args.hypothesis} nominal")
    ax.scatter([0.0], [0.0], s=190, color="#ffd21f", edgecolor="black", zorder=5, label="Sun")
    ax.text(
        -7.55,
        -7.45,
        format_element_box(kep, std, frac_e_gt_1),
        fontsize=15,
        va="bottom",
        ha="left",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78},
    )
    ax.text(
        7.55,
        -7.55,
        f"n={len(finite_samples)} DASST samples",
        ha="right",
        va="bottom",
        fontsize=11,
        color="0.45",
    )
    ax.set_xlabel("Heliocentric ecliptic x (AU)", fontsize=16)
    ax.set_ylabel("Heliocentric ecliptic y (AU)", fontsize=16)
    ax.set_xlim(-8.0, 8.0)
    ax.set_ylim(-8.0, 8.0)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.90", lw=1.0)
    ax.legend(loc="upper center", ncol=3, fontsize=12, framealpha=0.92)
    fig.savefig(args.output, dpi=220)
    print(args.output)


if __name__ == "__main__":
    main()
