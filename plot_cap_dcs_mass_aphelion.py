#!/usr/bin/env python3
"""Plot CAP/DCS orbital elements against initial mass and radius constraints."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_SUMMARY = Path(__file__).resolve().parent / "figs" / "cap_dcs_three_pulse_mass_summary.h5"
FIGURE_DIR = Path(__file__).resolve().parent / "figs"
DEFAULT_MASS_APHELION_OUTPUT = FIGURE_DIR / "cap_dcs_mass_vs_aphelion_bounds.png"
DEFAULT_RADIUS_APHELION_OUTPUT = FIGURE_DIR / "cap_dcs_radius_vs_aphelion_bounds.png"
DEFAULT_RADIUS_ECCENTRICITY_OUTPUT = FIGURE_DIR / "cap_dcs_radius_vs_eccentricity_bounds.png"


GROUPS = (
    (0, "DCS, late interval", "C0", "o"),
    (1, "CAP", "C1", "o"),
    (2, "DCS, early interval", "C0", "s"),
)


def plot_constraint_triptych(
    output: Path,
    constraint_values: tuple[np.ndarray, np.ndarray, np.ndarray],
    y: np.ndarray,
    orbit_valid: np.ndarray,
    passage_index: np.ndarray,
    maximum_at_bound: np.ndarray,
    upper_constrained: np.ndarray,
    xlabel: str,
    ylabel: str,
) -> None:
    fig, axes = plt.subplots(
        1, 3, figsize=(14.2, 4.6), sharex=True, sharey=True, constrained_layout=True
    )
    panels = (
        (constraint_values[0], "95% lower bound", np.zeros(len(y), dtype=bool)),
        (constraint_values[1], "Maximum likelihood", maximum_at_bound),
        (constraint_values[2], "95% upper bound", orbit_valid & ~upper_constrained),
    )
    for ax, (x, title, flagged) in zip(axes, panels, strict=True):
        valid = orbit_valid & np.isfinite(x) & (x > 0.0) & np.isfinite(y)
        for index, label, color, marker in GROUPS:
            selected = valid & (passage_index == index) & ~flagged
            ax.scatter(
                x[selected],
                y[selected],
                s=14,
                marker=marker,
                color=color,
                alpha=0.46,
                linewidths=0,
                label=label,
            )
        if np.any(valid & flagged):
            ax.scatter(
                x[valid & flagged],
                y[valid & flagged],
                s=20,
                marker="x",
                color="0.35",
                alpha=0.60,
                linewidths=0.7,
                label=(
                    "At fit bound"
                    if title == "Maximum likelihood"
                    else "Upper bound unconstrained"
                ),
            )
        ax.set_xscale("log")
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.grid(alpha=0.22, linewidth=0.5)

    axes[0].set_ylabel(ylabel)
    handles, labels = axes[-1].get_legend_handles_labels()
    unique = dict(zip(labels, handles, strict=True))
    axes[-1].legend(unique.values(), unique.keys(), frameon=False, fontsize=8.5)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240)
    plt.close(fig)
    print(output)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument(
        "--mass-aphelion-output", type=Path, default=DEFAULT_MASS_APHELION_OUTPUT
    )
    parser.add_argument(
        "--radius-aphelion-output", type=Path, default=DEFAULT_RADIUS_APHELION_OUTPUT
    )
    parser.add_argument(
        "--radius-eccentricity-output",
        type=Path,
        default=DEFAULT_RADIUS_ECCENTRICITY_OUTPUT,
    )
    args = parser.parse_args()

    with h5py.File(args.summary, "r") as handle:
        maximum_likelihood_mass = np.asarray(
            handle["maximum_likelihood_initial_mass_kg"], dtype=float
        )
        marginal_mass = np.asarray(handle["marginal_mass_quantiles_kg"], dtype=float)
        radius = np.asarray(handle["maximum_likelihood_initial_radius_um"], dtype=float)
        marginal_radius = np.asarray(
            handle["marginal_radius_quantiles_um"], dtype=float
        )
        upper_constrained = np.asarray(
            handle["radius_upper_limit_data_constrained"], dtype=bool
        )
        kepler = np.asarray(handle["kepler"], dtype=float)
        passage_index = np.asarray(handle["passage_index"], dtype=np.int8)
        status = np.asarray(handle["status"], dtype=np.int8)

    semimajor_axis = kepler[:, 0]
    eccentricity = kepler[:, 1]
    aphelion = semimajor_axis * (1.0 + eccentricity)
    orbit_valid = (
        (status == 1)
        & np.isfinite(aphelion)
        & (aphelion > 0.0)
    )
    maximum_at_bound = orbit_valid & (radius >= 9999.0)

    plot_constraint_triptych(
        args.mass_aphelion_output,
        (marginal_mass[:, 0], maximum_likelihood_mass, marginal_mass[:, 2]),
        aphelion,
        orbit_valid,
        passage_index,
        maximum_at_bound,
        upper_constrained,
        r"Initial mass, $m_0$ (kg)",
        r"Aphelion distance, $Q=a(1+e)$ (AU)",
    )
    radius_constraints = (
        marginal_radius[:, 0],
        radius,
        marginal_radius[:, 2],
    )
    plot_constraint_triptych(
        args.radius_aphelion_output,
        radius_constraints,
        aphelion,
        orbit_valid,
        passage_index,
        maximum_at_bound,
        upper_constrained,
        r"Initial radius, $r_0$ ($\mu$m)",
        r"Aphelion distance, $Q=a(1+e)$ (AU)",
    )
    plot_constraint_triptych(
        args.radius_eccentricity_output,
        radius_constraints,
        eccentricity,
        orbit_valid,
        passage_index,
        maximum_at_bound,
        upper_constrained,
        r"Initial radius, $r_0$ ($\mu$m)",
        r"Eccentricity, $e$",
    )
    print(
        f"valid={np.count_nonzero(orbit_valid)} "
        f"maximum_at_radius_bound={np.count_nonzero(maximum_at_bound)} "
        f"upper_bound_constrained={np.count_nonzero(orbit_valid & upper_constrained)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
