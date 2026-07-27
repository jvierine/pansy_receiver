#!/usr/bin/env python3
"""Compare 10-kyr Jovian dispersion of high- and low-inclination SER analogues."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import rebound


SER_MEAN = np.asarray(
    [
        3.44,   # a (AU)
        0.712,  # e
        91.22,  # i (deg)
        290.12, # Omega (deg)
        317.12, # omega (deg)
        42.45,  # true anomaly (deg)
    ],
    dtype=np.float64,
)
LOW_INCLINATION_DEG = 10.0
ELEMENT_COLUMNS = ("a_au", "e", "i_deg", "Omega_deg", "omega_deg", "q_au", "Q_au")

# Approximate J2000 osculating elements. The fixed initialization keeps the
# experiment deterministic and offline while retaining Jovian encounters.
JUPITER_MASS_MSUN = 9.543e-4
JUPITER_ELEMENTS = {
    "a": 5.2044,
    "e": 0.0489,
    "inc": 1.303,
    "Omega": 100.464,
    "omega": 273.867,
    "M": 20.020,
}


def wrap180(values):
    return (np.asarray(values, dtype=np.float64) + 180.0) % 360.0 - 180.0


def circular_std_deg(values, axis=None):
    """Circular standard deviation in degrees, robust to the 0/360 boundary."""
    angle = np.deg2rad(np.asarray(values, dtype=np.float64))
    resultant = np.abs(np.nanmean(np.exp(1j * angle), axis=axis))
    resultant = np.clip(resultant, np.finfo(np.float64).tiny, 1.0)
    return np.rad2deg(np.sqrt(-2.0 * np.log(resultant)))


def circular_mean_unwrapped_deg(values, axis=None):
    """Circular mean unwrapped along its remaining time axis."""
    angle = np.deg2rad(np.asarray(values, dtype=np.float64))
    mean_angle = np.angle(np.nanmean(np.exp(1j * angle), axis=axis))
    return np.rad2deg(np.unwrap(mean_angle, axis=-1))


def paired_initial_elements(
    count: int,
    seed: int,
    sigma_a_au: float,
    sigma_e: float,
    sigma_angle_deg: float,
) -> np.ndarray:
    """Return matched compact ensembles differing only in mean inclination."""
    rng = np.random.default_rng(seed)
    offsets = np.column_stack(
        (
            rng.normal(0.0, sigma_a_au, count),
            rng.normal(0.0, sigma_e, count),
            rng.normal(0.0, sigma_angle_deg, count),
            rng.normal(0.0, sigma_angle_deg, count),
            rng.normal(0.0, sigma_angle_deg, count),
            rng.normal(0.0, sigma_angle_deg, count),
        )
    )
    high = SER_MEAN[None, :] + offsets
    low = high.copy()
    low[:, 2] += LOW_INCLINATION_DEG - SER_MEAN[2]
    return np.stack((high, low), axis=0)


def make_simulation(initial_elements: np.ndarray) -> rebound.Simulation:
    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")
    sim.integrator = "ias15"
    sim.add(m=1.0, hash="Sun")
    sun = sim.particles["Sun"]
    sim.add(
        m=JUPITER_MASS_MSUN,
        a=JUPITER_ELEMENTS["a"],
        e=JUPITER_ELEMENTS["e"],
        inc=np.deg2rad(JUPITER_ELEMENTS["inc"]),
        Omega=np.deg2rad(JUPITER_ELEMENTS["Omega"]),
        omega=np.deg2rad(JUPITER_ELEMENTS["omega"]),
        M=np.deg2rad(JUPITER_ELEMENTS["M"]),
        primary=sun,
        hash="Jupiter",
    )
    sim.N_active = 2
    sim.move_to_com()
    sun = sim.particles["Sun"]

    for ensemble in initial_elements:
        for a, e, inc, node, argp, anomaly in ensemble:
            sim.add(
                m=0.0,
                a=float(a),
                e=float(e),
                inc=np.deg2rad(inc),
                Omega=np.deg2rad(node),
                omega=np.deg2rad(argp),
                f=np.deg2rad(anomaly),
                primary=sun,
            )
    return sim


def particle_elements(sim: rebound.Simulation, particle_index: int) -> np.ndarray:
    orbit = sim.particles[particle_index].orbit(primary=sim.particles["Sun"])
    a = float(orbit.a)
    e = float(orbit.e)
    return np.asarray(
        [
            a,
            e,
            np.rad2deg(orbit.inc),
            np.rad2deg(orbit.Omega) % 360.0,
            np.rad2deg(orbit.omega) % 360.0,
            a * (1.0 - e),
            a * (1.0 + e),
        ],
        dtype=np.float64,
    )


def integrate_ensembles(
    initial_elements: np.ndarray,
    duration_years: float,
    output_step_years: float,
) -> tuple[np.ndarray, np.ndarray]:
    count = initial_elements.shape[1]
    times = np.arange(0.0, duration_years + 0.5 * output_step_years, output_step_years)
    elements = np.full((2, count, len(times), len(ELEMENT_COLUMNS)), np.nan, dtype=np.float32)
    sim = make_simulation(initial_elements)

    for time_index, time_years in enumerate(times):
        sim.integrate(float(time_years))
        for ensemble_index in range(2):
            start = 2 + ensemble_index * count
            for particle_offset in range(count):
                values = particle_elements(sim, start + particle_offset)
                elements[ensemble_index, particle_offset, time_index] = values.astype(np.float32)
    return times, elements


def dispersion_series(elements: np.ndarray) -> dict[str, np.ndarray]:
    bound = (
        np.isfinite(elements[..., 0])
        & np.isfinite(elements[..., 1])
        & (elements[..., 0] > 0.0)
        & (elements[..., 1] < 1.0)
    )
    masked = np.where(bound[..., None], elements, np.nan)
    return {
        "bound_count": np.sum(bound, axis=1).astype(np.int32),
        "a_std_au": np.nanstd(masked[..., 0], axis=1),
        "e_std": np.nanstd(masked[..., 1], axis=1),
        "i_std_deg": np.nanstd(masked[..., 2], axis=1),
        "Omega_mean_unwrapped_deg": circular_mean_unwrapped_deg(masked[..., 3], axis=1),
        "Omega_circular_std_deg": circular_std_deg(masked[..., 3], axis=1),
        "omega_mean_unwrapped_deg": circular_mean_unwrapped_deg(masked[..., 4], axis=1),
        "omega_circular_std_deg": circular_std_deg(masked[..., 4], axis=1),
        "q_std_au": np.nanstd(masked[..., 5], axis=1),
    }


def save_results(
    output: Path,
    initial_elements: np.ndarray,
    times: np.ndarray,
    elements: np.ndarray,
    dispersion: dict[str, np.ndarray],
    args,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        handle.attrs["description"] = (
            "REBOUND Sun-Jupiter comparison of SER-like ensembles at 91.22 and 10 deg inclination"
        )
        handle.attrs["integrator"] = "IAS15"
        handle.attrs["dynamical_model"] = "Sun and Jupiter; massless meteoroids"
        handle.attrs["duration_years"] = np.float32(args.duration_years)
        handle.attrs["output_step_years"] = np.float32(args.output_step_years)
        handle.attrs["random_seed"] = np.int32(args.seed)
        handle.attrs["element_columns"] = ",".join(ELEMENT_COLUMNS)
        handle.attrs["ensemble_names"] = "SER i=91.22 deg,SER analogue i=10 deg"
        handle.attrs["sigma_a_au"] = np.float32(args.sigma_a_au)
        handle.attrs["sigma_e"] = np.float32(args.sigma_e)
        handle.attrs["sigma_angle_deg"] = np.float32(args.sigma_angle_deg)
        handle.create_dataset("time_years", data=times.astype(np.float32))
        handle.create_dataset(
            "initial_elements",
            data=initial_elements.astype(np.float32),
            compression="gzip",
            compression_opts=3,
        )
        handle.create_dataset(
            "osculating_elements",
            data=elements,
            compression="gzip",
            compression_opts=3,
        )
        group = handle.create_group("dispersion")
        for name, values in dispersion.items():
            dtype = np.int32 if name == "bound_count" else np.float32
            group.create_dataset(name, data=values.astype(dtype))


def plot_results(
    output: Path,
    times: np.ndarray,
    elements: np.ndarray,
    dispersion: dict[str, np.ndarray],
) -> None:
    colors = ("#1f77b4", "#d95f02")
    labels = (r"SER mean, $i=91.22^\circ$", r"Low-inclination analogue, $i=10^\circ$")
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 8.4), layout="constrained")
    fig.get_layout_engine().set(h_pad=0.08, hspace=0.12)

    for ensemble_index in range(2):
        omega_mean = dispersion["omega_mean_unwrapped_deg"][ensemble_index]
        omega_change = omega_mean - omega_mean[0]
        omega_std = dispersion["omega_circular_std_deg"][ensemble_index]
        axes[0, 0].plot(
            times,
            omega_change,
            color=colors[ensemble_index],
            label=labels[ensemble_index],
        )
        axes[0, 0].fill_between(
            times,
            omega_change - omega_std,
            omega_change + omega_std,
            color=colors[ensemble_index],
            alpha=0.12,
            linewidth=0,
        )
        axes[0, 1].plot(
            times,
            omega_std,
            color=colors[ensemble_index],
        )
        axes[0, 2].plot(
            times,
            dispersion["q_std_au"][ensemble_index],
            color=colors[ensemble_index],
        )

    axes[0, 0].set_ylabel(r"Mean $\Delta\omega$ (deg)")
    axes[0, 1].set_ylabel(r"Circular $\sigma_\omega$ (deg)")
    axes[0, 2].set_ylabel(r"$\sigma_q$ (AU)")
    axes[0, 1].set_yscale("log")
    axes[0, 2].set_yscale("log")
    axes[0, 0].legend(frameon=False, fontsize=9)
    axes[0, 1].axhline(8.89, color="0.25", lw=0.9, ls="--")
    axes[0, 1].text(
        0.98,
        8.89,
        r"Observed SER $\sigma_\omega=8.89^\circ$ ",
        transform=axes[0, 1].get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=8.5,
        color="0.25",
    )

    for ensemble_index in range(2):
        final = elements[ensemble_index, :, -1]
        bound = (
            np.isfinite(final[:, 0])
            & np.isfinite(final[:, 1])
            & (final[:, 0] > 0.0)
            & (final[:, 1] < 1.0)
        )
        axes[1, ensemble_index].scatter(
            wrap180(final[bound, 3] - SER_MEAN[3]),
            wrap180(final[bound, 4] - SER_MEAN[4]),
            s=12,
            color=colors[ensemble_index],
            alpha=0.7,
            linewidths=0,
        )
        axes[1, ensemble_index].set_xlim(-180.0, 180.0)
        axes[1, ensemble_index].set_ylim(-180.0, 180.0)
        axes[1, ensemble_index].set_xlabel(r"$\Delta\Omega$ (deg)")
        axes[1, ensemble_index].set_ylabel(r"$\Delta\omega$ (deg)")
        axes[1, ensemble_index].set_title(labels[ensemble_index])

    for ensemble_index in range(2):
        final = elements[ensemble_index, :, -1]
        bound = (
            np.isfinite(final[:, 0])
            & np.isfinite(final[:, 1])
            & (final[:, 0] > 0.0)
            & (final[:, 1] < 1.0)
        )
        axes[1, 2].scatter(
            final[bound, 0],
            final[bound, 1],
            s=12,
            color=colors[ensemble_index],
            alpha=0.65,
            linewidths=0,
            label=f"{labels[ensemble_index]} ({np.sum(bound)} bound)",
        )
    axes[1, 2].set_xlabel("$a$ (AU)")
    axes[1, 2].set_ylabel("$e$")
    axes[1, 2].legend(frameon=False, fontsize=8)

    for ax in axes[0]:
        ax.set_xlabel("Integration time (yr)")
    for ax in axes.flat:
        ax.grid(alpha=0.2)
    fig.suptitle(
        "Dispersion of compact SER-like ensembles under Jovian perturbations",
        fontsize=15,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duration-years", type=float, default=10_000.0)
    parser.add_argument("--output-step-years", type=float, default=25.0)
    parser.add_argument("--particles", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20260727)
    parser.add_argument("--sigma-a-au", type=float, default=0.01)
    parser.add_argument("--sigma-e", type=float, default=0.002)
    parser.add_argument("--sigma-angle-deg", type=float, default=0.2)
    parser.add_argument("--output", type=Path, default=Path("figs/ser_apsidal_dispersion_10000yr.h5"))
    parser.add_argument("--plot", type=Path, default=Path("figs/ser_apsidal_dispersion_10000yr.png"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    initial = paired_initial_elements(
        args.particles,
        args.seed,
        args.sigma_a_au,
        args.sigma_e,
        args.sigma_angle_deg,
    )
    times, elements = integrate_ensembles(initial, args.duration_years, args.output_step_years)
    dispersion = dispersion_series(elements)
    save_results(args.output, initial, times, elements, dispersion, args)
    plot_results(args.plot, times, elements, dispersion)
    for index, name in enumerate(("SER i=91.22 deg", "analogue i=10 deg")):
        print(
            f"{name}: bound={dispersion['bound_count'][index, -1]}/{args.particles}, "
            f"sigma_omega={dispersion['omega_circular_std_deg'][index, -1]:.2f} deg, "
            f"sigma_Omega={dispersion['Omega_circular_std_deg'][index, -1]:.2f} deg, "
            f"sigma_q={dispersion['q_std_au'][index, -1]:.3f} AU"
        )
    print(args.output)
    print(args.plot)


if __name__ == "__main__":
    main()
