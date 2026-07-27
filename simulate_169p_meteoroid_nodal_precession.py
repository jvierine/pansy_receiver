#!/usr/bin/env python3
"""Track nodal evolution of 100-micrometre meteoroids released by 169P/NEAT."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import rebound
import reboundx


SECONDS_PER_YEAR = 365.25 * 86400.0
AU_M = 149597870700.0
C_AU_PER_YEAR = 299792458.0 * SECONDS_PER_YEAR / AU_M
G_SI = 6.67430e-11
L_SUN_W = 3.828e26
M_SUN_KG = 1.98847e30

PLANETS = ("Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune")
SHOWER_SOLAR_LONGITUDES_DEG = {
    "CAP": 132.3,
    "OES": 294.2,
    "DCS": 312.53,
}


def radiation_beta(diameter_um: float, density_kg_m3: float, q_pr: float = 1.0) -> float:
    radius_m = 0.5 * diameter_um * 1e-6
    return 3.0 * L_SUN_W * q_pr / (
        16.0 * np.pi * 299792458.0 * G_SI * M_SUN_KG * density_kg_m3 * radius_m
    )


def wrap360(value):
    return np.asarray(value) % 360.0


def wrap180(value):
    return (np.asarray(value) + 180.0) % 360.0 - 180.0


def state_elements(position_au, velocity_au_yr, mu_au3_yr2):
    r_vec = np.asarray(position_au, dtype=np.float64)
    v_vec = np.asarray(velocity_au_yr, dtype=np.float64)
    r = np.linalg.norm(r_vec)
    v2 = np.dot(v_vec, v_vec)
    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)
    node_vec = np.cross([0.0, 0.0, 1.0], h_vec)
    node = np.linalg.norm(node_vec)
    e_vec = np.cross(v_vec, h_vec) / mu_au3_yr2 - r_vec / r
    e = np.linalg.norm(e_vec)
    energy = 0.5 * v2 - mu_au3_yr2 / r
    a = -mu_au3_yr2 / (2.0 * energy)
    inc = np.rad2deg(np.arccos(np.clip(h_vec[2] / h, -1.0, 1.0)))
    ascending_node = wrap360(np.rad2deg(np.arctan2(node_vec[1], node_vec[0]))).item()
    if node > 0.0 and e > 0.0:
        cos_argp = np.clip(np.dot(node_vec, e_vec) / (node * e), -1.0, 1.0)
        argp = np.arccos(cos_argp)
        if e_vec[2] < 0.0:
            argp = 2.0 * np.pi - argp
        argp_deg = np.rad2deg(argp)
    else:
        argp_deg = np.nan
    p = h * h / mu_au3_yr2
    cos_argp = np.cos(np.deg2rad(argp_deg))
    r_ascending = p / (1.0 + e * cos_argp)
    r_descending = p / (1.0 - e * cos_argp)
    return np.asarray(
        [a, e, inc, ascending_node, argp_deg, r_ascending, r_descending],
        dtype=np.float64,
    )


def heliocentric_state(sim: rebound.Simulation, particle_index: int):
    sun = sim.particles[0]
    particle = sim.particles[particle_index]
    position = np.asarray([particle.x - sun.x, particle.y - sun.y, particle.z - sun.z])
    velocity = np.asarray([particle.vx - sun.vx, particle.vy - sun.vy, particle.vz - sun.vz])
    return position, velocity


def make_present_simulation(epoch: str, step_years: float) -> rebound.Simulation:
    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")
    sim.integrator = "whfast"
    sim.dt = step_years
    sim.add("Sun", date=epoch)
    for body in PLANETS:
        sim.add(body, date=epoch)
    sim.add("169P", date=epoch, hash="169P")
    sim.N_active = 1 + len(PLANETS)
    sim.move_to_com()
    return sim


def add_radiation_forces(sim: rebound.Simulation):
    sim.integrator = "ias15"
    extras = reboundx.Extras(sim)
    force = extras.load_force("radiation_forces")
    force.params["c"] = C_AU_PER_YEAR
    extras.add_force(force)
    sim.particles[0].params["radiation_source"] = 1
    return extras


def release_particles(
    sim: rebound.Simulation,
    release_start_years: float,
    release_stop_years: float,
    beta: float,
    detection_step_years: float,
):
    parent_index = 1 + len(PLANETS)
    releases = []
    previous_radial_velocity = None
    while sim.t < release_stop_years:
        sim.integrate(min(sim.t + detection_step_years, release_stop_years))
        position, velocity = heliocentric_state(sim, parent_index)
        radial_velocity = float(np.dot(position, velocity) / np.linalg.norm(position))
        if previous_radial_velocity is not None and previous_radial_velocity < 0.0 <= radial_velocity:
            parent = sim.particles[parent_index]
            sim.add(
                x=parent.x,
                y=parent.y,
                z=parent.z,
                vx=parent.vx,
                vy=parent.vy,
                vz=parent.vz,
                m=0.0,
            )
            dust_index = len(sim.particles) - 1
            sim.particles[dust_index].params["beta"] = beta
            dust_position, dust_velocity = heliocentric_state(sim, dust_index)
            mu = sim.G * sim.particles[0].m * (1.0 - beta)
            elements = state_elements(dust_position, dust_velocity, mu)
            releases.append((dust_index, sim.t, elements))
        previous_radial_velocity = radial_velocity
    if not releases:
        raise RuntimeError(
            f"no 169P perihelia found between {release_start_years:g} and "
            f"{release_stop_years:g} years"
        )
    return releases


def run(args):
    beta = radiation_beta(args.diameter_um, args.density_kg_m3, args.q_pr)
    sim = make_present_simulation(args.epoch, args.step_days / 365.25)
    release_start = -abs(args.release_oldest_years)
    release_stop = -abs(args.release_youngest_years)
    sim.integrate(release_start)
    sim.dt = abs(sim.dt)
    extras = add_radiation_forces(sim)
    releases = release_particles(
        sim,
        release_start,
        release_stop,
        beta,
        args.step_days / 365.25,
    )
    sim.integrate(0.0)

    n = len(releases)
    release_epoch_years = np.empty(n)
    release_elements = np.empty((n, 7))
    present_elements = np.empty((n, 7))
    mu = sim.G * sim.particles[0].m * (1.0 - beta)
    for row, (particle_index, release_time, initial_elements) in enumerate(releases):
        position, velocity = heliocentric_state(sim, particle_index)
        release_epoch_years[row] = release_time
        release_elements[row] = initial_elements
        present_elements[row] = state_elements(position, velocity, mu)

    node_choice = np.argmin(np.abs(present_elements[:, 5:7] - 1.0), axis=1)
    node_distance = present_elements[np.arange(n), 5 + node_choice]
    node_solar_longitude = np.where(
        node_choice == 0,
        wrap360(present_elements[:, 3] - 180.0),
        wrap360(present_elements[:, 3]),
    )
    node_precession = wrap180(present_elements[:, 3] - release_elements[:, 3])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output, "w") as handle:
        handle.attrs["description"] = "Nodal precession of 169P/NEAT meteoroids"
        handle.attrs["present_epoch"] = args.epoch
        handle.attrs["integrator"] = (
            "REBOUND WHFast backward parent reconstruction; IAS15 forward "
            "integration with REBOUNDx radiation_forces"
        )
        handle.attrs["release_velocity_m_s"] = 0.0
        handle.attrs["diameter_um"] = args.diameter_um
        handle.attrs["density_kg_m3"] = args.density_kg_m3
        handle.attrs["q_pr"] = args.q_pr
        handle.attrs["radiation_beta"] = beta
        handle.attrs["step_days"] = args.step_days
        handle.create_dataset("release_epoch_years", data=release_epoch_years.astype(np.float32))
        handle.create_dataset("release_elements", data=release_elements.astype(np.float32))
        handle.create_dataset("present_elements", data=present_elements.astype(np.float32))
        handle.create_dataset("node_choice", data=node_choice.astype(np.int8))
        handle.create_dataset("node_distance_au", data=node_distance.astype(np.float32))
        handle.create_dataset(
            "node_solar_longitude_deg", data=node_solar_longitude.astype(np.float32)
        )
        handle.create_dataset("node_precession_deg", data=node_precession.astype(np.float32))
        handle["release_elements"].attrs["columns"] = "a,e,i,Omega,omega,r_asc,r_desc"
        handle["present_elements"].attrs["columns"] = "a,e,i,Omega,omega,r_asc,r_desc"

    plot_results(
        args.plot,
        release_epoch_years,
        node_precession,
        present_elements,
        node_distance,
        node_solar_longitude,
        beta,
    )
    _ = extras
    return n, beta


def plot_results(
    output,
    release_epoch_years,
    node_precession,
    present_elements,
    node_distance,
    node_solar_longitude,
    beta,
):
    output = Path(output)
    age = -release_epoch_years
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 9.0), sharex=True)
    axes[0].plot(age, node_precession, "o", ms=3.0)
    axes[0].set_ylabel(r"$\Delta\Omega$ (deg)")
    axes[1].plot(age, present_elements[:, 5], "o", ms=3.0, label="Ascending node")
    axes[1].plot(age, present_elements[:, 6], "o", ms=3.0, label="Descending node")
    axes[1].axhline(1.0, color="black", lw=0.8, ls="--")
    axes[1].set_ylabel("Present node distance (AU)")
    axes[1].legend(frameon=False)
    distance_mask = np.abs(node_distance - 1.0) <= 0.1
    axes[2].plot(age[distance_mask], node_solar_longitude[distance_mask], "o", ms=3.0)
    shower_colors = {"CAP": "C1", "OES": "C2", "DCS": "C3"}
    for name, longitude in SHOWER_SOLAR_LONGITUDES_DEG.items():
        axes[2].axhline(
            longitude,
            color=shower_colors[name],
            lw=0.9,
            ls="--",
            label=name,
        )
    axes[2].set_ylabel(r"Solar longitude at nearest node (deg)")
    axes[2].set_xlabel("Release age (yr)")
    axes[2].set_ylim(0.0, 360.0)
    axes[2].legend(frameon=False, ncol=3)
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.suptitle(rf"169P/NEAT, 100-$\mu$m diameter, $\beta={beta:.4f}$", y=0.985)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.07, top=0.94, hspace=0.07)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epoch", default="2026-01-01 00:00")
    parser.add_argument("--diameter-um", type=float, default=100.0)
    parser.add_argument("--density-kg-m3", type=float, default=3000.0)
    parser.add_argument("--q-pr", type=float, default=1.0)
    parser.add_argument("--release-oldest-years", type=float, default=5000.0)
    parser.add_argument("--release-youngest-years", type=float, default=4500.0)
    parser.add_argument("--step-days", type=float, default=3.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("figs/169p_100um_nodal_precession.h5"),
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=Path("figs/169p_100um_nodal_precession.png"),
    )
    return parser.parse_args()


if __name__ == "__main__":
    parsed = parse_args()
    count, particle_beta = run(parsed)
    print(f"released {count} particles; beta={particle_beta:.6f}")
    print(parsed.output)
    print(parsed.plot)
