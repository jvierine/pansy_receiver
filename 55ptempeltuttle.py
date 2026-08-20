#!/usr/bin/env python3
"""Animate secular Poynting-Robertson and Jupiter evolution of 55P dust.

The secular equations are Fitzpatrick (2012), equations 10.183--10.188:
https://farside.ph.utexas.edu/teaching/celestial/Celestial/node95.html

Jupiter's perturbation is averaged over both orbits and converted to element
rates with Fitzpatrick's Gauss equations I.55--I.58. A Jupiter-Hill-radius
softening suppresses individual close encounters, which are outside a secular
model. The evolution stops when the semimajor axis reaches 1 AU.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


AU_M = 149_597_870_700.0
YEAR_S = 365.25 * 86_400.0
C_M_S = 299_792_458.0
MU_SUN = 1.327_124_400_18e20
L_SUN_W = 3.846e26  # Value used by Fitzpatrick (2012).
MU_SUN_AU3_YR2 = 4.0 * np.pi**2
JUPITER_MASS_SOLAR = 9.543e-4
JUPITER_A_AU = 5.2044
JUPITER_E = 0.0489
JUPITER_I_DEG = 1.303
JUPITER_NODE_DEG = 100.56
JUPITER_ARG_PERI_DEG = 273.87
JUPITER_RING_SAMPLES = 40
COMET_RING_SAMPLES = 40
JUPITER_SOFTENING_AU = 0.35

# JPL SBDB osculating elements, retrieved 2026-08-20.
A0_AU = 10.338_338_229_757_7
E0 = 0.905_552_720_972_412
I_DEG = 162.486_575_379_434
NODE_DEG = 235.270_989_149_082
ARG_PERI_DEG = 172.500_273_682_805_9
M0_DEG = 4.978_339_684_688_16


def element_rotation(node: float, inclination: float, arg_peri: float) -> np.ndarray:
    """Return the perifocal-to-ecliptic rotation matrix."""
    def rz(angle: float) -> np.ndarray:
        return np.array(
            [
                [np.cos(angle), -np.sin(angle), 0.0],
                [np.sin(angle), np.cos(angle), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )

    rx = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(inclination), -np.sin(inclination)],
            [0.0, np.sin(inclination), np.cos(inclination)],
        ]
    )
    return rz(node) @ rx @ rz(arg_peri)


def solve_kepler(mean_anomaly: np.ndarray | float, eccentricity: float) -> np.ndarray | float:
    mean_anomaly = np.remainder(mean_anomaly, 2.0 * np.pi)
    eccentric_anomaly = np.array(mean_anomaly, dtype=np.float64, copy=True)
    for _ in range(12):
        eccentric_anomaly -= (
            eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly) - mean_anomaly
        ) / (1.0 - eccentricity * np.cos(eccentric_anomaly))
    if np.ndim(mean_anomaly) == 0:
        return float(eccentric_anomaly)
    return eccentric_anomaly


def sampled_orbit(
    a_au: float,
    eccentricity: float,
    inclination: float,
    node: float,
    arg_peri: float,
    mean_anomaly: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return positions and local orbital basis at uniform mean anomalies."""
    eccentric_anomaly = np.asarray(solve_kepler(mean_anomaly, eccentricity))
    true_anomaly = np.arctan2(
        np.sqrt(1.0 - eccentricity**2) * np.sin(eccentric_anomaly),
        np.cos(eccentric_anomaly) - eccentricity,
    )
    radius = a_au * (1.0 - eccentricity * np.cos(eccentric_anomaly))
    rotation = element_rotation(node, inclination, arg_peri)
    radial = (
        rotation
        @ np.vstack(
            [np.cos(true_anomaly), np.sin(true_anomaly), np.zeros_like(true_anomaly)]
        )
    ).T
    transverse = (
        rotation
        @ np.vstack(
            [-np.sin(true_anomaly), np.cos(true_anomaly), np.zeros_like(true_anomaly)]
        )
    ).T
    normal = rotation[:, 2]
    position = radius[:, None] * radial
    return position, radial, transverse, normal, true_anomaly


JUPITER_MEAN_ANOMALY = (np.arange(JUPITER_RING_SAMPLES) + 0.5) * (
    2.0 * np.pi / JUPITER_RING_SAMPLES
)
COMET_MEAN_ANOMALY = (np.arange(COMET_RING_SAMPLES) + 0.5) * (
    2.0 * np.pi / COMET_RING_SAMPLES
)
JUPITER_RING_POSITION = sampled_orbit(
    JUPITER_A_AU,
    JUPITER_E,
    np.deg2rad(JUPITER_I_DEG),
    np.deg2rad(JUPITER_NODE_DEG),
    np.deg2rad(JUPITER_ARG_PERI_DEG),
    JUPITER_MEAN_ANOMALY,
)[0]
JUPITER_INDIRECT_ACCELERATION = (
    MU_SUN_AU3_YR2
    * JUPITER_MASS_SOLAR
    * np.mean(
        JUPITER_RING_POSITION
        / np.linalg.norm(JUPITER_RING_POSITION, axis=1)[:, None] ** 3,
        axis=0,
    )
)


def jupiter_secular_rates(
    a_au: float,
    eccentricity: float,
    inclination: float,
    node: float,
    arg_peri: float,
) -> np.ndarray:
    """Return orbit-averaged de, di, dOmega, and domega per year."""
    position, radial, transverse, normal, true_anomaly = sampled_orbit(
        a_au,
        eccentricity,
        inclination,
        node,
        arg_peri,
        COMET_MEAN_ANOMALY,
    )
    separation = JUPITER_RING_POSITION[None, :, :] - position[:, None, :]
    softened_distance_cubed = (
        np.sum(separation**2, axis=2) + JUPITER_SOFTENING_AU**2
    ) ** 1.5
    acceleration = (
        MU_SUN_AU3_YR2
        * JUPITER_MASS_SOLAR
        * np.mean(separation / softened_distance_cubed[:, :, None], axis=1)
        - JUPITER_INDIRECT_ACCELERATION
    )
    radial_acceleration = np.sum(acceleration * radial, axis=1)
    transverse_acceleration = np.sum(acceleration * transverse, axis=1)
    normal_acceleration = acceleration @ normal

    p_au = a_au * (1.0 - eccentricity**2)
    h = np.sqrt(MU_SUN_AU3_YR2 * p_au)
    radius = np.linalg.norm(position, axis=1)
    argument_of_latitude = arg_peri + true_anomaly
    sin_i = np.sin(inclination)
    de_dt = (
        p_au * np.sin(true_anomaly) * radial_acceleration
        + (
            (p_au + radius) * np.cos(true_anomaly) + radius * eccentricity
        )
        * transverse_acceleration
    ) / h
    di_dt = radius * np.cos(argument_of_latitude) * normal_acceleration / h
    d_node_dt = (
        radius
        * np.sin(argument_of_latitude)
        * normal_acceleration
        / (h * sin_i)
    )
    d_arg_peri_dt = (
        -p_au * np.cos(true_anomaly) * radial_acceleration
        + (p_au + radius) * np.sin(true_anomaly) * transverse_acceleration
    ) / (h * eccentricity) - np.cos(inclination) * d_node_dt
    return np.array(
        [np.mean(de_dt), np.mean(di_dt), np.mean(d_node_dt), np.mean(d_arg_peri_dt)]
    )


def secular_rhs(
    _time_years: float,
    state: np.ndarray,
    mass_area: float,
    include_pr_drag: bool = True,
) -> np.ndarray:
    """Combine Fitzpatrick P-R drag with orbit-averaged Jupiter perturbations."""
    a_au, eccentricity, _mean_anomaly, inclination, node, arg_peri = state
    a_m = a_au * AU_M
    coefficient = L_SUN_W / (4.0 * np.pi * a_m**2 * mass_area * C_M_S**2)
    one_minus_e2 = 1.0 - eccentricity**2
    mean_motion = np.sqrt(MU_SUN / a_m**3)

    da_dt = (
        -a_m
        * coefficient
        * (2.0 + 3.0 * eccentricity**2)
        / one_minus_e2**1.5
    )
    de_dt = -coefficient * 5.0 * eccentricity / (2.0 * np.sqrt(one_minus_e2))
    d_m0_dt = -coefficient * 2.0 * C_M_S / (mean_motion * a_m)
    de_jupiter, di_dt, d_node_dt, d_arg_peri_dt = jupiter_secular_rates(
        a_au, eccentricity, inclination, node, arg_peri
    )
    pr_scale = 1.0 if include_pr_drag else 0.0
    return np.array(
        [
            pr_scale * da_dt * YEAR_S / AU_M,
            pr_scale * de_dt * YEAR_S + de_jupiter,
            (mean_motion + pr_scale * d_m0_dt) * YEAR_S,
            di_dt,
            d_node_dt,
            d_arg_peri_dt,
        ]
    )


def integrate_evolution(
    mass_area: float,
    include_pr_drag: bool = True,
    duration_years: float = 2_451_788.620_925_513,
) -> solve_ivp:
    """Integrate to 1 AU with P-R drag, or for a fixed Jupiter-only duration."""
    initial_state = np.array(
        [
            A0_AU,
            E0,
            np.deg2rad(M0_DEG),
            np.deg2rad(I_DEG),
            np.deg2rad(NODE_DEG),
            np.deg2rad(ARG_PERI_DEG),
        ]
    )

    def reaches_one_au(_time_s: float, state: np.ndarray) -> float:
        return state[0] - 1.0

    reaches_one_au.terminal = True
    reaches_one_au.direction = -1

    solution = solve_ivp(
        lambda t, y: secular_rhs(t, y, mass_area, include_pr_drag),
        (0.0, 3.0e6 if include_pr_drag else duration_years),
        initial_state,
        events=reaches_one_au if include_pr_drag else None,
        rtol=2e-6,
        atol=(1e-8, 1e-9, 1e-6, 1e-9, 1e-9, 1e-9),
        max_step=100.0,
    )
    if not solution.success:
        raise RuntimeError(f"secular integration failed: {solution.message}")
    if include_pr_drag and not solution.t_events[0].size:
        raise RuntimeError("P-R integration did not reach a = 1 AU")
    return solution


def resample_evolution(
    solution: solve_ivp,
    frames: int,
    include_pr_drag: bool,
) -> dict[str, np.ndarray]:
    """Resample in log(a) for P-R drag or uniformly in time for Jupiter only."""
    a_raw = solution.y[0]
    if include_pr_drag:
        a_au = np.geomspace(a_raw[0], 1.0, frames)
        time_years = np.interp(a_au[::-1], a_raw[::-1], solution.t[::-1])[::-1]
    else:
        time_years = np.linspace(solution.t[0], solution.t[-1], frames)
        a_au = np.interp(time_years, solution.t, a_raw)
    eccentricity = np.interp(time_years, solution.t, solution.y[1])
    mean_anomaly = np.remainder(
        np.interp(time_years, solution.t, solution.y[2]), 2.0 * np.pi
    )
    return {
        "time_years": time_years,
        "a_au": a_au,
        "eccentricity": eccentricity,
        "mean_anomaly_rad": mean_anomaly,
        "inclination_rad": np.interp(time_years, solution.t, solution.y[3]),
        "longitude_ascending_node_rad": np.interp(time_years, solution.t, solution.y[4]),
        "argument_perihelion_rad": np.interp(time_years, solution.t, solution.y[5]),
    }


def save_hdf5(
    path: Path,
    evolution: dict[str, np.ndarray],
    mass_area: float,
    include_pr_drag: bool,
) -> None:
    with h5py.File(path, "w") as h5:
        for key, values in evolution.items():
            h5.create_dataset(key, data=np.asarray(values, dtype=np.float32))
        h5.attrs["model"] = (
            "Fitzpatrick P-R equations 10.183-10.188 plus Jupiter-averaged "
            "Gauss equations I.55-I.58"
        )
        h5.attrs["source_url"] = (
            "https://farside.ph.utexas.edu/teaching/celestial/Celestial/node95.html"
        )
        h5.attrs["initial_body"] = "55P/Tempel-Tuttle"
        h5.attrs["mass_to_cross_section_kg_m2"] = mass_area
        h5.attrs["inclination_deg"] = I_DEG
        h5.attrs["longitude_ascending_node_deg"] = NODE_DEG
        h5.attrs["argument_perihelion_deg"] = ARG_PERI_DEG
        h5.attrs["jupiter_softening_au"] = JUPITER_SOFTENING_AU
        h5.attrs["secular_only"] = True
        h5.attrs["include_pr_drag"] = include_pr_drag


def orbit_xyz(
    a_au: float,
    eccentricity: float,
    inclination: float,
    node: float,
    arg_peri: float,
    samples: int = 700,
) -> np.ndarray:
    eccentric_anomaly = np.linspace(0.0, 2.0 * np.pi, samples)
    perifocal = np.vstack(
        (
            a_au * (np.cos(eccentric_anomaly) - eccentricity),
            a_au * np.sqrt(1.0 - eccentricity**2) * np.sin(eccentric_anomaly),
            np.zeros_like(eccentric_anomaly),
        )
    )
    return (element_rotation(node, inclination, arg_peri) @ perifocal).T


def particle_xyz(
    a_au: float,
    eccentricity: float,
    mean_anomaly: float,
    inclination: float,
    node: float,
    arg_peri: float,
) -> np.ndarray:
    eccentric_anomaly = solve_kepler(mean_anomaly, eccentricity)
    perifocal = np.array(
        [
            a_au * (np.cos(eccentric_anomaly) - eccentricity),
            a_au * np.sqrt(1.0 - eccentricity**2) * np.sin(eccentric_anomaly),
            0.0,
        ]
    )
    return element_rotation(node, inclination, arg_peri) @ perifocal


def camera_angles(frame: int, frame_count: int) -> tuple[float, float]:
    """Return a smooth, fixed-scale camera pan in azimuth and elevation."""
    phase = frame / max(frame_count - 1, 1)
    smooth_phase = phase * phase * (3.0 - 2.0 * phase)
    azimuth = np.deg2rad(-35.0 + 75.0 * smooth_phase)
    elevation = np.deg2rad(68.0 - 38.0 * np.sin(np.pi * phase) ** 2)
    return azimuth, elevation


def project_orthographic(
    points_xyz: np.ndarray,
    azimuth: float,
    elevation: float,
) -> np.ndarray:
    """Project 3D ecliptic coordinates with an orthographic camera."""
    points = np.atleast_2d(points_xyz)
    view_direction = np.array(
        [
            np.cos(elevation) * np.cos(azimuth),
            np.cos(elevation) * np.sin(azimuth),
            np.sin(elevation),
        ]
    )
    screen_right = np.array([-np.sin(azimuth), np.cos(azimuth), 0.0])
    screen_up = np.cross(view_direction, screen_right)
    return np.column_stack((points @ screen_right, points @ screen_up))


def render_animation(
    evolution: dict[str, np.ndarray],
    mass_area: float,
    include_pr_drag: bool,
    mp4_path: Path,
    poster_path: Path,
    fps: int,
) -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "text.color": "black",
            "axes.labelcolor": "black",
            "xtick.color": "0.2",
            "ytick.color": "0.2",
        }
    )
    fig, orbit_ax = plt.subplots(figsize=(8.0, 8.0), constrained_layout=True)
    for spine in orbit_ax.spines.values():
        spine.set_color("0.25")

    orbit_ax.set_aspect("equal")
    orbit_ax.set_xlabel("Projected x (AU)", fontsize=19)
    orbit_ax.set_ylabel("Projected y (AU)", fontsize=19)
    orbit_ax.tick_params(axis="both", labelsize=15)
    orbit_ax.grid(color="black", alpha=0.10, linewidth=0.5)

    theta = np.linspace(0.0, 2.0 * np.pi, 500)
    earth_xyz = np.column_stack(
        (np.cos(theta), np.sin(theta), np.zeros_like(theta))
    )
    jupiter_xyz = JUPITER_A_AU * earth_xyz
    earth_orbit, = orbit_ax.plot(
        [],
        [],
        color="#174a9c",
        lw=1.8,
        alpha=0.95,
        label="Earth orbit",
    )
    jupiter_orbit, = orbit_ax.plot(
        [],
        [],
        color="#c46b16",
        lw=1.6,
        alpha=0.90,
        label="Jupiter orbit",
    )
    orbit_ax.scatter([0.0], [0.0], s=110, color="#ffd51f", edgecolor="#fff3a0", zorder=8)

    original_xyz = orbit_xyz(
        A0_AU,
        E0,
        np.deg2rad(I_DEG),
        np.deg2rad(NODE_DEG),
        np.deg2rad(ARG_PERI_DEG),
    )
    original_orbit, = orbit_ax.plot(
        [],
        [],
        color="#d62728",
        lw=1.7,
        alpha=0.85,
        ls="--",
        label="Original 55P orbit",
    )
    evolution_label = (
        "Jupiter + P-R evolution" if include_pr_drag else "Jupiter-only evolution"
    )
    current_orbit, = orbit_ax.plot(
        [],
        [],
        color="#159947",
        lw=2.8,
        zorder=5,
        label=evolution_label,
    )
    particle = orbit_ax.scatter([], [], s=48, color="#ff3bd4", edgecolor="white", linewidth=0.7, zorder=9)
    orbit_ax.legend(loc="upper right", fontsize=13, framealpha=0.90)

    a_track = evolution["a_au"]
    e_track = evolution["eccentricity"]
    i_track = evolution["inclination_rad"]
    node_track = evolution["longitude_ascending_node_rad"]
    arg_peri_track = evolution["argument_perihelion_rad"]
    orbit_ax.set_title("55P/Tempel-Tuttle orbit evolution", fontsize=20, pad=12)
    view_radius = 1.10 * A0_AU * (1.0 + E0)
    orbit_ax.set_xlim(-view_radius, view_radius)
    orbit_ax.set_ylim(-view_radius, view_radius)
    time_text = orbit_ax.text(
        0.025,
        0.965,
        "",
        transform=orbit_ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        zorder=12,
    )

    def update(frame: int):
        a_au = float(a_track[frame])
        eccentricity = float(e_track[frame])
        elapsed_years = float(evolution["time_years"][frame])
        mean_anomaly = float(evolution["mean_anomaly_rad"][frame])
        inclination = float(i_track[frame])
        node = float(node_track[frame])
        arg_peri = float(arg_peri_track[frame])
        azimuth, elevation = camera_angles(frame, len(a_track))
        earth_xy = project_orthographic(earth_xyz, azimuth, elevation)
        jupiter_xy = project_orthographic(jupiter_xyz, azimuth, elevation)
        original_xy = project_orthographic(original_xyz, azimuth, elevation)
        current_xyz = orbit_xyz(
            a_au, eccentricity, inclination, node, arg_peri
        )
        current_xy = project_orthographic(current_xyz, azimuth, elevation)
        particle_xy = project_orthographic(
            particle_xyz(
                a_au,
                eccentricity,
                mean_anomaly,
                inclination,
                node,
                arg_peri,
            ),
            azimuth,
            elevation,
        )
        earth_orbit.set_data(earth_xy[:, 0], earth_xy[:, 1])
        jupiter_orbit.set_data(jupiter_xy[:, 0], jupiter_xy[:, 1])
        original_orbit.set_data(original_xy[:, 0], original_xy[:, 1])
        current_orbit.set_data(current_xy[:, 0], current_xy[:, 1])
        particle.set_offsets(particle_xy)

        time_text.set_text(
            f"{elapsed_years / 1000.0:,.1f} kyr\n"
            f"e = {eccentricity:.3f},  i = {np.rad2deg(inclination):.1f}$^\\circ$"
        )
        return (
            earth_orbit,
            jupiter_orbit,
            original_orbit,
            current_orbit,
            particle,
            time_text,
        )

    update(len(a_track) - 1)
    fig.savefig(poster_path, dpi=160)
    movie = animation.FuncAnimation(fig, update, frames=len(a_track), interval=1000 / fps, blit=False)
    writer = animation.FFMpegWriter(
        fps=fps,
        bitrate=3500,
        codec="libx264",
        extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"],
    )
    movie.save(mp4_path, writer=writer, dpi=100)
    plt.close(fig)


def make_gif(mp4_path: Path, gif_path: Path) -> None:
    palette = gif_path.with_suffix(".palette.png")
    subprocess.run(
        [
            "ffmpeg", "-y", "-i", str(mp4_path), "-vf",
            "fps=12,scale=960:-1:flags=lanczos,palettegen=max_colors=192",
            "-frames:v", "1", str(palette),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        [
            "ffmpeg", "-y", "-i", str(mp4_path), "-i", str(palette), "-lavfi",
            "fps=12,scale=960:-1:flags=lanczos[x];[x][1:v]paletteuse=dither=sierra2_4a",
            str(gif_path),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    palette.unlink(missing_ok=True)


def generate_products(
    output_prefix: str,
    mass_area: float,
    frames: int,
    fps: int,
    include_pr_drag: bool,
    duration_years: float,
) -> None:
    prefix = Path(output_prefix).expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    solution = integrate_evolution(
        mass_area,
        include_pr_drag=include_pr_drag,
        duration_years=duration_years,
    )
    evolution = resample_evolution(solution, frames, include_pr_drag)

    h5_path = prefix.with_suffix(".h5")
    mp4_path = prefix.with_suffix(".mp4")
    gif_path = prefix.with_suffix(".gif")
    poster_path = prefix.with_name(prefix.name + "_poster.png")
    save_hdf5(h5_path, evolution, mass_area, include_pr_drag)
    render_animation(
        evolution,
        mass_area,
        include_pr_drag,
        mp4_path,
        poster_path,
        fps,
    )
    make_gif(mp4_path, gif_path)

    final_e = evolution["eccentricity"][-1]
    final_time = evolution["time_years"][-1]
    model_name = "Jupiter + P-R" if include_pr_drag else "Jupiter only"
    print(
        f"{model_name}: t={final_time:,.1f} years, "
        f"a={evolution['a_au'][-1]:.6f} AU, e={final_e:.6f}"
    )
    print(f"Wrote {h5_path}")
    print(f"Wrote {mp4_path}")
    print(f"Wrote {gif_path}")
    print(f"Wrote {poster_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model",
        choices=("both", "jupiter", "jupiter-pr"),
        default="both",
    )
    parser.add_argument("--output-prefix", default=None)
    parser.add_argument("--mass-area", type=float, default=0.1, help="m/A in kg m^-2")
    parser.add_argument("--duration-years", type=float, default=2_451_788.620_925_513)
    parser.add_argument("--frames", type=int, default=240)
    parser.add_argument("--fps", type=int, default=24)
    args = parser.parse_args()

    products = []
    if args.model in ("both", "jupiter"):
        products.append(
            (
                args.output_prefix or "tempel_tuttle_jupiter_only",
                False,
            )
        )
    if args.model in ("both", "jupiter-pr"):
        products.append(
            (
                args.output_prefix or "tempel_tuttle_jupiter_pr_drag",
                True,
            )
        )
    if args.model == "both" and args.output_prefix:
        raise ValueError("--output-prefix can only be used with a single model")

    for output_prefix, include_pr_drag in products:
        generate_products(
            output_prefix,
            args.mass_area,
            args.frames,
            args.fps,
            include_pr_drag,
            args.duration_years,
        )


if __name__ == "__main__":
    main()
