#!/usr/bin/env python3
"""Animate orbit-averaged Poynting-Robertson evolution of 55P dust.

The secular equations are Fitzpatrick (2012), equations 10.183--10.188:
https://farside.ph.utexas.edu/teaching/celestial/Celestial/node95.html

The initial osculating elements are from the JPL Small-Body Database entry for
55P/Tempel-Tuttle. The evolution stops when the semimajor axis reaches 1 AU.
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

# JPL SBDB osculating elements, retrieved 2026-08-20.
A0_AU = 10.338_338_229_757_7
E0 = 0.905_552_720_972_412
I_DEG = 162.486_575_379_434
NODE_DEG = 235.270_989_149_082
ARG_PERI_DEG = 172.500_273_682_805_9
M0_DEG = 4.978_339_684_688_16


def secular_rhs(_time_s: float, state: np.ndarray, mass_area: float) -> np.ndarray:
    """Return da/dt, de/dt, and dM/dt from equations 10.183--10.185."""
    a_m, eccentricity, _mean_anomaly = state
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
    return np.array([da_dt, de_dt, mean_motion + d_m0_dt])


def integrate_evolution(mass_area: float) -> solve_ivp:
    """Integrate until the secular semimajor axis reaches 1 AU."""
    initial_state = np.array([A0_AU * AU_M, E0, np.deg2rad(M0_DEG)])

    def reaches_one_au(_time_s: float, state: np.ndarray) -> float:
        return state[0] - AU_M

    reaches_one_au.terminal = True
    reaches_one_au.direction = -1

    solution = solve_ivp(
        lambda t, y: secular_rhs(t, y, mass_area),
        (0.0, 2.0e6 * YEAR_S),
        initial_state,
        events=reaches_one_au,
        rtol=2e-10,
        atol=(1e-2, 1e-12, 1e-8),
        max_step=50.0 * YEAR_S,
    )
    if not solution.success or not solution.t_events[0].size:
        raise RuntimeError("P-R integration did not reach a = 1 AU")
    return solution


def resample_evolution(solution: solve_ivp, frames: int) -> dict[str, np.ndarray]:
    """Resample uniformly in log(a) so the contraction remains visible."""
    a_raw = solution.y[0] / AU_M
    a_au = np.geomspace(a_raw[0], 1.0, frames)
    time_s = np.interp(a_au[::-1], a_raw[::-1], solution.t[::-1])[::-1]
    eccentricity = np.interp(time_s, solution.t, solution.y[1])
    mean_anomaly = np.remainder(
        np.interp(time_s, solution.t, solution.y[2]), 2.0 * np.pi
    )
    return {
        "time_years": time_s / YEAR_S,
        "a_au": a_au,
        "eccentricity": eccentricity,
        "mean_anomaly_rad": mean_anomaly,
    }


def save_hdf5(path: Path, evolution: dict[str, np.ndarray], mass_area: float) -> None:
    with h5py.File(path, "w") as h5:
        for key, values in evolution.items():
            h5.create_dataset(key, data=np.asarray(values, dtype=np.float32))
        h5.attrs["model"] = "Fitzpatrick 2012 equations 10.183-10.188"
        h5.attrs["source_url"] = (
            "https://farside.ph.utexas.edu/teaching/celestial/Celestial/node95.html"
        )
        h5.attrs["initial_body"] = "55P/Tempel-Tuttle"
        h5.attrs["mass_to_cross_section_kg_m2"] = mass_area
        h5.attrs["inclination_deg"] = I_DEG
        h5.attrs["longitude_ascending_node_deg"] = NODE_DEG
        h5.attrs["argument_perihelion_deg"] = ARG_PERI_DEG


def rotation_matrix() -> np.ndarray:
    node = np.deg2rad(NODE_DEG)
    inclination = np.deg2rad(I_DEG)
    arg_peri = np.deg2rad(ARG_PERI_DEG)

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


ROTATION = rotation_matrix()


def orbit_xy(a_au: float, eccentricity: float, samples: int = 700) -> tuple[np.ndarray, np.ndarray]:
    eccentric_anomaly = np.linspace(0.0, 2.0 * np.pi, samples)
    perifocal = np.vstack(
        (
            a_au * (np.cos(eccentric_anomaly) - eccentricity),
            a_au * np.sqrt(1.0 - eccentricity**2) * np.sin(eccentric_anomaly),
            np.zeros_like(eccentric_anomaly),
        )
    )
    ecliptic = ROTATION @ perifocal
    return ecliptic[0], ecliptic[1]


def solve_kepler(mean_anomaly: float, eccentricity: float) -> float:
    mean_anomaly = np.remainder(mean_anomaly, 2.0 * np.pi)
    eccentric_anomaly = mean_anomaly
    for _ in range(12):
        eccentric_anomaly -= (
            eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly) - mean_anomaly
        ) / (1.0 - eccentricity * np.cos(eccentric_anomaly))
    return eccentric_anomaly


def particle_xy(a_au: float, eccentricity: float, mean_anomaly: float) -> tuple[float, float]:
    eccentric_anomaly = solve_kepler(mean_anomaly, eccentricity)
    perifocal = np.array(
        [
            a_au * (np.cos(eccentric_anomaly) - eccentricity),
            a_au * np.sqrt(1.0 - eccentricity**2) * np.sin(eccentric_anomaly),
            0.0,
        ]
    )
    ecliptic = ROTATION @ perifocal
    return float(ecliptic[0]), float(ecliptic[1])


def render_animation(
    evolution: dict[str, np.ndarray],
    mass_area: float,
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
    fig, orbit_ax = plt.subplots(figsize=(12.8, 7.2), constrained_layout=True)
    for spine in orbit_ax.spines.values():
        spine.set_color("0.25")

    orbit_ax.set_aspect("equal")
    orbit_ax.set_xlabel("Ecliptic x (AU)")
    orbit_ax.set_ylabel("Ecliptic y (AU)")
    orbit_ax.grid(color="black", alpha=0.10, linewidth=0.5)

    theta = np.linspace(0.0, 2.0 * np.pi, 500)
    orbit_ax.plot(np.cos(theta), np.sin(theta), color="#20dfe3", lw=1.4, alpha=0.85)
    orbit_ax.text(0.03, 0.90, "Earth orbit", color="#20dfe3", transform=orbit_ax.transAxes, va="top")
    orbit_ax.scatter([0.0], [0.0], s=110, color="#ffd51f", edgecolor="#fff3a0", zorder=8)

    original_x, original_y = orbit_xy(A0_AU, E0)
    orbit_ax.plot(original_x, original_y, color="black", lw=0.9, alpha=0.20, ls="--")
    current_orbit, = orbit_ax.plot([], [], color="#ffd21f", lw=2.5, zorder=5)
    particle = orbit_ax.scatter([], [], s=48, color="#ff3bd4", edgecolor="white", linewidth=0.7, zorder=9)

    a_track = evolution["a_au"]
    e_track = evolution["eccentricity"]
    orbit_ax.set_title("55P/Tempel-Tuttle dust: Poynting-Robertson contraction", fontsize=21, pad=12)
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
        x, y = orbit_xy(a_au, eccentricity)
        current_orbit.set_data(x, y)
        particle.set_offsets([particle_xy(a_au, eccentricity, mean_anomaly)])

        aphelion = a_au * (1.0 + eccentricity)
        view_radius = max(1.55, 1.10 * aphelion)
        orbit_ax.set_xlim(-view_radius, view_radius)
        orbit_ax.set_ylim(-view_radius, view_radius)
        time_text.set_text(f"{elapsed_years / 1000.0:,.1f} kyr")
        return current_orbit, particle, time_text

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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-prefix", default="tempel_tuttle_pr_drag")
    parser.add_argument("--mass-area", type=float, default=0.1, help="m/A in kg m^-2")
    parser.add_argument("--frames", type=int, default=240)
    parser.add_argument("--fps", type=int, default=24)
    args = parser.parse_args()

    prefix = Path(args.output_prefix).expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    solution = integrate_evolution(args.mass_area)
    evolution = resample_evolution(solution, args.frames)

    h5_path = prefix.with_suffix(".h5")
    mp4_path = prefix.with_suffix(".mp4")
    gif_path = prefix.with_suffix(".gif")
    poster_path = prefix.with_name(prefix.name + "_poster.png")
    save_hdf5(h5_path, evolution, args.mass_area)
    render_animation(evolution, args.mass_area, mp4_path, poster_path, args.fps)
    make_gif(mp4_path, gif_path)

    final_e = evolution["eccentricity"][-1]
    final_time = evolution["time_years"][-1]
    print(f"Reached a=1 AU after {final_time:,.1f} years with e={final_e:.6f}")
    print(f"Wrote {h5_path}")
    print(f"Wrote {mp4_path}")
    print(f"Wrote {gif_path}")
    print(f"Wrote {poster_path}")


if __name__ == "__main__":
    main()
