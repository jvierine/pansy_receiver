#!/usr/bin/env python3
"""Orbit determination from PANSY ballistic_fit metadata."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import digital_rf as drf
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import get_body_barycentric_posvel
from astropy.time import Time

import pansy_ballistic as pbal


AU_KM = 149_597_870.700
MU_SUN_KM3_S2 = 1.32712440018e11

ORBIT_WRITER_ARGS = {
    "subdirectory_cadence_seconds": 3600,
    "file_cadence_seconds": 60,
    "samples_per_second_numerator": 1000000,
    "samples_per_second_denominator": 1,
    "file_name": "orbit",
}


def metadata_writer(path: str | Path) -> drf.DigitalMetadataWriter:
    Path(path).mkdir(parents=True, exist_ok=True)
    return drf.DigitalMetadataWriter(
        str(path),
        ORBIT_WRITER_ARGS["subdirectory_cadence_seconds"],
        ORBIT_WRITER_ARGS["file_cadence_seconds"],
        ORBIT_WRITER_ARGS["samples_per_second_numerator"],
        ORBIT_WRITER_ARGS["samples_per_second_denominator"],
        ORBIT_WRITER_ARGS["file_name"],
    )


def angle_deg(value_rad: float) -> float:
    return float(np.rad2deg(value_rad) % 360.0)


def angle_diff_deg(a: np.ndarray, ref: float) -> np.ndarray:
    return np.rad2deg(np.angle(np.exp(1j * np.deg2rad(a)) * np.exp(-1j * np.deg2rad(ref))))


def heliocentric_state_from_gcrs(state_gcrs_m_mps: np.ndarray, epoch_unix: float) -> tuple[np.ndarray, np.ndarray]:
    """Approximate heliocentric ICRS state from a geocentric GCRS state."""
    obstime = Time(float(epoch_unix), format="unix", scale="utc")
    earth_pos, earth_vel = get_body_barycentric_posvel("earth", obstime)
    sun_pos, sun_vel = get_body_barycentric_posvel("sun", obstime)
    earth_helio_pos_km = (earth_pos.xyz - sun_pos.xyz).to_value(u.km)
    earth_helio_vel_km_s = (earth_vel.xyz - sun_vel.xyz).to_value(u.km / u.s)
    state = np.asarray(state_gcrs_m_mps, dtype=np.float64)
    return earth_helio_pos_km + state[:3] / 1e3, earth_helio_vel_km_s + state[3:] / 1e3


def kepler_from_state(r_km: np.ndarray, v_km_s: np.ndarray) -> np.ndarray:
    r = np.asarray(r_km, dtype=np.float64)
    v = np.asarray(v_km_s, dtype=np.float64)
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)
    h = np.cross(r, v)
    h_norm = np.linalg.norm(h)
    n_vec = np.cross([0.0, 0.0, 1.0], h)
    n_norm = np.linalg.norm(n_vec)
    e_vec = np.cross(v, h) / MU_SUN_KM3_S2 - r / r_norm
    e = np.linalg.norm(e_vec)
    energy = 0.5 * v_norm**2 - MU_SUN_KM3_S2 / r_norm
    a_km = -MU_SUN_KM3_S2 / (2.0 * energy)
    inc = angle_deg(np.arccos(np.clip(h[2] / h_norm, -1.0, 1.0)))
    raan = 0.0 if n_norm == 0 else angle_deg(np.arctan2(n_vec[1], n_vec[0]))
    if n_norm == 0 or e == 0:
        argp = 0.0
    else:
        argp = angle_deg(np.arctan2(np.dot(np.cross(n_vec, e_vec), h) / h_norm, np.dot(n_vec, e_vec)))
    if e == 0:
        nu = 0.0
    else:
        nu = angle_deg(np.arctan2(np.dot(np.cross(e_vec, r), h) / h_norm, np.dot(e_vec, r)))
    q_au = a_km * (1.0 - e) / AU_KM
    return np.array([a_km / AU_KM, e, inc, raan, argp, nu, q_au], dtype=np.float64)


def summarize_samples(samples: np.ndarray, nominal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if len(samples) == 0:
        return np.full_like(nominal, np.nan), np.full((len(nominal), len(nominal)), np.nan)
    centered = samples.copy()
    for col in (2, 3, 4, 5):
        centered[:, col] = nominal[col] + angle_diff_deg(centered[:, col], nominal[col])
    std = np.nanstd(centered, axis=0)
    cov = np.cov(centered.T)
    return std, cov


def orbit_from_ballistic_fit(fit: dict, n_samples: int = 200, seed: int = 0) -> dict:
    epoch = float(fit["above_atmosphere_epoch_unix"])
    state_gcrs = np.asarray(fit["above_atmosphere_state_gcrs_m_mps"], dtype=np.float64)
    r_km, v_km_s = heliocentric_state_from_gcrs(state_gcrs, epoch)
    kepler = kepler_from_state(r_km, v_km_s)

    rng = np.random.default_rng(seed)
    cov = np.asarray(fit["parameter_covariance"], dtype=np.float64)
    params = np.asarray(fit["params"], dtype=np.float64)
    samples = []
    state_samples = []
    if cov.shape == (7, 7) and np.all(np.isfinite(cov)):
        cov = 0.5 * (cov + cov.T)
        rho_of_alt_m, _meta = pbal.density_interpolator(float(fit["t0_unix"]))
        for sample_params in rng.multivariate_normal(params, cov, size=n_samples):
            sample_params[6] = np.clip(sample_params[6], -4.0, 6.0)
            try:
                rt, rp, rv, _beta = pbal.reverse_to_above_atmosphere(sample_params, rho_of_alt_m)
                sample_epoch = float(fit["t0_unix"]) + float(rt[-1])
                sample_state = pbal.enu_state_to_gcrs(np.concatenate([rp[-1], rv[-1]]), sample_epoch)
                sr, sv = heliocentric_state_from_gcrs(sample_state, sample_epoch)
                samples.append(kepler_from_state(sr, sv))
                state_samples.append(sample_state)
            except Exception:
                continue
    samples = np.asarray(samples, dtype=np.float64)
    state_samples = np.asarray(state_samples, dtype=np.float64)
    kepler_std, kepler_cov = summarize_samples(samples, kepler)

    return {
        "sample_idx": int(fit["sample_idx"]),
        "t0_unix": float(fit["t0_unix"]),
        "epoch_unix": epoch,
        "reference_frame": np.asarray("GCRS", dtype="S8"),
        "heliocentric_frame": np.asarray("ICRS", dtype="S8"),
        "above_atmosphere_state_gcrs_m_mps": state_gcrs,
        "heliocentric_state_km_kms": np.concatenate([r_km, v_km_s]),
        "kepler": kepler,
        "kepler_names": np.asarray(["a_au", "e", "i_deg", "raan_deg", "argp_deg", "nu_deg", "q_au"], dtype="S16"),
        "kepler_std": kepler_std,
        "kepler_covariance": kepler_cov,
        "kepler_samples": samples,
        "gcrs_state_samples_m_mps": state_samples,
        "n_uncertainty_samples": int(len(samples)),
        "source_ballistic_coefficient_kg_m2": float(fit["ballistic_coefficient_kg_m2"]),
        "speed_increase_to_above_atmosphere_mps": float(fit["speed_increase_to_above_atmosphere_mps"]),
    }


def write_orbit(writer: drf.DigitalMetadataWriter, orbit: dict) -> None:
    writer.write(int(orbit["sample_idx"]), orbit)


def plot_orbit(orbit: dict, output_dir: str | Path) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    kep = orbit["kepler"]
    samples = orbit["kepler_samples"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)

    ax = axes[0, 0]
    if len(samples):
        ax.scatter(samples[:, 0], samples[:, 1], s=8, alpha=0.25, label="covariance samples")
    ax.scatter([kep[0]], [kep[1]], color="black", label="nominal")
    ax.set_xlabel("a (AU)")
    ax.set_ylabel("e")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    if len(samples):
        ax.scatter(samples[:, 6], samples[:, 2], s=8, alpha=0.25)
    ax.scatter([kep[6]], [kep[2]], color="black")
    ax.set_xlabel("q (AU)")
    ax.set_ylabel("i (deg)")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    labels = ["a (AU)", "e", "i (deg)", "Omega (deg)", "omega (deg)", "nu (deg)", "q (AU)"]
    rows = ["element        nominal        1-sigma"]
    for label, value, sigma in zip(labels, kep, orbit["kepler_std"]):
        rows.append(f"{label:<12s} {value:11.4g} {sigma:12.4g}")
    ax.axis("off")
    ax.set_title("Element uncertainty summary")
    ax.text(0.02, 0.98, "\n".join(rows), va="top", family="monospace", fontsize=9)

    ax = axes[1, 1]
    text = (
        f"epoch: {dt.datetime.fromtimestamp(orbit['epoch_unix'], tz=dt.timezone.utc).isoformat()}\n"
        f"a = {kep[0]:.3g} AU\n"
        f"e = {kep[1]:.3f}\n"
        f"i = {kep[2]:.2f} deg\n"
        f"q = {kep[6]:.3f} AU\n"
        f"n samples = {orbit['n_uncertainty_samples']}\n"
        f"delta-v above atmosphere = {orbit['speed_increase_to_above_atmosphere_mps'] / 1e3:.2f} km/s"
    )
    ax.axis("off")
    ax.text(0.02, 0.98, text, va="top", family="monospace")

    path = output_dir / f"orbit_fit_{orbit['sample_idx']}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def read_first_ballistic(path: str | Path) -> tuple[int, dict]:
    dm = drf.DigitalMetadataReader(str(path))
    bounds = dm.get_bounds()
    data = dm.read(bounds[0], bounds[1])
    if not data:
        raise RuntimeError(f"No ballistic_fit records found in {path}")
    key = sorted(data.keys())[0]
    return int(key), data[key]


def main() -> int:
    parser = argparse.ArgumentParser(description="Determine orbit from PANSY ballistic_fit metadata.")
    parser.add_argument("--ballistic-dir", type=Path, default=Path("data/metadata/ballistic_fit"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/metadata/orbital_parameters"))
    parser.add_argument("--plot-dir", type=Path, default=Path("plots/orbital_parameters"))
    parser.add_argument("--n-samples", type=int, default=200)
    parser.add_argument("--reanalyze", action="store_true", help="Delete an existing orbital_parameters key before writing.")
    args = parser.parse_args()
    _key, fit = read_first_ballistic(args.ballistic_dir)
    orbit = orbit_from_ballistic_fit(fit, n_samples=args.n_samples)
    if args.reanalyze:
        pbal.delete_metadata_key(args.output_dir, orbit["sample_idx"], ORBIT_WRITER_ARGS["file_name"])
    writer = metadata_writer(args.output_dir)
    write_orbit(writer, orbit)
    plot_path = plot_orbit(orbit, args.plot_dir)
    print(f"wrote orbit_fit metadata for {orbit['sample_idx']}")
    print(f"plot: {plot_path}")
    print(f"kepler [a_au,e,i,raan,argp,nu,q_au] = {orbit['kepler']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
