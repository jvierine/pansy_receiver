#!/usr/bin/env python3
"""DASST orbit determination for PANSY candidate first-detection states.

Input is the HDF5 candidate-state product written by
plot_interferometric_disambiguation.py.  The states are fitted first-detection
GCRS states; DASST is used to remove terrestrial gravity/zenithal attraction and
report heliocentric mean-ecliptic orbital elements.
"""

from __future__ import annotations

import argparse
import fcntl
from pathlib import Path

import digital_rf as drf
import h5py
import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation, GCRS, ITRS
from astropy.time import Time

import dasst


AU_M = 149_597_870_700.0
DEFAULT_KERNEL = Path("/Users/jvi019/src/rebound_examples/data/de430.bsp")
ORBIT_METADATA_WRITER_ARGS = {
    "subdirectory_cadence_seconds": 3600,
    "file_cadence_seconds": 60,
    "samples_per_second_numerator": 1000000,
    "samples_per_second_denominator": 1,
    "file_name": "orbit",
}


def gcrs_to_itrs_state(state_gcrs_m_mps: np.ndarray, epoch: Time) -> np.ndarray:
    rep = CartesianRepresentation(
        state_gcrs_m_mps[0] * u.m,
        state_gcrs_m_mps[1] * u.m,
        state_gcrs_m_mps[2] * u.m,
        differentials=CartesianDifferential(
            state_gcrs_m_mps[3] * u.m / u.s,
            state_gcrs_m_mps[4] * u.m / u.s,
            state_gcrs_m_mps[5] * u.m / u.s,
        ),
    )
    itrs = GCRS(rep, obstime=epoch).transform_to(ITRS(obstime=epoch))
    pos = itrs.cartesian.without_differentials().xyz.to_value(u.m)
    vel = itrs.cartesian.differentials["s"].d_xyz.to_value(u.m / u.s)
    return np.concatenate([pos, vel])


def circular_std_deg(values_deg: np.ndarray) -> float:
    vals = np.deg2rad(np.asarray(values_deg, dtype=np.float64))
    vals = vals[np.isfinite(vals)]
    if len(vals) < 2:
        return np.nan
    center = np.angle(np.nanmean(np.exp(1j * vals)))
    delta = (vals - center + np.pi) % (2.0 * np.pi) - np.pi
    return float(np.rad2deg(np.nanstd(delta, ddof=1)))


def dasst_kepler_from_gcrs(states_gcrs_m_mps: np.ndarray, epoch_unix: float, kernel: Path) -> np.ndarray:
    epoch = Time(float(epoch_unix), format="unix", scale="utc")
    itrs_states = np.asarray([gcrs_to_itrs_state(s, epoch) for s in states_gcrs_m_mps], dtype=np.float64).T
    od = dasst.orbit_determination.rebound_od(
        itrs_states,
        epoch,
        str(kernel),
        kepler_out_frame=["HeliocentricMeanEcliptic"],
        radiant_out_frame=["GCRS", "GeocentricMeanEcliptic"],
        termination_check=True,
        dt=10.0,
        max_t=7 * 24 * 3600.0,
        progress_bar=False,
    )
    return np.asarray(od["kepler_HeliocentricMeanEcliptic"][:, -1, :], dtype=np.float64)


def reordered_elements(kepler: np.ndarray) -> np.ndarray:
    """Return a, e, i, node, argp, nu, q from DASST/pyorb order.

    DASST/pyorb order is a, e, i, omega, Omega, true anomaly.
    """
    a_au = kepler[0] / AU_M
    e = kepler[1]
    inc = kepler[2]
    argp = kepler[3]
    node = kepler[4]
    nu = kepler[5]
    q_au = abs(a_au * (1.0 - e))
    return np.array([a_au, e, inc, node, argp, nu, q_au], dtype=np.float64)


def summarize_samples(kepler_samples: np.ndarray) -> tuple[np.ndarray, float]:
    if kepler_samples.size == 0:
        return np.full(7, np.nan), np.nan
    elems = np.asarray([reordered_elements(kepler_samples[:, i]) for i in range(kepler_samples.shape[1])])
    std = np.array(
        [
            float(np.nanstd(elems[:, 0], ddof=1)) if len(elems) > 1 else np.nan,
            float(np.nanstd(elems[:, 1], ddof=1)) if len(elems) > 1 else np.nan,
            circular_std_deg(elems[:, 2]),
            circular_std_deg(elems[:, 3]),
            circular_std_deg(elems[:, 4]),
            circular_std_deg(elems[:, 5]),
            float(np.nanstd(elems[:, 6], ddof=1)) if len(elems) > 1 else np.nan,
        ],
        dtype=np.float64,
    )
    return std, float(np.nanmean(elems[:, 1] > 1.0))


def orbital_element_covariance(samples: np.ndarray, nominal: np.ndarray) -> np.ndarray:
    if len(samples) < 2:
        return np.full((7, 7), np.nan)
    centered = np.asarray(samples, dtype=np.float64).copy()
    for col in (2, 3, 4, 5):
        centered[:, col] = nominal[col] + np.rad2deg(
            np.angle(np.exp(1j * np.deg2rad(centered[:, col] - nominal[col])))
        )
    return np.cov(centered.T)


def metadata_writer(path: Path) -> drf.DigitalMetadataWriter:
    path.mkdir(parents=True, exist_ok=True)
    return drf.DigitalMetadataWriter(
        str(path),
        ORBIT_METADATA_WRITER_ARGS["subdirectory_cadence_seconds"],
        ORBIT_METADATA_WRITER_ARGS["file_cadence_seconds"],
        ORBIT_METADATA_WRITER_ARGS["samples_per_second_numerator"],
        ORBIT_METADATA_WRITER_ARGS["samples_per_second_denominator"],
        ORBIT_METADATA_WRITER_ARGS["file_name"],
    )


def write_orbit_metadata(output_dir: Path, sample_key: int, payload: dict) -> None:
    writer = metadata_writer(output_dir)
    lock_path = output_dir / ".orbit_metadata.lock"
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            writer.write(int(sample_key), payload)
        except ValueError as exc:
            if "name already exists" not in str(exc):
                raise
            writer = metadata_writer(output_dir)
            writer.write(int(sample_key), payload)
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run DASST orbits for PANSY candidate GCRS state samples.")
    parser.add_argument("state_h5", type=Path)
    parser.add_argument("--kernel", type=Path, default=DEFAULT_KERNEL)
    parser.add_argument("--output-h5", type=Path, default=None)
    parser.add_argument("--metadata-dir", type=Path, default=None)
    args = parser.parse_args()

    rows = []
    unc_rows = []
    sample_rows = {}
    with h5py.File(args.state_h5, "r") as h:
        labels = sorted(h.keys(), key=lambda k: int(h[k].attrs["combined_rank"]))
        for label in labels:
            grp = h[label]
            epoch = float(grp.attrs["epoch_unix"])
            nominal_state = grp["state_gcrs_m_mps"][()][None, :]
            nominal_kep = dasst_kepler_from_gcrs(nominal_state, epoch, args.kernel)[:, 0]
            elems = reordered_elements(nominal_kep)
            sample_states = grp["state_gcrs_samples_m_mps"][()]
            sample_states = sample_states[np.all(np.isfinite(sample_states), axis=1)]
            if len(sample_states) > 0:
                sample_kep = dasst_kepler_from_gcrs(sample_states, epoch, args.kernel)
                std, frac_e_gt_1 = summarize_samples(sample_kep)
                sample_elems = np.asarray([reordered_elements(sample_kep[:, i]) for i in range(sample_kep.shape[1])])
                elem_cov = orbital_element_covariance(sample_elems, elems)
            else:
                std = np.full(7, np.nan)
                frac_e_gt_1 = np.nan
                sample_elems = np.empty((0, 7), dtype=np.float64)
                elem_cov = np.full((7, 7), np.nan)
            attrs = {key: grp.attrs[key] for key in grp.attrs.keys()}
            rows.append((label, attrs, elems))
            unc_rows.append((label, len(sample_states), float(attrs["sigma_log10_beta"]), std, frac_e_gt_1, elem_cov))
            sample_rows[label] = sample_elems

    print(
        "dasst_candidate_orbit_columns "
        "hypothesis_id candidate_number combined_rank combined_score log10_beta_kg_m2 "
        "first_alt_km a_au e i_deg raan_deg argp_deg nu_deg q_au"
    )
    print(
        "dasst_candidate_orbit_uncertainty_columns "
        "hypothesis_id n_samples sigma_log10_beta sigma_a_au sigma_e sigma_i_deg "
        "sigma_raan_deg sigma_argp_deg sigma_nu_deg sigma_q_au frac_e_gt_1"
    )
    for label, attrs, elems in rows:
        print(
            "dasst_candidate_orbit",
            label,
            int(attrs["candidate_number"]),
            int(attrs["combined_rank"]),
            f"{float(attrs['combined_score']):.3f}",
            f"{float(attrs['log10_beta_kg_m2']):.3f}",
            f"{float(attrs['first_alt_km']):.3f}",
            f"{elems[0]:.6g}",
            f"{elems[1]:.6f}",
            f"{elems[2]:.3f}",
            f"{elems[3]:.3f}",
            f"{elems[4]:.3f}",
            f"{elems[5]:.3f}",
            f"{elems[6]:.6g}",
        )
    for label, n_samples, sigma_beta, std, frac_e_gt_1, _elem_cov in unc_rows:
        print(
            "dasst_candidate_orbit_uncertainty",
            label,
            n_samples,
            f"{sigma_beta:.3f}",
            f"{std[0]:.6g}",
            f"{std[1]:.6f}",
            f"{std[2]:.3f}",
            f"{std[3]:.3f}",
            f"{std[4]:.3f}",
            f"{std[5]:.3f}",
            f"{std[6]:.6g}",
            f"{frac_e_gt_1:.3f}",
        )

    if args.output_h5 is not None:
        args.output_h5.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(args.output_h5, "w") as h:
            h.attrs["source_program"] = "dasst_orbits_from_candidate_states.py"
            h.attrs["input_state_h5"] = str(args.state_h5)
            h.attrs["kernel"] = str(args.kernel)
            h.attrs["orbit_propagator"] = "dasst.orbit_determination.rebound_od"
            for (label, attrs, elems), (_label, n_samples, sigma_beta, std, frac_e_gt_1, elem_cov) in zip(rows, unc_rows):
                grp = h.create_group(label)
                for key, val in attrs.items():
                    grp.attrs[key] = val
                grp.attrs["n_uncertainty_samples"] = int(n_samples)
                grp.attrs["frac_e_gt_1"] = float(frac_e_gt_1)
                grp.attrs["sigma_log10_beta"] = float(sigma_beta)
                grp.create_dataset("kepler", data=elems)
                grp.create_dataset("kepler_names", data=np.asarray(["a_au", "e", "i_deg", "raan_deg", "argp_deg", "nu_deg", "q_au"], dtype="S16"))
                grp.create_dataset("kepler_std", data=std)
                grp.create_dataset("kepler_covariance", data=elem_cov)
                grp.create_dataset("kepler_samples", data=sample_rows[label])
        print(f"dasst_orbit_h5 {args.output_h5}")

    if args.metadata_dir is not None and rows:
        best_idx = int(np.argmin([int(attrs["combined_rank"]) for _label, attrs, _elems in rows]))
        label, attrs, elems = rows[best_idx]
        _ulabel, n_samples, sigma_beta, std, frac_e_gt_1, elem_cov = unc_rows[best_idx]
        with h5py.File(args.state_h5, "r") as h:
            state_grp = h[label]
            sample_key = int(round(float(state_grp.attrs["epoch_unix"]) * 1_000_000.0))
            payload = {
                "sample_idx": int(sample_key),
                "source_cut_sample_idx": int(round(float(h.attrs["sample_epoch_unix"]) * 1_000_000.0)),
                "selected_hypothesis": np.asarray(label, dtype="S8"),
                "candidate_number": int(attrs["candidate_number"]),
                "combined_rank": int(attrs["combined_rank"]),
                "combined_score": float(attrs["combined_score"]),
                "state_epoch": np.asarray("first_detection", dtype="S32"),
                "reference_frame": np.asarray("GCRS", dtype="S8"),
                "initial_state_gcrs_m_mps": state_grp["state_gcrs_m_mps"][()],
                "initial_state_samples_gcrs_m_mps": state_grp["state_gcrs_samples_m_mps"][()],
                "fit_parameter_names": np.asarray(
                    ["east_m", "north_m", "up_m", "ve_mps", "vn_mps", "vu_mps", "log10_beta_kg_m2"],
                    dtype="S32",
                ),
                "fit_parameters": state_grp["ballistic_params"][()],
                "fit_parameter_covariance": state_grp["ballistic_parameter_covariance"][()],
                "log10_beta_kg_m2": float(attrs["log10_beta_kg_m2"]),
                "sigma_log10_beta": float(sigma_beta),
                "orbit_frame": np.asarray("HeliocentricMeanEcliptic", dtype="S32"),
                "orbit_propagator": np.asarray("dasst.orbit_determination.rebound_od", dtype="S64"),
                "kepler_names": np.asarray(["a_au", "e", "i_deg", "raan_deg", "argp_deg", "nu_deg", "q_au"], dtype="S16"),
                "kepler": elems,
                "kepler_std": std,
                "kepler_covariance": elem_cov,
                "kepler_samples": sample_rows[label],
                "n_uncertainty_samples": int(n_samples),
                "frac_e_gt_1": float(frac_e_gt_1),
                "radiant_frame": np.asarray("GCRS", dtype="S8"),
                "radiant_ra_deg": float(np.nan),
                "radiant_dec_deg": float(np.nan),
                "radiant_speed_km_s": float(np.linalg.norm(state_grp["state_gcrs_m_mps"][()][3:]) / 1e3),
                "mass_density_assumption_g_cm3": float(3.0),
                "mass_estimate_kg": float(np.nan),
                "mass_estimate_note": np.asarray("placeholder; requires calibrated RCS/ablation model", dtype="S64"),
            }
        write_orbit_metadata(args.metadata_dir, sample_key, payload)
        print(f"orbit_metadata_write {args.metadata_dir} {sample_key}")


if __name__ == "__main__":
    main()
