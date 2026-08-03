#!/usr/bin/env python3
"""Run the established three-pulse mass fit over a fixed shower manifest."""

from __future__ import annotations

import argparse
import datetime as dt
import os
import subprocess
import sys
import traceback
from pathlib import Path

import h5py
import numpy as np

from fit_inter_pulse_acceleration_mass import fit_profile
from interferometer_alias_diagnostics import load_cut
from inter_pulse_phase_deceleration import (
    analyze_fit_observables,
    diagnostic_path,
    load_selected,
    write_fit_observables_h5,
)
from run_catalogue_mass_profiles import process_profile


REPO = Path(__file__).resolve().parent


def valid_h5(path: Path, required: str, schema: str | None = None) -> bool:
    if not path.is_file():
        return False
    try:
        with h5py.File(path, "r") as handle:
            if required not in handle:
                return False
            return schema is None or handle.attrs.get("schema", "") == schema
    except OSError:
        return False


def write_initial_fft(path: Path, sample_idx: int, observables: dict) -> None:
    decoded = observables["decoded"]
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["sample_idx"] = int(sample_idx)
        handle.attrs["schema"] = "pansy.three_pulse_fft_initialization.v1"
        handle.attrs["range_estimator"] = (
            "FFT-bandlimited ambiguity peak at 2x range interpolation and 16x Doppler padding"
        )
        handle.create_dataset("raw_idx", data=np.asarray(decoded["raw_idx"], dtype=np.int64))
        handle.create_dataset("range_km", data=np.asarray(decoded["range_km"], dtype=float))
        handle.create_dataset(
            "doppler_mps", data=np.asarray(decoded["coarse_doppler_mps"], dtype=float)
        )
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--worker-index",
        type=int,
        help="Worker index; defaults to the Open MPI or PMI rank.",
    )
    parser.add_argument(
        "--worker-count",
        type=int,
        help="Worker count; defaults to the Open MPI or PMI world size.",
    )
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    parser.add_argument("--max-events", type=int)
    args = parser.parse_args()
    worker_index = args.worker_index
    worker_count = args.worker_count
    if worker_index is None:
        for name in ("OMPI_COMM_WORLD_RANK", "PMI_RANK", "SLURM_PROCID"):
            if name in os.environ:
                worker_index = int(os.environ[name])
                break
    if worker_count is None:
        for name in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "SLURM_NTASKS"):
            if name in os.environ:
                worker_count = int(os.environ[name])
                break
    if worker_index is None or worker_count is None:
        parser.error(
            "worker index/count are required outside an MPI or scheduler launch"
        )
    if worker_count < 1 or not 0 <= worker_index < worker_count:
        parser.error("worker-index must satisfy 0 <= worker-index < worker-count")

    with h5py.File(args.manifest, "r") as handle:
        samples = np.asarray(handle["sample_idx"], dtype=np.int64)
    assigned = np.arange(worker_index, len(samples), worker_count)
    if args.max_events is not None:
        assigned = assigned[: args.max_events]

    baseline_dir = args.output_dir / "baseline_profiles"
    phase_dir = args.output_dir / "phase_profiles"
    initial_dir = args.output_dir / "initial_fft"
    fit_dir = args.output_dir / "three_pulse_fits"
    scratch_dir = args.output_dir / "scratch"
    log_dir = args.output_dir / "worker_logs"
    for directory in (
        baseline_dir,
        phase_dir,
        initial_dir,
        fit_dir,
        scratch_dir,
        log_dir,
    ):
        directory.mkdir(parents=True, exist_ok=True)

    log_path = log_dir / f"worker_{worker_index:03d}.tsv"
    if not log_path.exists():
        log_path.write_text("manifest_index\tsample_idx\tstatus\n")

    with log_path.open("a", buffering=1) as log:
        for manifest_index in assigned:
            sample_idx = int(samples[manifest_index])
            baseline = baseline_dir / f"mass_profile_{sample_idx}.h5"
            phase = phase_dir / f"mass_profile_phase_aware_{sample_idx}.h5"
            initial = initial_dir / f"highres_fft_i2_p16_{sample_idx}.h5"
            output = fit_dir / f"three_pulse_full_event_{sample_idx}.h5"
            observables_h5 = scratch_dir / f"fit_observables_{worker_index:03d}.h5"
            phase_partial = scratch_dir / f"phase_profile_{worker_index:03d}.h5"
            diagnostics = diagnostic_path(args.base, sample_idx)
            try:
                if valid_h5(
                    output,
                    "dynamics_refit/profile_radius_um",
                    "pansy.three_pulse_full_event.v1",
                ):
                    status = "already_complete"
                else:
                    if not diagnostics.is_file():
                        raise FileNotFoundError(f"missing diagnostics: {diagnostics}")
                    if not valid_h5(baseline, "result/free_best_mass_kg"):
                        process_profile(
                            diagnostics,
                            baseline,
                            grid_n=41,
                            radius_min_um=1.0,
                            radius_max_um=10000.0,
                            point_timeout_s=8,
                            max_nfev=45,
                        )

                    hypothesis = load_selected(diagnostics)
                    cut = load_cut(args.base / "metadata/cut", sample_idx)
                    observables = analyze_fit_observables(
                        cut, hypothesis, args.snr_threshold
                    )
                    write_fit_observables_h5(
                        sample_idx, observables, observables_h5
                    )
                    if not valid_h5(initial, "range_km"):
                        write_initial_fft(initial, sample_idx, observables)
                    if not valid_h5(phase, "result/best_radius_um"):
                        fit_profile(
                            diagnostics,
                            baseline,
                            observables_h5,
                            phase_partial,
                            None,
                            max_nfev=100,
                            global_starts=True,
                            adaptive_profile=True,
                            adaptive_spacing_dex=0.02,
                            adaptive_max_nfev=120,
                        )
                        os.replace(phase_partial, phase)

                    command = [
                        sys.executable,
                        str(REPO / "fit_full_event_three_pulse_complex_envelope.py"),
                        "--sample-idx",
                        str(sample_idx),
                        "--base",
                        str(args.base),
                        "--diagnostics-h5",
                        str(diagnostics),
                        "--initial-fit-h5",
                        str(initial),
                        "--prior-profile-h5",
                        str(phase),
                        "--output-dir",
                        str(fit_dir),
                        "--snr-threshold",
                        str(args.snr_threshold),
                    ]
                    subprocess.run(command, cwd=REPO, check=True)
                    if not valid_h5(
                        output,
                        "dynamics_refit/profile_radius_um",
                        "pansy.three_pulse_full_event.v1",
                    ):
                        raise RuntimeError("three-pulse fitter did not create a valid result")
                    for suffix in (".png",):
                        (fit_dir / f"three_pulse_full_event_{sample_idx}{suffix}").unlink(
                            missing_ok=True
                        )
                        (fit_dir / f"three_pulse_full_event_plot_{sample_idx}{suffix}").unlink(
                            missing_ok=True
                        )
                    status = "ok"
            except Exception as error:
                status = f"ERR {type(error).__name__}: {error}"
                traceback.print_exc()
            finally:
                observables_h5.unlink(missing_ok=True)
                phase_partial.unlink(missing_ok=True)
            timestamp = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
            print(
                timestamp,
                f"worker={worker_index:03d}",
                f"index={manifest_index}",
                f"sample={sample_idx}",
                status,
                flush=True,
            )
            log.write(f"{manifest_index}\t{sample_idx}\t{status}\n")

    (log_dir / f"worker_{worker_index:03d}.done").write_text(
        dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds") + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
