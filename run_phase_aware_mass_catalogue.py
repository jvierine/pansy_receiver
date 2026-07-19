#!/usr/bin/env python3
"""Refit a catalogue subset with beat-phase and interferometric phase-rate data."""

from __future__ import annotations

import argparse
import datetime as dt
import os
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


def completed_profile(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as handle:
            return "profile/chi2" in handle and "cross_phase_velocity" in handle
    except OSError:
        return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--baseline-profile-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--worker-index", type=int, required=True)
    parser.add_argument("--worker-count", type=int, required=True)
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    parser.add_argument("--max-events", type=int)
    args = parser.parse_args()
    if args.worker_count < 1 or not 0 <= args.worker_index < args.worker_count:
        parser.error("worker-index must satisfy 0 <= worker-index < worker-count")

    with h5py.File(args.manifest, "r") as handle:
        samples = np.asarray(handle["sample_idx"], dtype=np.int64)
    global_indices = np.arange(len(samples), dtype=int)[args.worker_index :: args.worker_count]
    if args.max_events is not None:
        global_indices = global_indices[: args.max_events]

    profiles = args.output_dir / "profiles"
    logs = args.output_dir / "worker_logs"
    scratch = args.output_dir / "scratch"
    for directory in (profiles, logs, scratch):
        directory.mkdir(parents=True, exist_ok=True)
    log_path = logs / f"worker_{args.worker_index:03d}.tsv"
    if not log_path.exists():
        log_path.write_text("global_index\tsample_idx\tstatus\n")

    with log_path.open("a") as log:
        for global_index in global_indices:
            sample_idx = int(samples[global_index])
            output = profiles / f"mass_profile_phase_aware_{sample_idx}.h5"
            if completed_profile(output):
                status = "already_complete"
                log.write(f"{global_index}\t{sample_idx}\t{status}\n")
                log.flush()
                continue
            baseline = args.baseline_profile_dir / f"mass_profile_{sample_idx}.h5"
            if not baseline.exists():
                status = "baseline_missing"
                log.write(f"{global_index}\t{sample_idx}\t{status}\n")
                log.flush()
                continue

            temporary_observables = scratch / f"fit_observables_{args.worker_index:03d}_{sample_idx}.h5"
            partial_output = scratch / f"mass_profile_{args.worker_index:03d}_{sample_idx}.partial.h5"
            try:
                diagnostics = diagnostic_path(args.base, sample_idx)
                hypothesis = load_selected(diagnostics)
                cut = load_cut(args.base / "metadata/cut", sample_idx)
                observables = analyze_fit_observables(cut, hypothesis, args.snr_threshold)
                write_fit_observables_h5(sample_idx, observables, temporary_observables)
                fit_profile(diagnostics, baseline, temporary_observables, partial_output, None)
                os.replace(partial_output, output)
                status = "ok"
            except Exception as exc:
                status = f"ERR {type(exc).__name__}: {exc}"
                traceback.print_exc()
            finally:
                temporary_observables.unlink(missing_ok=True)
                partial_output.unlink(missing_ok=True)
            print(
                dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                f"worker={args.worker_index:03d}",
                f"index={global_index}",
                f"sample={sample_idx}",
                status,
                flush=True,
            )
            log.write(f"{global_index}\t{sample_idx}\t{status}\n")
            log.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
