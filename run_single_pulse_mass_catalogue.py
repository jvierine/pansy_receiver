#!/usr/bin/env python3
"""Run the Figure-3 single-pulse raw-voltage mass fit over a catalogue manifest."""

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

from inter_pulse_phase_deceleration import diagnostic_path


def completed_profile(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as handle:
            return (
                handle.attrs.get("schema", "")
                == "pansy.catalogue_mass_profile.single_pulse_raw_voltage.v1"
                and "result/free_best_radius_um" in handle
                and "quality/path_length_km" in handle
            )
    except OSError:
        return False


def first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--initial-fit-dir", type=Path, required=True)
    parser.add_argument("--prior-profile-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--worker-index", type=int, required=True)
    parser.add_argument("--worker-count", type=int, required=True)
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    parser.add_argument("--max-events", type=int)
    parser.add_argument("--frequency-half-width-hz", type=float, default=1.5e4)
    parser.add_argument("--frequency-grid-size", type=int, default=2001)
    parser.add_argument("--keep-event-h5", action="store_true")
    args = parser.parse_args()
    if args.worker_count < 1 or not 0 <= args.worker_index < args.worker_count:
        parser.error("worker-index must satisfy 0 <= worker-index < worker-count")

    with h5py.File(args.manifest, "r") as handle:
        samples = np.asarray(handle["sample_idx"], dtype=np.int64)
    global_indices = np.arange(len(samples), dtype=int)[args.worker_index :: args.worker_count]
    if args.max_events is not None:
        global_indices = global_indices[: args.max_events]

    profiles = args.output_dir / "profiles"
    event_h5_dir = args.output_dir / "single_pulse_events"
    logs = args.output_dir / "worker_logs"
    for directory in (profiles, event_h5_dir, logs):
        directory.mkdir(parents=True, exist_ok=True)
    log_path = logs / f"worker_{args.worker_index:03d}.tsv"
    if not log_path.exists():
        log_path.write_text("global_index\tsample_idx\tstatus\n")

    script = Path(__file__).with_name("fit_full_event_single_pulse_raw_voltage.py")
    with log_path.open("a", buffering=1) as log:
        for global_index in global_indices:
            sample_idx = int(samples[global_index])
            output = profiles / f"mass_profile_{sample_idx}.h5"
            event_h5 = event_h5_dir / f"single_pulse_raw_voltage_{sample_idx}.h5"
            event_png = event_h5_dir / f"single_pulse_raw_voltage_plot_{sample_idx}.png"
            if completed_profile(output):
                status = "already_complete"
                log.write(f"{global_index}\t{sample_idx}\t{status}\n")
                continue

            diagnostics = diagnostic_path(args.base, sample_idx)
            initial_fit = first_existing(
                [
                    args.initial_fit_dir / f"fit_observables_catalogue_{sample_idx}.h5",
                    args.initial_fit_dir / f"highres_fft_i2_p16_{sample_idx}.h5",
                ]
            )
            prior_profile = first_existing(
                [
                    args.prior_profile_dir / f"mass_profile_phase_aware_{sample_idx}.h5",
                    args.prior_profile_dir / f"mass_profile_{sample_idx}.h5",
                ]
            )
            if not diagnostics.exists():
                status = "diagnostics_missing"
            elif initial_fit is None:
                status = "initial_fit_missing"
            elif prior_profile is None:
                status = "prior_profile_missing"
            else:
                command = [
                    sys.executable,
                    str(script),
                    "--sample-idx",
                    str(sample_idx),
                    "--base",
                    str(args.base),
                    "--diagnostics-h5",
                    str(diagnostics),
                    "--initial-fit-h5",
                    str(initial_fit),
                    "--prior-profile-h5",
                    str(prior_profile),
                    "--output-dir",
                    str(event_h5_dir),
                    "--mass-profile-output-dir",
                    str(profiles),
                    "--snr-threshold",
                    str(args.snr_threshold),
                    "--frequency-half-width-hz",
                    str(args.frequency_half_width_hz),
                    "--frequency-grid-size",
                    str(args.frequency_grid_size),
                    "--no-plot",
                ]
                try:
                    subprocess.run(command, check=True)
                    status = "ok"
                    if not args.keep_event_h5:
                        event_h5.unlink(missing_ok=True)
                        event_png.unlink(missing_ok=True)
                except Exception as exc:
                    status = f"ERR {type(exc).__name__}: {exc}"
                    traceback.print_exc()
            print(
                dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                f"worker={args.worker_index:03d}",
                f"index={global_index}",
                f"sample={sample_idx}",
                status,
                flush=True,
            )
            log.write(f"{global_index}\t{sample_idx}\t{status}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
