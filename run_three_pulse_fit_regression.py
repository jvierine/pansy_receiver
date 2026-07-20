#!/usr/bin/env python3
"""Run and score the fixed 50-event three-pulse trajectory-fit benchmark."""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

from three_pulse_fit_quality import QUALITY_SCALES


REPO = Path(__file__).resolve().parent
DEFAULT_SAMPLES = REPO / "three_pulse_regression_samples.txt"


def read_samples(path: Path) -> list[int]:
    return [
        int(line)
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def find_sample_file(directory: Path, sample_idx: int) -> Path:
    matches = sorted(directory.rglob(f"*{sample_idx}*.h5"))
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one HDF5 file for {sample_idx} below {directory}, found {len(matches)}"
        )
    return matches[0]


def catalogue_diagnostics_file(events_dir: Path, sample_idx: int) -> Path:
    event_date = dt.datetime.fromtimestamp(
        sample_idx / 1e6, tz=dt.timezone.utc
    ).date().isoformat()
    path = events_dir / event_date / f"pansy_disambiguation_diagnostics_{sample_idx}.h5"
    if not path.is_file():
        raise FileNotFoundError(f"catalogue diagnostics not found: {path}")
    return path


def fit_one(args: argparse.Namespace, sample_idx: int) -> Path:
    output = args.output_dir / "fits" / f"three_pulse_full_event_{sample_idx}.h5"
    if output.exists() and args.resume:
        return output
    diagnostics_h5 = (
        catalogue_diagnostics_file(args.events_dir, sample_idx)
        if args.events_dir is not None
        else find_sample_file(args.diagnostics_dir, sample_idx)
    )
    command = [
        args.python,
        str(REPO / "fit_full_event_three_pulse_complex_envelope.py"),
        "--sample-idx",
        str(sample_idx),
        "--base",
        str(args.base),
        "--diagnostics-h5",
        str(diagnostics_h5),
        "--initial-fit-h5",
        str(find_sample_file(args.initial_fit_dir, sample_idx)),
        "--prior-profile-h5",
        str(find_sample_file(args.prior_profile_dir, sample_idx)),
        "--output-dir",
        str(args.output_dir / "fits"),
    ]
    log = args.output_dir / "logs" / f"{sample_idx}.log"
    with log.open("w") as stream:
        subprocess.run(
            command,
            cwd=REPO,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return output


def read_quality(path: Path) -> dict[str, float]:
    with h5py.File(path, "r") as handle:
        group = handle["quality"]
        names = [*QUALITY_SCALES, "variance_weighted_score"]
        return {name: float(group.attrs[name]) for name in names}


def write_summary(output: Path, samples: list[int], rows: list[dict[str, float]]) -> None:
    names = [*QUALITY_SCALES, "variance_weighted_score"]
    temporary = output.with_suffix(output.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["schema"] = "pansy.three_pulse_fit_regression.v1"
        handle.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        handle.attrs["score_definition"] = "sum((RMS / reference_standard_deviation)**2)"
        handle.create_dataset("sample_idx", data=np.asarray(samples, dtype=np.int64))
        for name in names:
            values = np.asarray([row[name] for row in rows], dtype=float)
            handle.create_dataset(name, data=values)
            handle.attrs[f"median_{name}"] = float(np.median(values))
            handle.attrs[f"mean_{name}"] = float(np.mean(values))
            handle.attrs[f"p90_{name}"] = float(np.percentile(values, 90.0))
        for name, scale in QUALITY_SCALES.items():
            handle.attrs[f"reference_standard_deviation_{name}"] = scale
    temporary.replace(output)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path)
    parser.add_argument("--diagnostics-dir", type=Path)
    parser.add_argument(
        "--events-dir",
        type=Path,
        help="authoritative catalogue events tree; preferred over --diagnostics-dir",
    )
    parser.add_argument("--initial-fit-dir", type=Path)
    parser.add_argument("--prior-profile-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sample-list", type=Path, default=DEFAULT_SAMPLES)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    samples = read_samples(args.sample_list)
    if len(samples) != 50 or len(set(samples)) != 50:
        raise RuntimeError("the regression manifest must contain exactly 50 unique samples")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "fits").mkdir(exist_ok=True)
    (args.output_dir / "logs").mkdir(exist_ok=True)

    if not args.audit_only:
        required = (args.base, args.initial_fit_dir, args.prior_profile_dir)
        if any(path is None for path in required):
            parser.error("the input paths are required unless --audit-only is used")
        if args.events_dir is None and args.diagnostics_dir is None:
            parser.error("either --events-dir or --diagnostics-dir is required")
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(fit_one, args, sample): sample for sample in samples}
            for future in concurrent.futures.as_completed(futures):
                sample = futures[future]
                future.result()
                print(f"completed {sample}", flush=True)

    paths = [
        args.output_dir / "fits" / f"three_pulse_full_event_{sample}.h5"
        for sample in samples
    ]
    rows = [read_quality(path) for path in paths]
    write_summary(args.output_dir / "quality_summary.h5", samples, rows)
    scores = np.asarray([row["variance_weighted_score"] for row in rows])
    print(
        f"events={len(rows)} score_median={np.median(scores):.3f} "
        f"score_mean={np.mean(scores):.3f} score_p90={np.percentile(scores, 90):.3f}"
    )
    for name in QUALITY_SCALES:
        values = np.asarray([row[name] for row in rows])
        print(
            f"{name}: median={np.median(values):.6g} "
            f"mean={np.mean(values):.6g} p90={np.percentile(values, 90):.6g}"
        )
    (args.output_dir / "COMPLETE").write_text(
        dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds") + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
