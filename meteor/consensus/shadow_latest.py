"""Hourly bounded shadow sample, isolated from production metadata/checkpoints."""
import argparse
import fcntl
from pathlib import Path
import subprocess
import time

import h5py
import numpy as np


def sample_window(lower, upper, now):
    # Leave two seconds behind the newest batch to avoid a live writer.
    end = ((int(upper) - 2_000_000) // 2_000_000) * 2_000_000
    start = end - 60_000_000
    if start < lower or upper > now or now - upper > 3_600_000_000:
        return None
    return start, end


def prune(directory, keep=48, max_bytes=64 * 1024**2):
    files = sorted(directory.glob("hourly_*.h5"), reverse=True)
    total = 0
    for i, path in enumerate(files):
        total += path.stat().st_size
        if i >= keep or total > max_bytes:
            path.unlink()


def validate(path):
    points = members = tracks = 0
    max_r = max_v = 0.0
    with h5py.File(path, "r+") as data:
        for group in data["windows"].values():
            samples = group["sample_idx"][:]
            labels = group["track_id"][:]
            epochs = group["track_epoch_sample"][:]
            if np.any(labels < -1) or np.any(labels >= len(epochs)):
                raise ValueError("invalid membership label")
            points += len(samples)
            members += int(np.count_nonzero(labels >= 0))
            tracks += len(epochs)
            for k, epoch in enumerate(epochs):
                use = labels == k
                times = samples[use]
                if len(times) < 6 or len(np.unique(times)) != len(times):
                    raise ValueError("invalid track membership")
                dt = (times - epoch).astype(np.float64) / 1e6
                r0 = float(group["track_range_at_epoch_km"][k])
                v0 = float(group["track_radial_velocity_m_s"][k])
                a = float(group["track_radial_acceleration_m_s2"][k])
                dr = np.abs(group["range_km"][:][use] - r0 - (v0*dt + .5*a*dt**2)/1000)
                dv = np.abs(group["radial_velocity_m_s"][:][use] - v0 - a*dt)
                if not np.all(np.isfinite(dr)) or not np.all(np.isfinite(dv)):
                    raise ValueError("nonfinite residual")
                max_r = max(max_r, float(dr.max()))
                max_v = max(max_v, float(dv.max()))
        # Allow float32 storage rounding, not a relaxation of fitting gates.
        if max_r > .5001 or max_v > 3000.1:
            raise ValueError("stored model violates fit tolerance")
        data.attrs["candidate_points"] = np.int32(points)
        data.attrs["assigned_points"] = np.int32(members)
        data.attrs["track_fragments"] = np.int32(tracks)
        data.attrs["max_range_residual_km"] = np.float32(max_r)
        data.attrs["max_velocity_residual_m_s"] = np.float32(max_v)
        data.attrs["membership_checks_passed"] = np.int8(1)
    return points, members, tracks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / ".hourly.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        # Only our own bounded diagnostic products are eligible for cleanup.
        for partial in args.output.glob("hourly_*.partial"):
            partial.unlink()
        prune(args.output)
        import digital_rf
        reader = digital_rf.DigitalMetadataReader(str(args.metadata))
        if reader.get_samples_per_second() != 1e6:
            raise ValueError("expected 1 MHz metadata sample indices")
        window = sample_window(*reader.get_bounds(), time.time_ns() // 1000)
        if window is None:
            print("No recent complete shadow window; no production changes", flush=True)
            return
        start, end = window
        destination = args.output / f"hourly_{start}.h5"
        if destination.exists():
            print("Window already sampled", flush=True)
            return
        staging = args.output / f"hourly_{start}.partial"
        import sys
        began = time.monotonic()
        try:
            subprocess.run([sys.executable, str(Path(__file__).with_name("shadow_batch.py")),
                            str(args.metadata), str(start), str(end), str(staging)],
                           check=True, timeout=240)
            counts = validate(staging)
            with h5py.File(staging, "r+") as data:
                data.attrs["elapsed_s"] = np.float32(time.monotonic() - began)
            staging.replace(destination)
        finally:
            staging.unlink(missing_ok=True)
        prune(args.output)
        print(f"Shadow checks passed: candidates={counts[0]}, assigned={counts[1]}, "
              f"track_fragments={counts[2]}, elapsed={time.monotonic()-began:.3f}s", flush=True)


if __name__ == "__main__":
    main()
