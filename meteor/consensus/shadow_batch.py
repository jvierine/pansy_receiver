"""Finite shadow replay: native reading/fitting, one HDF5 result, no live writes."""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("metadata", type=Path)
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    p.add_argument("output", type=Path)
    p.add_argument("--binary", type=Path, default=Path(__file__).parent/"target/release/pansy-meteor-consensus")
    args = p.parse_args()
    if not 0 < args.end-args.start <= 300_000_000:
        p.error("request an explicit interval no longer than 300 s (1 MHz indices)")
    if args.output.exists():
        p.error("output already exists")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_suffix(f".{os.getpid()}.partial")
    try:
        with tempfile.TemporaryDirectory(prefix="pansy-shadow-") as tmp, h5py.File(partial, "w") as out:
            out.attrs["schema_version"] = np.int32(1)
            out.attrs["shadow_only"] = np.int8(1)
            out.attrs["window_seconds"] = np.float32(2)
            out.attrs["seed_snr_db"] = np.float32(10*np.log10(7))
            out.attrs["grow_snr_db"] = np.float32(10*np.log10(2))
            out.attrs["tx_power_min"] = np.float32(1e9)
            out.attrs["range_tolerance_km"] = np.float32(0.5)
            out.attrs["radial_velocity_tolerance_m_s"] = np.float32(3000)
            out.attrs["min_points"] = np.int32(6)
            out.attrs["min_seed_points"] = np.int32(4)
            out.attrs["min_duration_s"] = np.float32(.06)
            out.attrs["max_gap_s"] = np.float32(.15)
            out.attrs["trials"] = np.int32(512)
            windows = out.create_group("windows")
            for start in range(args.start, args.end, 2_000_000):
                result = Path(tmp)/f"{start}.h5"
                end = min(start+2_000_000, args.end)
                subprocess.run([str(args.binary), str(args.metadata), str(start), str(end), "--output", str(result)], check=True, timeout=30)
                group = windows.create_group(str(start))
                with h5py.File(result) as source:
                    for key in source:
                        source.copy(key, group)
            out.flush()
        os.replace(partial, args.output)
    finally:
        partial.unlink(missing_ok=True)
    print(f"Wrote shadow results: {args.output}", flush=True)


if __name__ == "__main__":
    main()
