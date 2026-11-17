"""Audit all declared raw IQ words; keep file padding and markers out of statistics."""

import argparse
from pathlib import Path
import h5py
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("raw", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(a.output, "w") as out:
        for path in sorted(a.raw.glob("20*.raw")):
            with path.open("rb") as f:
                header = f.read(8192)
            n = int(np.frombuffer(header, "<u4", 1, 80)[0])
            g = int(np.frombuffer(header, "<u4", 1, 240)[0])
            c = int(np.frombuffer(header, "<u2", 1, 426)[0])
            dtype = np.dtype(
                [("start", "<u4", (2,)), ("iq", "<i4", (g, c, 2)), ("end", "<u4", (2,))]
            )
            raw = np.memmap(path, dtype=dtype, mode="r", offset=8192, shape=(n,))
            lo = np.full((c, 2), 2**31, np.int64)
            hi = np.full((c, 2), -(2**31), np.int64)
            outside = np.zeros((c, 2), np.uint64)
            bits = np.zeros((c, 2), np.uint32)
            counts = np.zeros((c, 2), np.uint64)
            for start in range(0, n, 128):
                iq = np.asarray(raw["iq"][start : start + 128])
                lo = np.minimum(lo, np.min(iq, axis=(0, 1)))
                hi = np.maximum(hi, np.max(iq, axis=(0, 1)))
                outside += np.sum(
                    (iq < -32768) | (iq > 32767), axis=(0, 1), dtype=np.uint64
                )
                bits |= np.bitwise_or.reduce(
                    np.bitwise_or.reduce(
                        abs(iq.astype(np.int64)).astype(np.uint32), axis=0
                    ),
                    axis=0,
                )
                counts += iq.shape[0] * g
            group = out.create_group(path.name)
            for k, v in dict(
                minimum=lo,
                maximum=hi,
                outside_int16=outside,
                absolute_bitwise_or=bits,
                word_count=counts,
            ).items():
                group[k] = v
            print(
                path.name,
                "min/max",
                lo.min(),
                hi.max(),
                "out of int16",
                outside.sum(),
                "OR",
                hex(int(np.bitwise_or.reduce(bits.ravel()))),
                flush=True,
            )


if __name__ == "__main__":
    main()
