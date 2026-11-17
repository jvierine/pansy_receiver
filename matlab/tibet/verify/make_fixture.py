"""Generate deterministic binary radar test data, not substitute observations."""

import argparse
from pathlib import Path
import struct
import h5py
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("header_cut", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with h5py.File(args.header_cut) as h:
        header = bytearray(h["raw_header"][()].tobytes())
    n = 512
    struct.pack_into("<I", header, 80, n)
    struct.pack_into("<H", header, 42, 1)
    rng = np.random.default_rng(5617)
    a = np.repeat([1, -1, 1, 1, 1, -1, -1, -1], 4)
    b = np.repeat([1, -1, 1, 1, -1, 1, 1, 1], 4)
    wavelength = 299792458 / (struct.unpack_from("<I", header, 84)[0] * 120e6 / 2**32)
    frequency = 80000 / wavelength
    phases = np.exp(1j * np.linspace(-1.2, 1.4, 15))
    gain = np.array([0.7, 0.6, 3.0, 0.8, 0.7] + [0.8] * 10)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as f:
        f.write(header)
        for pulse in range(n):
            z = rng.normal(0, 1500, (640, 15)) + 1j * rng.normal(0, 1500, (640, 15))
            code = a if pulse % 2 == 0 else b
            common = np.exp(1j * 0.02 * pulse)
            z[7:39, 2] += 12000 * common * code
            if 40 <= pulse <= 250:
                gate = 355 - int(round(40000 * pulse * 0.00128 / (299792458 / 1e6)))
                envelope = 12000 * np.exp(-0.5 * ((pulse - 145) / 45) ** 2)
                echo = code * np.exp(
                    2j * np.pi * frequency * (pulse * 0.00128 + np.arange(32) / 500000)
                )
                z[gate : gate + 32] += envelope * common * echo[:, None] * phases * gain
            iq = np.stack([np.rint(z.real), np.rint(z.imag)], axis=-1).astype("<i4")
            f.write(struct.pack("<II", 0xAABBCCDD, pulse + 1))
            f.write(iq.tobytes())
            f.write(struct.pack("<II", 0xAABBCCDD, pulse + 1))
    print(args.output)


if __name__ == "__main__":
    main()
