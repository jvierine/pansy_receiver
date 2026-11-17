"""Independent direct-sum checks of native paired-code and missing-data spectra."""

import argparse
from pathlib import Path
import h5py
import numpy as np


def verify(directory, raw_directory):
    decode_errors, acf_errors = [], []
    for path in sorted(directory.glob("2026*.h5")):
        with h5py.File(path) as h:
            raw = np.memmap(
                raw_directory / (path.stem + ".raw"),
                dtype="<i4",
                mode="r",
                offset=8192,
                shape=(6144, 19204),
            )
            values = h["decoded_iq"][()]
            decoded = values[..., 0].astype(float) + 1j * values[..., 1]
            mask = h["mask_reason"][()] == 0
            for pair in [0, 511, 1537, 3071]:
                for r in [0, 8, 31, 75, 140]:
                    gate = 39 + 4 * r
                    expected = np.zeros(5, complex)
                    for pulse in [2 * pair, 2 * pair + 1]:
                        iq = raw[pulse, 2:-2].reshape(640, 15, 2)
                        z = iq[..., 0].astype(float) + 1j * iq[..., 1]
                        tx = z[7:39, 2]
                        expected += (
                            np.sum(z[gate : gate + 32, :5] * tx[:, None].conj(), axis=0)
                            / np.linalg.norm(tx)
                            / np.sqrt(2)
                        )
                    error = np.max(abs(expected - decoded[pair, :, r])) / max(
                        np.max(abs(expected)), 1.0
                    )
                    decode_errors.append(error)
                    assert error < 2e-5, (path, pair, r, error)
            for c in [0, 2, 4]:
                spectrum = h[f"fine_channel_{c}/spectrum_clean"][()]
                counts = h[f"fine_channel_{c}/valid_pair_count"][()]
                recovered = np.fft.ifft(np.fft.ifftshift(spectrum, axes=-1), axis=-1)[
                    :, :1024
                ]
                recovered /= 1 - np.arange(1024)[None, :] / 1024
                for r in [8, 31, 75]:
                    z = decoded[:, c, r]
                    good = mask[:, r]
                    z = z - z[good].mean()
                    scale = np.mean(abs(z[good]) ** 2)
                    for lag in [0, 1, 127, 511, 1023]:
                        valid = good[lag:] & good[: len(z) - lag]
                        assert counts[r, lag] == valid.sum()
                        direct = np.mean((z[lag:] * z[: len(z) - lag].conj())[valid])
                        error = abs(direct - recovered[r, lag]) / max(scale, 1.0)
                        acf_errors.append(error)
                        assert error < 2e-4, (path, c, r, lag, error)
    with h5py.File(directory / "verification.h5", "w") as h:
        h["decoded_pair_relative_errors"] = decode_errors
        h["acf_relative_errors"] = acf_errors
        h.attrs["scope"] = (
            "Independent direct voltage sums, native FFT ACF against direct valid-pair sums, and exact lag counts. Not validation of complete interference rejection or atmospheric origin."
        )
    print(
        "PASS: max decoded relative error",
        max(decode_errors),
        "; max ACF relative error",
        max(acf_errors),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("raw_directory", type=Path)
    args = parser.parse_args()
    verify(args.directory, args.raw_directory)
