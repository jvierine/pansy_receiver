"""Combine independently reduced files, retaining lag counts and signed spectra."""

import argparse
from pathlib import Path
import h5py
import numpy as np


def combine(directory):
    paths = sorted(directory.glob("2026*.h5"))
    assert paths
    with h5py.File(paths[0]) as h:
        ranges = h["range_km"][()]
        fine = "fine_doppler_hz" in h
        frequency_key = "fine_doppler_hz" if fine else "doppler_hz"
        prefix = "fine_" if fine else ""
        frequencies = h[frequency_key][()]
    before, clean, count, means = [], [], [], []
    fractions = []
    for path in paths:
        with h5py.File(path) as h:
            np.testing.assert_array_equal(h["range_km"][()], ranges)
            np.testing.assert_array_equal(h[frequency_key][()], frequencies)
            before.append(
                [h[f"{prefix}channel_{c}/spectrum_before"][()] for c in range(5)]
            )
            clean.append(
                [h[f"{prefix}channel_{c}/spectrum_clean"][()] for c in range(5)]
            )
            count.append(h[f"{prefix}channel_2/valid_pair_count"][()])
            means.append(h["channel_2/removed_coherent_mean_iq"][()])
            mask = h["mask_reason"][()]
            fractions.append(
                [np.mean(mask != 0)] + [np.mean((mask & b) != 0) for b in [1, 2, 4, 8]]
            )
    # Equal-duration integrations, each individually pair-count normalized.
    before = np.mean(before, axis=0)
    clean = np.nanmean(clean, axis=0)
    count = np.sum(count, axis=0)
    # A local spectral baseline, not a calibrated receiver temperature.
    noise = np.nanmedian(clean[:, :, abs(frequencies) > 100], axis=-1)
    noise = np.maximum(noise, 1.0)
    output = directory / "summary.h5"
    with h5py.File(output, "w") as h:
        h["range_km"] = ranges
        h["doppler_hz"] = frequencies
        h["spectrum_before"] = before
        h["spectrum_clean"] = clean
        h["spectral_noise_baseline"] = noise
        h["before_db"] = 10 * np.log10(np.maximum(before / noise[:, :, None], 1e-3))
        h["clean_db"] = 10 * np.log10(np.maximum(clean / noise[:, :, None], 1e-3))
        h["valid_pair_count"] = count
        h["mask_fraction_per_file"] = fractions
        h["removed_coherent_mean_iq"] = means
        h.attrs["mask_bits"] = (
            "1 bad pulse/reference; 2 impulsive voltage; 4 known meteor; 8 transient target"
        )
        h.attrs["noise"] = (
            "median signed spectral power at |f|>100 Hz, separately for each range/channel; not Tsys calibration"
        )
        h.attrs["limitations"] = (
            "No proof of complete RFI/satellite rejection. Stationary mean removed; |f|<0.5 Hz not interpreted. Nominal ranges; even-IPP aliases unresolved. Five full-array receive looks with zenith transmission, not five transmit beams."
        )
        h.attrs["integration"] = (
            f"five independent 7.86432-s records; gaps are never bridged; {count.shape[1]}-lag Bartlett window, {len(frequencies)} frequency samples"
        )
    print(output)
    print(
        "Mean mask fractions [all,bad,impulse,meteor,transient]:",
        np.mean(fractions, axis=0),
    )
    print(
        "Lag-pair retention at range 20/50/100 km:",
        [count[np.argmin(abs(ranges - r)), [0, 127]].tolist() for r in [20, 50, 100]],
    )
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    combine(parser.parse_args().directory)
