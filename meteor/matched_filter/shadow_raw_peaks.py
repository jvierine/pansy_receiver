"""Explicit 20-pulse raw replay with independent complex-FFT verification."""
import argparse
from pathlib import Path
import sys
import time

import digital_rf
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from meteor_rust_filter import RustMatchedFilter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("raw")
    parser.add_argument("mf")
    parser.add_argument("sample", type=int, help="Exact existing MF batch key")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output exists")
    records = digital_rf.DigitalMetadataReader(args.mf).read(args.sample, args.sample,
                                                          ["tx_idxs", "tx_pwr"])
    record = records[args.sample]
    indices = np.asarray(record["tx_idxs"], dtype=np.int64)
    if not np.array_equal(indices, args.sample + np.arange(20)*1600):
        raise ValueError("not an aligned mesosphere batch")
    if not np.all(np.isfinite(record["tx_pwr"])) or not np.all(np.asarray(record["tx_pwr"]) > 1e9):
        raise ValueError("failed TX gate")
    reader = digital_rf.DigitalRFReader(args.raw)
    with h5py.File(ROOT/"data/mesocal.h5") as cal:
        pwr = cal["pwr"][:]
    scale = np.real(np.sqrt(pwr[0])/np.sqrt(pwr[:7]))
    rx = np.array([reader.read_vector_1d(args.sample, 32000, f"ch{ch:03}")*scale[ch]
                   for ch in range(7)], dtype=np.complex64)
    tx = np.asarray(reader.read_vector_1d(args.sample, 32000, "ch007"), dtype=np.complex64)
    gates = np.arange(int(np.floor(60/(299792458/2e9))),
                      int(np.floor(170/(299792458/2e9)))+1, 2)
    filt = RustMatchedFilter(gates, 132, 8, 256)
    partial = args.output.with_suffix(".partial")
    started = time.monotonic()
    count = 0
    try:
        with h5py.File(partial, "x") as data:
            data.attrs["shadow_only"] = np.int8(1)
            data.attrs["phase_calibrated"] = np.int8(0)
            data.attrs["decimation"] = np.int32(8)
            data.attrs["nfft"] = np.int32(256)
            data.attrs["sample_rate_hz"] = np.float32(1e6)
            data.attrs["receiver_amplitude_scale"] = scale.astype(np.float32)
            data.attrs["threshold_excess_over_median"] = np.float32(7)
            for i, sample in enumerate(indices):
                voltage = rx[:, i*1600:(i+1)*1600]
                code = tx[i*1600:i*1600+132]
                if np.sum(abs(code)**2) <= 1e9:
                    raise ValueError("raw TX gate failed")
                peaks = filt.peaks(voltage, code, threshold=7)
                for k, (gate, bin_) in enumerate(zip(peaks["range_sample"], peaks["doppler_bin"])):
                    gate, bin_ = int(gate), int(bin_)
                    product = (voltage[:, gate:gate+132]*code.conj())[:, :128].reshape(7, 16, 8).sum(axis=2)
                    expected = np.fft.fftshift(np.fft.fft(product, n=256, axis=1), axes=1)[:, bin_]
                    error = np.max(abs(peaks["receiver_amplitudes"][k]-expected))
                    if error > max(1e-3, 3e-6*np.max(abs(expected))):
                        raise ValueError("complex matched-filter mismatch")
                group = data.create_group(str(sample))
                group["sample_idx"] = np.int64(sample)
                group["range_sample"] = peaks["range_sample"].astype(np.int32)
                group["doppler_bin"] = peaks["doppler_bin"].astype(np.int32)
                group["excess_over_median"] = peaks["excess_over_median"]
                group["receiver_amplitudes"] = peaks["receiver_amplitudes"]
                group["median_power"] = np.float32(peaks["median_power"])
                count += len(peaks["range_sample"])
            data.attrs["complex_reference_passed"] = np.int8(1)
            data.attrs["elapsed_s"] = np.float32(time.monotonic()-started)
        partial.replace(args.output)
    finally:
        filt.close()
        partial.unlink(missing_ok=True)
    print(f"PASS: 20 station pulses, {count} candidate peaks, complex amplitudes verified; "
          f"elapsed {time.monotonic()-started:.3f}s", flush=True)


if __name__ == "__main__":
    main()
