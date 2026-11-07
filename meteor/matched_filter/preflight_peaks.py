"""Check the peak ABI against independent NumPy complex matched filtering."""
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from meteor_rust_filter import RustMatchedFilter


def main():
    rng = np.random.default_rng(319)
    gates = np.arange(400, 1134, 2, dtype=np.uintp)
    durations = []
    for decimation in (1, 8):
        filt = RustMatchedFilter(gates, 132, decimation, 256)
        for _ in range(4):
            tx = np.exp(1j*np.pi*rng.integers(0, 2, 132)).astype(np.complex64)
            rx = (.1*(rng.normal(size=(7, 1600)) + 1j*rng.normal(size=(7, 1600)))).astype(np.complex64)
            phase = np.exp(1j*np.arange(7)*.43).astype(np.complex64)
            for gate, frequency in ((600, .013), (920, -.021)):
                rx[:, gate:gate+132] += 3*phase[:, None]*tx*np.exp(2j*np.pi*frequency*np.arange(132))
            t = time.perf_counter()
            peaks = filt.peaks(rx, tx, threshold=7, max_peaks=8)
            durations.append(time.perf_counter()-t)
            assert all(gate in peaks["range_sample"] for gate in (600, 920))
            power, median = filt.process(rx, tx)
            np.testing.assert_allclose(peaks["median_power"], median)
            for i, (gate, bin_) in enumerate(zip(peaks["range_sample"], peaks["doppler_bin"])):
                gate, bin_ = int(gate), int(bin_)
                product = rx[:, gate:gate+132]*tx.conj()
                product = product[:, :132//decimation*decimation].reshape(7, -1, decimation).sum(axis=2)
                expected = np.fft.fftshift(np.fft.fft(product, n=256, axis=1), axes=1)[:, bin_]
                np.testing.assert_allclose(peaks["receiver_amplitudes"][i], expected, rtol=2e-5, atol=2e-4)
                cell = power[np.flatnonzero(gates == gate)[0], bin_]
                np.testing.assert_allclose(np.sum(abs(peaks["receiver_amplitudes"][i])**2), cell, rtol=2e-6)
                np.testing.assert_allclose(peaks["excess_over_median"][i], cell/median-1, rtol=2e-6)
        for bad in (0, 257, -1):
            try:
                filt.peaks(rx, tx, max_peaks=bad)
            except ValueError:
                pass
            else:
                raise AssertionError("invalid capacity accepted")
        filt.close()
        try:
            filt.peaks(rx, tx)
        except RuntimeError:
            pass
        else:
            raise AssertionError("closed filter accepted")
    print(f"PASS: 8 two-target pulses, 7 complex receiver amplitudes per peak, "
          f"native and legacy decimation; mean {np.mean(durations)*1000:.3f} ms/pulse")


if __name__ == "__main__":
    main()
