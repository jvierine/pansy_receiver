"""Optional native kernel; metadata and candidate selection remain unchanged."""
import ctypes as ct
import os
from pathlib import Path
import sys
import threading

import numpy as np


class RustMatchedFilter:
    def __init__(self, gates, txlen, decimation, nfft):
        suffix = "dylib" if sys.platform == "darwin" else "so"
        default = Path(__file__).parent / "meteor/matched_filter/target/release" / ("libpansy_meteor_matched_filter." + suffix)
        self.lib = ct.CDLL(os.environ.get("PANSY_MF_RUST_LIBRARY", str(default)))
        self.lib.pansy_mf_new.argtypes = [ct.c_void_p] + [ct.c_size_t]*4
        self.lib.pansy_mf_new.restype = ct.c_void_p
        self.lib.pansy_mf_process.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_size_t, ct.c_size_t, ct.c_void_p, ct.c_size_t, ct.c_void_p, ct.c_size_t, ct.c_void_p]
        self.lib.pansy_mf_process.restype = ct.c_int
        self.lib.pansy_mf_free.argtypes = [ct.c_void_p]
        self.lib.pansy_mf_free.restype = None
        gates = np.ascontiguousarray(gates, dtype=np.uintp)
        self.shape = (len(gates), int(nfft))
        self.txlen = int(txlen)
        self.lock = threading.Lock()
        self.handle = self.lib.pansy_mf_new(gates.ctypes.data, len(gates), self.txlen, int(decimation), int(nfft))
        if not self.handle:
            raise ValueError("Rust matched-filter initialization failed")

    def process(self, rx, tx):
        rx = np.ascontiguousarray(rx, dtype=np.complex64)
        tx = np.ascontiguousarray(tx, dtype=np.complex64)
        if rx.ndim != 2 or tx.shape != (self.txlen,):
            raise ValueError("Invalid RX/TX dimensions")
        power = np.empty(self.shape, dtype=np.float32)
        median = ct.c_float()
        with self.lock:
            if not self.handle:
                raise RuntimeError("Filter is closed")
            status = self.lib.pansy_mf_process(self.handle, rx.ctypes.data, *rx.shape, tx.ctypes.data, len(tx), power.ctypes.data, power.size, ct.byref(median))
        if status:
            raise RuntimeError(f"Rust matched-filter rejected pulse: status {status}")
        return power, median.value

    def close(self):
        with self.lock:
            if self.handle:
                self.lib.pansy_mf_free(self.handle)
                self.handle = None

    def peaks(self, rx, tx, threshold=2.0, max_peaks=8, range_guard=2, doppler_guard=2):
        """Diagnostic candidates with uncalibrated complex receiver amplitudes.

        Threshold is linear excess over the legacy median, NOT calibrated SNR.
        No phase correction or interferometric alias selection is performed.
        """
        for value in (max_peaks, range_guard, doppler_guard):
            if not isinstance(value, (int, np.integer)) or value < 0:
                raise ValueError("Peak capacity and guards must be nonnegative integers")
        if not 1 <= max_peaks <= 256 or not np.isfinite(threshold) or threshold < 0:
            raise ValueError("Invalid peak capacity or threshold")
        rx = np.ascontiguousarray(rx, dtype=np.complex64)
        tx = np.ascontiguousarray(tx, dtype=np.complex64)
        if rx.ndim != 2 or rx.shape[0] == 0 or tx.shape != (self.txlen,):
            raise ValueError("Invalid RX/TX dimensions")
        gates = np.empty(max_peaks, dtype=np.uintp)
        bins = np.empty(max_peaks, dtype=np.uintp)
        excess = np.empty(max_peaks, dtype=np.float32)
        amplitudes = np.empty((max_peaks, rx.shape[0]), dtype=np.complex64)
        count, median = ct.c_size_t(), ct.c_float()
        with self.lock:
            if not self.handle:
                raise RuntimeError("Filter is closed")
            # Resolve lazily so the legacy API works with older deployed libraries.
            function = self.lib.pansy_mf_peaks
            function.argtypes = ([ct.c_void_p, ct.c_void_p, ct.c_size_t, ct.c_size_t,
                                  ct.c_void_p, ct.c_size_t, ct.c_float] +
                                 [ct.c_size_t]*3 + [ct.c_void_p]*6)
            function.restype = ct.c_int
            status = function(self.handle, rx.ctypes.data, *rx.shape, tx.ctypes.data,
                              len(tx), threshold, max_peaks, range_guard, doppler_guard,
                              gates.ctypes.data, bins.ctypes.data, excess.ctypes.data,
                              amplitudes.ctypes.data, ct.byref(count), ct.byref(median))
        if status:
            raise RuntimeError(f"Rust peak extraction rejected pulse: status {status}")
        n = count.value
        return dict(range_sample=gates[:n].copy(), doppler_bin=bins[:n].copy(),
                    excess_over_median=excess[:n].copy(),
                    receiver_amplitudes=amplitudes[:n].copy(), median_power=median.value)

    def __del__(self):
        if getattr(self, "handle", None):
            self.close()
