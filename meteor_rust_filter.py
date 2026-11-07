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

    def __del__(self):
        if getattr(self, "handle", None):
            self.close()
