"""Native RANSAC with short-lag receiver-interferometry compatibility."""
import ctypes as ct
import os
from pathlib import Path
import sys
import time

import digital_rf
import h5py
import numpy as np
from meteor_rust_filter import RustMatchedFilter

ROOT = Path(__file__).resolve().parent


class PhaseGrouper:
    def __init__(self, raw_root="/media/archive"):
        suffix = "dylib" if sys.platform == "darwin" else "so"
        library = ROOT/"meteor/consensus/target/release"/f"libpansy_meteor_consensus.{suffix}"
        self.lib = ct.CDLL(os.environ.get("PANSY_CONSENSUS_LIBRARY", str(library)))
        self.lib.pansy_consensus_phase.argtypes = [ct.c_void_p]*5 + [ct.c_size_t, ct.c_double, ct.c_double, ct.c_void_p]
        self.lib.pansy_consensus_phase.restype = ct.c_int
        self.raw_root = raw_root
        self.filters = {}
        with h5py.File(ROOT/"data/mesocal.h5") as cal:
            pwr = cal["pwr"][:]
        self.scale = np.real(np.sqrt(pwr[0])/np.sqrt(pwr[:7])).astype(np.float32)
        self.correction = np.empty((5,7), dtype=np.complex64)
        for beam in range(5):
            with h5py.File(ROOT/f"data/phases{beam}.h5") as cal:
                # Match the phase sign in pansy_interferometry.image_points1.
                self.correction[beam] = np.exp(1j*cal["phasecal"][:])
        if not np.all(np.isfinite(self.scale)) or not np.all(np.isfinite(self.correction)):
            raise ValueError("invalid receiver calibration")

    def signatures(self, samples, ranges, dopplers, beams):
        samples = np.asarray(samples, dtype=np.int64)
        beams = np.asarray(beams, dtype=np.int32)
        if np.any((beams < 0) | (beams >= 20)):
            raise ValueError("MF beam_pos_idx must identify position within 20-pulse batch")
        starts = samples-beams.astype(np.int64)*1600
        reader = digital_rf.DigitalRFReader(self.raw_root)
        result = np.empty((len(samples),7), dtype=np.complex64)
        gates = np.rint(np.asarray(ranges)/(299792458/2e9)).astype(int)
        frequencies = np.fft.fftshift(np.fft.fftfreq(256, d=8e-6))
        velocities = frequencies*299792458/(2*47e6)
        for start in np.unique(starts):
            ids = np.flatnonzero(starts == start)
            rx = np.array([reader.read_vector_1d(int(start),32000,f"ch{ch:03}")*self.scale[ch]
                           for ch in range(7)],dtype=np.complex64)
            tx = np.asarray(reader.read_vector_1d(int(start),32000,"ch007"),dtype=np.complex64)
            for i in ids:
                offset = int(beams[i])*1600
                gate = int(gates[i])
                if gate not in self.filters:
                    self.filters[gate] = RustMatchedFilter([gate],132,8,256)
                code = tx[offset:offset+132]
                if not np.all(np.isfinite(code)) or np.sum(abs(code)**2) <= 1e9:
                    raise ValueError("raw TX pulse failed quality gate")
                power, _, amplitudes = self.filters[gate].process(
                    rx[:,offset:offset+1600],code,return_amplitudes=True)
                measured = velocities[np.argmax(power[0])]
                if abs(measured-dopplers[i]) > 1.:
                    raise ValueError("raw Doppler peak disagrees with MF metadata")
                result[i] = amplitudes*self.correction[int(beams[i])%5]
        return result

    def group(self, samples, ranges, dopplers, snr, beams):
        if len(samples) > 16384:
            raise ValueError("phase-grouping window too large; never silently truncate")
        if not len(samples):
            return np.empty(0,dtype=np.int32), np.empty((0,7),dtype=np.complex64)
        start = time.monotonic()
        arrays = [np.ascontiguousarray(samples,dtype=np.uint64),
                  np.ascontiguousarray(ranges,dtype=np.float32),
                  np.ascontiguousarray(dopplers,dtype=np.float32),
                  np.ascontiguousarray(snr,dtype=np.float32)]
        if any(a.shape != (len(samples),) for a in arrays) or not all(np.all(np.isfinite(a)) for a in arrays):
            raise ValueError("invalid phase-grouping observations")
        signatures = self.signatures(*arrays[:3], beams)
        labels = np.full(len(samples),-1,dtype=np.int32)
        count = self.lib.pansy_consensus_phase(
            *[a.ctypes.data for a in arrays],signatures.ctypes.data,len(samples),
            .65,.05,labels.ctypes.data)
        if count < 0:
            raise RuntimeError(f"native interferometric consensus failed ({count})")
        print(f"Rust interferometric RANSAC: points={len(samples)} tracks={count} "
              f"assigned={np.count_nonzero(labels>=0)} elapsed={time.monotonic()-start:.3f}s",flush=True)
        return labels, signatures

    def close(self):
        for filt in self.filters.values():
            filt.close()
        self.filters.clear()
