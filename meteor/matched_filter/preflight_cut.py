"""Compare the actual cut-filter class, without importing its MPI/startup code."""
import argparse
import ast
import itertools
import os
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import numpy as np
import scipy.constants
import scipy.fftpack

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def compare(tx, rx):
    tree = ast.parse((ROOT/"process_cut_meteor.py").read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "range_doppler_search")
    scope = dict(n=np, c=scipy.constants, pc=SimpleNamespace(freq=47e6),
                 fft=scipy.fftpack.fft, os=os, itertools=itertools)
    exec(compile(ast.Module(body=[cls], type_ignores=[]), "process_cut_meteor.py", "exec"), scope)
    elapsed = np.zeros(2)
    count = 0
    for interp in (1, 2):
        filters = []
        for backend in ("python", "rust"):
            os.environ["PANSY_CUT_MF_BACKEND"] = backend
            filters.append(scope["range_doppler_search"](tx.shape[1], rx.shape[2], rx.shape[1], interp))
        try:
            for i in range(len(tx)):
                outputs = []
                for j, filt in enumerate(filters):
                    start = time.perf_counter()
                    outputs.append(filt.mf(tx[i], rx[i]))
                    elapsed[j] += time.perf_counter()-start
                a,b = outputs
                assert np.max(abs(a[0]-b[0])) <= 4e-6*np.max(a[0])
                np.testing.assert_allclose(a[1], b[1], rtol=4e-6)
                np.testing.assert_array_equal(a[2], b[2])
                np.testing.assert_allclose(a[3], b[3], rtol=4e-6)
                assert a[5] == b[5]
                assert np.max(abs(a[4]-b[4])) <= 6e-6*np.max(abs(a[4]))
                # Relative phase directly, without phase unwrapping.
                strong = abs(a[4]) > .01*np.max(abs(a[4]))
                assert np.max(abs(np.angle(a[4][strong]*b[4][strong].conj()))) < 1e-4
                count += 1
        finally:
            filters[1].native.close()
    print(f"PASS: {count} cut pulse/interpolation cases; power, profiles, Doppler, median, "
          f"range peak, all 21 complex cross-products. Python {elapsed[0]/count*1000:.3f} "
          f"ms/pulse, Rust {elapsed[1]/count*1000:.3f} ms/pulse", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cut-root")
    args = parser.parse_args()
    rng = np.random.default_rng(733)
    tx = (rng.normal(size=(8,132))+1j*rng.normal(size=(8,132))).astype(np.complex64)
    rx = (rng.normal(size=(8,7,190))+1j*rng.normal(size=(8,7,190))).astype(np.complex64)
    for i in range(8):
        rx[i,:,24:156] += 4*tx[i].conj()*np.exp(2j*np.pi*.015*np.arange(132))
    compare(tx, rx)
    if args.cut_root:
        import digital_rf
        import h5py
        reader = digital_rf.DigitalMetadataReader(args.cut_root)
        upper = reader.get_bounds()[1]
        candidates = reader.read(upper-60_000_000, upper, ["c_snr"])
        key = max(candidates, key=lambda k: np.max(candidates[k]["c_snr"]))
        record = reader.read(key, key)[key]
        tx = np.asarray(record["ztx_pulses_re"], dtype=np.complex64) + 1j*np.asarray(record["ztx_pulses_im"], dtype=np.float32)
        rx = np.asarray(record["zrx_echoes_re"], dtype=np.complex64) + 1j*np.asarray(record["zrx_echoes_im"], dtype=np.float32)
        with h5py.File(ROOT/"data/mesocal.h5") as calibration:
            pwr = calibration["pwr"][:]
        rx *= np.real(np.sqrt(pwr[0])/np.sqrt(pwr[:7]))[None,:,None]
        indices = np.unique(np.linspace(0, len(tx)-1, min(16,len(tx)), dtype=int))
        print(f"Recorded cut sample {key}, testing {len(indices)} pulses", flush=True)
        compare(tx[indices], rx[indices])


if __name__ == "__main__":
    main()
