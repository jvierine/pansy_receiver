"""Read-only check of the production adapter, independent of MPI and metadata."""
import ast
import os
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import numpy as np
import scipy.constants
import pyfftw

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root))
tree = ast.parse((root / "quick_search_meteor.py").read_text())
cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "range_doppler_search")
scope = dict(n=np, c=scipy.constants, pc=SimpleNamespace(freq=47e6),
             fft=pyfftw.interfaces.scipy_fftpack.fft, os=os, log0=print)
exec(compile(ast.Module(body=[cls], type_ignores=[]), "quick_search_meteor.py", "exec"), scope)
gates = np.arange(int(np.floor(60/(scipy.constants.c/2e9))), int(np.floor(170/(scipy.constants.c/2e9)))+1, 2)
os.environ["PANSY_MF_BACKEND"] = "python"
python = scope["range_doppler_search"](rg=gates)
os.environ["PANSY_MF_BACKEND"] = "rust"
rust = scope["range_doppler_search"](rg=gates)
rng = np.random.default_rng(4109)
times = np.zeros(2)
for pulse in range(20):
    tx = (rng.normal(size=132)+1j*rng.normal(size=132)).astype(np.complex64)
    rx = (rng.normal(size=(7,1600))+1j*rng.normal(size=(7,1600))).astype(np.complex64)
    rx[:,600:732] += 5*tx*np.exp(2j*np.pi*0.023*np.arange(132))
    outputs = []
    for i, filt in enumerate([python, rust]):
        start = time.perf_counter()
        outputs.append(filt.mf(rx, tx))
        times[i] += time.perf_counter()-start
    a, b = outputs
    assert np.max(np.abs(a[0]-b[0])) < np.max(a[0])*3e-6
    np.testing.assert_allclose(a[1], b[1], rtol=3e-6)
    np.testing.assert_array_equal(a[2], b[2])
    assert abs(a[3]/b[3]-1) < 3e-6
rust.native.close()
print(f"PASS: 20 seven-receiver pulses, full 60-170 km search; Python {times[0]/20*1000:.3f} ms/pulse; Rust {times[1]/20*1000:.3f} ms/pulse")
