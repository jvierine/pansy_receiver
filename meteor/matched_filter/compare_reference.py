"""Offline Rust parity against the actual legacy class and independent complex FFTs.

Uses recorded cut voltages when --cut-root is supplied. No production writes.
"""
import argparse
import ast
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace

import h5py
import numpy as np
import scipy.constants
from scipy.fft import fft

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cut-root", type=Path)
    args = parser.parse_args()
    tree = ast.parse((ROOT / "quick_search_meteor.py").read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "range_doppler_search")
    scope = {"n": np, "c": scipy.constants, "pc": SimpleNamespace(freq=47e6), "fft": fft, "os": os, "log0": lambda _: None}
    exec(compile(ast.Module(body=[cls], type_ignores=[]), "quick_search_meteor.py", "exec"), scope)
    rng = np.random.default_rng(3107)
    tx = rng.choice([-1, 1], 132).astype(np.complex64)
    rx = (rng.normal(size=(7, 400)) + 1j * rng.normal(size=(7, 400))).astype(np.complex64)*0.2
    for gate, freq in [(20, 23/256), (180, -37/256)]:
        phase = np.exp(1j*(np.arange(7)[:, None]*0.37 + 2*np.pi*freq*np.arange(132)[None, :]))
        rx[:, gate:gate+132] += 4*phase*tx
    cases = [("two_targets", rx, tx, np.arange(0, 250, 2))]
    if args.cut_root:
        path = next(iter(sorted(args.cut_root.glob("*/*.h5"))))
        with h5py.File(path) as f:
            g = f[next(iter(f))]
            for pulse in range(min(3, g["zrx_echoes_re"].shape[0])):
                rx = (g["zrx_echoes_re"][pulse].astype(np.float32) + 1j*g["zrx_echoes_im"][pulse]).astype(np.complex64)
                tx = (g["ztx_pulses_re"][pulse].astype(np.float32) + 1j*g["ztx_pulses_im"][pulse]).astype(np.complex64)
                cases.append((f"recorded_{pulse}", rx, tx, np.arange(0, rx.shape[1]-len(tx)+1, 2)))
    with tempfile.TemporaryDirectory(prefix="pansy-mf-") as tmp:
        path = Path(tmp)/"reference.h5"
        with h5py.File(path, "w") as f:
            for name, rx, tx, gates in cases:
                for dec in (1, 8):
                    os.environ["PANSY_MF_BACKEND"] = "python"
                    ref = scope["range_doppler_search"](txlen=len(tx), rg=gates, fdec=dec, fftlen=256)
                    power, _, _, median = ref.mf(rx, tx)
                    start = time.perf_counter()
                    for _ in range(20):
                        ref.mf(rx, tx)
                    print(f"{name}_d{dec}: Python warm {(time.perf_counter()-start)*1000/20:.3f} ms/pulse", flush=True)
                    os.environ["PANSY_MF_BACKEND"] = "rust"
                    native = scope["range_doppler_search"](txlen=len(tx), rg=gates, fdec=dec, fftlen=256)
                    actual, profile, velocity, actual_median = native.mf(rx, tx)
                    assert np.max(np.abs(actual-power)) < 3e-6*np.max(power)
                    assert abs(actual_median/median-1) < 3e-6
                    np.testing.assert_allclose(profile, power.max(axis=1), rtol=3e-6)
                    np.testing.assert_array_equal(velocity, ref.dopv[np.argmax(power, axis=1)])
                    native.native.close()
                    decoded = rx[:, gates[:, None]+np.arange(len(tx))] * tx.conj()
                    if dec > 1:
                        width = len(tx)//dec
                        decoded = decoded[:, :, :width*dec].reshape(len(rx), len(gates), width, dec).sum(axis=-1)
                    spectra = np.fft.fftshift(fft(decoded, 256, axis=-1), axes=-1)
                    g = f.create_group(f"{name}_d{dec}")
                    for key, value in dict(rx_re=rx.real, rx_im=rx.imag, tx_re=tx.real, tx_im=tx.imag, gates=gates.astype(np.uint64), power=power, complex_re=spectra.real, complex_im=spectra.imag).items():
                        g[key] = value
                    g.attrs["median"] = np.float32(median)
                    g.attrs["decimation"] = np.uint64(dec)
        subprocess.run(["cargo", "test", "--release", "--manifest-path", str(Path(__file__).with_name("Cargo.toml")), "python_reference", "--", "--ignored", "--nocapture"], env={**os.environ, "MF_REFERENCE": str(path)}, check=True)


if __name__ == "__main__":
    main()
