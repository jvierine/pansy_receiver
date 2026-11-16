"""Numpy interface to the reusable Rust kernels; no acquisition/service side effects."""

import ctypes
from pathlib import Path
import sys
import numpy as np

LIB = (
    Path(__file__).resolve().parents[1]
    / "meteor/analysis/target/release"
    / (
        "libmeteor_analysis.dylib"
        if sys.platform == "darwin"
        else "libmeteor_analysis.so"
    )
)
_library = None


def library():
    global _library
    if _library is None:
        _library = ctypes.CDLL(str(LIB))
        p = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
        _library.meteor_matched.argtypes = [
            p,
            ctypes.c_size_t,
            ctypes.c_size_t,
            p,
            ctypes.c_size_t,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            p,
        ]
        _library.meteor_matched.restype = ctypes.c_int
        _library.meteor_propagate.argtypes = [
            p,
            p,
            ctypes.c_size_t,
            p,
            p,
            ctypes.c_size_t,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
            p,
        ]
        _library.meteor_propagate.restype = ctypes.c_int
        i = np.ctypeslib.ndpointer(dtype=np.int64, flags="C_CONTIGUOUS")
        b = np.ctypeslib.ndpointer(dtype=np.uint8, flags="C_CONTIGUOUS")
        _library.meteor_array_response.argtypes = [
            p,
            i,
            b,
            ctypes.c_size_t,
            ctypes.c_size_t,
            p,
            ctypes.c_size_t,
            ctypes.c_double,
            p,
        ]
        _library.meteor_array_response.restype = ctypes.c_int
    return _library


def matched(rx, tx, start, fd, fs):
    rx = np.ascontiguousarray(rx, np.complex128)
    tx = np.ascontiguousarray(tx, np.complex128)
    out = np.empty(rx.shape[0], np.complex128)
    status = library().meteor_matched(
        rx.view(np.float64),
        *rx.shape,
        tx.view(np.float64),
        len(tx),
        float(start),
        float(fd),
        float(fs),
        out.view(np.float64),
    )
    if status:
        raise ValueError(f"native matched-filter status {status}")
    return out


def array_response(xy, membership, active, uv, wavelength):
    xy = np.ascontiguousarray(xy, np.float64)
    membership = np.ascontiguousarray(membership, np.int64)
    active = np.ascontiguousarray(active, np.uint8)
    uv = np.ascontiguousarray(uv, np.float64)
    out = np.empty((len(uv), membership.shape[1]), np.complex128)
    status = library().meteor_array_response(
        xy,
        membership,
        active,
        len(xy),
        membership.shape[1],
        uv,
        len(uv),
        wavelength,
        out.view(np.float64),
    )
    if status:
        raise ValueError(f"native array response status {status}")
    return out


def propagate(
    state,
    radius,
    times,
    altitudes,
    densities,
    station_alt,
    material=3000.0,
    sigma=1e-8,
    dt=0.0005,
    gravity=True,
):
    initial = np.ascontiguousarray(np.r_[state, radius], np.float64)
    times = np.ascontiguousarray(times, np.float64)
    altitudes = np.ascontiguousarray(altitudes, np.float64)
    densities = np.ascontiguousarray(densities, np.float64)
    out = np.empty((len(times), 7), np.float64)
    status = library().meteor_propagate(
        initial,
        times,
        len(times),
        altitudes,
        densities,
        len(altitudes),
        station_alt,
        material,
        sigma,
        dt,
        int(gravity),
        out,
    )
    if status:
        raise ValueError(f"native dynamics status {status}")
    return out
