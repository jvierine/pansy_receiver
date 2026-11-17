"""Lossless compressed int16 voltage arrays, with legacy int32 compatibility."""

import numpy as np
from scipy.io import loadmat
from scipy.io.matlab import mat_struct


def load_mat_cut(path):
    """Convert MATLAB column-major cut axes into the native HDF5 axis order."""

    def convert(value, key=""):
        if isinstance(value, mat_struct):
            return {k: convert(getattr(value, k), k) for k in value._fieldnames}
        if isinstance(value, np.ndarray) and value.dtype == object and value.size == 1:
            return convert(value.flat[0], key)
        value = np.asarray(value)
        if value.dtype.kind in "US":
            return np.frombuffer("".join(value.ravel()).encode(), np.uint8)
        if key.endswith("_overflow_words_i16"):
            return value.reshape(2, -1).T
        if value.ndim > 1:
            value = value.transpose(tuple(reversed(range(value.ndim))))
            if value.ndim <= 2 and 1 in value.shape:
                value = value.reshape(-1)
        return value.reshape(()) if value.size == 1 and value.ndim <= 2 else value

    return {
        k: convert(v, k)
        for k, v in loadmat(path, struct_as_record=False).items()
        if not k.startswith("__")
    }


def read_iq(group, name):
    key = name + "_iq_i16"
    if key not in group:
        return group[name + "_iq_i32"][()]
    values = group[key][()].astype(np.int32)
    if name + "_overflow_index" not in group:
        assert name + "_overflow_index_empty" in group
        return values
    indices = group[name + "_overflow_index"][()].astype(np.intp)
    words = group[name + "_overflow_words_i16"][()].astype(np.int32)
    assert words.shape == (len(indices), 2)
    assert np.all(indices < values.size) and np.all(indices >= 0)
    values.ravel()[indices] = (words[:, 1] << 16) | (words[:, 0] & 65535)
    return values
