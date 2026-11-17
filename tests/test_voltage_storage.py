"""Storage contracts shared by native Tibet cuts and the MATLAB adapter."""

import h5py
import numpy as np
from scipy.io import savemat
from meteor_analysis.voltage import read_iq, load_mat_cut


def fixture():
    original = np.array(
        [-2147483648, -32769, -32768, -1, 0, 32767, 32768, 2147483647], np.int32
    ).reshape(1, 1, 4, 2)
    compact = np.clip(original, -32768, 32767).astype(np.int16)
    index = np.flatnonzero(original.ravel() != compact.ravel()).astype(np.uint64)
    words = np.stack(
        [
            (original.ravel()[index] & 65535).astype(np.uint16).view(np.int16),
            (original.ravel()[index] >> 16).astype(np.int16),
        ],
        axis=-1,
    )
    return original, dict(
        echo_iq_i16=compact, echo_overflow_index=index, echo_overflow_words_i16=words
    )


def test_hdf_and_mat_iq_are_identical(tmp_path):
    original, fields = fixture()
    with h5py.File(tmp_path / "cut.h5", "w") as h:
        for k, v in fields.items():
            h.create_dataset(k, data=v, compression="gzip", shuffle=True)
        np.testing.assert_array_equal(read_iq(h, "echo"), original)
    mat_fields = {
        k: (v.transpose(tuple(reversed(range(v.ndim)))) if v.ndim > 1 else v)
        for k, v in fields.items()
    }
    savemat(
        tmp_path / "cut.mat", mat_fields, do_compression=True, long_field_names=True
    )
    np.testing.assert_array_equal(
        read_iq(load_mat_cut(tmp_path / "cut.mat"), "echo"), original
    )


def test_empty_overflow_and_legacy():
    small = np.arange(8, dtype=np.int16).reshape(1, 1, 4, 2)
    np.testing.assert_array_equal(
        read_iq(dict(echo_iq_i16=small, echo_overflow_index_empty=np.uint8(1)), "echo"),
        small,
    )
    np.testing.assert_array_equal(
        read_iq(dict(echo_iq_i32=small.astype(np.int32)), "echo"), small
    )
