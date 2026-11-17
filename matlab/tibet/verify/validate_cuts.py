"""Prepare real cut fixtures and normalize MATLAB cuts for native Rust replay."""

import argparse
from pathlib import Path
import h5py
import numpy as np
from scipy.io import loadmat, savemat
from scipy.io.matlab import mat_struct
from meteor_analysis.voltage import read_iq


def prepare(source, output):
    output.mkdir(parents=True, exist_ok=False)

    def contents(group):
        return {
            k: (
                contents(v)
                if isinstance(v, h5py.Group)
                else (
                    v[()].transpose(tuple(reversed(range(v.ndim))))
                    if v.ndim > 1
                    else v[()]
                )
            )
            for k, v in group.items()
        }

    for path in sorted((source / "cuts").glob("*.h5")):
        with h5py.File(path) as h:
            savemat(
                output / (path.stem + ".mat"),
                contents(h),
                do_compression=True,
                long_field_names=True,
            )


def unpack(value):
    if isinstance(value, mat_struct):
        return {key: unpack(getattr(value, key)) for key in value._fieldnames}
    if isinstance(value, np.ndarray) and value.dtype == object and value.size == 1:
        return unpack(value.flat[0])
    return value


def normalize(source, output, reference=None):
    output.mkdir(parents=True, exist_ok=False)
    components = 0
    for path in sorted(source.glob("event_*.mat")):
        record = {
            k: unpack(v)
            for k, v in loadmat(path, struct_as_record=False).items()
            if not k.startswith("__")
        }

        def write(group, fields):
            for key, value in fields.items():
                if isinstance(value, dict):
                    write(group.create_group(key), value)
                    continue
                value = np.asarray(value)
                if value.dtype.kind in "US":
                    value = np.frombuffer("".join(value.ravel()).encode(), np.uint8)
                elif key.endswith("_overflow_words_i16"):
                    value = value.reshape(2, -1).T
                elif value.ndim > 1:
                    value = value.transpose(tuple(reversed(range(value.ndim))))
                    if value.ndim <= 2 and 1 in value.shape:
                        value = value.reshape(-1)
                if value.size == 1 and value.ndim <= 2:
                    value = value.reshape(())
                kw = (
                    dict(compression="gzip", compression_opts=4, shuffle=True)
                    if value.ndim and value.size
                    else {}
                )
                group.create_dataset(key, data=value, **kw)

        with h5py.File(output / (path.stem + ".h5"), "w") as h:
            write(h, record)
            if reference:
                with h5py.File(reference / "cuts" / (path.stem + ".h5")) as old:
                    for name in ["zrx_echoes", "ztx_history"]:
                        a = read_iq(h, name)
                        b = read_iq(old, name)
                        np.testing.assert_array_equal(a, b)
                        components += a.size
    print(
        "Normalized",
        len(list(source.glob("event_*.mat"))),
        "MAT files;",
        components,
        "exact real IQ components checked",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["prepare", "normalize"])
    p.add_argument("source", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--reference", type=Path)
    a = p.parse_args()
    if a.mode == "prepare":
        prepare(a.source, a.output)
    else:
        normalize(a.source, a.output, a.reference)
