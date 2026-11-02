"""Write the downloadable HDF5 for one changed JSON packet; never access the archive."""
import json
import hashlib
import os
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np


def packet_digest(item):
    return hashlib.sha256(json.dumps(
        item, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode()).hexdigest()


def file_signature(path):
    stat = path.stat()
    return [stat.st_ino, stat.st_size, stat.st_mtime_ns]


def write_sidecar(packet_path, output):
    packet = json.loads(Path(packet_path).read_text())
    output = Path(output)
    temporary = output.with_suffix(f".h5.{os.getpid()}.partial")
    index_path = output.with_suffix(".h5.index.json")
    fingerprints = {}
    try:
        cached = json.loads(index_path.read_text())
        if cached["signature"] == file_signature(output):
            fingerprints = cached["packets"]
    except (OSError, ValueError, KeyError):
        pass
    try:
        if output.exists():
            shutil.copyfile(output, temporary)
        with h5py.File(temporary, "r+" if temporary.exists() else "w") as h5:
            if not fingerprints and packet.get("selected_day"):
                # The published JSON is an exact identity index for its HDF5
                # revision. This avoids reopening hundreds of unchanged groups
                # when adopting a sidecar produced by the previous writer.
                previous = output.parent / "days" / f"{packet['selected_day']}.json"
                try:
                    old_packet = json.loads(previous.read_text())
                    if old_packet["revision"] == h5.attrs.get("revision"):
                        fingerprints = {p["key"]: packet_digest(p)
                                        for p in old_packet["frame_packets"]}
                except (OSError, ValueError, KeyError):
                    pass
            h5.attrs.update(
                schema="pansy_rvd_plot_v2", revision=packet["revision"],
                generator="pmse/write_gui_sidecar.py", snr_definition="excess signal/noise",
                moment_censoring_db=3,
            )
            for name, value in packet["geometry"].items():
                if name in h5 and np.array_equal(h5[name][()], value):
                    continue
                if name in h5:
                    del h5[name]
                h5[name] = value
            frames = h5.require_group("frames")
            keep = {item["key"] for item in packet["frame_packets"]}
            for key in set(frames) - keep:
                del frames[key]
            for item in packet["frame_packets"]:
                digest = packet_digest(item)
                if fingerprints.get(item["key"]) == digest:
                    continue
                fingerprints[item["key"]] = digest
                if item["key"] in frames:
                    old = frames[item["key"]]
                    if old.attrs.get("packet_sha256") == digest:
                        continue
                    # Legacy sidecars use the native source_revision as their
                    # measurement identity. Seed their fingerprint once without
                    # recreating unchanged datasets and groups.
                    expected = {k for k, v in item.items() if k != "key" and v is not None}
                    previous = set(old) | (set(old.attrs) - {"packet_sha256"})
                    if ("packet_sha256" not in old.attrs and expected == previous
                            and old.attrs.get("source_revision") == item.get("source_revision")
                            and bool(old["spectra_ready"][()] if "spectra_ready" in old else True)
                            == bool(item.get("spectra_ready", True))):
                        old.attrs["packet_sha256"] = digest
                        continue
                    del frames[item["key"]]
                group = frames.create_group(item["key"])
                group.attrs["packet_sha256"] = digest
                for name, value in item.items():
                    if name == "key" or value is None:
                        continue
                    if isinstance(value, str):
                        group.attrs[name] = value
                    else:
                        group[name] = value
        temporary.replace(output)
        index_temporary = index_path.with_suffix(f".json.{os.getpid()}.partial")
        index_temporary.write_text(json.dumps({
            "signature": file_signature(output),
            "packets": {k: v for k, v in fingerprints.items() if k in keep},
        }, separators=(",", ":")))
        index_temporary.replace(index_path)
    finally:
        temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    write_sidecar(*sys.argv[1:])
