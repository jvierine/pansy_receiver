"""Write the downloadable HDF5 for one changed JSON packet; never access the archive."""
import json
import os
import sys
from pathlib import Path

import h5py


def write_sidecar(packet_path, output):
    packet = json.loads(Path(packet_path).read_text())
    output = Path(output)
    temporary = output.with_suffix(f".h5.{os.getpid()}.partial")
    try:
        with h5py.File(temporary, "w") as h5:
            h5.attrs.update(
                schema="pansy_rvd_plot_v2", revision=packet["revision"],
                generator="pmse/write_gui_sidecar.py", snr_definition="excess signal/noise",
                moment_censoring_db=3,
            )
            for name, value in packet["geometry"].items():
                h5[name] = value
            frames = h5.create_group("frames")
            for item in packet["frame_packets"]:
                group = frames.create_group(item["key"])
                for name, value in item.items():
                    if name == "key" or value is None:
                        continue
                    if isinstance(value, str):
                        group.attrs[name] = value
                    else:
                        group[name] = value
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    write_sidecar(*sys.argv[1:])
