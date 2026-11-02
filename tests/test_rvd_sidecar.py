import json

import h5py
import numpy as np

from pmse.write_gui_sidecar import write_sidecar


def test_updates_only_changed_frames(tmp_path, monkeypatch):
    source = tmp_path / "day.json"
    target = tmp_path / "day.h5"
    packet = {"revision": "first", "geometry": {"range_km": [80.0, 81.0]},
              "frame_packets": [{"key": "1", "source_revision": "a", "snr_db": [5.0]},
                                {"key": "2", "source_revision": "b", "snr_db": [6.0]}]}
    source.write_text(json.dumps(packet))
    write_sidecar(source, target)
    with h5py.File(target) as h5:
        address = h5py.h5o.get_info(h5["frames/1/snr_db"].id).addr
    packet["revision"] = "second"
    packet["frame_packets"][1]["snr_db"] = [8.0]
    source.write_text(json.dumps(packet))
    original_getitem = h5py.Group.__getitem__

    def checked_getitem(group, name):
        if group.name == "/frames" and name == "1":
            raise AssertionError("unchanged frame metadata must not be reopened")
        return original_getitem(group, name)

    with monkeypatch.context() as patch:
        patch.setattr(h5py.Group, "__getitem__", checked_getitem)
        write_sidecar(source, target)
    with h5py.File(target) as h5:
        assert h5.attrs["revision"] == "second"
        assert h5py.h5o.get_info(h5["frames/1/snr_db"].id).addr == address
        np.testing.assert_equal(h5["frames/2/snr_db"][()], [8.0])


def test_removes_deleted_frames_and_keeps_previous_file_on_error(tmp_path):
    source = tmp_path / "day.json"
    target = tmp_path / "day.h5"
    packet = {"revision": "first", "geometry": {}, "frame_packets": [
        {"key": "1", "snr_db": [5.0]}, {"key": "2", "snr_db": [6.0]}]}
    source.write_text(json.dumps(packet))
    write_sidecar(source, target)
    packet["frame_packets"].pop()
    packet["revision"] = "second"
    source.write_text(json.dumps(packet))
    write_sidecar(source, target)
    with h5py.File(target) as h5:
        assert list(h5["frames"]) == ["1"]
