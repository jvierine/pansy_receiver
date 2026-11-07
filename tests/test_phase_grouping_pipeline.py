import os
import unittest
from unittest.mock import patch
import numpy as np
import cluster_mf
from types import SimpleNamespace
from meteor_progress import read_progress, write_progress


def test_grouping_advances_through_completed_quiet_interval(tmp_path,monkeypatch):
    monkeypatch.setenv("HOME",str(tmp_path))
    monkeypatch.setenv("PANSY_GROUPING_BACKEND","rust_phase")
    state=tmp_path/".local/state/pansy_receiver"
    state.mkdir(parents=True)
    checkpoint=state/"cluster_phase_v1.h5"
    write_progress(checkpoint,90_000_000)
    det=tmp_path/"detections"
    det.mkdir()
    with patch.object(cluster_mf.pc,"detections_metadata_dir",str(det)), \
         patch.object(cluster_mf,"metadata_bounds",return_value=(object(),(10_000_000,20_000_000))), \
         patch.object(cluster_mf,"read_completed_frontier",return_value=100_000_000), \
         patch.object(cluster_mf.drf,"DigitalMetadataWriter",return_value=object()), \
         patch.object(cluster_mf.drf,"DigitalRFReader",return_value=SimpleNamespace(get_bounds=lambda ch:(0,200_000_000))), \
         patch.object(cluster_mf,"read_mf_output",return_value=tuple(np.array([]) for _ in range(6))), \
         patch.object(cluster_mf,"find_clusters"):
        cluster_mf.analyze_until_now()
    assert read_progress(checkpoint)==100_000_000


class PhasePipelineTests(unittest.TestCase):
    def run_detector(self, writer):
        t=np.arange(10,dtype=np.int64)*8000+1_746_000_000_000_000
        with patch.dict(os.environ,{"PANSY_GROUPING_BACKEND":"rust_phase"}), \
             patch.object(cluster_mf,"phase_clusters",return_value=[{"idx":np.arange(10),"xhat":np.array([100.,-20.,0.])}]):
            cluster_mf.find_clusters(np.full(10,1e10),t,100-20*np.arange(10)*.008,
                                     np.full(10,-20000),np.full(10,20.),np.arange(10),writer)

    def test_metadata_contract(self):
        class Writer:
            def write(self, keys, data):
                self.keys,self.data=keys,data
        writer=Writer()
        self.run_detector(writer)
        self.assertEqual(set(writer.data),{"xhat","tx_idx","range","doppler","snr","beam"})
        self.assertEqual(writer.keys,[1_746_000_000_000_000])
        self.assertEqual(len(writer.data["tx_idx"]),10)

    def test_write_failure_is_not_swallowed(self):
        class Writer:
            def write(self,*args): raise OSError("disk write failed")
        with self.assertRaisesRegex(OSError,"disk write failed"):
            self.run_detector(Writer())

    def test_replayed_existing_record_is_idempotent(self):
        class Writer:
            def write(self,*args): raise ValueError("sample already in data")
        with patch.object(cluster_mf,"write_complete_record",return_value=False):
            self.run_detector(Writer())


if __name__ == "__main__":
    unittest.main()
