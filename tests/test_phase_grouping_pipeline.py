import os
import unittest
from unittest.mock import patch
import numpy as np
import cluster_mf


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
        self.run_detector(Writer())


if __name__ == "__main__":
    unittest.main()
