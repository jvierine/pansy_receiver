import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from monitor.publish import atomic_json, finite, metadata_rows, packet


class ReceiverMonitorTests(unittest.TestCase):
    def test_packet_is_read_only_24h_and_preserves_complex_raw_components(self):
        now=1793386800.
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            directory=root/'simple_meteor_fit'/'2026-10-30T12-00-00'
            directory.mkdir(parents=True)
            path=directory/'fit@1793361600.h5'
            with h5py.File(path,'w') as h:
                for sample in [int((now-90000)*1e6),int((now-100)*1e6)]:
                    g=h.create_group(str(sample))
                    for name,value in dict(v0=[3.,4.,0.],r0=[0.,0.,90.],eclon=20.,eclat=-30.,
                        txidx=np.array([sample,sample+1000]),ew=[1.,2.],ns=[3.,4.],up=[90.,89.],snr=[10.,100.]).items():g[name]=value
            before=path.read_bytes()

            class RawReader:
                def __init__(self, path):pass
                def get_bounds(self,ch):return int((now-86400)*1e6),int((now-1)*1e6)
                def read_vector_1d(self,start,count,ch):return np.arange(count)+2j*np.arange(count)

            data=packet(root,root/'raw',now=now,reader_factory=RawReader)
            self.assertEqual(len(data['fits']),1)
            self.assertEqual(data['fits'][0][1:],[5.,90.,20.,-30.])
            self.assertEqual(len(data['raw']['channels']),8)
            self.assertEqual(data['raw']['channels'][0]['real'][:3],[0.,8.,16.])
            self.assertEqual(data['raw']['channels'][0]['imag'][:3],[0.,16.,32.])
            self.assertEqual(data['meteor']['snr_db'],[10.,20.])
            self.assertEqual(path.read_bytes(),before)
            json.dumps(data,allow_nan=False)

    def test_missing_raw_does_not_prevent_telemetry_and_nan_becomes_null(self):
        with tempfile.TemporaryDirectory() as temporary:
            with patch('monitor.publish.latest_sample',return_value=None):
                def failed(path):raise OSError('raw temporarily absent')
                data=packet(temporary,temporary,now=1793386800.,reader_factory=failed)
            self.assertIsNone(data['raw'])
            self.assertIn('temporarily absent',data['errors']['raw'])
            self.assertEqual(finite([1.,np.nan,np.inf]),[1.,None,None])
            out=Path(temporary)/'out.json';atomic_json(out,data)
            self.assertEqual(json.loads(out.read_text()),data)
            self.assertFalse(out.with_suffix('.json.partial').exists())

    def test_page_uses_existing_rvd_components_not_historical_images(self):
        root=Path(__file__).resolve().parents[1]
        html=(root/'web/index.html').read_text()
        self.assertNotIn('.gif',html)
        for legacy in ['status.png','latest_hist.png','latest_radiants.png','latest_meteor.png','raw.png','processing.png']:
            self.assertNotIn(legacy,html)
        source=(root/'web/monitor.js').read_text()
        self.assertIn("from '/rvd/lib/interactive_plot.js",source)
        self.assertIn("from '/rvd/lib/interactive_series.js",source)
        self.assertIn("'/rvd/data/manifest.json'",source)
        self.assertIn("'/rvd/data/sky-noise/status.json'",source)
        self.assertNotIn('48 hours',html)


if __name__=='__main__':unittest.main()
