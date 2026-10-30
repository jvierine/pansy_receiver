import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from monitor.publish import atomic_json, finite, metadata_rows, packet


class ReceiverMonitorTests(unittest.TestCase):
    def test_saved_event_field_is_bounded_and_upstream_trace_matches_exact_times(self):
        from monitor.event import save_event
        now=1793386800.;key=int((now-10)*1e6)
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);metadata=root/'metadata'
            fit=metadata/'simple_meteor_fit'/'2026-10-30T19-00-00';fit.mkdir(parents=True)
            times=np.array([key,key+1000],dtype=np.int64)
            with h5py.File(fit/'fit.h5','w') as h:
                g=h.create_group(str(key))
                for name,value in dict(v0=[1.,2.,3.],r0=[0.,0.,90.],eclon=20.,eclat=-30.,
                    txidx=times,ew=[0.,1.],ns=[0.,1.],up=[90.,89.],snr=[10.,20.],std=[.1,.2,.3],mfs=[21.,18.]).items():g[name]=value
            cut=metadata/'cut'/'2026-10-30T19-00-00';cut.mkdir(parents=True)
            with h5py.File(cut/'cut.h5','w') as h:
                g=h.create_group(str(key))
                for name,value in dict(c_tx_idx=times,c_range_km=[90.,89.],c_snr=[10.,20.],c_doppler=[-1000.,-2000.]).items():g[name]=value
            event=root/'.local/state/pansy-monitor/latest-event.h5'
            save_event(key,np.arange(600)/1000+key/1e6,np.arange(513)/10+50,np.ones((513,600)),[-1000.,-2000.],event)
            before=event.read_bytes()
            with patch('monitor.publish.Path.home',return_value=root):
                def missing(_):raise OSError('No raw fixture')
                data=packet(metadata,root,now,reader_factory=missing)
            m=data['meteor'];field=m['range_time']
            self.assertLessEqual(len(field['seconds']),256);self.assertLessEqual(len(field['range_km']),128)
            self.assertEqual(m['doppler_km_s'],[-1.,-2.]);self.assertEqual(m['detection_trace']['doppler_km_s'],[-1.,-2.])
            self.assertEqual(before,event.read_bytes(),'Publishing must not change the scientific event arrays')

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
        self.assertIn("from './lib/interactive_series.js",source)
        self.assertNotIn('/rvd/api/frame',source)
        self.assertIn("'/rvd/data/manifest.json'",source)
        self.assertIn("'/rvd/data/sky-noise/status.json'",source)
        self.assertNotIn('48 hours',html)


if __name__=='__main__':unittest.main()
