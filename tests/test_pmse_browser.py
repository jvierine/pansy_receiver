import json
from pathlib import Path
import struct
import tempfile
import unittest
import h5py
import numpy as np
from pmse.serve_browser import existing_packet


class BrowserTransportTests(unittest.TestCase):
    def test_selected_range_all_channels_and_bounded_looks(self):
        key=1792433051921664
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);directory=root/'2026-10-19T18-00-00';directory.mkdir()
            values=(np.arange(150*7*3).reshape(150,7,3)*(1+2j)).astype('complex64')
            with h5py.File(directory/f'pmse@{key//1_000_000//600*600}.h5','w') as h:
                g=h.create_group(str(key))
                fields=dict(complete=1,sparse_index=[[0,0,0],[1,1,1],[4,2,1]],slant_range_km=[84.,85.],
                            spectra=values,look_bounds_sample=np.column_stack([key+np.arange(150)*2_000_000,key+(np.arange(150)+1)*2_000_000]),
                            mode_start_sample=key,mode_stop_sample_exclusive=key+300_000_000,mode_duration_seconds=300.,
                            radial_velocity_mps=[-99.,0.,99.],noise_power=np.ones((7,5,2)),schema_version='test',
                            config_json=json.dumps(dict(ipp_samples=1600,sample_rate_hz=1e6,radar_frequency_hz=47e6,complementary_pairs=True)))
                for name,value in fields.items():g[name]=value
            raw=existing_packet(root,key,85.,120.)
            length=struct.unpack('<I',raw[:4])[0]
            header=json.loads(raw[4:4+length]);offset=4+(length+3)//4*4
            data=np.frombuffer(raw,dtype='<c8',offset=offset).reshape(header['looks'],7,-1)
            self.assertEqual(header['looks'],60)
            self.assertEqual(header['stop_us']-header['start_us'],120_000_000)
            self.assertEqual(header['index'],[[1,1,1],[4,2,1]])
            np.testing.assert_array_equal(data,values[60:120,:,1:])
            self.assertAlmostEqual(header['nyquist_mps'],99.66504587765957)


if __name__=='__main__':unittest.main()
