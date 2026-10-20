import json
from pathlib import Path
import struct
import tempfile
import io
import tarfile
import unittest
from unittest.mock import patch
import h5py
import numpy as np
from pmse.serve_browser import existing_packet,day_packet,day_archive_files,stream_day_archive
from pmse.sparse_record import pending_modes,bounded_modes
from pmse.publish_gui import MetadataCache


class BrowserTransportTests(unittest.TestCase):
    def test_day_download_preserves_native_files_and_excludes_other_days(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary).resolve()
            (root/'dmd_properties.h5').write_bytes(b'native properties')
            files=[]
            for day in ['2026-10-18','2026-10-19']:
                directory=root/f'{day}T18-00-00';directory.mkdir()
                path=directory/'pmse@1792432800.h5'
                with h5py.File(path,'w') as h:
                    h['spectra']=np.array([1+2j,3-4j],dtype='complex64')
                if day=='2026-10-19':files.append(path)
            selected=day_archive_files(root,'2026-10-19')
            self.assertEqual(selected,[root/'dmd_properties.h5']+files)
            output=io.BytesIO();stream_day_archive(output,root,selected);output.seek(0)
            with tarfile.open(fileobj=output) as archive:
                self.assertEqual(archive.getnames(),[str(p.relative_to(root)) for p in selected])
                for path in selected:
                    self.assertEqual(archive.extractfile(str(path.relative_to(root))).read(),path.read_bytes())
            with self.assertRaises(FileNotFoundError):day_archive_files(root,'2026-10-20')
            for invalid in ['../../2026-10-19','2026-1-1','2026-10-19\r\nInjected: yes']:
                with self.assertRaises(ValueError):day_archive_files(root,invalid)

    def test_one_worker_prioritizes_newest_then_backfills(self):
        choices=[(i,i+1) for i in range(10)]
        self.assertEqual(pending_modes(choices,{9},3),[(8,9),(7,8),(0,1)])
        self.assertEqual(pending_modes(choices,set(range(10)),3),[])
        self.assertEqual(bounded_modes([(0,300_000_000-32_000)]),[(0,119_968_000),(120_000_000,239_968_000),(240_000_000,299_968_000)])

    def test_day_selection_does_not_read_spectra(self):
        key=1792433051921664
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);directory=root/'2026-10-19T18-00-00';directory.mkdir()
            with h5py.File(directory/'pmse@1792432800.h5','w') as h:
                g=h.create_group(str(key))
                for name,value in dict(complete=1,mode_start_sample=key,mode_stop_sample_exclusive=key+80_000_000,
                    mode_duration_seconds=80.,coverage_seconds=80.,n_coherent_looks=40,quality='ok',
                    schema_version='pansy_pmse_sparse16_v2',sparse_index=[[0,0,0]],detector_snr_db=[4.]).items():g[name]=value
            m=day_packet(root,'2026-10-19',{})
            self.assertEqual(m['frame_count'],1)
            self.assertEqual(m['latest_mode_us'],key+80_000_000)
            self.assertEqual(m['first_mode_us'],key)
            self.assertEqual(m['axis_stop_us']-m['axis_start_us'],86400_000_000)
            self.assertEqual(day_packet(root,'2026-10-18',{})['frame_count'],0)
            self.assertIsNone(m['frame_packets'][0]['look_dominance'])
            cache=MetadataCache(root/'cache')
            cached=day_packet(root,'2026-10-19',{},cache)
            # A subsequent page refresh (even a fresh publisher process) must
            # use local metadata and never reopen unchanged archive HDF5.
            with patch('pmse.serve_browser.h5py.File',side_effect=AssertionError('Archive reopened')):
                again=day_packet(root,'2026-10-19',{},MetadataCache(root/'cache'))
            self.assertEqual(again,cached)
            with h5py.File(directory/'pmse@1792432800.h5','a') as h:
                h[str(key)]['detector_snr_db'][...]=[9.]
            updated=day_packet(root,'2026-10-19',{})
            self.assertNotEqual(m['revision'],updated['revision'])
            self.assertEqual(updated['frame_packets'][0]['snr_db'],[9.])
            self.assertEqual(day_packet(root,'2026-10-19',{},cache)['frame_packets'][0]['snr_db'],[9.])

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
