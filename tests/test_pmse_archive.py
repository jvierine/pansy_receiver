import tempfile
from pathlib import Path
import unittest
import h5py
from pmse.web_retention import expire,sha256
from pmse.archive_inventory import inventory
from pmse.snapshot_archive import snapshot
from pmse.replicate_sparse import catchup_files,changed_files,historical_turn
import numpy as np


class ArchiveTests(unittest.TestCase):
    def test_same_size_changed_native_chunk_is_forced_to_web(self):
        import os
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);relative='2026-10-20T00-00-00/pmse@600.h5'
            path=root/relative;path.parent.mkdir();path.write_bytes(b'Hann')
            receipt=root/'receipt.h5';st=path.stat()
            with h5py.File(receipt,'w') as h:
                g=h.require_group('files/'+relative)
                g.attrs.update(mtime_ns=st.st_mtime_ns,bytes=st.st_size)
            self.assertEqual(changed_files(root,receipt,[relative]),[])
            path.write_bytes(b'Rect')
            os.utime(path,ns=(st.st_atime_ns,st.st_mtime_ns+1_000_000))
            self.assertEqual(changed_files(root,receipt,[relative]),[relative])

    def test_snapshot_compacts_without_changing_scientific_values(self):
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/'source';source.mkdir();dest=Path(directory)/'copy'
            path=source/'pmse@600.h5'
            with h5py.File(path,'w') as h:
                h.attrs['provenance']='test';h['abandoned']=np.ones(100000)
                h['keep']=np.arange(1000,dtype='complex64');h['keep'].attrs['units']='V'
                del h['abandoned']
            self.assertEqual(snapshot(source,dest,reserve_gib=0),1)
            self.assertLess((dest/path.name).stat().st_size,path.stat().st_size/2)
            with h5py.File(path) as src,h5py.File(dest/path.name) as dst:
                np.testing.assert_array_equal(src['keep'][()],dst['keep'][()])
                self.assertEqual(dict(src.attrs),dict(dst.attrs))
                self.assertEqual(dict(src['keep'].attrs),dict(dst['keep'].attrs))
            self.assertEqual(snapshot(source,dest,reserve_gib=0),0)
            self.assertEqual(snapshot(source,dest,reserve_gib=0,force=True),1)

    def test_recent_chunks_precede_historical_backfill(self):
        with tempfile.TemporaryDirectory() as directory:
            files=[(f'2026-10-19T20-00-00/pmse@{t}.h5',100,1.) for t in [20000,19400,18800,600]]
            recent,cold=catchup_files(files,directory,20000)
            self.assertEqual(recent,[files[1][0]])
            self.assertEqual(cold,[files[-1][0]])

    def test_oldest_transfer_and_processing_frontier_both_get_turns(self):
        historical=['newest-backfilled.h5','middle.h5','oldest.h5']
        self.assertEqual(historical_turn(historical,0),['oldest.h5'])
        self.assertEqual(historical_turn(historical,1),['newest-backfilled.h5'])
        self.assertEqual(historical_turn(historical,2),['oldest.h5'])
        advanced=['even-newer.h5']+historical
        # A large partially transferred frontier chunk must finish, rather
        # than chasing a new file on every turn and exposing no complete data.
        self.assertEqual(historical_turn(advanced,3,'newest-backfilled.h5'),['newest-backfilled.h5'])
        self.assertEqual(historical_turn(advanced,4,'newest-backfilled.h5'),['oldest.h5'])
        self.assertEqual(historical_turn(advanced,5,'already-complete.h5'),['even-newer.h5'])
        self.assertEqual(historical_turn([],3),[])
        actual=['2026-10-19T18-00-00/pmse@1792433400.h5',
                '2026-10-14T19-00-00/pmse@1792006200.h5',
                '2026-10-12T05-00-00/pmse@1791782400.h5']
        self.assertEqual(historical_turn(actual,1,processed_frontier_us=1792006419897755),[actual[1]])
        self.assertEqual(historical_turn(actual,2,processed_frontier_us=1792006419897755),[actual[2]])
        with tempfile.TemporaryDirectory() as directory:
            files=[(f'2026-10-19T20-00-00/pmse@{t}.h5',100,1.) for t in [50000,49400,19400,600]]
            recent,cold=catchup_files(files,directory,50000)
            self.assertEqual(recent,[files[1][0]])
            self.assertEqual(cold,[files[2][0],files[3][0]])

    def test_expiration_requires_matching_retained_source_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)/'secondary';hour=root/'1970-01-01T00-00-00';hour.mkdir(parents=True)
            a=hour/'pmse@600.h5';b=hour/'pmse@1200.h5';a.write_bytes(b'preserved');b.write_bytes(b'mismatch')
            receipt=Path(directory)/'receipt.h5'
            with h5py.File(receipt,'w') as h:
                g=h.require_group('files').require_group(str(a.relative_to(root)));g.attrs['sha256']=sha256(a)
                g=h['files'].require_group(str(b.relative_to(root)));g.attrs['sha256']='incorrect'
            result=expire(root,receipt,now=1_000_000,dry_run=True)
            self.assertEqual(result['expired_files'],[str(a.relative_to(root))])
            self.assertEqual(result['unverified_files_preserved'],[str(b.relative_to(root))])
            self.assertTrue(a.exists())
            with self.assertRaises(ValueError):expire(root,receipt,now=1_000_000,dry_run=False)

    def test_inventory_ignores_incomplete_records_and_no_spectra_reads(self):
        with tempfile.TemporaryDirectory() as directory:
            hour=Path(directory)/'2026-10-19T00-00-00';hour.mkdir()
            path=hour/'pmse@1000.h5'
            with h5py.File(path,'w') as h:
                for key,complete in [(1_000_000_000,1),(1_100_000_000,0)]:
                    g=h.create_group(str(key));g['complete']=complete;g['mode_stop_sample_exclusive']=key+86_000_000
            result=inventory(directory)
            self.assertEqual(result['records'],1);self.assertEqual(result['files'],1)
            self.assertEqual(result['observed_span_seconds'],86)
            self.assertIsNone(result['projected_bytes_per_day'])
            self.assertIsNone(result['complex_sample_retained_fraction'])

    def test_retention_weights_looks_and_counts_empty_records_without_reading_spectra(self):
        with tempfile.TemporaryDirectory() as directory:
            hour=Path(directory)/'2026-10-19T00-00-00';hour.mkdir()
            path=hour/'pmse@1000.h5'
            with h5py.File(path,'w') as h:
                for key,looks,cells,complete in [(1_000_000_000,2,10,1),(1_100_000_000,8,0,1),(1_200_000_000,20,100,0)]:
                    g=h.create_group(str(key));g['complete']=complete;g['mode_stop_sample_exclusive']=key+86_000_000
                    g['n_tested_cells']=100
                    g.create_dataset('noise_power',shape=(7,5,10),dtype='float32')
                    g.create_dataset('radial_velocity_mps',shape=(4,),dtype='float64')
                    g.create_dataset('slant_range_km',shape=(10,),dtype='float64')
                    # An unavailable external dataset ensures inventory only reads its shape.
                    g.create_dataset('spectra',shape=(looks,7,cells),dtype='complex64',
                                     external=[('unavailable-spectra.bin',0,h5py.h5f.UNLIMITED)])
            for _ in range(2):  # Cached and fresh inventories must agree.
                result=inventory(directory)
                self.assertEqual(result['retention_records'],2)
                self.assertEqual(result['retained_complex_samples'],140)
                self.assertEqual(result['tested_complex_samples'],7000)
                self.assertAlmostEqual(result['complex_sample_retained_fraction'],.02)
                self.assertEqual(result['dense_complex_samples'],14000)
                self.assertEqual(result['dense_retention_records'],2)
                self.assertAlmostEqual(result['dense_complex_sample_retained_fraction'],.01)
                self.assertEqual(result['hourly_dense_retained_fraction'],[.01])
            with h5py.File(path,'a') as h:
                g=h.create_group('1300000000');g['complete']=1;g['mode_stop_sample_exclusive']=1_400_000_000
            self.assertIsNone(inventory(directory)['complex_sample_retained_fraction'])
            self.assertIsNone(inventory(directory)['dense_complex_sample_retained_fraction'])
            self.assertEqual(inventory(directory)['hourly_dense_retained_fraction'],[None])

    def test_hourly_retention_uses_sample_totals_and_leaves_missing_hours_blank(self):
        with tempfile.TemporaryDirectory() as directory:
            hour=Path(directory)/'1970-01-01T00-00-00';hour.mkdir()
            for seconds,looks,cells in [(600,2,10),(1200,8,0),(7200,2,40)]:
                with h5py.File(hour/f'pmse@{seconds}.h5','w') as h:
                    g=h.create_group(str(seconds*1_000_000));g['complete']=1
                    g['mode_stop_sample_exclusive']=seconds*1_000_000+86_000_000
                    g['n_tested_cells']=100
                    g.create_dataset('noise_power',shape=(7,5,10),dtype='float32')
                    g.create_dataset('radial_velocity_mps',shape=(4,),dtype='float64')
                    g.create_dataset('slant_range_km',shape=(10,),dtype='float64')
                    g.create_dataset('spectra',shape=(looks,7,cells),dtype='complex64')
            for _ in range(2):
                result=inventory(directory)
                self.assertEqual(result['hourly_unix_s'],[0,3600,7200])
                self.assertEqual(result['hourly_retained_complex_samples'],[140,0,560])
                self.assertEqual(result['hourly_dense_complex_samples'],[14000,0,2800])
                self.assertEqual(result['hourly_dense_retained_fraction'],[.01,None,.2])


if __name__=='__main__':unittest.main()
