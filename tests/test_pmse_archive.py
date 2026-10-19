import tempfile
from pathlib import Path
import unittest
import h5py
from pmse.web_retention import expire,sha256
from pmse.archive_inventory import inventory
from pmse.snapshot_archive import snapshot
from pmse.replicate_sparse import catchup_files
import numpy as np


class ArchiveTests(unittest.TestCase):
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


if __name__=='__main__':unittest.main()
