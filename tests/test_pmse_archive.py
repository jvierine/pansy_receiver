import tempfile
from pathlib import Path
import unittest
import h5py
from pmse.web_retention import expire,sha256
from pmse.archive_inventory import inventory


class ArchiveTests(unittest.TestCase):
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
