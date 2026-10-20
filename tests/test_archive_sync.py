import hashlib
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from pmse.archive_sync import pending_files, update_receipt, pull_once, publish_once
from pmse.replicate_sparse import rsync
from types import SimpleNamespace


class ArchiveSyncTests(unittest.TestCase):
    def make_archive(self, root):
        path = root / '2026-10-20T12-00-00/pmse@1792497600.h5'
        path.parent.mkdir(parents=True)
        with h5py.File(path, 'w') as h:
            h['spectra'] = np.array([1+2j, 3-4j], dtype='complex64')
        (root / 'dmd_properties.h5').write_bytes(b'properties')
        return path

    def test_receipt_detects_same_size_reprocessed_snapshot(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / 'archive'
            path = self.make_archive(root)
            receipt = Path(temporary) / 'receipt.h5'
            update_receipt(root, receipt, [path])
            self.assertEqual(pending_files(root, receipt), [])
            size = path.stat().st_size
            old_time = path.stat().st_mtime_ns
            with h5py.File(path, 'r+') as h:
                h['spectra'][0] = 8+9j
            os.utime(path, ns=(old_time+1_000_000_000, old_time+1_000_000_000))
            self.assertEqual(path.stat().st_size, size)
            self.assertEqual(pending_files(root, receipt), [path])
            update_receipt(root, receipt, [path])
            with h5py.File(receipt) as h:
                self.assertEqual(h['files/'+str(path.relative_to(root))].attrs['sha256'],
                                 hashlib.sha256(path.read_bytes()).hexdigest())

    @unittest.skipUnless(shutil.which('rsync'), 'rsync is required')
    def test_checksum_pull_preserves_native_bytes_after_equal_size_and_mtime_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / 'source'
            destination = Path(temporary) / 'destination'
            destination.mkdir()
            path = self.make_archive(source)
            rsync(str(source)+'/', str(destination)+'/', checksum=True)
            before = path.stat()
            with h5py.File(path, 'r+') as h:
                h['spectra'][0] = 11-12j
            os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
            self.assertEqual(path.stat().st_size, before.st_size)
            rsync(str(source)+'/', str(destination)+'/', checksum=True)
            self.assertEqual(path.read_bytes(), (destination/path.relative_to(source)).read_bytes())

    def test_partial_pull_never_acknowledges_a_complete_archive(self):
        with tempfile.TemporaryDirectory() as temporary:
            a = SimpleNamespace(host='j@juha.no', archive_root='/mnt/shovel/rvd/pmse_sparse16',
                                root=Path(temporary), receipt=Path(temporary)/'receipt.h5',
                                metrics=Path(temporary)/'metrics.h5', bandwidth=1024)
            with patch('pmse.archive_sync.check_archive'), patch('pmse.archive_sync.save_metrics'), \
                 patch('pmse.archive_sync.rsync', return_value={'partial':True}) as transfer:
                pull_once(a)
            self.assertEqual(transfer.call_count, 1)
            self.assertTrue(transfer.call_args.kwargs['checksum'])
            self.assertTrue(transfer.call_args.args[0].startswith('j@juha.no:'))
            self.assertFalse(a.receipt.exists())

    def test_partial_source_publication_does_not_mark_snapshot_delivered(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)/'snapshot'
            path = self.make_archive(root)
            a = SimpleNamespace(source=root, root=root, host='j@juha.no',
                                archive_root='/mnt/shovel/rvd/pmse_sparse16',
                                receipt=Path(temporary)/'receipt.h5',
                                metrics=Path(temporary)/'metrics.h5', bandwidth=10)
            with patch('pmse.snapshot_archive.snapshot'), patch('pmse.processing_status.status', return_value={}), \
                 patch('pmse.archive_sync.check_archive'), patch('pmse.archive_sync.save_metrics'), \
                 patch('pmse.archive_sync.rsync', return_value={'partial':True}):
                publish_once(a, 1)
            self.assertEqual(pending_files(root, a.receipt), [path])

    def test_service_roles_have_no_reverse_forwarding(self):
        root = Path(__file__).resolve().parents[1]/'pmse/systemd'
        self.assertFalse((root/'pansy-rvd-browser-tunnel.service').exists())
        source = (root/'pansy-rvd-source-publish.service').read_text()
        pull = (root/'pansy-pmse-sparse-backup.service').read_text()
        api = (root/'pansy-rvd-browser-api.service').read_text()
        self.assertIn('pmse.archive_sync publish', source)
        self.assertIn('--bandwidth 10', source)
        self.assertIn('pmse.archive_sync pull', pull)
        self.assertIn('/mnt/shovel/rvd/pmse_sparse16', api)
        for text in [source, pull, api]:
            self.assertNotIn(' -R ', text)
            self.assertNotIn('RemoteForward', text)


if __name__ == '__main__':
    unittest.main()
