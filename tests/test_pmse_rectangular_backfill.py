"""Native replacement safety: recovery, atomicity, idempotence and ordering."""
from dataclasses import asdict
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from pmse.rectangular_backfill import native_path, pending, promote, rectangular_record
from pmse.sparse_record import Config, read_complete, write_record


class RectangularBackfillTests(unittest.TestCase):
    def test_worker_services_have_separate_single_core_budgets(self):
        root=Path(__file__).resolve().parents[1]/'pmse/systemd'
        live=(root/'pansy-pmse-sparse.service').read_text()
        backfill=(root/'pansy-pmse-rectangular-backfill.service').read_text()
        self.assertIn('CPUAffinity=0',live)
        self.assertIn('--lookback-hours 3',live)
        self.assertIn('CPUAffinity=1',backfill)
        self.assertIn('--replace-hann',backfill)
        for text in (live,backfill):
            self.assertIn('CPUQuota=100%',text)
            self.assertIn('OPENBLAS_NUM_THREADS=1',text)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name) / 'live'
        self.backup = Path(self.temporary.name) / 'recovery'
        self.key = 1_792_435_200_000_000

    def record(self, key, rectangular):
        window = np.ones(64) if rectangular else np.hanning(64)
        window = (window / np.linalg.norm(window)).astype('float32')
        cfg = asdict(Config())
        if not rectangular:
            del cfg['window_kind']
        return dict(config_json=json.dumps(cfg),look_window=np.tile(window,(4,1)),
                    mode_start_sample=key,mode_stop_sample_exclusive=key+8_000_000,
                    coverage_seconds=8.,quality='ok',range_sample=np.arange(4),
                    doppler_frequency_hz=np.arange(64),slant_range_km=np.arange(4)+80.,
                    sparse_index=np.c_[np.zeros(32,int),np.arange(32),np.zeros(32,int)],
                    spectra=np.full((4,7,32),3+4j if rectangular else 5-6j,dtype='complex64'),
                    noise_power=np.ones((7,5,4),dtype='float32'),
                    noise_power_standard_error=np.full((7,5,4),.1,dtype='float32'))

    def test_replacement_preserves_other_groups_and_recovery_and_is_idempotent(self):
        old = self.record(self.key,False)
        write_record(self.root,self.key,old)
        other = self.key+10_000_000
        write_record(self.root,other,self.record(other,False))
        new = self.record(self.key,True)
        self.assertEqual(promote(self.root,self.key,new,self.backup,0),'replaced_hann')
        self.assertTrue(rectangular_record(self.root,self.key))
        np.testing.assert_array_equal(read_complete(self.root,self.key)['spectra'],new['spectra'])
        np.testing.assert_array_equal(read_complete(self.root,other)['spectra'],self.record(other,False)['spectra'])
        recovery = next(self.backup.glob('*/*.h5'))
        with h5py.File(recovery) as h:
            np.testing.assert_array_equal(h[str(self.key)]['spectra'][()],old['spectra'])
        self.assertEqual(promote(self.root,self.key,new,self.backup,0),'already_rectangular')
        self.assertEqual(len(list(self.backup.glob('*/*.h5'))),1)
        self.assertEqual(pending(self.root,[(other,other+1),(self.key,self.key+1)]),[(other,other+1)])

    def test_missing_record_is_created_with_real_native_properties(self):
        new = self.record(self.key,True)
        self.assertEqual(promote(self.root,self.key,new,self.backup,0),'created')
        self.assertTrue((self.root/'dmd_properties.h5').exists())
        self.assertTrue(rectangular_record(self.root,self.key))

    def test_failure_before_atomic_replace_leaves_old_readable_and_backup_durable(self):
        import os
        old = self.record(self.key,False)
        write_record(self.root,self.key,old)
        original_replace = os.replace
        def interrupted(source, destination):
            if Path(destination)==native_path(self.root,self.key):
                raise OSError('injected interrupted publication')
            return original_replace(source,destination)
        with patch('pmse.rectangular_backfill.os.replace',side_effect=interrupted):
            with self.assertRaisesRegex(OSError,'injected'):
                promote(self.root,self.key,self.record(self.key,True),self.backup,0)
        np.testing.assert_array_equal(read_complete(self.root,self.key)['spectra'],old['spectra'])
        self.assertFalse(rectangular_record(self.root,self.key))
        self.assertEqual(len(list(self.backup.glob('*/*.h5'))),1)
        self.assertEqual(promote(self.root,self.key,self.record(self.key,True),self.backup,0),'replaced_hann')

    def test_refuse_unsafe_replacements_and_preserve_original(self):
        old = self.record(self.key,False)
        write_record(self.root,self.key,old)
        for change in [dict(coverage_seconds=1.),dict(mode_stop_sample_exclusive=self.key+1),
                       dict(quality='insufficient_looks'),dict(range_sample=np.arange(4)+1)]:
            new = dict(self.record(self.key,True),**change)
            with self.assertRaises((ValueError,AssertionError)):
                promote(self.root,self.key,new,self.backup,0)
            np.testing.assert_array_equal(read_complete(self.root,self.key)['spectra'],old['spectra'])
        with self.assertRaises(ValueError):
            promote(self.root,self.key,self.record(self.key,True),self.root/'backup',0)
        with self.assertRaises(RuntimeError):
            promote(self.root,self.key,self.record(self.key,True),self.backup,1e12)


if __name__=='__main__':
    unittest.main()
