import tempfile
from pathlib import Path
import unittest
import h5py
import numpy as np
from monitor.event import save_event

class CurrentEventTests(unittest.TestCase):
    def test_only_one_event_is_retained_and_older_completion_cannot_replace_it(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'latest-event.h5'
            data=np.arange(12,dtype=np.float32).reshape(3,4)
            save_event(20,[1,2,3,4],[80,81,82],data,[1,2],path)
            save_event(10,[1,2,3,4],[80,81,82],data+100,[3,4],path)
            with h5py.File(path) as h:
                self.assertEqual(int(h['key'][()]),20)
                np.testing.assert_array_equal(h['power_snr'][()],data)
            save_event(30,[1,2,3,4],[80,81,82],data+1,[1,2],path)
            with h5py.File(path) as h:self.assertEqual(int(h['key'][()]),30)
            self.assertEqual(len(list(Path(directory).glob('*.h5'))),1)
            self.assertFalse(path.with_suffix('.h5.partial').exists())

if __name__=='__main__':unittest.main()
