import importlib.util
from pathlib import Path
import tempfile
import unittest
import h5py
import numpy as np

spec = importlib.util.spec_from_file_location("shadow_latest", Path(__file__).parents[1]/"shadow_latest.py")
shadow = importlib.util.module_from_spec(spec)
spec.loader.exec_module(shadow)


class ShadowTests(unittest.TestCase):
    def test_validation_rejects_bad_model(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"sample.h5"
            with h5py.File(path, "w") as data:
                group = data.create_group("windows/0")
                group["sample_idx"] = np.arange(6, dtype=np.int64)*20000
                group["track_id"] = np.zeros(6, dtype=np.int32)
                group["track_epoch_sample"] = np.array([0], dtype=np.int64)
                group["range_km"] = np.full(6, 100, dtype=np.float32)
                group["radial_velocity_m_s"] = np.zeros(6, dtype=np.float32)
                group["track_range_at_epoch_km"] = np.array([100], dtype=np.float32)
                group["track_radial_velocity_m_s"] = np.array([0], dtype=np.float32)
                group["track_radial_acceleration_m_s2"] = np.array([0], dtype=np.float32)
            self.assertEqual(shadow.validate(path), (6, 6, 1))
            with h5py.File(path, "r+") as data:
                data["windows/0/track_range_at_epoch_km"][0] = 101
            with self.assertRaisesRegex(ValueError, "violates"):
                shadow.validate(path)

    def test_window(self):
        self.assertEqual(shadow.sample_window(0, 100_100_000, 101_000_000), (38_000_000, 98_000_000))
        self.assertIsNone(shadow.sample_window(50_000_000, 100_100_000, 101_000_000))
        self.assertIsNone(shadow.sample_window(0, 100_100_000, 4_000_000_000))
        self.assertIsNone(shadow.sample_window(0, 100_100_000, 90_000_000))

    def test_retention_only_owned_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for i in range(5):
                (root/f"hourly_{i:03}.h5").write_bytes(b"1234")
            (root/"replay_reference.h5").write_bytes(b"keep")
            shadow.prune(root, keep=3, max_bytes=8)
            self.assertEqual(sorted(p.name for p in root.iterdir()),
                             ["hourly_003.h5", "hourly_004.h5", "replay_reference.h5"])


if __name__ == "__main__":
    unittest.main()
