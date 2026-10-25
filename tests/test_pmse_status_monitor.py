import ast
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import digital_rf as drf
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pmse.status_monitor import MONITOR_FIELDS, read_sparse_status, plot_sparse_status, zenith_cells


class SparseStatusTests(unittest.TestCase):
    def row(self, stop, index=((0, 1, 1), (0, 2, 1), (1, 1, 0)), snr=(4., 9., 20.), complete=1):
        return dict(complete=complete, mode_stop_sample_exclusive=stop,
                    sparse_index=np.asarray(index, dtype=np.int32).reshape(-1, 3),
                    detector_snr_db=np.asarray(snr),
                    nominal_altitude_km=np.tile([80., 81., 82.], (5, 1)),
                    spectra=np.ones((2, 7, len(snr)), dtype=np.complex64))

    def test_real_metadata_ignores_incomplete_tail_and_counts_noise_only_modes(self):
        with tempfile.TemporaryDirectory() as temporary:
            writer = drf.DigitalMetadataWriter(temporary, 3600, 600, 1_000_000, 1, "pmse")
            key = 1792920000000000
            writer.write([key], [self.row(key + 80_000_000)])
            writer.write([key + 100_000_000], [self.row(key + 180_000_000, index=(), snr=())])
            writer.write([key + 200_000_000], [self.row(key + 280_000_000, complete=0)])
            # A crashed writer may not yet have created the other datasets.
            for path in Path(temporary).glob("*/pmse@*.h5"):
                with h5py.File(path, "a") as handle:
                    if str(key + 200_000_000) in handle:
                        del handle[str(key + 200_000_000)]["sparse_index"]
            with patch("pmse.status_monitor.h5py.File", wraps=h5py.File) as opener:
                bounds, records = read_sparse_status(temporary)
            self.assertTrue(opener.call_args_list)
            self.assertTrue(all(call.args[1] == "r" for call in opener.call_args_list))
            self.assertEqual(bounds, [key, key + 180_000_000])
            self.assertEqual(len(records), 2)
            self.assertEqual(set(records[0][2]), set(MONITOR_FIELDS))
            self.assertNotIn("spectra", records[0][2])
            self.assertEqual(len(zenith_cells(records[-1][2])[2]), 0)

    def test_unreadable_input_is_never_deleted_or_repaired(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "2026-10-25T12-00-00"
            directory.mkdir()
            path = directory / "pmse@1792929600.h5"
            content = b"interrupted HDF5 transfer"
            path.write_bytes(content)
            with self.assertWarns(UserWarning):
                self.assertEqual(read_sparse_status(temporary), ([-1, -1], []))
            self.assertEqual(path.read_bytes(), content)

    def test_zenith_doppler_peak_preserves_missing_cells_and_snr_units(self):
        low, high, snr = zenith_cells(self.row(80_000_000))
        np.testing.assert_array_equal(low, [80.5])
        np.testing.assert_array_equal(high, [81.5])
        np.testing.assert_array_equal(snr, [9.])

    def test_plot_does_not_fill_gaps_or_old_data(self):
        end = 1792920000.
        key = int((end - 1000) * 1e6)
        row = self.row(key + 80_000_000)
        fig, ax = plt.subplots()
        try:
            plot_sparse_status(fig, ax, [(key - 2 * 86400_000_000, key - 86400_000_000, row),
                                         (key, key + 80_000_000, row),
                                         (key + 200_000_000, key + 280_000_000, row)], end)
            self.assertEqual(len(ax.collections), 1)
            self.assertEqual(len(ax.collections[0].get_paths()), 2)
            np.testing.assert_array_equal(ax.collections[0].get_array(), [9., 9.])
            fig.canvas.draw()
        finally:
            plt.close(fig)

    def test_all_noise_plot_is_valid_without_fabricated_signal(self):
        key = 1792920000000000
        fig, ax = plt.subplots()
        try:
            plot_sparse_status(fig, ax, [(key, key + 80_000_000, self.row(key + 80_000_000, (), ()))],
                               key / 1e6 + 100)
            self.assertEqual(len(ax.collections), 0)
            self.assertIn("No retained", ax.texts[0].get_text())
            fig.canvas.draw()
        finally:
            plt.close(fig)

    def test_live_monitor_no_longer_reads_xc_and_installer_requires_opt_in(self):
        root = Path(__file__).resolve().parents[1]
        source = (root / "status_plot.py").read_text()
        tree = ast.parse(source)
        live = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "plot_status")
        live_source = ast.get_source_segment(source, live)
        self.assertNotIn("get_xc(", live_source)
        self.assertNotIn("xc_metadata_dir", live_source)
        self.assertIn('"Spectra"', live_source)
        services = next(node.value for node in tree.body if isinstance(node, ast.Assign)
                        and any(isinstance(t, ast.Name) and t.id == "PANSY_SERVICE_NAMES" for t in node.targets))
        self.assertIn("pansy-pmse-sparse.service", ast.literal_eval(services))
        self.assertNotIn("pansy-meso-xc.service", ast.literal_eval(services))
        installer = (root / "receiver/install_user_service.sh").read_text()
        setup = installer[installer.index("UNITS=("):installer.index('mkdir -p "$CONFIG_DIR"')]
        for opt_in, expected in [("0", False), ("1", True)]:
            result = subprocess.run(["bash", "-c", f'PANSY_ENABLE_LEGACY_XC={opt_in}\n' + setup +
                                     '\nprintf "%s\\n" "${UNITS[@]}"'],
                                    check=True, capture_output=True, text=True)
            self.assertEqual("pansy-meso-xc.service" in result.stdout.splitlines(), expected)


if __name__ == "__main__":
    unittest.main()
