"""Real Digital RF -> production Rust processor -> native archive -> voltages.

No mocked readers, detector, writer, DSP, or inverse transform. The voltage
oracle uses direct measured-TX dot products and explicit DFT matrices, not
the production decoder/window/FFT. Requires digital-rf and the built Rust DSP.
"""
from dataclasses import replace
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import digital_rf as drf
import h5py
import numpy as np

import pansy_modes
from pmse.reconstruct_echoes import reconstruct
from pmse.sparse_record import Config, axes, process_mode, read_complete, write_record


class DecodedVoltageEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="pansy-voltage-e2e-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name)
        cls.cfg = replace(Config(), nfft=64, complementary_pairs=True,
                          velocity_limit_mps=None, altitude_min_km=80.,
                          altitude_max_km=84., range_step_samples=8, min_looks=4)
        cfg = cls.cfg
        gates = axes(cfg)[0]
        cls.key = 1_792_435_200_000_000
        cls.cycle_counts = [32, 32, 32, 16]  # Last look must be zero-padded.
        cycles = sum(cls.cycle_counts)
        samples = cycles * 20 * cfg.ipp_samples
        raw = cls.root / "raw"
        tx_metadata = cls.root / "tx"
        tx_metadata.mkdir()
        raw.mkdir()
        rng = np.random.default_rng(20261020)
        codes = np.repeat(np.asarray(pansy_modes.get_m_mode()["codes"][:2]), 8, axis=1)
        family = np.tile(np.repeat([0, 1, 1, 0], 5), cycles)
        p = np.arange(cycles * 20)
        measured = np.zeros((len(p), cfg.ipp_samples), dtype="complex64")
        measured[:, :128] = (codes[family] *
            ((.8 + .2 * np.cos(p * .13)) * np.exp(1j * p * .17))[:, None])
        # Keep a metadata marker for this cycle, but reject its zero TX in DSP.
        bad_cycle = cls.cycle_counts[0] + 5
        measured[20 * bad_cycle:20 * (bad_cycle + 1)] = 0
        energy = np.sum(abs(measured[:, :cfg.tx_samples].astype("complex128"))**2, axis=1)
        markers = cls.key + np.arange(cycles, dtype="int64") * 20 * cfg.ipp_samples
        writer = drf.DigitalMetadataWriter(str(tx_metadata), 3600, 600, 1_000_000, 1, "tx")
        # A missing marker must be recovered without inventing voltage samples.
        writer.write(np.delete(markers, 9).tolist(), {"id": 1})
        del writer

        def write_channel(channel, values):
            directory = raw / f"ch{channel:03d}"
            directory.mkdir()
            writer = drf.DigitalRFWriter(str(directory), np.complex64, 3600, 1000,
                cls.key, 1_000_000, 1, marching_periods=False)
            writer.rf_write(np.ascontiguousarray(values.reshape(samples)))
            writer.close()

        write_channel(7, measured)
        cls.reference = [np.zeros((7, 5, cfg.nfft, len(gates)), dtype="complex128")
                         for _ in cls.cycle_counts]
        cls.windows = []
        offset = 0
        for cycle_count in cls.cycle_counts:
            valid = np.ones(cycle_count, dtype=bool)
            if offset <= bad_cycle < offset + cycle_count:
                valid[bad_cycle - offset] = False
            window = np.repeat(valid, 2).astype(float)
            window /= np.sqrt(np.sum(window**2))
            cls.windows.append(np.pad(window, (0, cfg.nfft - len(window))))
            offset += cycle_count
        frequency = np.arange(cfg.nfft) - cfg.nfft // 2
        # Forward matrix [frequency,time]; the inverse is its conjugate / F.
        cls.dft = np.exp(-2j * np.pi * frequency[:, None] *
                         np.arange(cfg.nfft)[None, :] / cfg.nfft)
        for channel in range(7):
            rx = ((rng.normal(size=measured.shape) + 1j * rng.normal(size=measured.shape))
                  * .3 / np.sqrt(2)).astype("complex64")
            beam = p % 5
            pulse_per_beam = p // 5
            # Distinct beam Dopplers, channel phases, ranges, amplitudes, looks.
            tones = np.array([-19, -7, 0, 11, 24])
            for r, scale in [(1, 1.), (len(gates) - 1, .65)]:
                amplitude = (20 * scale * (1 + .09 * channel) *
                    np.exp(1j * (.31 * channel + .23 * beam + .41 * r +
                                2 * np.pi * tones[beam] * pulse_per_beam / (2 * cfg.nfft))))
                amplitude *= 1 + .15 * np.sin(pulse_per_beam * .07)
                gain = np.divide(amplitude, np.sqrt(energy),
                                 out=np.zeros_like(amplitude), where=energy > 0)
                lag = gates[r]
                rx[:, lag:lag + cfg.tx_samples] += (gain[:, None] *
                    measured[:, :cfg.tx_samples]).astype("complex64")
            write_channel(channel, rx)
            # Independent, direct linear matched filtering at physical lags.
            kernel = np.divide(measured[:, :cfg.tx_samples].astype("complex128"),
                np.sqrt(energy)[:, None], out=np.zeros((len(p), cfg.tx_samples), complex),
                where=energy[:, None] > 0)
            decoded = np.stack([np.sum(rx[:, lag:lag + cfg.tx_samples] * kernel.conj(), axis=1)
                                for lag in gates], axis=-1)
            offset = 0
            for look, cycle_count in enumerate(cls.cycle_counts):
                pulses = decoded[20 * offset:20 * (offset + cycle_count)]
                # Independently assign chronological sweeps to their five beams.
                ordered = np.stack([pulses[b::5] for b in range(5)])
                pair = (ordered[:, ::2] + ordered[:, 1::2]) / np.sqrt(2.)
                cls.reference[look][channel, :, :pair.shape[1], :] = (
                    pair * cls.windows[look][:pair.shape[1]][None, :, None])
                offset += cycle_count
        cls.reference = np.stack(cls.reference)
        cls.dense_reference = np.einsum("ft,lcbtr->lcbfr", cls.dft, cls.reference)
        cls.produced = process_mode(raw, tx_metadata, cls.key, int(markers[-1]), cfg, backend="rust")
        cls.archive = cls.root / "archive"
        if not write_record(cls.archive, cls.key, cls.produced):
            raise AssertionError("First native archive write must publish a new record")
        cls.saved = read_complete(cls.archive, cls.key)

    def test_native_archive_and_actual_production_coverage(self):
        saved = self.saved
        self.assertEqual(int(saved["complete"]), 1)
        self.assertIn("Rust DSP", saved["generator"])
        self.assertEqual(json.loads(saved["skipped_looks_json"]), [])
        self.assertEqual(int(saved["n_coherent_looks"]), 4)
        self.assertEqual(json.loads(saved["config_json"])["window_kind"], "rectangular")
        self.assertIn("rectangular", saved["normalization"])
        self.assertEqual(int(saved["n_recovered_cycle_markers"]), 1)
        self.assertEqual(int(np.count_nonzero(saved["tx_cycle_verified"] == 0)), 1)
        np.testing.assert_array_equal(saved["scheduled_pulses_per_beam_per_look"], [128, 128, 128, 64])
        np.testing.assert_allclose(saved["look_window"], self.windows, atol=3e-8)
        # No artificial endpoint zeros, including the short final look.
        self.assertTrue(np.all(saved["look_window"][:, 0] > 0))
        for look, cycles in enumerate(self.cycle_counts):
            self.assertGreater(saved["look_window"][look, 2 * cycles - 1], 0)
        np.testing.assert_array_equal(saved["spectra"], self.produced["spectra"])
        np.testing.assert_array_equal(saved["sparse_index"], self.produced["sparse_index"])
        self.assertEqual(saved["spectra"].dtype, np.dtype("complex64"))
        self.assertFalse(write_record(self.archive, self.key, self.produced))
        files = list(self.archive.glob("*T*/pmse@*.h5"))
        self.assertEqual(len(files), 1)
        with h5py.File(files[0]) as h:
            self.assertEqual(h[str(self.key)]["spectra"].compression, "gzip")

    def test_complex_voltage_against_independent_decoded_reference(self):
        index = self.saved["sparse_index"]
        self.assertGreater(len(index), 0, "Empty archives cannot prove recovery")
        self.assertLess(len(index), 5 * self.cfg.nfft * len(axes(self.cfg)[0]))
        self.assertEqual(set(index[:, 0]), set(range(5)))
        self.assertGreater(len(set(index[:, 2])), 1)
        bins = index[:, 1] - self.cfg.nfft // 2
        self.assertTrue(np.any(bins < 0) and np.any(bins > 0))
        # Check the saved complex coefficients against direct decoding + DFT,
        # not merely FFT(IFFT(saved)): that would miss forward-path mistakes.
        oracle = self.dense_reference[:, :, index[:, 0], index[:, 1], index[:, 2]]
        np.testing.assert_allclose(self.saved["spectra"], oracle, rtol=3e-5, atol=8e-5)
        errors = []
        for beam in range(5):
            mask = np.zeros((self.cfg.nfft, len(axes(self.cfg)[0])), bool)
            selected = index[index[:, 0] == beam]
            mask[selected[:, 1], selected[:, 2]] = True
            expected = np.einsum("ft,lcfr->lctr", self.dft.conj(),
                self.dense_reference[:, :, beam] * mask[None, None]) / self.cfg.nfft
            actual = reconstruct(self.saved, beam)
            np.testing.assert_allclose(actual, expected, rtol=4e-5, atol=1e-5)
            relative_error = np.linalg.norm(actual - expected) / np.linalg.norm(expected)
            self.assertLess(relative_error, 2e-6)
            errors.append(float(relative_error))
            # Explicit amplitude and phase contracts, excluding near-zero phases.
            np.testing.assert_allclose(abs(actual), abs(expected), rtol=4e-5, atol=1e-5)
            strong = abs(expected) > .01 * abs(expected).max()
            self.assertLess(float(abs(np.angle(actual[strong] * expected[strong].conj())).max()), 2e-5)
            self.assertGreater(np.linalg.norm(actual - expected.conj()), .1 * np.linalg.norm(expected))
            # Constant normalization is reversible, even with gaps/short looks.
            pair_voltage = reconstruct(self.saved, beam, restore_pair_scale=True)
            pair_expected = expected / np.max(self.windows, axis=1)[:, None, None, None]
            np.testing.assert_allclose(pair_voltage, pair_expected, rtol=4e-5, atol=8e-5)
        print("\nDecoded voltage E2E: retained cells=" + str(len(index)) +
              ", relative complex L2 errors by beam=" + str(errors))

    def test_public_reconstruction_command_reads_native_archive(self):
        output = self.root / "decoded-voltage.h5"
        subprocess.run([sys.executable, "-m", "pmse.reconstruct_echoes", str(self.archive),
                        "--key", str(self.key), "--beam", "3", "--output", str(output)],
                       cwd=Path(__file__).resolve().parents[1], check=True, capture_output=True)
        with h5py.File(output) as h:
            np.testing.assert_array_equal(h["decoded_echo"][()], reconstruct(self.saved, 3))
            np.testing.assert_array_equal(h["look_window"][()], self.saved["look_window"])
            self.assertEqual(h.attrs["axis_order"], "coherent_look,rx_channel,pair_time,range")
        subprocess.run([sys.executable, "-m", "pmse.reconstruct_echoes", str(self.archive),
                        "--key", str(self.key), "--beam", "3", "--pair-voltage", "--output", str(output)],
                       cwd=Path(__file__).resolve().parents[1], check=True, capture_output=True)
        with h5py.File(output) as h:
            np.testing.assert_array_equal(h["decoded_echo"][()], reconstruct(self.saved, 3, restore_pair_scale=True))
            self.assertEqual(h.attrs["amplitude_domain"], "decoded complementary pair voltage")


if __name__ == "__main__":
    unittest.main()
