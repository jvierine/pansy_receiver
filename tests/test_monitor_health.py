import tempfile
import unittest
from pathlib import Path
from monitor.health import phase_health


class PhaseHealthTest(unittest.TestCase):
    def read(self, text):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'watchdog.log'
            path.write_text(text)
            return phase_health(path)

    def test_new_pulse_skip_does_not_refresh_measurement_age(self):
        result = self.read('2026-10-31T05:00:00+00:00 raw phase check ok: raw TX phase differs from reference by 0.6 deg on 7\n'
                           '2026-10-31T05:10:00+00:00 raw phase check skipped: no new raw TX pulse\n')
        self.assertTrue(result['ok'])
        self.assertEqual(result['max_deviation_deg'], .6)
        self.assertEqual(result['checked_unix_s'], 1793422800)

    def test_failure_overrides_good(self):
        result = self.read('2026-10-31T05:00:00+00:00 raw phase check ok: fine\n'
                           '2026-10-31T05:01:00+00:00 bad raw TX phase detected during restart cooldown: raw TX phase differs from reference by 90.0 deg on 2\n')
        self.assertFalse(result['ok'])
        self.assertEqual(result['max_deviation_deg'], 90)

    def test_failed_measurement_invalidates_good(self):
        result = self.read('2026-10-31T05:00:00+00:00 raw phase check ok: fine\n'
                           '2026-10-31T05:01:00+00:00 raw phase check skipped: only 2 valid raw TX phase channels\n')
        self.assertIsNone(result['ok'])

    def test_missing_and_malformed_are_unknown(self):
        self.assertIsNone(phase_health('/nonexistent/watchdog.log')['ok'])
        self.assertIsNone(self.read('malformed raw phase check ok: nope\n')['ok'])

    def test_tail_is_bounded(self):
        result = self.read('x' * 100000 + '\n2026-10-31T05:00:00+00:00 raw phase check ok: good\n')
        self.assertTrue(result['ok'])


if __name__ == '__main__':
    unittest.main()
