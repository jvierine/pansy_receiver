"""Regression: metadata-only monitoring must not require raw-RF/orbit tools."""
import subprocess
import sys
import unittest
from pathlib import Path


class MonitorImports(unittest.TestCase):
    def test_monitor_imports_without_raw_rf_orbit_runtime(self):
        code = '''
import importlib.abc
import sys
class BlockOptional(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'digital_rf', 'dasst', 'pansy_ballistic'}:
            raise AssertionError('Unnecessary monitor dependency: ' + fullname)
sys.meta_path.insert(0, BlockOptional())
import publish_radiant_monitor
import plot_fitted_radiant_distribution as plots
assert list(plots.wrap180([0, 181, 359])) == [0, -179, -1]
print('metadata-only publisher imports without raw-RF/orbit packages')
'''
        subprocess.run([sys.executable, '-c', code],
                       cwd=Path(__file__).resolve().parents[1], check=True)


if __name__ == '__main__':
    unittest.main()
