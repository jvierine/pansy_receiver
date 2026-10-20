"""Render real native auxiliary plots with wide scientific-notation ticks."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest


class AuxiliaryLayoutTests(unittest.TestCase):
    def test_high_sky_temperature_and_tiny_retention_fit_without_disabling_checks(self):
        renderer=Path(__file__).resolve().parents[1]/'pmse/plotting/target/release/pansy-pmse-plots'
        self.assertTrue(renderer.is_file(),'Build the native plotting release gate first')
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);spec=root/'input.json';out=root/'axes'
            spec.write_text(json.dumps(dict(first_mode_us=0,latest_mode_us=7_200_000_000,
                time_tick_labels=['00:00','01:00','02:00'],noise_model_data=dict(
                    times_h=[0.,1.,2.],time_tick_labels=['00:00','01:00','02:00'],
                    sky_temperature_k=[[9000.+3000.*t+100.*b for b in range(5)] for t in range(3)],
                    noise_power=[[[100.+10.*t+b+c for c in range(7)] for b in range(5)] for t in range(3)]),
                transfer=dict(archive=dict(hourly_unix_s=[0.,3600.],hourly_bytes=[1e6,2e6],
                    hour_tick_labels=['10-12 00:00','10-12 01:00','10-12 02:00'],
                    hourly_dense_retained_fraction=[0.,1e-7],latest_mode_us=7_200_000_000)))))
            result=subprocess.run([str(renderer),str(spec),str(out)],capture_output=True,text=True,timeout=120)
            self.assertEqual(result.returncode,0,result.stderr)
            for width in [640,960]:
                for name in [f'noise-module-0-{width}.svg',f'archive-growth-{width}.svg']:
                    self.assertIn('<svg',(out/name).read_text())
            # Metadata-only publication must not require an RTI time axis,
            # noise arrays, or render any of the unrelated static frames.
            archive_only=root/'archive-only'
            manifest=json.loads(spec.read_text())
            spec.write_text(json.dumps(dict(transfer=manifest['transfer'])))
            result=subprocess.run([str(renderer),'--archive-only',str(spec),str(archive_only)],
                                  capture_output=True,text=True,timeout=120)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(sorted(p.name for p in archive_only.glob('*.svg')),
                             ['archive-growth-640.svg','archive-growth-960.svg'])
            self.assertFalse((archive_only/'layouts.json').exists())
            noise_only=root/'noise-only'
            spec.write_text(json.dumps(dict(noise_model_data=manifest['noise_model_data'])))
            result=subprocess.run([str(renderer),'--noise-only',str(spec),str(noise_only)],
                                  capture_output=True,text=True,timeout=120)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(len(list(noise_only.glob('noise-module-*.svg'))),14)
            self.assertFalse((noise_only/'layouts.json').exists())
            self.assertFalse(list(noise_only.glob('archive-growth-*.svg')))


if __name__=='__main__':unittest.main()
