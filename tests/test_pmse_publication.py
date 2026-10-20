"""Current native data must publish even when optional plot rendering fails."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from pmse.publish_rvd import export
from pmse.processing_status import parse_worker


class PublicationTests(unittest.TestCase):
    def test_plot_error_or_timeout_cannot_freeze_manifest_or_date_listing(self):
        for failure in [subprocess.CalledProcessError(1,['renderer'],stderr='Label would be cropped'),
                        subprocess.TimeoutExpired(['renderer'],120),OSError('renderer unavailable')]:
            with self.subTest(failure=type(failure).__name__),tempfile.TemporaryDirectory() as temporary:
                base=Path(temporary);root=base/'metadata'/'sparse';output=base/'web';key=1792454400000000
                directory=root/'2026-10-20T00-00-00';directory.mkdir(parents=True)
                # A historical day must remain discoverable even with zero local
                # plotting success, and without loading its voltage spectra.
                historical=root/'2026-10-14T00-00-00';historical.mkdir()
                with h5py.File(historical/'pmse@1791936000.h5','w'):pass
                with h5py.File(directory/'pmse@1792454400.h5','w') as h:
                    g=h.create_group(str(key))
                    fields=dict(complete=1,mode_start_sample=key,mode_stop_sample_exclusive=key+86_000_000,
                        mode_duration_seconds=86.,coverage_seconds=85.,n_coherent_looks=4,quality='ok',
                        schema_version='pansy_pmse_sparse16_v2',sparse_index=[[0,1,0]],
                        spectra=np.ones((4,7,1),dtype='complex64'),noise_power=np.ones((7,5,2)),
                        detector_snr_db=[4.],n_recovered_cycle_markers=0,null_power_relative_standard_deviation=.5,
                        nominal_altitude_km=np.full((5,2),85.),slant_range_km=[85.,86.],
                        radial_velocity_mps=[-99.,0.,99.],beam_azimuth_deg=[0.,0.,90.,180.,270.],
                        beam_elevation_deg=[90.,80.,80.,80.,80.])
                    for name,value in fields.items():g[name]=value
                status=dict(generated_unix_s=1.,backfill=dict(current_index=250,total_modes=2374))
                (root/'processing-status.json').write_text(json.dumps(status))
                renderer=base/'renderer';renderer.touch()
                with patch('pmse.publish_rvd.subprocess.run',side_effect=failure):
                    export(root,output,renderer=renderer)
                manifest=json.loads((output/'manifest.json').read_text())
                self.assertEqual(manifest['plot_generation_status'],'delayed')
                self.assertEqual(manifest['available_days'],['2026-10-14','2026-10-20'])
                self.assertEqual(manifest['latest_mode_us'],key+86_000_000)
                self.assertEqual(manifest['source_processing'],status)
                with h5py.File(output/Path(manifest['sidecar']).name) as h:
                    self.assertIn(str(key),h['frames'])

    def test_source_worker_status_is_small_and_contains_no_access_details(self):
        b=parse_worker(['Rectangular backfill 254/2374: 1791992495895760..1791992582199772',
                       '{"start":1791992277895729,"result":"created","selected":4,"seconds":86.9}'],True)
        self.assertEqual(b['current_index'],254)
        self.assertEqual(b['total_modes'],2374)
        self.assertEqual(b['last_written_start_us'],1791992277895729)
        live=parse_worker(['{"start":1792479598928320,"seconds":85.856012,"processing_seconds":29.3}'])
        self.assertEqual(live['last_written_stop_us'],1792479684784332)
        self.assertNotIn('seconds',live)
        self.assertTrue(parse_worker(['Rectangular backfill complete for indexed retained raw buffer'],True)['complete'])


if __name__=='__main__':unittest.main()
