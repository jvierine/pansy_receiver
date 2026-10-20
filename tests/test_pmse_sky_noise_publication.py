"""Sky-noise diagnostics publish independently, without opening signal spectra."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from pmse.publish_sky_noise import fit_specification,publish


class SkyNoisePublicationTests(unittest.TestCase):
    def fit_file(self,path):
        with h5py.File(path,'w') as h:
            h['times_s']=[0.,3600.,7200.]
            h['sky_temperature_k']=np.ones((3,5))*10000
            h['noise_power']=np.ones((3,5,7))*100
            h.attrs.update(generated_unix_s=7201.,fit_unix_s=7200.,
                fits_json=json.dumps([dict(status='insufficient coverage or sky modulation',receiver_temp_k=None)]*7))

    @staticmethod
    def renderer(command,**kwargs):
        if command[0]=='rsync':return
        assert command[1]=='--noise-only'
        directory=Path(command[-1]);directory.mkdir(exist_ok=True)
        for ch in range(7):
            for width in (640,960):(directory/f'noise-module-{ch}-{width}.svg').write_text('<svg/>')

    def test_native_noise_only_publication_and_status_transferred_last(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);fit=root/'fit.h5';out=root/'plots';self.fit_file(fit)
            with patch('pmse.publish_sky_noise.update') as update,patch(
                    'pmse.publish_sky_noise.subprocess.run',side_effect=self.renderer) as run:
                status=publish(root/'native',root/'cache.h5',fit,out,'web:sky-noise/',renderer='/renderer')
            update.assert_called_once_with(root/'native',root/'cache.h5',fit)
            self.assertEqual(len(status['noise_model_status']),7)
            self.assertEqual(status['noise_model_observation_stop_unix_s'],7200.)
            self.assertEqual(status['noise_plot']['status'],'ready')
            self.assertEqual(len(status['noise_plot']['files']),7)
            self.assertEqual(json.loads((out/'status.json').read_text()),status)
            self.assertEqual(run.call_args_list[-1].args[0][-2:], [str(out/'status.json'),'web:sky-noise/'])
            self.assertIn('--exclude=*.json',run.call_args_list[-2].args[0])
            spec=json.loads((out/'plot-spec.json').read_text())
            self.assertEqual(set(spec),{'noise_model_data'})
            self.assertEqual(spec['noise_model_data']['times_h'],[0.,1.,2.])

    def test_failed_render_preserves_previous_published_status(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);fit=root/'fit.h5';self.fit_file(fit)
            out=root/'plots';out.mkdir();(out/'status.json').write_text('{"previous":true}')
            with patch('pmse.publish_sky_noise.update'),patch('pmse.publish_sky_noise.subprocess.run',
                    side_effect=subprocess.CalledProcessError(1,['renderer'])):
                with self.assertRaises(subprocess.CalledProcessError):
                    publish(root,root/'cache.h5',fit,out,renderer='/renderer')
            self.assertEqual(json.loads((out/'status.json').read_text()),{'previous':True})

    def test_invalid_model_does_not_claim_ready(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'fit.h5';self.fit_file(path)
            with h5py.File(path,'a') as h:h['noise_power'][0,0,0]=0
            with self.assertRaises(ValueError):fit_specification(path)


if __name__=='__main__':unittest.main()
