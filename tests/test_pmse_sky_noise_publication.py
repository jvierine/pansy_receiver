"""Sky-noise arrays publish independently, without signal spectra or SVG builds."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from pmse.publish_sky_noise import fit_specification,publish
from pmse.sky_noise import update as update_fit


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

    def test_numeric_noise_publication_never_builds_or_transfers_plot_images(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);fit=root/'fit.h5';out=root/'plots';self.fit_file(fit)
            def numeric_transfer(command,**kwargs):self.assertEqual(command[0],'rsync')
            with patch('pmse.publish_sky_noise.update') as update,patch(
                    'pmse.publish_sky_noise.subprocess.run',side_effect=numeric_transfer) as run:
                status=publish(root/'native',root/'cache.h5',fit,out,'web:sky-noise/',renderer='/renderer')
            update.assert_called_once_with(root/'native',root/'cache.h5',fit)
            self.assertEqual(len(status['noise_model_status']),7)
            self.assertEqual(status['noise_model_observation_stop_unix_s'],7200.)
            self.assertEqual(status['noise_plot']['status'],'ready')
            self.assertEqual(status['noise_plot']['format'],'rustplotlib-series-data-v1')
            self.assertEqual(status['schema'],'pansy_sky_noise_v2')
            self.assertEqual(len(status['noise_model_data']['noise_power'][0][0]),7)
            self.assertEqual(status['noise_model_origin_unix_s'],0.)
            self.assertEqual(json.loads((out/'status.json').read_text()),status)
            self.assertEqual(run.call_args_list[-1].args[0][-2:], [str(out/'status.json'),'web:sky-noise/'])
            self.assertEqual(run.call_count,1)
            self.assertFalse(list(out.glob('*.svg')))
            self.assertEqual(status['noise_model_data']['times_h'],[0.,1.,2.])

    def test_failed_fit_preserves_previous_published_status(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);fit=root/'fit.h5';self.fit_file(fit)
            out=root/'plots';out.mkdir();(out/'status.json').write_text('{"previous":true}')
            with patch('pmse.publish_sky_noise.update',side_effect=RuntimeError('Model unavailable')),patch('pmse.publish_sky_noise.subprocess.run') as run:
                with self.assertRaises(RuntimeError):
                    publish(root,root/'cache.h5',fit,out,renderer='/renderer')
                run.assert_not_called()
            self.assertEqual(json.loads((out/'status.json').read_text()),{'previous':True})

    def test_invalid_model_does_not_claim_ready(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'fit.h5';self.fit_file(path)
            with h5py.File(path,'a') as h:h['noise_power'][0,0,0]=0
            with self.assertRaises(ValueError):fit_specification(path)

    def test_zenith_fits_use_only_beam_zero_and_keep_existing_rvd_fits(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);hour=root/'2026-10-30T12-00-00';hour.mkdir()
            with h5py.File(hour/'pmse@1793361600.h5','w') as h:
                g=h.create_group('1793361600000000')
                power=np.arange(1,8)[:,None,None]*np.ones((7,5,3))
                power[:,1:]*=1000
                for name,value in dict(complete=1,mode_start_sample=1793361600000000,
                    mode_stop_sample_exclusive=1793361680000000,noise_power=power,
                    noise_power_standard_error=power*.1,altitude_valid=np.ones((5,3),bool)).items():g[name]=value
            sky=np.array([[10000.,20000.,30000.,40000.,50000.]])
            with patch('pmse.sky_noise.time.time',return_value=1793361700.),patch(
                'pmse.sky_noise.interpolate_cache',return_value=sky),patch('pmse.sky_noise.fit_module',
                return_value=dict(status='test',receiver_temp_k=300,power_per_kelvin=2)) as fit:
                update_fit(root,root/'cache.h5',root/'fit.h5')
            self.assertEqual(fit.call_count,14)
            for ch,call in enumerate(fit.call_args_list[7:]):
                np.testing.assert_array_equal(call.args[0],sky[:,0])
                np.testing.assert_array_equal(call.args[1],[ch+1])
            data,status=fit_specification(root/'fit.h5')
            self.assertEqual(len(status['noise_zenith_model_status']),7)
            self.assertEqual(len(status['noise_model_status']),7)


if __name__=='__main__':unittest.main()
