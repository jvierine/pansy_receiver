"""Execute the Python recipes printed in the UiT download memo."""
import contextlib
import io
import json
from pathlib import Path
import re
import tempfile
import unittest

import h5py
import numpy as np
from pmse.sparse_record import coherent_window


class DownloadMemoExamples(unittest.TestCase):
    def test_inverse_and_noise_examples_preserve_the_retained_signal(self):
        self.check_examples('rectangular')

    def test_historical_hann_examples_remain_readable_but_cannot_be_unwindowed(self):
        self.check_examples('hann')

    def check_examples(self,kind):
        text=(Path(__file__).resolve().parents[1]/'memos/013_sparse_spectra_download.tex').read_text()
        snippets=re.findall(r'\\begin\{lstlisting\}\[language=Python\]\n(.*?)\\end\{lstlisting\}',text,re.S)
        self.assertEqual(len(snippets),4)
        rng=np.random.default_rng(1994);K=256;F=128
        index=np.array([[0,10,0],[0,84,0],[1,64,1]])
        spectra=(rng.normal(size=(K,7,3))+1j*rng.normal(size=(K,7,3))).astype('complex64')
        channel_noise=np.linspace(1.,4.,7)
        window=coherent_window(F) if kind=='rectangular' else np.hanning(F).astype('float32')
        window/=np.linalg.norm(window)
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)/'extracted-day';directory=root/'2026-10-19T19-00-00'
            directory.mkdir(parents=True);path=directory/'pmse@1792438200.h5'
            with h5py.File(path,'w') as h:
                g=h.create_group('1792438490922442')
                for name,value in dict(complete=1,quality='ok',spectra=spectra,sparse_index=index,
                    doppler_frequency_hz=np.arange(F),slant_range_km=[85.,86.],
                    config_json=json.dumps(dict(sample_rate_hz=1e6,window_kind=kind)),effective_ipp_seconds=.016,
                    look_pair_time0_sample_by_beam=np.zeros((K,5),dtype='int64'),
                    scheduled_pulses_per_beam_per_look=np.full(K,256,dtype='int32'),
                    look_window=np.tile(window,(K,1)),
                    noise_power=np.broadcast_to(channel_noise[:,None,None],(7,5,2))).items():g[name]=value
                h.create_group('1792438490922443')['complete']=0
            source=snippets[0].replace('root = Path("extracted-day")',f'root = Path({str(root)!r})')
            scope={}
            with contextlib.redirect_stdout(io.StringIO()):exec(compile(source,'memo-inverse','exec'),scope)
            expected=np.zeros((K,7,F),dtype='complex64');expected[:,:,index[:2,1]]=spectra[:,:,:2]
            np.testing.assert_array_equal(scope['X'],expected)
            np.testing.assert_allclose(np.fft.fftshift(np.fft.fft(scope['q'],axis=2),axes=2),expected,atol=1e-6)
            exec(compile(snippets[2],'memo-noise-model','exec'),scope)
            with h5py.File(path,'r') as h:
                g=h['1792438490922442'];scope['g']=g
                if kind=='hann':
                    with self.assertRaisesRegex(ValueError,'historical Hann'):
                        exec(compile(snippets[1],'memo-pair-voltage','exec'),scope)
                else:
                    exec(compile(snippets[1],'memo-pair-voltage','exec'),scope)
                    np.testing.assert_allclose(scope['pair_voltage'],scope['q']*np.sqrt(F),rtol=1e-6)
                kept=np.arange(F)%3==0
                noise,full,retained=scope['synthetic_noise'](g,0,0,kept,seed=2026)
                again=scope['synthetic_noise'](g,0,0,kept,seed=2026)
                np.testing.assert_array_equal(noise,again[0])
                # Correct per-bin spectral power under the unit-L2 window.
                np.testing.assert_allclose(np.mean(abs(noise)**2,axis=(0,2)),channel_noise,rtol=.06)
                # Exact realization-level Parseval energy, including 1/F.
                np.testing.assert_allclose(np.sum(abs(retained)**2,axis=2),
                    np.sum(abs(noise[:,:,kept])**2,axis=2)/F,rtol=1e-12,atol=1e-12)
                np.testing.assert_allclose(np.mean(np.sum(abs(retained)**2,axis=2),axis=0),
                    channel_noise*kept.sum()/F,rtol=.08)
                if kind=='hann':self.assertTrue(np.all(full[:,:,[0,-1]]==0))
                else:self.assertTrue(np.all(abs(full[:,:,[0,-1]])>0))
                exec(compile(snippets[3],'memo-noise-filling','exec'),scope)
                back=np.fft.fftshift(np.fft.fft(scope['q_filled'],axis=2),axes=2)
                np.testing.assert_allclose(back[:,:,scope['kept']],expected[:,:,scope['kept']],rtol=1e-6,atol=1e-6)


if __name__=='__main__':unittest.main()
