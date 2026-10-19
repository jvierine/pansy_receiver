import unittest
from dataclasses import replace
import numpy as np
from pmse.sparse_record import Config,axes,coherent_window,coherent_spectrum,complementary_spectrum,sparsify
from pmse.rust_dsp import RustDSP
from pmse.reconstruct_echoes import reconstruct


class ComplementaryTests(unittest.TestCase):
    def setUp(self):
        self.cfg=replace(Config(),nfft=128,complementary_pairs=True,velocity_limit_mps=None)
        self.dsp=RustDSP();self.rng=np.random.default_rng(7881)

    def test_full_nyquist_band_and_integration(self):
        f=axes(self.cfg)[4];v=axes(self.cfg)[5]
        self.assertAlmostEqual(f[0],-31.25)
        self.assertAlmostEqual(v[0],-99.66504587765958,places=5)
        self.assertEqual(len(f),128)
        self.assertAlmostEqual(v[1]-v[0],1.55726634,places=6)

    def test_rust_matches_pair_reference_and_unit_noise_gain(self):
        x=(self.rng.normal(size=(7,5,256,32))+1j*self.rng.normal(size=(7,5,256,32))).astype('c8')/np.sqrt(2)
        valid=np.ones(256,bool);valid[16:24]=False
        a=complementary_spectrum(x,self.cfg,valid);b=self.dsp.complementary_spectrum(x,self.cfg,valid)
        np.testing.assert_allclose(a,b,rtol=3e-5,atol=2e-6)
        self.assertAlmostEqual(float(np.mean(abs(a)**2)),1.,delta=.03)

    def test_both_high_doppler_signs_have_known_two_pulse_response(self):
        f=axes(self.cfg)[4];v=axes(self.cfg)[5]
        for target in [-95,95]:
            k=np.argmin(abs(v-target));t=np.arange(256)*.008
            tone=np.exp(2j*np.pi*f[k]*t)
            x=np.broadcast_to(tone[None,None,:,None],(7,5,256,1)).astype('c8')
            out=self.dsp.complementary_spectrum(x,self.cfg)
            self.assertEqual(np.argmax(abs(out[0,0,:,0])),k)
            expected=np.sqrt(2)*coherent_window(128).sum()*abs(np.cos(np.pi*f[k]*.008))
            self.assertAlmostEqual(float(abs(out[0,0,k,0])),expected,places=4)

    def test_stationary_complementary_range_sidelobes_cancel(self):
        import pansy_modes
        codes=np.repeat(np.asarray(pansy_modes.get_m_mode()['codes'][:2]),8,axis=1)
        response=np.array([np.correlate(c,c,'full')[127::8]/np.sqrt(128) for c in codes])
        family=np.tile([0,1,1,0],64)
        f=0
        pulse=response[family]*np.exp(2j*np.pi*f*np.arange(256)*.008)[:,None]
        x=np.broadcast_to(pulse[None,None],(7,5,256,pulse.shape[1])).astype('c8')
        out=self.dsp.complementary_spectrum(x,self.cfg)
        self.assertGreater(abs(out[0,0,64,0]),100)
        np.testing.assert_allclose(out[0,0,64,1:],0,atol=3e-5)

    def test_inverse_fft_recovers_retained_windowed_echo(self):
        n=128;time=np.arange(n);wave=np.exp(2j*np.pi*9*time/n)
        spec=np.fft.fftshift(np.fft.fft(wave))
        record=dict(sparse_index=np.array([[0,73,0]]),spectra=np.broadcast_to(spec[73],(2,7,1)),
            doppler_frequency_hz=np.zeros(n),slant_range_km=np.ones(1))
        echo=reconstruct(record,0)
        np.testing.assert_allclose(echo[0,0,:,0],wave,atol=2e-7)


if __name__=='__main__':unittest.main()
