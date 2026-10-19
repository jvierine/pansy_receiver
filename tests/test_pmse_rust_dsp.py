"""Rust/reference parity; build pmse/processor with cargo before running."""
import unittest
from pathlib import Path
import numpy as np
from pmse.sparse_record import Config, axes, decode_pulses, coherent_spectrum
from pmse.rust_dsp import RustDSP


class RustParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Deliberately not skipped: the Rust release gate must fail if unbuilt.
        cls.dsp=RustDSP()

    def test_random_measured_tx_selected_lags_and_power(self):
        cfg=Config(); rng=np.random.default_rng(733)
        x=(rng.normal(size=(80,1600))+1j*rng.normal(size=(80,1600))).astype('c8')
        tx=(rng.normal(size=x.shape)+1j*rng.normal(size=x.shape)).astype('c8')
        ranges=axes(cfg)[0]
        reference,power=decode_pulses(x,tx,cfg)
        rust,rust_power=self.dsp.decode_selected(x,tx,cfg,ranges)
        np.testing.assert_allclose(rust,reference[:,ranges],rtol=2e-5,atol=2e-6)
        np.testing.assert_allclose(rust_power,power,rtol=1e-12)
        scaled,_=self.dsp.decode_selected(x,3*tx,cfg,ranges)
        np.testing.assert_allclose(rust,scaled,rtol=2e-5,atol=2e-6)

    def test_masked_padded_full_doppler_and_channel_phase(self):
        cfg=Config(); rng=np.random.default_rng(12)
        for count in (12,252,256):
            x=(rng.normal(size=(7,5,count,8))+1j*rng.normal(size=(7,5,count,8))).astype('c8')
            valid=np.ones(count,dtype=bool); valid[6:9]=False
            a=coherent_spectrum(x,cfg,valid)
            b=self.dsp.coherent_spectrum(x,cfg,valid)
            np.testing.assert_allclose(a,b,rtol=3e-5,atol=2e-6)

    def test_rejects_zero_template_and_invalid_shapes(self):
        cfg=Config(); x=np.ones((4,1600),dtype='c8')
        with self.assertRaises(ValueError): self.dsp.decode_selected(x,x*0,cfg,[500])
        with self.assertRaises(ValueError): self.dsp.decode_selected(x,x,cfg,[-1])
        with self.assertRaises(ValueError): self.dsp.decode_selected(x,x[:3],cfg,[500])


if __name__=='__main__': unittest.main()
