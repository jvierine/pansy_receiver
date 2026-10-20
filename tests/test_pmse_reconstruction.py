"""Reversibility boundaries of the retained complex spectrum, not raw ADC."""
import unittest
import numpy as np
from pmse.reconstruct_echoes import reconstruct
from pmse.sparse_record import coherent_window


class RetainedReconstructionTests(unittest.TestCase):
    def test_complex_multichannel_multilook_roundtrip_and_range_separation(self):
        rng=np.random.default_rng(6401)
        index=np.array([[0,3,0],[0,27,0],[0,64,1],[1,73,0],[4,127,1]])
        z=(rng.normal(size=(3,7,len(index)))+1j*rng.normal(size=(3,7,len(index)))).astype('complex64')
        record=dict(sparse_index=index,spectra=z,doppler_frequency_hz=np.zeros(128),slant_range_km=[85.,86.])
        for beam in range(5):
            q=reconstruct(record,beam)
            back=np.fft.fftshift(np.fft.fft(q,axis=2),axes=2)
            expected=np.zeros((3,7,128,2),dtype='complex64')
            cells=np.flatnonzero(index[:,0]==beam)
            expected[:,:,index[cells,1],index[cells,2]]=z[:,:,cells]
            np.testing.assert_allclose(back,expected,atol=5e-7,rtol=1e-6)

    def test_full_spectrum_recovers_windowed_pair_amplitudes(self):
        rng=np.random.default_rng(6402)
        pulse=(rng.normal(size=(2,7,256))+1j*rng.normal(size=(2,7,256))).astype('complex64')
        paired=(pulse[:,:,::2]+pulse[:,:,1::2])/np.sqrt(2.)
        expected=paired*coherent_window(128)[None,None,:]
        z=np.fft.fftshift(np.fft.fft(expected,axis=2),axes=2)
        record=dict(sparse_index=np.c_[np.zeros(128,int),np.arange(128),np.zeros(128,int)],
            spectra=z,doppler_frequency_hz=np.zeros(128),slant_range_km=[85.])
        np.testing.assert_allclose(reconstruct(record,0)[:,:,:,0],expected,atol=2e-7,rtol=2e-6)

    def test_omitted_bin_is_not_recoverable(self):
        full=np.zeros(128,complex);full[73]=2+3j
        alternate=full.copy();alternate[17]=9-4j
        record=dict(sparse_index=np.array([[0,73,0]]),spectra=np.full((2,7,1),full[73]),
            doppler_frequency_hz=np.zeros(128),slant_range_km=[85.])
        q=reconstruct(record,0)[0,0,:,0]
        np.testing.assert_allclose(q,np.fft.ifft(np.fft.ifftshift(full)),atol=1e-8)
        self.assertGreater(np.linalg.norm(q-np.fft.ifft(np.fft.ifftshift(alternate))),0.5)
        self.assertEqual(full[73],alternate[73])

    def test_window_zero_and_pair_sum_are_not_invertible(self):
        a=np.ones(128,complex);b=np.ones(128,complex)*2
        different_a=a+3j;different_b=b-3j
        np.testing.assert_array_equal(a+b,different_a+different_b)
        window=coherent_window(128)
        source=np.ones(128,complex);changed=source.copy();changed[[0,-1]]=123+45j
        np.testing.assert_array_equal(window*source,window*changed)


if __name__=='__main__':unittest.main()
