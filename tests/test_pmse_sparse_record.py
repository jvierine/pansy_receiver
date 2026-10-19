"""Numerical contracts for the sparse PMSE archive and detection statistics."""
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pmse.sparse_record import (Config, C, axes, decode_pulses, coherent_spectrum,
                                sparsify, write_record, read_complete, contiguous_looks,mode_is_closed,pulse_indices,
                                recover_cycle_grid,valid_transmit_cycles,coherent_window)


class SparseRecordTests(unittest.TestCase):
    def setUp(self):
        self.cfg=replace(Config(),nfft=64,altitude_min_km=80,altitude_max_km=82)
        self.rng=np.random.default_rng(913)

    def noise(self,looks=32,correlated=False):
        nr=len(axes(self.cfg)[0])
        shape=(looks,1 if correlated else 7,5,self.cfg.nfft,nr)
        z=(self.rng.normal(size=shape)+1j*self.rng.normal(size=shape))/np.sqrt(2)
        return np.repeat(z,7,axis=1).astype('c8') if correlated else z.astype('c8')

    def test_matched_filter_tx_power_invariance_and_relative_phase(self):
        cfg=Config(); n=cfg.ipp_samples
        tx=self.rng.choice([-1,1],cfg.tx_samples).astype('c8')*np.exp(.4j)
        padded=np.zeros(n,dtype='c8'); padded[:cfg.tx_samples]=tx
        rx=np.zeros((7,n),dtype='c8'); phase=np.arange(7)*.3
        for ch in range(7):
            rx[ch,540:540+cfg.tx_samples]=tx*np.exp(1j*phase[ch])
        reference=np.repeat(padded[None],7,axis=0)
        decoded,power=decode_pulses(rx,reference,cfg)
        scaled,_=decode_pulses(rx,3*reference,cfg)
        np.testing.assert_allclose(decoded,scaled,atol=5e-6)
        np.testing.assert_allclose(np.angle(decoded[:,540]),phase,atol=2e-6)
        np.testing.assert_allclose(abs(decoded[:,540]),np.sqrt(cfg.tx_samples),rtol=2e-6)
        self.assertTrue(np.all(np.argmax(abs(decoded),axis=1)==540))

    def test_unit_noise_gain_and_linear_no_wrap(self):
        cfg=Config(); trials=500
        tx=(self.rng.normal(size=(trials,cfg.ipp_samples))+
            1j*self.rng.normal(size=(trials,cfg.ipp_samples))).astype('c8')
        rx=(self.rng.normal(size=tx.shape)+1j*self.rng.normal(size=tx.shape)).astype('c8')/np.sqrt(2)
        decoded,_=decode_pulses(rx,tx,cfg)
        self.assertAlmostEqual(float(np.mean(abs(decoded[:,400:800])**2)),1.,delta=.03)

    def test_velocity_and_altitude_axes(self):
        cfg=Config(); ri,r,h,valid,f,v=axes(cfg)
        np.testing.assert_allclose(f,2*cfg.radar_frequency_hz*v/C)
        self.assertTrue(np.all((h[valid]>=50)&(h[valid]<=120)))
        self.assertGreater(r[-1],120)
        self.assertGreater(max(v),190)

    def test_coherent_fft_detects_both_extreme_velocities(self):
        cfg=Config(); _,_,_,_,f,v=axes(cfg)
        for target in [-90,90]:
            index=np.argmin(abs(v-target))
            t=np.arange(cfg.nfft)*5*cfg.ipp_samples/cfg.sample_rate_hz
            decoded=np.exp(2j*np.pi*f[index]*t)[None,None,:,None]
            spec=coherent_spectrum(decoded,cfg)
            self.assertEqual(int(np.argmax(abs(spec[0,0,:,0]))),int(index))

    def test_signal_phase_survives_incoherent_averaging(self):
        z=self.noise(); _,_,_,valid,_,v=axes(self.cfg)
        beam=0; r=int(np.flatnonzero(valid[beam])[0]); f=int(np.argmin(abs(v-75)))
        phases=np.arange(7)*.25
        signal=20*np.exp(1j*phases)
        for look in range(len(z)):
            z[look,:,beam,f,r]=signal*(-1)**look
        record=sparsify(z,self.cfg)
        rows=record['sparse_index']
        hit=np.flatnonzero(np.all(rows==[beam,f,r],axis=1))
        self.assertEqual(len(hit),1)
        values=record['spectra'][:,:,hit[0]]
        np.testing.assert_allclose(abs(values),20,rtol=1e-6)
        self.assertLess(float(abs(values.mean(axis=0)).max()),1e-5)
        cross=(values[:,:,None]*values[:,None,:].conj()).mean(axis=0)
        np.testing.assert_allclose(np.angle(cross[1:,0]),phases[1:],atol=1e-6)
        np.testing.assert_allclose(record['mean_power'][:,hit[0]],(abs(values)**2).mean(axis=0),rtol=1e-6)

    def test_null_independent_and_fully_correlated_receivers(self):
        for correlated in [False,True]:
            for trial in range(8):
                record=sparsify(self.noise(correlated=correlated),self.cfg)
                self.assertEqual(len(record['sparse_index']),0)
                np.testing.assert_allclose(np.median(record['noise_power']),1,atol=.12)

    def test_insufficient_looks_retained_as_quality_record(self):
        z=self.noise(looks=4)
        z[:,:,0,32,:]+=30
        result=sparsify(z,self.cfg)
        self.assertEqual(result['quality'],'insufficient_looks')
        self.assertEqual(len(result['sparse_index']),0)

    def test_gap_does_not_create_fake_coherent_sequence(self):
        cfg=Config(); keys=np.r_[np.arange(40)*32000,5_000_000+np.arange(50)*32000]
        looks=list(contiguous_looks(keys,cfg))
        self.assertEqual([len(x) for x in looks],[40,50])

    def test_boundary_closure_handles_schedule_without_inter_mode_metadata(self):
        last=10_000_000
        self.assertTrue(mode_is_closed(last,last+20_000_000,{}))
        self.assertFalse(mode_is_closed(last,last+1_000_000,{}))
        self.assertFalse(mode_is_closed(last,last+20_000_000,{last+32000:{'id':1}}))
        self.assertTrue(mode_is_closed(last,last+20_000_000,{last+100000:{'id':0}}))

    def test_clock_tracking_does_not_fragment_looks_or_misalign_cycles(self):
        cfg=Config()
        keys=np.arange(128,dtype=np.int64)*32000+np.arange(128)//30
        groups=list(contiguous_looks(keys,cfg))
        self.assertEqual([len(g) for g in groups],[64,64])
        for group in groups:
            idx=pulse_indices(group,cfg)
            np.testing.assert_array_equal(idx[::20,0],group-group[0])
            self.assertEqual(int(idx[-1,-1])+1,int(group[-1]-group[0])+32000)

    def test_missing_metadata_cycles_are_recovered_without_voltage_interpolation(self):
        cfg=Config(); original=1_700_000_000_000_000+np.arange(128)*32000+np.arange(128)//30
        markers=np.delete(original,[4,5,8,20,25,50,101])
        markers=np.r_[markers,markers[-1]+28800]  # transition false marker
        keys,observed=recover_cycle_grid(markers,cfg)
        self.assertEqual(len(keys),128)
        self.assertEqual(int((~observed).sum()),7)
        np.testing.assert_allclose(keys-original,0,atol=1)

    def test_measured_tx_verification_rejects_wrong_codes_and_zero_pulses(self):
        cfg=Config()
        template=self.rng.choice([-1,1],(20,cfg.ipp_samples)).astype('c8')
        reference=np.tile(template,(3,1))
        reference[20:40]=np.roll(template,1,axis=0)
        reference[40]=0
        valid,similarity=valid_transmit_cycles(reference,template,cfg)
        np.testing.assert_array_equal(valid,[True,False,False])
        self.assertGreater(similarity[0],.99)

    def test_gaps_have_zero_weight_and_unit_noise_gain(self):
        cfg=Config(); valid=np.ones(cfg.nfft,dtype=bool); valid[48:52]=False
        window=coherent_window(cfg.nfft,valid)
        self.assertTrue(np.all(window[48:52]==0))
        self.assertAlmostEqual(float((window**2).sum()),1.,places=6)

    def test_native_digital_metadata_roundtrip_empty_and_signal(self):
        z=self.noise(); ri,_,_,valid,_,_=axes(self.cfg)
        for has_signal in [False,True]:
            if has_signal:
                z[:,:,0,32,np.flatnonzero(valid[0])[0]]+=20
            record=sparsify(z,self.cfg)
            with tempfile.TemporaryDirectory() as root:
                key=1792430659921323
                self.assertTrue(write_record(root,key,record))
                self.assertFalse(write_record(root,key,record))
                saved=read_complete(root,key)
                np.testing.assert_array_equal(saved['spectra'],record['spectra'])
                np.testing.assert_array_equal(saved['sparse_index'],record['sparse_index'])
                self.assertEqual(int(saved['complete']),1)


if __name__=='__main__':
    unittest.main()
