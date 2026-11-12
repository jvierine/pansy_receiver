import numpy as np
import pytest
from meteor_doppler import refine_doppler,doppler_outliers
from meteor_fit import robust_fit_meteor


def pulse(seed=1,frequency_hz=4137.2,noise=.15,tx_offset_hz=0.0):
    rng=np.random.default_rng(seed);n=132;t=np.arange(n)/1e6
    tx=rng.choice([-1.,1.],n).astype(complex)*np.exp(2j*np.pi*tx_offset_hz*t)
    tx[:4]=0;tx[-4:]=0
    amplitude=(.5+rng.random(7))*np.exp(1j*rng.uniform(-np.pi,np.pi,7))
    rx=noise*(rng.normal(size=(7,260))+1j*rng.normal(size=(7,260)))/np.sqrt(2)
    rx[:,40:40+n]+=amplitude[:,None]*tx[None,:]*np.exp(2j*np.pi*frequency_hz*t)
    conversion=299792458/(2*47e6);bin_hz=125000/512
    coarse=np.round((frequency_hz+2*tx_offset_hz)/bin_hz)*bin_hz*conversion
    return tx,rx,coarse,frequency_hz*conversion


@pytest.mark.parametrize('f',[-38123.4,-4137.2,0.0,4137.2,38123.4])
def test_off_bin_complex_voltage_frequency_and_sign(f):
    tx,rx,coarse,truth=pulse(frequency_hz=f,noise=1e-5)
    result=refine_doppler(tx,rx,40,coarse)
    assert result['valid']
    assert result['velocity_mps']==pytest.approx(truth,abs=.05)
    assert result['sigma_mps']>0


def test_corrects_measured_tx_phase_ramp_instead_of_doubling_it():
    tx,rx,coarse,truth=pulse(noise=1e-5,tx_offset_hz=150)
    result=refine_doppler(tx,rx,40,coarse)
    assert result['valid'];assert abs(result['velocity_mps']-truth)<.05
    assert abs(coarse-truth)>700


def test_precision_beats_fft_quantization_and_uncertainties_have_coverage():
    errors=[];old=[];standardized=[]
    for seed in range(80):
        tx,rx,coarse,truth=pulse(seed,frequency_hz=4000+seed*7.3,noise=.25)
        result=refine_doppler(tx,rx,40,coarse)
        assert result['valid'];errors.append(result['velocity_mps']-truth);old.append(coarse-truth)
        standardized.append(errors[-1]/result['sigma_mps'])
    assert np.sqrt(np.mean(np.square(errors)))<.6*np.sqrt(np.mean(np.square(old)))
    assert .5<np.std(standardized)<1.5
    assert np.mean(abs(np.asarray(standardized))<=2)>.85


def test_invalid_and_boundary_doppler_is_not_silently_coarse():
    tx,rx,coarse,_=pulse()
    for invalid in [refine_doppler(tx,rx,-1,coarse),refine_doppler(tx*0,rx,40,coarse),
                    refine_doppler(tx,rx,40,coarse+20000)]:
        assert not invalid['valid'];assert np.isnan(invalid['velocity_mps'])


def test_reference_noise_increases_common_error_floor():
    tx,rx,coarse,_=pulse(noise=.01)
    perfect=refine_doppler(tx,rx,40,coarse)
    noisy_tx=tx.copy();noisy_tx[:4]=.02+.01j;noisy_tx[-4:]=-.02-.01j
    noisy=refine_doppler(noisy_tx,rx,40,coarse)
    assert noisy['sigma_mps']>perfect['sigma_mps']


def track(seed=1):
    rng=np.random.default_rng(seed);t=1794440000+np.linspace(-.2,.2,80)
    p0=np.array([1.,2.,100.]);v=np.array([30.,-10.,-25.])
    points=p0+(t-t.mean())[:,None]*v+rng.normal(size=(80,3))*[.05,.05,.02]
    return t,points,p0,v


def test_joint_3d_rejection_handles_endpoint_and_clustered_outliers():
    t,points,p0,v=track();bad=np.r_[0:9,30:37,78:80]
    points[bad]+=[5,-8,4]
    fit=robust_fit_meteor(t,*points.T)
    assert not fit['inlier'][bad].any()
    assert fit['inlier'].sum()>=59
    assert np.linalg.norm(fit['v0']-v)<.15
    assert np.linalg.norm(fit['r0']-p0)<.03
    assert fit['model'].shape==points.shape
    assert fit['residual'].shape==points.shape


def test_clean_track_does_not_lose_most_observations():
    t,points,p0,v=track();fit=robust_fit_meteor(t,*points.T)
    assert fit['inlier'].sum()>=78
    assert np.linalg.norm(fit['v0']-v)<.15


def test_nonfinite_and_degenerate_trajectory():
    t,points,_,_=track();points[12,1]=np.nan
    fit=robust_fit_meteor(t,*points.T);assert not fit['inlier'][12]
    with pytest.raises(ValueError):robust_fit_meteor(t*0,*points.T)
    with pytest.raises(ValueError):robust_fit_meteor(t[:2],*points[:2].T)


def test_doppler_mask_is_independent_and_preserves_curvature():
    t=np.linspace(-.1,.1,80);true=-12000+500*t+40000*t*t
    measured=true+np.random.default_rng(3).normal(0,50,len(t));measured[32]+=5000
    keep,model=doppler_outliers(t,measured,np.full(len(t),50))
    assert not keep[32];assert keep.sum()>=77
    assert np.sqrt(np.mean((model[keep]-true[keep])**2))<30


def test_doppler_incoherence_cannot_inflate_its_own_acceptance_scale():
    t=np.linspace(-.1,.1,50)
    measured=np.random.default_rng(71).uniform(-20000,20000,len(t))
    keep,model=doppler_outliers(t,measured,np.full(len(t),100))
    assert not keep.any()
    assert np.isnan(model).all()


def test_doppler_consensus_rejects_cluster_and_high_leverage_endpoints():
    t=np.linspace(-.2,.2,80);true=-20000+200*t+10000*t*t
    measured=true+np.random.default_rng(2).normal(0,80,len(t))
    bad=np.r_[0:8,35:45,78:80];measured[bad]+=15000
    keep,model=doppler_outliers(t,measured,np.full(len(t),80))
    assert not keep[bad].any();assert keep.sum()>=58


def test_appending_v2_preserves_complete_legacy_rows_and_registers_fields(tmp_path,monkeypatch):
    import h5py
    import process_cut_meteor as processor
    from meteor_metadata import StrictDigitalMetadataReader
    root=tmp_path/'fit';root.mkdir()
    writer=processor.simple_fit_writer(str(root))
    monkeypatch.setattr(processor.pc,'simple_fit_metadata_dir',str(root))
    monkeypatch.setattr(processor,'SIMPLE_FIT_WRITE_LOCK',str(tmp_path/'lock'))
    legacy={k:np.array([1.,2.,3.]) for k in processor.LEGACY_SIMPLE_FIT_FIELDS}
    assert processor.write_simple_fit_metadata(writer,1_700_000_000_000_000,legacy)
    payload={**legacy,'simple_fit_version':2,'doppler_mps':np.array([30.,40.,50.]),'fit_inlier':np.array([True,False,True])}
    assert not processor.write_simple_fit_metadata(writer,1_700_000_000_000_000,payload),'Do not upgrade/overwrite an existing complete legacy fit'
    assert processor.write_simple_fit_metadata(writer,1_700_000_001_000_000,payload)
    reader=StrictDigitalMetadataReader(str(root))
    assert 'doppler_mps' in reader.get_fields()
    records=reader.read(*reader.get_bounds())
    assert 'simple_fit_version' not in records[1_700_000_000_000_000]
    assert records[1_700_000_001_000_000]['doppler_mps'].tolist()==[30.,40.,50.]
