"""Sparse seven-channel PMSE/PMWE spectra; scientific provenance stays in metadata.

One Digital Metadata record per closed mesomode. Sparse cells share a common
mask across channels AND looks, retaining voltage phase and unbiased cross
products at every selected cell. Incoherent averaging is applied to power.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time

import numpy as np
from scipy import fft, special, stats

C = 299792458.0


@dataclass(frozen=True)
class Config:
    sample_rate_hz: float = 1e6
    radar_frequency_hz: float = 47e6
    ipp_samples: int = 1600
    tx_samples: int = 132
    nfft: int = 256
    range_step_samples: int = 4
    altitude_min_km: float = 50.0
    altitude_max_km: float = 120.0
    velocity_limit_mps: float | None = 90.0
    complementary_pairs: bool = False
    snr_threshold_db: float = 3.0
    false_alarm_per_record: float = 1e-3
    min_looks: int = 16
    beam_elevation_deg: tuple = (90., 80., 80., 80., 80.)

    def validate(self):
        if self.nfft < 8 or self.nfft % 4 or self.ipp_samples <= self.tx_samples:
            raise ValueError('nfft must be a multiple of four >=8; pulse must fit IPP')
        nyquist = self.sample_rate_hz/(10*self.ipp_samples*(2 if self.complementary_pairs else 1))*C/(2*self.radar_frequency_hz)
        if self.velocity_limit_mps is not None and not 0 < self.velocity_limit_mps < nyquist:
            raise ValueError('Requested velocity exceeds safe five-beam Nyquist coverage')
        if not 0 < self.altitude_min_km < self.altitude_max_km or self.range_step_samples < 1:
            raise ValueError('Invalid altitude bounds or range spacing')
        if self.min_looks < 2 or not 0 < self.false_alarm_per_record < 1:
            raise ValueError('Invalid integration or false-alarm configuration')


def axes(cfg):
    cfg.validate()
    dr = C/(2*cfg.sample_rate_hz)/1e3
    top = cfg.altitude_max_km/min(np.sin(np.deg2rad(cfg.beam_elevation_deg)))
    ri = np.arange(int(np.ceil(cfg.altitude_min_km/dr)),int(np.floor(top/dr))+1,
                   cfg.range_step_samples,dtype=np.int32)
    if ri[-1]+cfg.tx_samples > cfg.ipp_samples:
        raise ValueError('Upper altitude extends beyond unambiguous decoded range')
    ranges = ri*dr
    height = np.sin(np.deg2rad(cfg.beam_elevation_deg))[:,None]*ranges
    valid = (height >= cfg.altitude_min_km)&(height <= cfg.altitude_max_km)
    freq = fft.fftshift(fft.fftfreq(cfg.nfft,d=5*cfg.ipp_samples*(2 if cfg.complementary_pairs else 1)/cfg.sample_rate_hz))
    velocity = freq*C/(2*cfg.radar_frequency_hz)
    return ri,ranges,height,valid,freq,velocity


def decode_pulses(rx,tx,cfg):
    """Linear matched filter, per-pulse measured TX normalized to unit energy.

    Equivalent to first setting TX mean pulse power to one and dividing the
    correlator by sqrt(tx_samples). This fixes white receiver-noise gain at one.
    TX complex phase is retained, removing the common transmitted pulse phase.
    """
    tx = np.asarray(tx,dtype=np.complex64)[...,:cfg.tx_samples]
    energy = np.sum(np.abs(tx.astype(np.complex128))**2,axis=-1)
    if not np.all(np.isfinite(energy)) or np.any(energy <= 0):
        raise ValueError('Nonfinite/zero transmit pulse; cannot decode this look')
    kernel = tx/np.sqrt(energy)[...,None]
    n = fft.next_fast_len(cfg.ipp_samples+cfg.tx_samples-1)
    decoded = fft.ifft(fft.fft(rx,n=n,axis=-1,workers=1)*
                       np.conj(fft.fft(kernel,n=n,axis=-1,workers=1)),axis=-1,workers=1)
    return decoded[...,:cfg.ipp_samples].astype(np.complex64),energy/cfg.tx_samples


def coherent_window(count,pulse_valid=None):
    window=np.hanning(count).astype(np.float32)
    if pulse_valid is not None:
        window*=np.asarray(pulse_valid,dtype=np.float32)
    norm=np.linalg.norm(window)
    if norm<=0:
        raise ValueError('No usable coherent window support')
    return window/norm


def coherent_spectrum(decoded,cfg,pulse_valid=None):
    """Transform [channel,beam,pulse,range]; no pulse presumming or notch removal."""
    count = decoded.shape[2]
    if count < 4 or count > cfg.nfft:
        raise ValueError('A coherent look needs 4..nfft pulses per beam')
    window = coherent_window(count,pulse_valid)
    return fft.fftshift(fft.fft(decoded*window[None,None,:,None],n=cfg.nfft,
                               axis=2,workers=1),axes=2).astype(np.complex64)


def complementary_spectrum(decoded,cfg,pulse_valid=None):
    """Add adjacent complementary decoded returns, /sqrt(2) for unit noise gain.

    Each pair is timestamped at its midpoint and sampled every 16 ms. The known
    two-pulse response is cos(pi*f*8 ms) in voltage; do not mistake its high-
    Doppler sensitivity loss or residual range sidelobes for an unbiased PSD.
    """
    if decoded.shape[2]%2:raise ValueError('Incomplete complementary pulse pair')
    valid=None if pulse_valid is None else np.asarray(pulse_valid).reshape(-1,2).all(axis=1)
    paired=(decoded[:,:,0::2,:]+decoded[:,:,1::2,:])/np.sqrt(np.float32(2))
    return coherent_spectrum(paired,cfg,valid)


def complementary_template_families(template):
    """Verify paired A/B families; measured TX is still the decoder kernel."""
    import pansy_modes
    codes=np.repeat(np.asarray(pansy_modes.get_m_mode()['codes'][:2]),8,axis=1)
    autocorrelation=sum(np.correlate(c,c,'full') for c in codes[:,::8])
    if np.any(np.delete(autocorrelation,len(autocorrelation)//2)):
        raise ValueError('Configured code families are not complementary')
    families=[];similarity=[]
    for pulse in template:
        scores=np.array([[abs(np.vdot(code,pulse[s:s+128]))/
            max(np.linalg.norm(code)*np.linalg.norm(pulse[s:s+128]),1e-30) for s in range(5)] for code in codes])
        family,shift=np.unravel_index(np.argmax(scores),scores.shape)
        families.append(family);similarity.append(scores[family,shift])
    family=np.asarray(families).reshape(4,5).T
    if min(similarity)<.7 or np.any(family[:,0::2]==family[:,1::2]):
        raise ValueError('Recorded TX does not verify as adjacent complementary A/B pairs')
    return family,np.asarray(similarity).reshape(4,5).T


def sparsify(spectra,cfg,noise_reference=None):
    """Select from ensemble power; never threshold channels/looks separately."""
    ri,ranges,height,valid,freq,velocity = axes(cfg)
    z = np.asarray(spectra,dtype=np.complex64)
    if z.ndim != 5 or z.shape[1:] != (7,5,cfg.nfft,len(ri)) or not np.isfinite(z).all():
        raise ValueError('Expected finite [look,7,5,nfft,range] complex spectra')
    k = len(z)
    if k < 2:
        raise ValueError('Need at least two coherent looks for variance')
    power = np.abs(z.astype(np.complex128))**2
    mean = power.mean(axis=0)
    variance = power.var(axis=0,ddof=1)
    # Use full-band noise guards OUTSIDE the saved +/-90 m/s interval. Median
    # correction is exact for a Gamma(k,1/k) white-noise averaged periodogram.
    bin_halfwidth = abs(velocity[1]-velocity[0])/2
    archived_band = np.ones(len(velocity),dtype=bool) if cfg.velocity_limit_mps is None else np.abs(velocity)<=cfg.velocity_limit_mps+bin_halfwidth
    if noise_reference is None:
        if cfg.complementary_pairs:raise ValueError('Full-band paired data require independent pre-pair noise guards')
        noise_mean,noise_freq,noise_velocity=mean,freq,velocity
        guard=(~archived_band)&(np.abs(freq)<.95*max(abs(freq)))
    else:
        noise_mean,noise_freq,noise_velocity=noise_reference
        nyquist=cfg.sample_rate_hz/(20*cfg.ipp_samples)*C/(2*cfg.radar_frequency_hz)
        guard=(abs(noise_velocity)>nyquist)&(abs(noise_freq)<.95*max(abs(noise_freq)))
    if guard.sum() < 16:
        raise ValueError('Insufficient independent noise guard bandwidth')
    median_bias = special.gammaincinv(k,.5)/k
    noise = np.median(noise_mean[:,:,guard,:],axis=2)/median_bias
    if not np.isfinite(noise).all() or np.any(noise <= 0):
        raise ValueError('Invalid decoded-noise estimate')
    # Hann adjacent spectral bins correlate. Half the guard count is a
    # conservative effective count for the median-uncertainty approximation.
    guard_neff = max(1.,guard.sum()/2.)
    noise_rel_se = 1.253314137/np.sqrt(k*guard_neff)
    normalized = power/noise[None,:,:,None,:]
    detection_looks = normalized.mean(axis=1)
    detection = detection_looks.mean(axis=0)
    detection_se = detection_looks.std(axis=0,ddof=1)/np.sqrt(k)
    signal_snr = np.maximum(detection-1,0)
    inside = archived_band
    tested = valid[:,None,:]&inside[None,:,None]
    ntests = int(tested.sum())
    # Conservative shape k rather than 7*k: receiver channels may share sky
    # noise. Finite median-noise uncertainty is approximated by the F denominator.
    noise_shape = 1/noise_rel_se**2
    pvalue = stats.f.sf(detection,2*k,2*noise_shape)
    selected = tested&(signal_snr>=10**(cfg.snr_threshold_db/10))
    selected &= pvalue<=cfg.false_alarm_per_record/max(ntests,1)
    if k < cfg.min_looks:
        selected[:] = False
    beam,fi,rj = np.nonzero(selected)
    values = z[:,:,beam,fi,rj]  # [look,channel,cell]
    out = dict(
        schema_version='pansy_pmse_sparse16_v2' if cfg.complementary_pairs else 'pansy_pmse_sparse_v1',config_json=json.dumps(asdict(cfg)),
        sparse_index=np.column_stack([beam,fi,rj]).astype(np.int32),
        spectra=values,mean_power=mean[:,beam,fi,rj].astype(np.float32),
        power_variance=variance[:,beam,fi,rj].astype(np.float32),
        power_standard_error=np.sqrt(variance[:,beam,fi,rj]/k).astype(np.float32),
        noise_power=noise.astype(np.float32),
        noise_power_standard_error=(noise*noise_rel_se).astype(np.float32),
        detector_snr_db=(10*np.log10(signal_snr[beam,fi,rj])).astype(np.float32),
        detector_standard_error=detection_se[beam,fi,rj].astype(np.float32),
        null_probability=pvalue[beam,fi,rj],n_coherent_looks=k,
        noise_gamma_shape_conservative=k,noise_guard_effective_bins=guard_neff,
        null_power_relative_standard_deviation=1/np.sqrt(k),
        null_power_degrees_of_freedom_conservative=2*k,
        n_tested_cells=ntests,altitude_valid=valid,range_sample=ri,
        slant_range_km=ranges,nominal_altitude_km=height,
        doppler_frequency_hz=freq,radial_velocity_mps=velocity,
        noise_guard_mask=guard,beam_azimuth_deg=np.array([0,0,90,180,270.]),
        beam_elevation_deg=np.array(cfg.beam_elevation_deg),
        rx_channels=np.array([f'ch{i:03d}'.encode() for i in range(7)]),
        high_velocity_candidate=(abs(velocity[fi])>=12.5),
        snr_definition='10log10((ensemble_total_power/noise)-1); whitened channel-power mean',
        spectra_axis_order='coherent_look,rx_channel,sparse_cell',
        sparse_index_columns='beam,full_fft_bin,range_axis_index',
        normalization='measured TX unit pulse energy; Hann window unit L2; unnormalized FFT',
        detection_statistics='F tail approximation; channel correlation conservatively bounded; Bonferroni per record',
        quality='insufficient_looks' if k<cfg.min_looks else 'ok',
        generator='pmse/sparse_record.py',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    if cfg.complementary_pairs:
        out.update(effective_ipp_seconds=10*cfg.ipp_samples/cfg.sample_rate_hz,
            noise_doppler_frequency_hz=noise_freq,noise_radial_velocity_mps=noise_velocity,
            noise_estimator_domain='pre-pair 8 ms full-band guard periodograms outside paired Nyquist; unit noise gain',
            pair_normalization='(decoded_first+decoded_second)/sqrt(2); unit white-noise gain',
            pair_phase_reference='midpoint of two pulses separated by 8 ms',
            pair_voltage_doppler_response=np.cos(np.pi*freq*5*cfg.ipp_samples/cfg.sample_rate_hz),
            normalization='measured TX unit pulse energy; complementary sum /sqrt(2); pair Hann unit L2; unnormalized FFT')
    return out


def write_record(root,key,record):
    """Digital Metadata writes serialized; complete marker is written last.

    Readers accept only complete=1. An incomplete group can be repaired by this
    writer after interruption; completed records are immutable and idempotent.
    """
    import digital_rf as drf
    import h5py
    root = Path(root); root.mkdir(parents=True,exist_ok=True)
    with (root/'.writer.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        sec = int(key)//1_000_000
        hour = datetime.fromtimestamp(sec//3600*3600,timezone.utc).strftime('%Y-%m-%dT%H-%M-%S')
        path = root/hour/f'pmse@{sec//600*600}.h5'
        if path.exists():
            with h5py.File(path,'a') as h:
                if str(key) in h:
                    if int(h[str(key)].get('complete',np.array(0))[()])==1:
                        return False
                    del h[str(key)]  # our interrupted, explicitly incomplete record
        writer = drf.DigitalMetadataWriter(str(root),3600,600,1_000_000,1,'pmse')
        if writer._fields is not None and set(writer._fields)!=set(record)|{'complete'}:
            raise ValueError('Digital Metadata schema changed: use a new versioned channel')
        writer.write([int(key)],[dict(record,complete=np.int8(0))])
        with h5py.File(path,'a') as h:
            g = h[str(key)]
            # Digital Metadata exposes no compression option. Repack only this
            # private incomplete group's arrays before publishing complete=1.
            for name in list(g):
                ds=g[name]
                if isinstance(ds,h5py.Dataset) and ds.ndim and ds.size>256:
                    value=ds[()]; del g[name]
                    g.create_dataset(name,data=value,compression='gzip',compression_opts=1,shuffle=True)
            h.flush(); g['complete'][...]=1; h.flush()
        return True


def read_complete(root,key):
    import digital_rf as drf
    data=drf.DigitalMetadataReader(str(root)).read(int(key),int(key))[int(key)]
    if int(data['complete'])!=1:
        raise ValueError('Incomplete sparse record')
    return data


def completed_keys(root,start,stop):
    import digital_rf as drf
    try:
        rows=drf.DigitalMetadataReader(str(root)).read(int(start),int(stop),'complete')
        return {int(k) for k,v in rows.items() if int(v)==1}
    except (IOError,OSError,ValueError):
        return set()


def channel_covariance(record):
    """Recover all 7x7 complex correlations without ever archiving 28 products."""
    z=np.asarray(record['spectra'],dtype=np.complex128)
    return np.einsum('lcn,ldn->ncd',z,z.conj(),optimize=True)/len(z)


def instrument_snapshot():
    """Embed geometry and calibration inputs; stored voltages remain uncalibrated."""
    import h5py
    repo=Path(__file__).resolve().parents[1]
    groups={}
    path=repo/'cfg/antpos.csv'
    for line in path.read_text().splitlines()[1:]:
        fields=line.split(',')
        groups.setdefault(fields[1],[]).append([float(x) for x in fields[4:7]])
    connections=[line.split(',')[0] for line in (repo/'cfg/connections.txt').read_text().splitlines()]
    names=[name for name in connections if name!='RFTX' and name in groups]
    positions=np.array([np.mean(groups[name],axis=0) for name in names])
    if positions.shape!=(7,3):
        raise ValueError('Instrument geometry does not map to seven receiver modules')
    phases=[]; hashes=[]
    for beam in range(5):
        cal=repo/f'data/phases{beam}.h5'
        if cal.exists():
            with h5py.File(cal) as h:
                phases.append(np.asarray(h['phasecal'][()],dtype=float))
            hashes.append(hashlib.sha256(cal.read_bytes()).hexdigest())
        else:
            phases.append(np.full(7,np.nan)); hashes.append('missing')
    return dict(antenna_positions_m=positions,module_names=np.array(names,dtype='S16'),
        phasecal_snapshot_rad=np.asarray(phases),phasecal_sha256=np.array(hashes,dtype='S64'),
        geometry_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        phasecal_applied=np.int8(0),
        phasecal_status='snapshot only; historical sign/validity requires independent verification')


def contiguous_looks(keys,cfg):
    keys=np.asarray(sorted(keys),dtype=np.int64)
    # TX clock tracking occasionally advances a cycle marker by one sample.
    # Keep it coherent but align each cycle to its measured absolute sample.
    # Actual missing cycles still split the coherent sequence.
    groups=np.split(keys,np.flatnonzero(abs(np.diff(keys)-20*cfg.ipp_samples)>2)+1)
    cycles=cfg.nfft//(2 if cfg.complementary_pairs else 4)
    for group in groups:
        for start in range(0,len(group),cycles):
            yield group[start:start+cycles]


def pulse_indices(group,cfg):
    return ((group-group[0])[:,None,None]+
            np.arange(20)[None,:,None]*cfg.ipp_samples+
            np.arange(cfg.ipp_samples)[None,None,:]).reshape(-1,cfg.ipp_samples)


def recover_cycle_grid(keys,cfg):
    """Interpolate clock *timing*, never received voltages, between TX markers.

    Every recovered cycle must subsequently pass measured-TX signature checks.
    Reject markers displaced by a substantial fraction of an IPP (transition
    detections can otherwise rotate beam assignments).
    """
    keys=np.asarray(sorted(set(keys)),dtype=np.int64)
    cycle=20*cfg.ipp_samples
    number=np.rint((keys-keys[0])/cycle).astype(np.int64)
    residual=keys-keys[0]-number*cycle
    good=abs(residual)<cfg.ipp_samples/4
    number,residual=number[good],residual[good]
    number,unique=np.unique(number,return_index=True); residual=residual[unique]
    grid=np.arange(number[-1]+1,dtype=np.int64)
    offset=np.rint(np.interp(grid,number,residual)).astype(np.int64)
    recovered=keys[0]+grid*cycle+offset
    observed=np.isin(grid,number)
    return recovered,observed


def valid_transmit_cycles(reference,template,cfg):
    """Verify 20-pulse code/beam order using a measured, power-normalized template."""
    pulses=np.asarray(reference[:,:cfg.tx_samples],dtype=np.complex128)
    template=np.asarray(template[:,:cfg.tx_samples],dtype=np.complex128)
    energy=np.sum(abs(pulses)**2,axis=1)
    te=np.sum(abs(template)**2,axis=1)
    normalized=np.divide(pulses,np.sqrt(energy)[:,None],out=np.zeros_like(pulses),
                         where=energy[:,None]>0).reshape(-1,20,cfg.tx_samples)
    target=np.divide(template,np.sqrt(te)[:,None],out=np.zeros_like(template),where=te[:,None]>0)
    # The marker tracker rounds fractional-sample clock drift. Permit small
    # pulse-envelope shifts for verification, while decoding each actual pulse.
    similarities=[]
    for shift in range(-2,3):
        moved=np.zeros_like(target)
        if shift>=0:
            moved[:,shift:]=target[:,:cfg.tx_samples-shift]
        else:
            moved[:,:shift]=target[:,-shift:]
        moved/=np.maximum(np.linalg.norm(moved,axis=1)[:,None],1e-30)
        similarities.append(abs(np.sum(normalized*moved.conj()[None],axis=2)))
    similarity=np.max(similarities,axis=0)
    valid=(similarity.mean(axis=1)>=.8)&(similarity.min(axis=1)>=.65)
    valid &= np.all(np.isfinite(energy.reshape(-1,20))&(energy.reshape(-1,20)>0),axis=1)
    return valid,similarity.mean(axis=1)


def process_mode(raw_root,tx_root,start,last,cfg,backend='python'):
    import digital_rf as drf
    if backend not in ('python','rust'):
        raise ValueError('Unknown numerical backend')
    if backend=='rust':
        from .rust_dsp import RustDSP
        dsp=RustDSP()
    reader=drf.DigitalRFReader(str(raw_root))
    tx=drf.DigitalMetadataReader(str(tx_root))
    records=tx.read(int(start),int(last))
    markers=sorted(int(k) for k,v in records.items() if int(v['id'])==1)
    keys,observed=recover_cycle_grid(markers,cfg)
    # Learn the pulse-code sequence from a measured central marker, away from
    # start/stop transitions. This is a verifier, never a decoding substitute.
    anchor=markers[len(markers)//2]
    template=reader.read_vector_1d(int(anchor),20*cfg.ipp_samples,'ch007').astype(np.complex64).reshape(20,-1)
    if cfg.complementary_pairs:families,family_similarity=complementary_template_families(template)
    ri,*_=axes(cfg)
    spectra=[]; look_bounds=[]; counts=[]; tx_mean=[]; skipped=[]; gains=[]; scheduled=[];windows=[]
    noise_sum=None;noise_cfg=replace(cfg,complementary_pairs=False,nfft=cfg.nfft*2,velocity_limit_mps=90.)
    verified=np.zeros(len(keys),dtype=np.int8)
    used=np.zeros(len(keys),dtype=np.int8)
    similarities=np.full(len(keys),np.nan,dtype=np.float32)
    for group in contiguous_looks(keys,cfg):
        pulse_count=4*len(group)
        try:
            count=int(group[-1]-group[0])+20*cfg.ipp_samples
            indices=pulse_indices(group,cfg)
            ref=reader.read_vector_1d(int(group[0]),count,'ch007').astype(np.complex64)
            ref=ref[indices]
            valid_cycle,similarity=valid_transmit_cycles(ref,template,cfg)
            group_index=np.searchsorted(keys,group)
            verified[group_index]=valid_cycle
            similarities[group_index]=similarity
            valid_pulse=np.repeat(valid_cycle,20)
            valid_per_beam=np.repeat(valid_cycle,4)
            if np.count_nonzero(valid_per_beam)<4:
                raise ValueError('Measured TX rejects all cycles in coherent look')
            decoded=np.zeros((7,5,pulse_count,len(ri)),dtype=np.complex64)
            for channel in range(7):
                rx=reader.read_vector_1d(int(group[0]),count,f'ch{channel:03d}').astype(np.complex64)
                if backend=='rust':
                    selected,pulse_power=dsp.decode_selected(rx[indices[valid_pulse]],ref[valid_pulse],cfg,ri)
                else:
                    matched,pulse_power=decode_pulses(rx[indices[valid_pulse]],ref[valid_pulse],cfg)
                    selected=matched[:,ri]
                full=np.zeros((len(indices),len(ri)),dtype=np.complex64)
                full[valid_pulse]=selected
                # pulse order: cycle, four repeated five-beam sweeps, beam
                decoded[channel]=full.reshape(-1,4,5,len(ri)).transpose(2,0,1,3).reshape(5,pulse_count,len(ri))
            if cfg.complementary_pairs:
                transform=dsp.complementary_spectrum if backend=='rust' else complementary_spectrum
                paired=transform(decoded,cfg,valid_per_beam)
                noise_transform=dsp.coherent_spectrum if backend=='rust' else coherent_spectrum
                raw_noise=noise_transform(decoded,noise_cfg,valid_per_beam)
                contribution=abs(raw_noise.astype(np.complex128))**2
                noise_sum=contribution if noise_sum is None else noise_sum+contribution
                spectra.append(paired)
                window=coherent_window(pulse_count//2,valid_per_beam.reshape(-1,2).all(axis=1))
                padded=np.zeros(cfg.nfft,dtype=np.float32);padded[:len(window)]=window;windows.append(padded)
            else:
                transform=dsp.coherent_spectrum if backend=='rust' else coherent_spectrum
                spectra.append(transform(decoded,cfg,valid_per_beam))
            used[group_index]=valid_cycle
            look_bounds.append([group[0],group[-1]+20*cfg.ipp_samples])
            counts.append(int(valid_per_beam.sum()))
            scheduled.append(pulse_count)
            gains.append(float(2*window.sum()**2) if cfg.complementary_pairs else float(coherent_window(pulse_count,valid_per_beam).sum()**2))
            tx_mean.append([pulse_power.mean(),pulse_power.min(),pulse_power.max()])
        except (OSError,IOError,ValueError) as error:
            skipped.append((int(group[0]),str(error)))
    if len(spectra)<2:
        raise ValueError(f'Insufficient usable looks: {len(spectra)}; skipped={skipped[:3]}')
    noise_reference=None
    if cfg.complementary_pairs:noise_reference=(noise_sum/len(spectra),axes(noise_cfg)[4],axes(noise_cfg)[5])
    out=sparsify(np.stack(spectra),cfg,noise_reference)
    if cfg.complementary_pairs:
        out.update(complementary_code_family=families,complementary_code_family_similarity=family_similarity,
            look_window=np.asarray(windows),pairs_per_beam_per_look=np.asarray(counts,dtype=np.int32)//2,
            look_pair_time0_sample_by_beam=np.asarray(look_bounds,dtype=np.int64)[:,0,None]+np.arange(5)[None,:]*cfg.ipp_samples+5*cfg.ipp_samples//2,
            reconstruction='zero-fill omitted Fourier bins, ifftshift, inverse FFT per channel/look/range; returns retained windowed complementary pair echoes, not original ADC')
    if backend=='rust':
        out['generator']='pmse/processor (Rust DSP); pmse/sparse_record.py (Digital RF I/O and detection statistics)'
        out['code_sha256']=hashlib.sha256(Path(__file__).read_bytes()+
            (Path(__file__).parent/'processor/src/lib.rs').read_bytes()).hexdigest()
    out.update(instrument_snapshot())
    out.update(mode_start_sample=int(start),mode_stop_sample_exclusive=int(last)+20*cfg.ipp_samples,
        mode_duration_seconds=(int(last)+20*cfg.ipp_samples-int(start))/cfg.sample_rate_hz,
        look_bounds_sample=np.asarray(look_bounds,dtype=np.int64),
        pulses_per_beam_per_look=np.asarray(counts,dtype=np.int32),
        tx_mean_min_max_power=np.asarray(tx_mean),
        coverage_seconds=sum(counts)*5*cfg.ipp_samples/cfg.sample_rate_hz,
        skipped_looks_json=json.dumps(skipped),raw_source=str(raw_root),tx_source=str(tx_root),
        per_pulse_tx_normalization='reference / sqrt(sum(abs(reference[:132])**2))',
        interferometer_calibration='not applied; original receiver channel order preserved',
        radial_velocity_sign='algebraic f_D=2*f_radar*v/c; physical up/down sign requires verification',
        physical_height_convention='nominal beam-center height above radar; actual sidelobe angle requires imaging',
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).parent,text=True).strip())
    out['coverage_fraction']=out['coverage_seconds']/out['mode_duration_seconds']
    out['n_tx_cycles']=len(keys)
    out['n_tx_metadata_markers']=len(markers)
    out['n_recovered_cycle_markers']=int(np.count_nonzero(~observed))
    out['tx_cycle_start_sample']=keys
    out['tx_cycle_marker_observed']=observed
    out['tx_cycle_similarity']=similarities
    out['tx_cycle_verified']=verified
    out['rx_cycle_used']=used
    out['tx_reference_template']=template[:,:cfg.tx_samples]
    out['total_pulses_per_beam']=int(sum(counts))
    out['look_coherent_seconds']=np.asarray(scheduled)*5*cfg.ipp_samples/cfg.sample_rate_hz
    out['scheduled_pulses_per_beam_per_look']=np.asarray(scheduled,dtype=np.int32)
    gains=np.asarray(gains)
    out['look_coherent_signal_power_gain']=gains
    out['look_equivalent_noise_bandwidth_hz']=cfg.sample_rate_hz/(5*cfg.ipp_samples*gains)
    if out['quality']=='ok' and out['coverage_fraction']<.999:
        out['quality']='partial_coverage'
    return out


def mode_is_closed(last,tx_latest,following_records,ipp_samples=1600):
    if tx_latest < last+5_000_000:
        return False
    meso=sorted(int(k) for k,v in following_records.items() if int(v['id'])==1)
    return not meso or meso[0]-last>20*ipp_samples+ipp_samples


def closed_modes(mode_root,tx_root,raw_root,lookback_hours=3):
    import digital_rf as drf
    tx=drf.DigitalMetadataReader(str(tx_root)); tb=tx.get_bounds()
    mode=drf.DigitalMetadataReader(str(mode_root)); mb=mode.get_bounds()
    raw=drf.DigitalRFReader(str(raw_root))
    oldest=max(raw.get_bounds(f'ch{i:03d}')[0] for i in range(8))
    beginning=oldest if lookback_hours<=0 else max(oldest,mb[1]-int(lookback_hours*3600e6))
    modes=mode.read(beginning,mb[1])
    unique={}
    for value in modes.values():
        start,last=int(value['start']),int(value['end'])
        if last<=start or start<oldest or last+20*1600>=tb[1]:
            continue
        after=tx.read(last+1,min(last+5_000_000,tb[1]))
        # Old boundary metadata can publish partial/open records. Do not use
        # these until a genuinely observed gap exceeds a full 20-pulse cycle.
        # TX metadata contains mesomode markers only during some schedules;
        # no rows in a five-second gap is meaningful once TX has advanced.
        if not mode_is_closed(last,tb[1],after):
            continue
        unique[start]=max(last,unique.get(start,last))
    # The legacy mode-boundary producer keeps a ~20 min lag. Find newly closed
    # modes directly in its TX-metadata tail for the realtime archive. A 0.5 s
    # split matches the current boundary producer and tolerates missed markers.
    tail_begin=max(oldest,mb[1]+1) if lookback_hours<=0 else max(oldest,mb[1]+1,tb[1]-int(lookback_hours*3600e6))
    tail=tx.read(tail_begin,tb[1],'id')
    tail_keys=np.array([int(k) for k,v in tail.items() if int(v)==1],dtype=np.int64)
    if len(tail_keys):
        for group in np.split(tail_keys,np.flatnonzero(np.diff(tail_keys)>500_000)+1):
            if tb[1]-group[-1]>=5_000_000 and len(group)>1:
                unique[int(group[0])]=int(group[-1])
    return sorted(unique.items())


def pending_modes(choices,done,limit=0):
    """One sequential queue: two freshest records then the oldest pending one."""
    remaining=[item for item in choices if item[0] not in done]
    if limit<=0:return remaining
    newest=remaining[-min(2,limit):][::-1]
    return (newest+[item for item in remaining if item not in newest])[:limit]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw-root',type=Path,default=Path('/media/archive'))
    p.add_argument('--tx-root',type=Path,default=Path('/media/analysis/metadata/tx'))
    p.add_argument('--mode-root',type=Path,default=Path('/media/analysis/metadata/mesomode'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--latest',type=int,default=1)
    p.add_argument('--backend',choices=('python','rust'),default='python',
                   help='Rust executes pulse decoding and Doppler FFTs; reference I/O/statistics are shared')
    p.add_argument('--legacy-8ms',action='store_true',help='Historical v1 only; default is full-band complementary pairs at 16 ms')
    p.add_argument('--reserve-free-gib',type=float,default=200)
    p.add_argument('--max-archive-gib',type=float,default=500)
    p.add_argument('--lookback-hours',type=float,default=3)
    p.add_argument('--max-per-cycle',type=int,default=0,help='Bound each scheduling batch; freshest two first, then backfill')
    p.add_argument('--loop',action='store_true')
    p.add_argument('--poll-seconds',type=float,default=60)
    p.add_argument('--start',type=int); p.add_argument('--last',type=int)
    a=p.parse_args(); cfg=Config() if a.legacy_8ms else replace(Config(),nfft=128,complementary_pairs=True,velocity_limit_mps=None)
    if a.latest<0 or a.reserve_free_gib<0 or a.max_archive_gib<=0:
        p.error('latest must be >=0 (zero=all); storage budgets must be positive')
    a.output.mkdir(parents=True,exist_ok=True)
    choices_cache=None;choices_checked=0.
    while True:
        free=shutil.disk_usage(a.output).free/2**30
        used=sum(path.stat().st_size for path in a.output.rglob('*.h5'))/2**30
        if free<a.reserve_free_gib or used>=a.max_archive_gib:
            message=f'Archive paused: free={free:.1f} GiB, archive={used:.1f} GiB; storage budget reached'
            print(message,flush=True)
            if not a.loop:
                raise RuntimeError(message)
            time.sleep(max(a.poll_seconds,10))
            continue
        if choices_cache is None or time.monotonic()-choices_checked>=60:
            choices_cache=[(a.start,a.last)] if a.start is not None and a.last is not None else closed_modes(a.mode_root,a.tx_root,a.raw_root,a.lookback_hours)
            choices_checked=time.monotonic()
        choices=choices_cache
        done=completed_keys(a.output,choices[0][0],choices[-1][0]) if choices else set()
        batch=pending_modes(choices[-a.latest:] if a.latest else choices,done,a.max_per_cycle)
        made_progress=False
        for start,last in batch:
            if start in done:
                continue
            if shutil.disk_usage(a.output).free/2**30<a.reserve_free_gib:
                break
            print(f'Processing closed mesomode {start}..{last}',flush=True)
            begun=time.monotonic()
            try:
                out=process_mode(a.raw_root,a.tx_root,start,last,cfg,a.backend)
                out['processing_wall_seconds']=time.monotonic()-begun
                import sys
                out['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
                write_record(a.output,start,out)
                made_progress=True
                print(json.dumps(dict(start=start,seconds=out['mode_duration_seconds'],looks=out['n_coherent_looks'],
                    selected=len(out['sparse_index']),processing_seconds=out['processing_wall_seconds'],
                    coverage_seconds=out['coverage_seconds'],quality=out['quality'])),flush=True)
            except Exception as error:
                import traceback
                traceback.print_exc()
                if not a.loop:
                    raise
        if not a.loop:
            break
        if made_progress:continue
        time.sleep(max(a.poll_seconds,10))


if __name__=='__main__':
    main()
