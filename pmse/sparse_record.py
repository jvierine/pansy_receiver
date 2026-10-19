"""Sparse seven-channel PMSE/PMWE spectra; scientific provenance stays in metadata.

One Digital Metadata record per closed mesomode. Sparse cells share a common
mask across channels AND looks, retaining voltage phase and unbiased cross
products at every selected cell. Incoherent averaging is applied to power.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
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
    velocity_limit_mps: float = 90.0
    snr_threshold_db: float = 3.0
    false_alarm_per_record: float = 1e-3
    min_looks: int = 16
    beam_elevation_deg: tuple = (90., 80., 80., 80., 80.)

    def validate(self):
        if self.nfft < 8 or self.nfft % 4 or self.ipp_samples <= self.tx_samples:
            raise ValueError('nfft must be a multiple of four >=8; pulse must fit IPP')
        nyquist = self.sample_rate_hz/(10*self.ipp_samples)*C/(2*self.radar_frequency_hz)
        if not 0 < self.velocity_limit_mps < .9*nyquist:
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
    freq = fft.fftshift(fft.fftfreq(cfg.nfft,d=5*cfg.ipp_samples/cfg.sample_rate_hz))
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


def coherent_spectrum(decoded,cfg):
    """Transform [channel,beam,pulse,range]; no pulse presumming or notch removal."""
    count = decoded.shape[2]
    if count < 4 or count > cfg.nfft:
        raise ValueError('A coherent look needs 4..nfft pulses per beam')
    window = np.hanning(count).astype(np.float32)
    window /= np.linalg.norm(window)
    return fft.fftshift(fft.fft(decoded*window[None,None,:,None],n=cfg.nfft,
                               axis=2,workers=1),axes=2).astype(np.complex64)


def sparsify(spectra,cfg):
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
    archived_band = np.abs(velocity)<=cfg.velocity_limit_mps+bin_halfwidth
    guard = (~archived_band)&(np.abs(freq)<.95*max(abs(freq)))
    if guard.sum() < 16:
        raise ValueError('Insufficient independent noise guard bandwidth')
    median_bias = special.gammaincinv(k,.5)/k
    noise = np.median(mean[:,:,guard,:],axis=2)/median_bias
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
        schema_version='pansy_pmse_sparse_v1',config_json=json.dumps(asdict(cfg)),
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
    groups=np.split(keys,np.flatnonzero(np.diff(keys)!=20*cfg.ipp_samples)+1)
    cycles=cfg.nfft//4
    for group in groups:
        for start in range(0,len(group),cycles):
            yield group[start:start+cycles]


def process_mode(raw_root,tx_root,start,last,cfg):
    import digital_rf as drf
    reader=drf.DigitalRFReader(str(raw_root))
    tx=drf.DigitalMetadataReader(str(tx_root))
    records=tx.read(int(start),int(last))
    keys=sorted(int(k) for k,v in records.items() if int(v['id'])==1)
    ri,*_=axes(cfg)
    spectra=[]; look_bounds=[]; counts=[]; tx_mean=[]; skipped=[]
    for group in contiguous_looks(keys,cfg):
        pulse_count=4*len(group)
        try:
            count=20*cfg.ipp_samples*len(group)
            ref=reader.read_vector_1d(int(group[0]),count,'ch007').astype(np.complex64)
            ref=ref.reshape(-1,cfg.ipp_samples)
            decoded=np.empty((7,5,pulse_count,len(ri)),dtype=np.complex64)
            for channel in range(7):
                rx=reader.read_vector_1d(int(group[0]),count,f'ch{channel:03d}').astype(np.complex64)
                matched,pulse_power=decode_pulses(rx.reshape(-1,cfg.ipp_samples),ref,cfg)
                # pulse order: cycle, four repeated five-beam sweeps, beam
                decoded[channel]=matched[:,ri].reshape(-1,4,5,len(ri)).transpose(2,0,1,3).reshape(5,pulse_count,len(ri))
            spectra.append(coherent_spectrum(decoded,cfg))
            look_bounds.append([group[0],group[-1]+20*cfg.ipp_samples])
            counts.append(pulse_count)
            tx_mean.append([pulse_power.mean(),pulse_power.min(),pulse_power.max()])
        except (OSError,IOError,ValueError) as error:
            skipped.append((int(group[0]),str(error)))
    if len(spectra)<2:
        raise ValueError(f'Insufficient usable looks: {len(spectra)}; skipped={skipped[:3]}')
    out=sparsify(np.stack(spectra),cfg)
    out.update(instrument_snapshot())
    out.update(mode_start_sample=int(start),mode_stop_sample_exclusive=int(last)+20*cfg.ipp_samples,
        mode_duration_seconds=(int(last)+20*cfg.ipp_samples-int(start))/cfg.sample_rate_hz,
        look_bounds_sample=np.asarray(look_bounds,dtype=np.int64),
        pulses_per_beam_per_look=np.asarray(counts,dtype=np.int32),
        tx_mean_min_max_power=np.asarray(tx_mean),
        coverage_seconds=sum(np.diff(np.asarray(look_bounds),axis=1).ravel())/cfg.sample_rate_hz,
        skipped_looks_json=json.dumps(skipped),raw_source=str(raw_root),tx_source=str(tx_root),
        per_pulse_tx_normalization='reference / sqrt(sum(abs(reference[:132])**2))',
        interferometer_calibration='not applied; original receiver channel order preserved',
        radial_velocity_sign='algebraic f_D=2*f_radar*v/c; physical up/down sign requires verification',
        physical_height_convention='nominal beam-center height above radar; actual sidelobe angle requires imaging',
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).parent,text=True).strip())
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
    modes=mode.read(max(oldest,mb[1]-int(lookback_hours*3600e6)),mb[1])
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
    return sorted(unique.items())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw-root',type=Path,default=Path('/media/archive'))
    p.add_argument('--tx-root',type=Path,default=Path('/media/analysis/metadata/tx'))
    p.add_argument('--mode-root',type=Path,default=Path('/media/analysis/metadata/mesomode'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--latest',type=int,default=1)
    p.add_argument('--lookback-hours',type=float,default=3)
    p.add_argument('--loop',action='store_true')
    p.add_argument('--poll-seconds',type=float,default=60)
    p.add_argument('--start',type=int); p.add_argument('--last',type=int)
    a=p.parse_args(); cfg=Config()
    while True:
        choices=[(a.start,a.last)] if a.start is not None and a.last is not None else closed_modes(a.mode_root,a.tx_root,a.raw_root,a.lookback_hours)
        for start,last in choices[-a.latest:]:
            if a.output.exists():
                try:
                    read_complete(a.output,start)
                    continue
                except (IOError,OSError,KeyError,ValueError):
                    pass
            begun=time.monotonic()
            try:
                out=process_mode(a.raw_root,a.tx_root,start,last,cfg)
                out['processing_wall_seconds']=time.monotonic()-begun
                out['peak_rss_platform_units']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                write_record(a.output,start,out)
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
        time.sleep(max(a.poll_seconds,10))


if __name__=='__main__':
    main()
