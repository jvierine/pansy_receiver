"""Five-beam HSV imaging on a common horizontal plane (height above radar).

Beam-specific 18-mode inversions are evaluated at r=sqrt(h*h+x*x+y*y).
Angular spectral power is a preliminary density proxy: squared inversion
magnitude, normalized over each beam domain, times excess spectral power/noise.
Cosine-taper overlap weights are display weights, not measured TX gains.
"""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import h5py
import numpy as np
from image_five_beams import C, F_RADAR
from render_five_beam_hsv import render


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--start',default='2025-01-30T04:00:00')
    p.add_argument('--end',default='2025-01-30T04:30:00')
    p.add_argument('--height-km',type=float,default=85.290954301)
    p.add_argument('--grid-size',type=int,default=101)
    p.add_argument('--modes',type=int,default=18)
    p.add_argument('--doppler-max',type=float,default=9.568)
    p.add_argument('--width-max',type=float,default=6.379)
    p.add_argument('--dynamic-range-db',type=float,default=30)
    p.add_argument('--output',type=Path,default=Path('figures/five_beam_hsv'))
    a=p.parse_args()
    start=datetime.fromisoformat(a.start).replace(tzinfo=timezone.utc).timestamp()*1e6
    end=datetime.fromisoformat(a.end).replace(tzinfo=timezone.utc).timestamp()*1e6
    records=[]
    for path in sorted(a.root.glob('*/xc@*.h5')):
        # File blocks are five minutes long; avoid opening other days/hours.
        stamp=int(path.stem.split('@')[1])*1e6
        if stamp+300e6 < start or stamp > end:
            continue
        with h5py.File(path) as h:
            records.extend((int(k),path,k) for k in h if k.isdigit() and start<=int(k)<end)
    records.sort()
    if not records:
        p.error('No metadata records in requested interval')
    groups={}
    for line in (a.reference/'cfg/antpos.csv').read_text().splitlines()[1:]:
        fields=line.split(',')
        groups.setdefault(fields[1],[]).append([float(x) for x in fields[4:7]])
    connections=[line.split(',')[0] for line in (a.reference/'cfg/connections.txt').read_text().splitlines()]
    positions=np.array([np.mean(groups[x],axis=0) for x in connections if x!='RFTX'])
    phases=[]
    hashes=[]
    for beam in range(5):
        path=a.reference/f'data/phases{beam}.h5'
        with h5py.File(path) as h:
            phases.append(h['phasecal'][()])
        hashes.append(hashlib.sha256(path.read_bytes()).hexdigest())
    phases=np.asarray(phases)
    hgt=a.height_km
    extent=hgt*np.tan(np.deg2rad(20))
    axis=np.linspace(-extent,extent,a.grid_size)
    xx,yy=np.meshgrid(axis,axis)
    rr=np.sqrt(hgt**2+xx**2+yy**2)
    directions=np.column_stack([xx.ravel()/rr.ravel(),yy.ravel()/rr.ravel(),-hgt/rr.ravel()])
    az=np.deg2rad([0,0,90,180,270])
    el=np.deg2rad([90,80,80,80,80])
    centers=np.column_stack([np.cos(el)*np.sin(az),np.cos(el)*np.cos(az),-np.sin(el)])
    angles=np.arccos(np.clip(directions@centers.T,-1,1))
    radius=np.deg2rad(9)
    weights=np.where(angles<radius,np.cos(np.pi*angles/(2*radius))**2,0).T
    masks=weights>0
    weight_sum=weights.sum(axis=0)
    cover=weight_sum>0
    inverses=[]
    singular=[]
    with h5py.File(records[0][1]) as h:
        pairs=h[records[0][2]]['ch_pairs'][()]
    assert np.array_equal(pairs[:7],np.column_stack([np.arange(7)]*2))
    baseline=positions[pairs[:,1]]-positions[pairs[:,0]]
    for mask in masks:
        A=np.exp(1j*2*np.pi*F_RADAR/C*(baseline@directions[mask].T))
        U,s,Vh=np.linalg.svd(A,full_matrices=False)
        keep=min(a.modes,len(s))
        inverses.append((Vh[:keep].conj().T/(s[:keep]+1))@U[:,:keep].conj().T)
        singular.append(s)
    cell_area=(axis[1]-axis[0])**2
    power_frames=[]
    mean_frames=[]
    width_frames=[]
    beam_moments=[]
    a.output.parent.mkdir(parents=True,exist_ok=True)
    partial=a.output.with_suffix('.h5.partial')
    sidecar=h5py.File(partial,'w')
    spectra=sidecar.create_group('spectra')
    spectra.attrs['axis_order']='time, covered_pixel, doppler'
    spectra.attrs['units']='relative excess spectral power/noise per km2 (proxy)'
    spectra.attrs['pixel_order']='C-order flattened north,east grid; see flat_pixel_index'
    spectral_datasets=[]
    for beam,mask in enumerate(masks):
        bg=spectra.create_group(f'beam{beam}')
        bg['flat_pixel_index']=np.flatnonzero(mask)
    combined=spectra.create_group('combined')
    combined['flat_pixel_index']=np.flatnonzero(cover)
    for frame_index,(stamp,path,key) in enumerate(records):
        with h5py.File(path) as h:
            g=h[key]
            assert np.array_equal(g['ch_pairs'][()],pairs)
            xc=g['xc_arr'][()].astype(np.complex128)
            ranges=np.arange(g['r0'][()],g['r1'][()],g['rdec'][()])*C/2e9
            f=np.fft.fftshift(np.fft.fftfreq(int(g['n_fft'][()]),
                 d=5*int(g['ipp'][()])/1e6))[int(g['f0'][()]):int(g['f1'][()])]
        velocity=f*C/(2*F_RADAR)
        if frame_index==0:
            sidecar['velocity_mps']=velocity
            sidecar['doppler_frequency_hz']=f
            sidecar['source_range_km']=ranges
            for beam,mask in enumerate(masks):
                spectral_datasets.append(spectra[f'beam{beam}'].create_dataset(
                    'power',shape=(len(records),mask.sum(),len(f)),dtype='f4',
                    chunks=(1,min(256,mask.sum()),len(f)),compression='gzip',shuffle=True))
            combined_ds=combined.create_dataset('power',shape=(len(records),cover.sum(),len(f)),
                dtype='f4',chunks=(1,min(256,cover.sum()),len(f)),compression='gzip',shuffle=True)
        else:
            assert np.array_equal(sidecar['velocity_mps'][()],velocity)
            assert np.array_equal(sidecar['source_range_km'][()],ranges)
        if not (rr[cover.reshape(rr.shape)].min()>=ranges.min() and
                rr[cover.reshape(rr.shape)].max()<=ranges.max()):
            raise ValueError('Altitude plane extends outside archived range gates')
        power=abs(xc[:7])
        outer=abs(f)>=.75*max(abs(f))
        noise=np.median(power[:,:,outer,:],axis=(2,3))
        if not np.all(noise>0):
            raise ValueError('Nonpositive noise estimate')
        blend=np.zeros((directions.shape[0],len(f)))
        this_beam=[]
        for beam,mask in enumerate(masks):
            factor=np.exp(1j*(phases[beam,pairs[:,0]]-phases[beam,pairs[:,1]]))
            signal=np.maximum(power[:,beam]-noise[:,beam,None,None],noise[:,beam,None,None])
            coherence=xc[:,beam]*factor[:,None,None]/np.sqrt(signal[pairs[:,0]]*signal[pairs[:,1]])
            recon=abs(inverses[beam]@coherence.reshape(len(pairs),-1))**2
            recon=recon.reshape(mask.sum(),len(f),len(ranges))
            norm=recon.sum(axis=0)*cell_area
            excess=np.maximum(power[:,beam].sum(axis=0)/noise[:,beam].sum()-1,0)
            recon=np.divide(recon,norm[None],out=np.zeros_like(recon),where=norm[None]>0)*excess[None]
            target=rr.ravel()[mask]
            hi=np.clip(np.searchsorted(ranges,target),1,len(ranges)-1)
            lo=hi-1
            frac=(target-ranges[lo])/(ranges[hi]-ranges[lo])
            rows=np.arange(mask.sum())
            plane=(1-frac[:,None])*recon[rows,:,lo]+frac[:,None]*recon[rows,:,hi]
            spectral_datasets[beam][frame_index]=plane
            blend[mask]+=plane*weights[beam,mask,None]
            bm=np.full((3,directions.shape[0]),np.nan)
            total=plane.sum(axis=1)
            mean=np.divide(plane@velocity,total,out=np.zeros_like(total),where=total>0)
            variance=np.divide(plane@(velocity**2),total,out=np.zeros_like(total),where=total>0)-mean**2
            bm[:,mask]=np.array([total,mean,np.sqrt(np.maximum(variance,0))])
            this_beam.append(bm.reshape(3,*xx.shape))
        blend=np.divide(blend,weight_sum[:,None],out=np.zeros_like(blend),where=weight_sum[:,None]>0)
        combined_ds[frame_index]=blend[cover]
        total=blend.sum(axis=1)
        mean=np.divide(blend@velocity,total,out=np.zeros_like(total),where=total>0)
        variance=np.divide(blend@(velocity**2),total,out=np.zeros_like(total),where=total>0)-mean**2
        assert np.isfinite(total).all() and np.isfinite(mean).all()
        power_frames.append(total.reshape(xx.shape))
        mean_frames.append(mean.reshape(xx.shape))
        width_frames.append(np.sqrt(np.maximum(variance,0)).reshape(xx.shape))
        beam_moments.append(this_beam)
        print(datetime.fromtimestamp(stamp/1e6,timezone.utc).isoformat(),flush=True)
        sidecar.flush()
    sidecar.close()
    intensity=np.asarray(power_frames)
    mean=np.asarray(mean_frames)
    width=np.asarray(width_frames)
    reference=intensity.max()
    db=10*np.log10(np.maximum(intensity/reference,1e-12))
    with h5py.File(partial,'a') as h:
        for key,value in dict(time_unix_us=[r[0] for r in records],east_km=axis,north_km=axis,
             height_km=hgt,slant_range_km=rr,beam_weights=weights.reshape(5,*xx.shape),
             relative_power_db=db,power_density_proxy=intensity,centroid_mps=mean,width_mps=width,
             beam_moments=beam_moments,phasecal_rad=phases,singular_values=singular,
             antenna_positions_m=positions,beam_az_deg=np.rad2deg(az),beam_el_deg=np.rad2deg(el)).items():
            h.create_dataset(key,data=value,compression='gzip' if np.ndim(value)>0 else None)
        h.attrs.update(generator='pmse/animate_five_beams.py',convention='positive meteor signs',
            calibration_sha256=hashes,retained_modes=a.modes,doppler_limit_mps=a.doppler_max,
            width_limit_mps=a.width_max,dynamic_range_db=a.dynamic_range_db,
            spatial_blend='cosine squared taper to 9 deg from each TX center; not TX gain',
            power_definition='normalized squared TSVD magnitude times excess spectral power/noise per km2',
            temporal_interpolation='none; one frame per recorded integration',
            source_files=[str(r[1]) for r in records],source_groups=[r[2] for r in records])
        h['coverage_mask']=cover.reshape(xx.shape)
        h['direction_uvw']=directions.reshape(*xx.shape,3)
        h['ch_pairs']=pairs
        h.attrs.update(schema_version='pansy-five-beam-hsv-v1',geometry='constant_altitude_plane',
            height_reference='above radar; flat local plane, not sea-level altitude',
            beam_moments_axis_order='time,beam,moment,north,east',
            beam_moments_names='power_density_proxy,centroid_mps,width_mps',
            power_reference=reference,radar_frequency_hz=F_RADAR,speed_of_light_mps=C)
    partial.replace(a.output.with_suffix('.h5'))
    render(a.output.with_suffix('.h5'),a.output,doppler_max=a.doppler_max,
           width_max=a.width_max,db_min=-a.dynamic_range_db)
    print('Saved',a.output,'frames',len(records),flush=True)


if __name__=='__main__':
    main()
