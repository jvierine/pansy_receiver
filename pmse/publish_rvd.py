"""Export sparse archive arrays for the Rust/WASM viewer; HDF5 is the data record.

JSON frame packets are a small browser transport. Each immutable frame is
sent once; manifest and HDF5 plot sidecar are published atomically. No plots
are rendered here and no raw-voltage data are sent to the browser.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess
import time
import h5py
import numpy as np


def packet(group,key):
    index=group['sparse_index'][()]
    z=group['spectra'][()].astype(np.complex128)
    noise=group['noise_power'][()]
    if len(index):
        power=(abs(z)**2/noise[:,index[:,0],index[:,2]][None]).mean(axis=1)
        signal=np.maximum(power-1,0)
        dominance=signal.max(axis=0)/np.maximum(signal.sum(axis=0),1e-30)
        persistence=np.sum(power>=2,axis=0)
    else:
        dominance=np.empty(0); persistence=np.empty(0,dtype=int)
    return dict(key=str(key),start_us=int(group['mode_start_sample'][()]),
        stop_us=int(group['mode_stop_sample_exclusive'][()]),
        duration_s=float(group['mode_duration_seconds'][()]),
        coverage_s=float(group['coverage_seconds'][()]),looks=int(group['n_coherent_looks'][()]),
        quality=group['quality'][()].decode(),
        index=index.tolist(),snr_db=group['detector_snr_db'][()].astype(float).tolist(),
        look_dominance=dominance.tolist(),look_persistence=persistence.tolist(),
        noise_power=noise.astype(float).tolist(),
        recovered_cycles=int(group['n_recovered_cycle_markers'][()]),
        null_relative_std=float(group['null_power_relative_standard_deviation'][()]))


def atomic_json(path,value):
    temporary=path.with_suffix(path.suffix+'.partial')
    temporary.write_text(json.dumps(value,separators=(',',':'),allow_nan=False))
    temporary.replace(path)


def export(root,output,hours=48):
    output.mkdir(parents=True,exist_ok=True); (output/'frames').mkdir(exist_ok=True)
    begin=time.time()-hours*3600
    date_cut=datetime.fromtimestamp(begin-3600,timezone.utc).strftime('%Y-%m-%dT%H')
    packets=[]; geometry=None
    for directory in sorted(root.glob('20*')):
        if directory.name[:13]<date_cut:
            continue
        for file in sorted(directory.glob('pmse@*.h5')):
            try:
                with h5py.File(file,'r') as h:
                    for name in sorted(h):
                        if not name.isdigit() or int(name)/1e6<begin:
                            continue
                        g=h[name]
                        if 'complete' not in g or int(g['complete'][()])!=1:
                            continue
                        if geometry is None:
                            geometry=dict(altitude_km=g['nominal_altitude_km'][()].tolist(),
                                slant_range_km=g['slant_range_km'][()].tolist(),
                                velocity_mps=g['radial_velocity_mps'][()].tolist(),
                                beam_azimuth_deg=g['beam_azimuth_deg'][()].tolist(),
                                beam_elevation_deg=g['beam_elevation_deg'][()].tolist())
                        path=output/'frames'/f'{name}.json'
                        if path.exists():
                            item=json.loads(path.read_text())
                        else:
                            item=packet(g,int(name)); atomic_json(path,item)
                        packets.append(item)
            except (OSError,KeyError,ValueError) as error:
                print(f'Skipping incomplete transfer {file}: {error}',flush=True)
    packets.sort(key=lambda p:p['start_us'])
    if not packets:
        raise ValueError('No complete sparse records in viewer window')
    latest=packets[-1]
    identity=hashlib.sha256(','.join(p['key'] for p in packets).encode()).hexdigest()[:16]
    manifest=dict(schema='pansy_rvd_plot_v1',revision=identity,
        generated_utc=datetime.now(timezone.utc).isoformat(),
        first_mode_us=packets[0]['start_us'],latest_mode_us=latest['stop_us'],
        frame_count=len(packets),geometry=geometry,
        frames=[{k:item[k] for k in ('key','start_us','stop_us','duration_s','coverage_s','looks','quality')}
                for item in packets],
        snr_threshold_db=3,altitude_limits_km=[50,120],velocity_limits_mps=[-90.321448,90.321448],
        sidecar='../archive/plot-data/plot-data-'+datetime.now(timezone.utc).strftime('%Y-%m-%d')+'.h5',
        display_fill='nearest mesomode for inter-mode gaps <=10 min; observations unchanged',
        moment_definition='threshold-censored moments of retained excess spectral power',
        generator='pmse/publish_rvd.py')
    transfer={}
    metrics=root.parent.parent/'rvd-transfer.h5'
    if metrics.exists():
        with h5py.File(metrics,'r') as h:
            for g in h.get('runs',{}).values():
                transfer[str(g.attrs['leg'])]={k:float(g[k][()]) for k in g}
    from .archive_inventory import inventory
    transfer['archive']=inventory(root)
    transfer['archive'].update(destination='/mnt/shovel/rvd/pmse_sparse',web_retention_days=7)
    manifest['transfer']=transfer
    cache=root.parent.parent/'sky-model-gsm2008-47MHz.h5'
    fit_path=root.parent.parent/'rvd-noise-fit.h5'
    if cache.exists():
        if not fit_path.exists() or time.time()-fit_path.stat().st_mtime>3600:
            from .sky_noise import update
            update(root,cache,fit_path)
        with h5py.File(fit_path,'r') as h:
            times=h['times_s'][()];first_noise=float(times[0])
            manifest['noise_model_data']=dict(times_h=((times-first_noise)/3600).tolist(),
                time_tick_labels=[datetime.fromtimestamp(t,timezone.utc).strftime('%H:%M')
                                  for t in [times[0],(times[0]+times[-1])/2,times[-1]]],
                sky_temperature_k=h['sky_temperature_k'][()].tolist(),noise_power=h['noise_power'][()].tolist(),
                fits=json.loads(h.attrs['fits_json']))
            manifest['noise_model_status']=manifest['noise_model_data']['fits']
            manifest['noise_model_updated_unix_s']=float(h.attrs['generated_unix_s'])
    first=packets[0]['start_us']; last=latest['stop_us']
    manifest['time_tick_labels']=[datetime.fromtimestamp((first+(last-first)*q)/1e6,timezone.utc).strftime('%H:%M:%S')
                                  for q in (0,.5,1)]
    renderer=Path(__file__).parent/'plotting/target/release/pansy-pmse-plots'
    if renderer.exists():
        specification=output/'plot-spec.json'
        atomic_json(specification,manifest)
        subprocess.run([str(renderer),str(specification),str(output/'axes')],check=True,timeout=120)
        manifest['plot_layouts']='axes/layouts.json'
    sidecar=output/Path(manifest['sidecar']).name
    write_sidecar=not sidecar.exists()
    if not write_sidecar:
        with h5py.File(sidecar,'r') as h:write_sidecar=h.attrs.get('revision','')!=identity
    if write_sidecar:
        tmp=sidecar.with_suffix('.h5.partial')
        with h5py.File(tmp,'w') as h:
            h.attrs.update(schema=manifest['schema'],generator='pmse/publish_rvd.py',source=str(root),
                snr_definition='excess signal/noise',moment_censoring_db=3,revision=identity)
            for key,value in geometry.items(): h[key]=value
            fg=h.create_group('frames')
            for item in packets:
                g=fg.create_group(item['key'])
                for key,value in item.items():
                    if key=='key': continue
                    if key=='quality': g.attrs[key]=value
                    else:
                        array=np.asarray(value)
                        g.create_dataset(key,data=array,compression='gzip' if array.ndim and array.size else None)
        tmp.replace(sidecar)
    atomic_json(output/'manifest.json',manifest)
    print(f'Exported {len(packets)} modes, latest {latest["key"]}',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/mnt/data/juha/pansy/metadata/pmse_sparse'))
    p.add_argument('--output',type=Path,default=Path('/mnt/data/juha/pansy/rvd-web/data'))
    p.add_argument('--rsync-dest'); p.add_argument('--hours',type=float,default=48)
    p.add_argument('--plot-archive-dest',default='j@juha.no:/mnt/shovel/rvd/plot-data/')
    p.add_argument('--loop',action='store_true')
    a=p.parse_args()
    while True:
        try:
            export(a.root,a.output,a.hours)
            if a.rsync_dest:
                # Transfer data packets and sidecar first, manifest last.
                subprocess.run(['rsync','-rt','--exclude=manifest.json','--exclude=*.partial','--exclude=*.h5',
                                '--exclude=.latex-cache/',
                                str(a.output)+'/',a.rsync_dest],check=True)
                manifest=json.loads((a.output/'manifest.json').read_text())
                sidecar=a.output/Path(manifest['sidecar']).name
                subprocess.run(['rsync','-rt','--no-perms','--no-owner','--no-group',str(sidecar),a.plot_archive_dest],check=True)
                subprocess.run(['rsync','-rt',str(a.output/'manifest.json'),a.rsync_dest],check=True)
        except Exception:
            import traceback; traceback.print_exc()
            if not a.loop: raise
        if not a.loop: break
        time.sleep(60)


if __name__=='__main__':
    main()
