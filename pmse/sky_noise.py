"""Cached paper-model sky convolution and inexpensive rolling noise fits.

Numerical work only: figures are rendered by the Rust rustplotlib executable.
Heavy sky-map work belongs on Revontuli, never the acquisition computer.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import time
import h5py
import numpy as np

SIDEREAL_DAY_S=86164.0905


def build_cache(path,day,n_az=180,n_el=60):
    from astropy.time import Time
    import pansy_config as pc
    from paper_plot_noise import sky_temperature_matrix
    import importlib.metadata
    start=datetime.strptime(day,'%Y-%m-%d').replace(tzinfo=timezone.utc).timestamp()
    times=start+np.linspace(0,SIDEREAL_DAY_S,289)
    temperature=sky_temperature_matrix(times,gain_model='rx',n_az=n_az,n_el=n_el,
        min_elevation_deg=.5,rx_channel=None,freq_mhz=pc.freq/1e6,
        include_element_pattern=True,element_pattern_blend=.1)
    lst=np.asarray(Time(times,format='unix').sidereal_time('apparent',longitude=pc.lon).hour)
    order=np.argsort(lst); path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix('.h5.partial')
    with h5py.File(temporary,'w') as h:
        h['times_s']=times;h['lst_hour']=lst[order];h['sky_temperature_k']=temperature[order]
        h.attrs.update(model='PyGDSM GSM2008',pygdsm_version=importlib.metadata.version('pygdsm'),
            frequency_hz=pc.freq,latitude_deg=pc.lat,longitude_deg=pc.lon,
            receive_pattern='one-way single-module; 90% canonical crossed-Yagi + 10% isotropic element',
            element_pattern_blend=.1,n_az=n_az,n_el=n_el,reference_day=day,
            created_unix_s=time.time(),generator='pmse/sky_noise.py',
            geometry_sha256=hashlib.sha256(Path('cfg/antpos.csv').read_bytes()).hexdigest())
    temporary.replace(path)


def interpolate_cache(path,times_s):
    from astropy.time import Time
    with h5py.File(path,'r') as h:
        lst=Time(times_s,format='unix').sidereal_time('apparent',longitude=float(h.attrs['longitude_deg'])).hour
        x=h['lst_hour'][()];y=h['sky_temperature_k'][()]
        return np.column_stack([np.interp(lst,x,y[:,beam],period=24) for beam in range(5)])


def fit_module(tsky,power,times_s):
    # Same bounded relative-power objective used for the paper. No figure call.
    from fit_cut_noise_pygdsm import fit_receiver_temperature
    x=np.asarray(tsky).ravel();y=np.asarray(power).ravel()
    valid=np.isfinite(x)&np.isfinite(y)&(y>0)
    span=np.ptp(times_s) if len(times_s) else 0
    if len(times_s)<50 or span<12*3600 or not valid.any() or np.ptp(x[valid])/np.mean(x[valid])<.1:
        return dict(status='insufficient coverage or sky modulation',receiver_temp_k=None)
    use=valid.copy()
    for _ in range(3):
        fit=fit_receiver_temperature(x[use],y[use],np.ones(use.sum()),min_trec_k=273,max_trec_k=20000)
        residual=(y-float(fit['slope'])*(x+float(fit['receiver_temp_k'])))/y
        center=np.median(residual[use]);mad=1.4826*np.median(abs(residual[use]-center))
        use=valid&(abs(residual-center)<=max(6*mad,.10))
        if use.sum()<50:return dict(status='too few uncontaminated measurements',receiver_temp_k=None)
    trec=float(fit['receiver_temp_k']);gain=float(fit['slope'])
    rms=float(np.sqrt(np.mean(residual[use]**2)))
    status='fit at physical bound' if trec<273.01 or trec>19999.99 else 'provisional colored-noise calibration'
    return dict(status=status,receiver_temp_k=trec,power_per_kelvin=gain,fractional_rms=rms,
                used=int(use.sum()),rejected=int(valid.sum()-use.sum()))


def update(root,cache,output,hours=72):
    from .export_noise import noise_rows
    rows=[];latest_key=0;cutoff=(time.time()-hours*3600)*1e6
    for path in sorted(Path(root).glob('*T*/pmse@*.h5')):
        with h5py.File(path,'r') as h:
            for key,g in h.items():
                if not key.isdecimal() or int(key)<cutoff or int(g.get('complete',np.array(0))[()])!=1:continue
                fields=['noise_power','noise_power_standard_error','altitude_valid','complete','mode_start_sample','mode_stop_sample_exclusive']
                rows.append(noise_rows({k:g[k][()] for k in fields},int(key)))
                latest_key=max(latest_key,int(key))
    if not rows:raise ValueError('No complete noise records')
    rows.sort(key=lambda r:r['times_s'][0]);times=np.asarray([r['times_s'][0] for r in rows]);noise=np.stack([r['noise_power'] for r in rows])
    sky=interpolate_cache(cache,times)
    out=Path(output);fit_time=0;fits=None
    if out.exists():
        with h5py.File(out,'r') as h:
            fit_time=float(h.attrs.get('fit_unix_s',0))
            if time.time()-fit_time<3600:fits=json.loads(h.attrs['fits_json'])
    if fits is None:
        fits=[fit_module(sky,noise[:,:,ch],times) for ch in range(7)];fit_time=time.time()
    out.parent.mkdir(parents=True,exist_ok=True);tmp=out.with_suffix('.h5.partial')
    with h5py.File(tmp,'w') as h:
        h['times_s']=times;h['noise_power']=noise;h['sky_temperature_k']=sky
        h.attrs.update(model_cache=str(cache),generator='pmse/sky_noise.py',fits_json=json.dumps(fits),
            model='PyGDSM GSM2008; paper one-way receive pattern',
            fit_definition='P=G*(T_sky+T_receiver); relative-power objective; T_receiver>=273 K',
            caution='Guard-bin colored-noise response and interference may bias absolute receiver temperature',
            generated_unix_s=time.time(),fit_unix_s=fit_time,latest_mode_start_sample=latest_key)
    tmp.replace(out);print(json.dumps(fits),flush=True)
    return fits


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    b=sub.add_parser('cache');b.add_argument('--output',type=Path,required=True);b.add_argument('--day',required=True)
    f=sub.add_parser('fit');f.add_argument('--root',type=Path,required=True);f.add_argument('--cache',type=Path,required=True)
    f.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.command=='cache':build_cache(a.output,a.day)
    else:update(a.root,a.cache,a.output)
