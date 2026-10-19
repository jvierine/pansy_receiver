"""Uncensored per-module noise measurements for the existing sky-noise fitter.

The science archive contains a dense noise estimate even when zero signal cells
survive. Never derive system noise from the threshold-selected signal voltages.
"""
import argparse
from datetime import datetime,timezone
from pathlib import Path
import h5py
import numpy as np


def noise_rows(record,key):
    if int(record.get('complete',0))!=1:
        raise ValueError('Only complete records are usable')
    noise=np.asarray(record['noise_power'],dtype=float)
    se=np.asarray(record['noise_power_standard_error'],dtype=float)
    valid=np.asarray(record['altitude_valid'],dtype=bool)
    if noise.ndim!=3 or noise.shape[:2]!=(7,5) or se.shape!=noise.shape or valid.shape!=noise.shape[1:]:
        raise ValueError('Invalid dense noise axes')
    power=[]; errors=[]
    for beam in range(5):
        gate=valid[beam]&np.all(np.isfinite(noise[:,beam])&(noise[:,beam]>0),axis=0)
        if not gate.any(): raise ValueError('No valid dense noise range gates')
        power.append(np.median(noise[:,beam,gate],axis=1))
        # Do not claim sqrt(n_range) precision: adjacent decoded gates correlate.
        errors.append(np.median(se[:,beam,gate],axis=1))
    start=int(record.get('mode_start_sample',key))
    stop=int(record.get('mode_stop_sample_exclusive',key))
    return dict(times_s=np.full(5,(start+stop)/2e6),beam_id=np.arange(5,dtype=np.int16),
                noise_power=np.asarray(power),noise_power_standard_error=np.asarray(errors))


def export(root,day,output):
    start=int(datetime.strptime(day,'%Y-%m-%d').replace(tzinfo=timezone.utc).timestamp()*1e6)
    stop=start+86400_000000; rows=[]
    for path in sorted(Path(root).glob('*T*/*.h5')):
        with h5py.File(path,'r') as h:
            for key in sorted(h,key=lambda x:int(x) if x.isdecimal() else -1):
                if not key.isdecimal() or not start<=int(key)<stop: continue
                g=h[key]
                if int(g.get('complete',np.array(0))[()])!=1: continue
                fields=['noise_power','noise_power_standard_error','altitude_valid','complete',
                        'mode_start_sample','mode_stop_sample_exclusive']
                rows.append(noise_rows({k:g[k][()] for k in fields},int(key)))
    if not rows: raise ValueError('No complete sparse records on requested UTC day')
    output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
    with h5py.File(output,'w') as h:
        h.attrs.update(day=day,source=str(root),per_module=1,n_modules=7,guard_samples=0,
            power_estimator='median range of bias-corrected guard-bin mean periodograms',
            calibration_h5='',normalization='unit-energy measured TX; unit-L2 Hann; raw receiver voltage power',
            description='Five beam rows per full mode; dense noise retained independently of detection. '
                        'White-noise gain is one. Colored-noise/filter response needs validation before Kelvin-fit reuse.',
            temporal_sampling='mode midpoints, not individual pulses')
        for k in rows[0]: h.create_dataset(k,data=np.concatenate([r[k] for r in rows]),compression='gzip',shuffle=True)
    return sum(len(r['times_s']) for r in rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--metadata',type=Path,required=True); p.add_argument('--day',required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); print(export(a.metadata,a.day,a.output))
