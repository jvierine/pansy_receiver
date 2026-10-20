"""Count native records without reading the (potentially large) spectra arrays."""
from pathlib import Path
from datetime import datetime,timezone
import h5py
import numpy as np

_cache={}


def inventory(root):
    files=list(Path(root).glob('*T*/pmse@*.h5')); count=0; first=None; stop=None; total=0;hour_bytes={}
    retained_samples=0;tested_samples=0;retention_records=0;dense_samples=0;dense_records=0
    hour_retention={};days={}
    for path in files:
        st=path.stat();total+=st.st_size;state=(st.st_mtime_ns,st.st_size)
        hour=int(path.stem.split('@')[1])//3600*3600;hour_bytes[hour]=hour_bytes.get(hour,0)+st.st_size
        entry=_cache.get(str(path))
        if entry is None or entry[0]!=state:
            keys=[];ends=[];kept=0;tested=0;measured=0;dense=0;dense_measured=0;per_hour={}
            try:
                with h5py.File(path,'r') as h:
                    for key,g in h.items():
                        if not key.isdecimal() or int(g.get('complete',np.array(0))[()])!=1:continue
                        keys.append(int(key));ends.append(int(g['mode_stop_sample_exclusive'][()]))
                        bucket=per_hour.setdefault(int(key)//3_600_000_000*3600,[0,0,0,0])
                        bucket[2]+=1
                        if 'spectra' in g and 'n_tested_cells' in g:
                            # Inspect dimensions only; never load complex spectra.
                            looks,channels,cells=g['spectra'].shape
                            kept+=looks*channels*cells
                            bucket[0]+=looks*channels*cells
                            tested+=looks*channels*int(g['n_tested_cells'][()])
                            measured+=1
                            if all(name in g for name in ('noise_power','radial_velocity_mps','slant_range_km')):
                                beams=g['noise_power'].shape[1]
                                possible=looks*channels*beams*g['radial_velocity_mps'].size*g['slant_range_km'].size
                                dense+=possible;bucket[1]+=possible;bucket[3]+=1
                                dense_measured+=1
            except OSError:continue
            entry=(state,len(keys),min(keys,default=None),max(ends,default=None),kept,tested,measured,dense,dense_measured,per_hour)
            _cache[str(path)]=entry
        count+=entry[1]
        day=path.parent.name[:10]
        d=days.setdefault(day,dict(records=0,bytes=0,first_mode_us=None,latest_mode_us=None))
        d['records']+=entry[1];d['bytes']+=st.st_size
        if entry[2] is not None:
            d['first_mode_us']=min(d['first_mode_us'] or entry[2],entry[2])
            d['latest_mode_us']=max(d['latest_mode_us'] or entry[3],entry[3])
        retained_samples+=entry[4];tested_samples+=entry[5];retention_records+=entry[6]
        dense_samples+=entry[7];dense_records+=entry[8]
        for t,values in entry[9].items():
            bucket=hour_retention.setdefault(t,[0,0,0,0])
            for i,value in enumerate(values):bucket[i]+=value
        if entry[2] is not None:first=entry[2] if first is None else min(first,entry[2]);stop=entry[3] if stop is None else max(stop,entry[3])
    properties=Path(root)/'dmd_properties.h5'
    if properties.exists():
        overhead=properties.stat().st_size;total+=overhead
        if hour_bytes:hour_bytes[min(hour_bytes)]+=overhead
    span=(stop-first)/1e6 if first is not None else 0
    hourly=sorted(hour_bytes.keys()|hour_retention.keys())
    if hourly:
        hourly=list(range(hourly[0],hourly[-1]+3600,3600))
    labels=[]
    if hourly and stop:
        labels=[datetime.fromtimestamp(t,timezone.utc).strftime('%m-%d %H:%M') for t in [hourly[0],(hourly[0]+stop/1e6)/2,stop/1e6]]
    return dict(bytes=total,files=len(files),records=count,first_mode_us=first,latest_mode_us=stop,days=days,
        retained_complex_samples=retained_samples,tested_complex_samples=tested_samples,
        retention_records=retention_records,
        complex_sample_retained_fraction=retained_samples/tested_samples if tested_samples and retention_records==count else None,
        dense_complex_samples=dense_samples,dense_retention_records=dense_records,
        dense_complex_sample_retained_fraction=retained_samples/dense_samples if dense_samples and dense_records==count else None,
        observed_span_seconds=span,
        hourly_unix_s=hourly,hourly_bytes=[hour_bytes.get(h,0) for h in hourly],hour_tick_labels=labels,
        hourly_retained_complex_samples=[hour_retention.get(h,[0,0])[0] for h in hourly],
        hourly_dense_complex_samples=[hour_retention.get(h,[0,0])[1] for h in hourly],
        hourly_dense_retained_fraction=[v[0]/v[1] if v[1] and v[2]==v[3] else None
            for v in (hour_retention.get(h,[0,0,0,0]) for h in hourly)],
        projected_bytes_per_day=total/span*86400 if count>=3 and span>=600 else None,
        projection_status='provisional until chronological backfill is complete')
