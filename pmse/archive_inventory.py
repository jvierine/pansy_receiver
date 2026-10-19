"""Count native records without reading the (potentially large) spectra arrays."""
from pathlib import Path
from datetime import datetime,timezone
import h5py
import numpy as np

_cache={}


def inventory(root):
    files=list(Path(root).glob('*T*/pmse@*.h5')); count=0; first=None; stop=None; total=0;hour_bytes={}
    retained_samples=0;tested_samples=0;retention_records=0
    for path in files:
        st=path.stat();total+=st.st_size;state=(st.st_mtime_ns,st.st_size)
        hour=int(path.stem.split('@')[1])//3600*3600;hour_bytes[hour]=hour_bytes.get(hour,0)+st.st_size
        entry=_cache.get(str(path))
        if entry is None or entry[0]!=state:
            keys=[];ends=[];kept=0;tested=0;measured=0
            try:
                with h5py.File(path,'r') as h:
                    for key,g in h.items():
                        if not key.isdecimal() or int(g.get('complete',np.array(0))[()])!=1:continue
                        keys.append(int(key));ends.append(int(g['mode_stop_sample_exclusive'][()]))
                        if 'spectra' in g and 'n_tested_cells' in g:
                            # Inspect dimensions only; never load complex spectra.
                            looks,channels,cells=g['spectra'].shape
                            kept+=looks*channels*cells
                            tested+=looks*channels*int(g['n_tested_cells'][()])
                            measured+=1
            except OSError:continue
            entry=(state,len(keys),min(keys,default=None),max(ends,default=None),kept,tested,measured)
            _cache[str(path)]=entry
        count+=entry[1]
        retained_samples+=entry[4];tested_samples+=entry[5];retention_records+=entry[6]
        if entry[2] is not None:first=entry[2] if first is None else min(first,entry[2]);stop=entry[3] if stop is None else max(stop,entry[3])
    properties=Path(root)/'dmd_properties.h5'
    if properties.exists():
        overhead=properties.stat().st_size;total+=overhead
        if hour_bytes:hour_bytes[min(hour_bytes)]+=overhead
    span=(stop-first)/1e6 if first is not None else 0
    hourly=sorted(hour_bytes)
    if hourly:
        hourly=list(range(hourly[0],hourly[-1]+3600,3600))
    labels=[]
    if hourly and stop:
        labels=[datetime.fromtimestamp(t,timezone.utc).strftime('%m-%d %H:%M') for t in [hourly[0],(hourly[0]+stop/1e6)/2,stop/1e6]]
    return dict(bytes=total,files=len(files),records=count,first_mode_us=first,latest_mode_us=stop,
        retained_complex_samples=retained_samples,tested_complex_samples=tested_samples,
        retention_records=retention_records,
        complex_sample_retained_fraction=retained_samples/tested_samples if tested_samples and retention_records==count else None,
        observed_span_seconds=span,
        hourly_unix_s=hourly,hourly_bytes=[hour_bytes.get(h,0) for h in hourly],hour_tick_labels=labels,
        projected_bytes_per_day=total/span*86400 if count>=3 and span>=600 else None,
        projection_status='provisional until chronological backfill is complete')
