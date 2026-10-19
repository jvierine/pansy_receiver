"""Count native records without reading the (potentially large) spectra arrays."""
from pathlib import Path
import h5py
import numpy as np

_cache={}


def inventory(root):
    files=list(Path(root).glob('*T*/pmse@*.h5')); count=0; first=None; stop=None; total=0
    for path in files:
        st=path.stat();total+=st.st_size;state=(st.st_mtime_ns,st.st_size)
        entry=_cache.get(str(path))
        if entry is None or entry[0]!=state:
            keys=[];ends=[]
            try:
                with h5py.File(path,'r') as h:
                    for key,g in h.items():
                        if not key.isdecimal() or int(g.get('complete',np.array(0))[()])!=1:continue
                        keys.append(int(key));ends.append(int(g['mode_stop_sample_exclusive'][()]))
            except OSError:continue
            entry=(state,len(keys),min(keys,default=None),max(ends,default=None))
            _cache[str(path)]=entry
        count+=entry[1]
        if entry[2] is not None:first=entry[2] if first is None else min(first,entry[2]);stop=entry[3] if stop is None else max(stop,entry[3])
    properties=Path(root)/'dmd_properties.h5'
    if properties.exists():total+=properties.stat().st_size
    span=(stop-first)/1e6 if first is not None else 0
    return dict(bytes=total,files=len(files),records=count,first_mode_us=first,latest_mode_us=stop,
        observed_span_seconds=span,
        projected_bytes_per_day=total/span*86400 if count>=12 and span>=3600 else None,
        projection_status='provisional until chronological backfill is complete')
