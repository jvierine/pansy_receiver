"""Keep one already-computed meteor diagnostic, never an accumulating product."""
import fcntl
import os
from pathlib import Path
import h5py
import numpy as np

def save_event(key, times, range_km, power, doppler_mps, path=None):
    """Copy existing processor results only: no decoding, FFT or refitting."""
    target=Path(path or Path.home()/'.local/state/pansy-monitor/latest-event.h5')
    target.parent.mkdir(parents=True,exist_ok=True)
    with open(str(target)+'.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if target.exists():
            with h5py.File(target,'r') as old:
                if int(old['key'][()])>=int(key):return
        temporary=target.with_suffix('.h5.partial')
        with h5py.File(temporary,'w') as h:
            h['key']=int(key)
            h['time_unix_s']=np.asarray(times,dtype=np.float64)
            h['range_km']=np.asarray(range_km,dtype=np.float32)
            h.create_dataset('power_snr',data=np.asarray(power,dtype=np.float32),compression='gzip')
            h['doppler_mps']=np.asarray(doppler_mps,dtype=np.float32)
            h.attrs['provenance']='Existing range-Doppler matched-filter result; no extra signal processing'
        os.replace(temporary,target)
