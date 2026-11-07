"""Re-evaluate saved phase observations without raw-data reads or service changes."""
import argparse
from pathlib import Path
import sys
import h5py
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from meteor_phase_grouping import PhaseGrouper

p=argparse.ArgumentParser()
p.add_argument("input")
p.add_argument("output")
a=p.parse_args()
with h5py.File(a.input) as h:
    all_t=h["sample_idx"][:]
    t=h["phase_sample_idx"][:].astype(np.uint64)
    idx=np.searchsorted(all_t,t)
    np.testing.assert_array_equal(all_t[idx],t)
    r=np.ascontiguousarray(h["range_km"][:][idx],dtype=np.float32)
    v=np.ascontiguousarray(h["doppler_m_s"][:][idx],dtype=np.float32)
    s=np.ascontiguousarray(10**(h["snr_db"][:][idx]/10),dtype=np.float32)
    amp=np.ascontiguousarray(h["receiver_amplitudes"][:],dtype=np.complex64)
    g=PhaseGrouper()
    labels=np.full(len(t),-1,dtype=np.int32)
    count=g.lib.pansy_consensus_phase(t.ctypes.data,r.ctypes.data,v.ctypes.data,s.ctypes.data,
                                     amp.ctypes.data,len(t),.65,.05,labels.ctypes.data)
    if count<0: raise RuntimeError(f"native grouping failed: {count}")
    full=np.full(len(all_t),-1,dtype=np.int32)
    full[idx]=labels
    old=h["legacy_track_id"][:]
    with h5py.File(a.output,"x") as out:
        for key in h:
            if key != "track_id": h.copy(key,out)
        out["track_id"]=full
        out.attrs["phase_coherence_min"]=np.float32(.65)
        out.attrs["max_phase_link_lag_s"]=np.float32(.05)
    print(f"Legacy {len(np.unique(old[old>=0]))} tracks, {np.count_nonzero(old>=0)} points; "
          f"native {count} tracks, {np.count_nonzero(full>=0)} points")
