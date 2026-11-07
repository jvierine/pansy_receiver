"""Regression for an inspected two-track station replay, without raw-data I/O."""
import argparse
from pathlib import Path
import sys
import h5py
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from meteor_phase_grouping import PhaseGrouper

p=argparse.ArgumentParser()
p.add_argument("input")
a=p.parse_args()
with h5py.File(a.input) as h:
    t=h["sample_idx"][:]
    samples=h["phase_sample_idx"][:].astype(np.uint64)
    idx=np.searchsorted(t,samples)
    np.testing.assert_array_equal(t[idx],samples)
    r=np.ascontiguousarray(h["range_km"][:][idx],dtype=np.float32)
    v=np.ascontiguousarray(h["doppler_m_s"][:][idx],dtype=np.float32)
    s=np.ascontiguousarray(10**(h["snr_db"][:][idx]/10),dtype=np.float32)
    amp=np.ascontiguousarray(h["receiver_amplitudes"][:],dtype=np.complex64)
    legacy=h["legacy_track_id"][:][idx]
g=PhaseGrouper()
labels=np.full(len(samples),-1,dtype=np.int32)
count=g.lib.pansy_consensus_phase(samples.ctypes.data,r.ctypes.data,v.ctypes.data,s.ctypes.data,
                                 amp.ctypes.data,len(samples),.65,.05,labels.ctypes.data)
assert count==2, f"expected two intact meteors, got {count}"
for old in np.unique(legacy[legacy>=0]):
    ids=labels[legacy==old]
    ids=ids[ids>=0]
    counts=np.bincount(ids,minlength=count)
    assert counts.max() >= .9*np.count_nonzero(legacy==old), "fragmented or lost meteor"
for new in range(count):
    assigned=legacy[labels==new]
    assert len(np.unique(assigned[assigned>=0]))==1, "combined separate meteors"
print("PASS: both inspected station meteors remain intact and separate with >=90% membership retention")
