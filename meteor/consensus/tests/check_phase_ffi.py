"""Exercise the production ctypes ABI with two indistinguishable radial tracks."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from meteor_phase_grouping import PhaseGrouper

g = PhaseGrouper()
n=60
t=np.arange(n,dtype=np.uint64)*1600+1_746_000_000_000_000
r=np.ascontiguousarray(100-20*np.arange(n)*.0016,dtype=np.float32)
v=np.full(n,-20000,dtype=np.float32)
s=np.full(n,20,dtype=np.float32)
phi=np.arange(n)[:,None]*1.2+np.where(np.arange(n)%2==0,.1,1.5)[:,None]*np.arange(7)
a=np.ascontiguousarray(np.exp(1j*phi),dtype=np.complex64)
labels=np.full(n,-1,dtype=np.int32)
count=g.lib.pansy_consensus_phase(t.ctypes.data,r.ctypes.data,v.ctypes.data,s.ctypes.data,
                                 a.ctypes.data,n,.65,.017,labels.ctypes.data)
assert count==2
assert np.all(labels>=0)
assert len(np.unique(labels[::2]))==1 and len(np.unique(labels[1::2]))==1
assert labels[0]!=labels[1]
a[0,0]=np.nan
assert g.lib.pansy_consensus_phase(t.ctypes.data,r.ctypes.data,v.ctypes.data,s.ctypes.data,
                                  a.ctypes.data,n,.65,.017,labels.ctypes.data)==-2
g.close()
print("PASS: native phase-grouping ABI separates two radial look-alikes and rejects missing phase")
