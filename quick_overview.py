import numpy as n
import pansy_detect as pd
import matplotlib.pyplot as plt
import pansy_modes as pm
import scipy.signal.windows as sw
import scipy.constants as c
import digital_rf as drf
import os
import stuffr
import time
import scipy.fftpack as fp

d=drf.DigitalRFReader("/media/archive")
b=d.get_bounds("ch000")
i0=b[0]+10000000
n_windows=1800

# 17th of Jan, 2025
ut0=1737072000

n_days=int(n.ceil((b[1]/1e6-ut0)/(24*3600)))
#n_days=int(n.ceil((b[1]-i0)/(3600*1000000)))
dt=int(n.round(24*3600*1000000/n_windows))

window_len=1600

for di in range(n_days):
    ut00=ut0+di*24*3600
    fname="overview-%d.png"%(ut00)
    if os.path.exists(fname):
        print("already exists")
        continue
    if (b[1]-ut00*1000000) < (24*3600*1000000):
        print("not enough to fill day")
        continue
    if ut00*1000000 < i0:
        print("day starts before samples. skipping")
        continue
    
    S=n.zeros([n_windows,window_len],dtype=n.float64)
    tv=[]
    i00=ut00*1000000
    for i in range(n_windows):
        print(i)
        read_idx=i00+dt*i
        try:
            z=d.read_vector_c81d(read_idx,window_len,"ch007")
            za=n.abs(z)
            za=n.convolve(za,n.repeat(1/8,8),mode="same")
            S[i,:]=za
        except:
            print("can't read")
        tv.append(n.datetime64(int(read_idx/1e6), 's'))
        
    plt.figure(figsize=(8*2,6.5*2))
    plt.pcolormesh(tv,n.arange(1600),S.T)
    plt.title(stuffr.unix2datestr(ut00))
    plt.xlabel("Time")
    plt.ylabel(r"IPP ($\mu$s)")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
