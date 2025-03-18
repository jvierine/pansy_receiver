import pansy_config as pc
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
import itertools
import glob
import h5py 

dm=drf.DigitalMetadataReader(pc.mf_isr_metadata_dir)


b=dm.get_bounds()

start_idx=b[1]-8*3600*1000000
dt=60000000

n_blocks=int(n.floor(b[1]-start_idx)/dt)

for i in range(n_blocks):
    data=dm.read(start_idx+i*dt,start_idx+i*dt+dt)
    snr=[]
    rangek=[]
    dop=[]
    t=[]
    for k in data.keys():
        snr.append(data[k]["max_snr"])
        rangek.append(data[k]["max_range"])
        dop.append(data[k]["max_dopvel"])
        t.append(data[k]["tx_idxs"])
    t=n.array(t)
    snr=n.array(snr)
    rangek=n.array(rangek)
    dop=n.array(dop)
    plt.scatter(t,rangek,c=dop,cmap="turbo",vmin=-70e3,vmax=20e3)
    plt.colorbar()
    plt.show()