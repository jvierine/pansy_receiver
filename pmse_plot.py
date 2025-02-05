import numpy as n
import matplotlib.pyplot as plt
import digital_rf as drf
import os
import stuffr
import time
import pansy_config as pc



mf_metadata_dir=pc.mf_metadata_dir
md = drf.DigitalMetadataReader(mf_metadata_dir)
b = md.get_bounds()

dt=15000000

L=3600*1000000
n_windows=int(n.floor(L/dt))

rvec=n.arange(0,150)
n_r=len(rvec)
rti=n.zeros([n_windows,n_r],dtype=n.float32)
n_avg=n.zeros([n_windows,n_r],dtype=n.float32)

for i in range(n_windows):
    print(i)
    i0=b[0]+i*dt
    i1=b[0]+i*dt+dt
    data_dict = md.read(i0, i1, ("tx_pwr","max_range","tx_idxs","max_dopvel","max_snr","beam_pos_idx"))
    for k in data_dict.keys():
        ridx=int(n.round(data_dict[k]["max_range"][0]))
#        print(ridx)
        rti[i,ridx]+=data_dict[k]["max_snr"][0]
        n_avg[i,ridx]+=1.0
rti=rti/(n_avg+1.0)
plt.pcolormesh(rti.T,vmin=0,vmax=10)
plt.colorbar()
plt.show()
