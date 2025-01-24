
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


mf_metadata_dir = "/media/archive/metadata/mf"
dm_mf = drf.DigitalMetadataReader(mf_metadata_dir)
db_mf = dm_mf.get_bounds()

dt=10000000
#n_min=int(n.floor((db_mf[1]-db_mf[0])/dt))
n_min=5
#i0=db_mf[0]
start_idx=db_mf[1]-5*60*1000000
for i in range(n_min):
    i0=start_idx+i*dt
    i1=start_idx+i*dt+dt
    txpa=n.array([],dtype=n.float32)
    txidxa=n.array([],dtype=n.uint64)
    rnga=n.array([],dtype=n.float32)
    dopa=n.array([],dtype=n.float32)
    snra=n.array([],dtype=n.float32)
    beam=n.array([],dtype=n.int32)            
    data_dict = dm_mf.read(i0, i1, ("tx_pwr","max_range","tx_idxs","max_dopvel","max_snr","beam_pos_idx"))

    if len(data_dict.keys()) == 0:
        print("no data")
        continue
    for k in data_dict.keys():
        data=data_dict[k]
        
        txp=data["tx_pwr"]
        if n.min(txp) > 1e10:
            txpa=n.concatenate((txpa,txp))
            txidxa=n.concatenate((txidxa,data["tx_idxs"]))
            rnga=n.concatenate((rnga,data["max_range"]))
            dopa=n.concatenate((dopa,data["max_dopvel"]))
            snra=n.concatenate((snra,data["max_snr"]))
            beam=n.concatenate((beam,data["beam_pos_idx"]))            
        else:
            print("low txpower. skipping")
    gidx=n.where(snra>7)[0]
    plt.subplot(311)
    plt.scatter((txidxa[gidx]-n.min(txidxa[gidx]))/1e6,rnga[gidx],s=1,c=10.0*n.log10(snra[gidx]),vmin=13,vmax=30)
    cb=plt.colorbar()
    plt.title("%s"%(stuffr.unix2datestr(i0/1e6)))
    cb.set_label("SNR (dB)")
    plt.ylabel("Range (km)")
    plt.subplot(312)
    plt.scatter((txidxa[gidx]-n.min(txidxa[gidx]))/1e6,dopa[gidx]/1e3,s=1,c=10.0*n.log10(snra[gidx]),vmin=13,vmax=30)
    cb=plt.colorbar()
    cb.set_label("SNR (dB)")
    plt.ylim([-100,10])
    plt.ylabel("Doppler velocity (km/s)")
    plt.subplot(313)
    plt.scatter((txidxa[gidx]-n.min(txidxa[gidx]))/1e6,10.0*n.log10(snra[gidx]),s=1,c=beam[gidx]%5,vmin=0,vmax=4,cmap="turbo")
    cb=plt.colorbar()
    cb.set_label("beam position")
    plt.ylabel("SNR (dB)")
    plt.xlabel("Time (s)")
    plt.ylim([7,50])
    plt.tight_layout()
    plt.show()
        
#        print(data.keys())
        
    
