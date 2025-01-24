
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

dt=60000000
n_min=int(n.floor((db_mf[1]-db_mf[0])/dt))
for i in range(n_min):
    i0=db_mf[0]+i*dt
    i1=db_mf[0]+i*dt+dt
    txpa=n.array([],dtype=n.float32)
    txidxa=n.array([],dtype=n.uint64)
    rnga=n.array([],dtype=n.float32)
    dopa=n.array([],dtype=n.float32)        
    data_dict = dm_mf.read(i0, i1, ("tx_pwr","max_range","tx_idxs","max_dopvel"))
    for k in data_dict.keys():
        data=data_dict[k]
        
        txp=data["tx_pwr"]
        if n.min(txp) > 1e10:
            txpa=n.concatenate((txpa,txp))
            txidxa=n.concatenate((txidxa,data["tx_idxs"]))
            rnga=n.concatenate((rnga,data["max_range"]))
            dopa=n.concatenate((dopa,data["max_dopvel"]))            
        else:
            print("low txpower. skipping")
    plt.subplot(211)
    plt.plot(txidxa,rnga,".")
    plt.subplot(212)
    plt.plot(txidxa,dopa,".")
    plt.show()
        
#        print(data.keys())
        
    
