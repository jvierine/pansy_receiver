
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
    data_dict = dm_mf.read(i0, i1, ("tx_pwr","max_range","tx_idxs"))
    for k in data_dict.keys():
        data=data_dict[k]
        print(data.keys())
        
    
