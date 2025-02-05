import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
import digital_rf as drf
import os
import stuffr
import time
import pansy_config as pc

dm = drf.DigitalMetadataReader(pc.detections_metadata_dir)
bm = dm.get_bounds()

dmf = drf.DigitalMetadataReader(pc.mf_metadata_dir)
bmf = dmf.get_bounds()

d=drf.DigitalRFReader("/media/archive/")
b=d.get_bounds("ch007")

dt=60
sr=1000000

n_block=int(n.floor((bm[1]-bm[0])/dt))

for i in range(n_block):
    i0=bm[0]+i*dt*sr
    i1=bm[0]+i*dt*sr + dt*sr
    
    #odata_dict["xhat"]=det["xhat"]
    #odata_dict["tx_idx"]=txidxa[det["idx"]]
    #odata_dict["range"]=rnga[det["idx"]]

    data_dict = dm.read(i0, i1, ("tx_idx","xhat","range","doppler","snr","beam"))

    for k in data_dict.keys():
        data=data_dict[k]
        print(data)
        xhat=data["xhat"]
        tx_idx=data["tx_idx"]
        range_km=data["range"]
        snr=data["snr"]
        beam_idx=data["beam"]
        doppler_ms=data["doppler"]
        plt.subplot(221)
        plt.plot(tx_idx/1e6,range_km,".")
        plt.subplot(222)
        plt.plot(tx_idx/1e6,doppler_ms/1e3,".")
        plt.subplot(223)
        plt.scatter(tx_idx/1e6,10.0*n.log10(snr),s=1,c=beam_idx,cmap="turbo")
        cb=plt.colorbar()
        cb.set_label("Beam index")
        plt.show()
    


