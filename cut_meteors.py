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

        
        n_ipp=len(tx_idx)
        RTI=n.zeros([n_ipp,1600],dtype=n.float32)
        chs=["ch000","ch002","ch003","ch004","ch005","ch006"]
        for ipp in range(n_ipp):
            z=d.read_vector_c81d(tx_idx[ipp],1600,"ch000")
            RTI[ipp,:]=n.abs(z)**2.0

        
        plt.subplot(221)
        plt.plot(tx_idx/1e6,range_km,".")
        plt.ylabel("Range (km)")
        plt.xlabel("Time (unix)")        
        plt.subplot(222)
        plt.plot(tx_idx/1e6,doppler_ms/1e3,".")
        plt.ylabel("Doppler (km/s)")
        plt.xlabel("Time (unix)")
        plt.subplot(223)
        plt.scatter(tx_idx/1e6,10.0*n.log10(snr),s=1,c=n.mod(beam_idx,5),cmap="berlin",vmin=0,vmax=4)
        plt.ylabel("SNR (dB)")
        plt.xlabel("Time (unix)")        
        
        cb=plt.colorbar(location="top")
        cb.set_label("Beam index")

        plt.subplot(224)
        plt.pcolormesh(tx_idx/1e6,n.arange(1600),10.0*n.log10(RTI.T),cmap="plasma")
        plt.colorbar(location="top")
        plt.xlabel("Time (unix)")                
        plt.colorbar()
        plt.tight_layout()
        plt.show()
    


