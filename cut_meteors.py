import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
import digital_rf as drf
import os
import stuffr
import time
import pansy_config as pc
import cluster_mf as cmf

# meteor detections
dm = drf.DigitalMetadataReader(pc.detections_metadata_dir)
bm = dm.get_bounds()

# match function outputs
dmf = drf.DigitalMetadataReader(pc.mf_metadata_dir)
bmf = dmf.get_bounds()

# raw voltage
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

    # for each detection
    for k in data_dict.keys():
        data=data_dict[k]
        print(data)
        xhat=data["xhat"]
        tx_idx=data["tx_idx"]
        range_km=data["range"]
        snr=data["snr"]
        beam_idx=data["beam"]
        doppler_ms=data["doppler"]

        idxidx=n.argsort(tx_idx)
        tx_idx=tx_idx[idxidx]
        snr=snr[idxidx]
        beam_idx=beam_idx[idxidx]
        range_km=range_km[idxidx]
        doppler_ms=doppler_ms[idxidx]
        
        n_ipp=len(tx_idx)
        RTI=n.zeros([n_ipp,256],dtype=n.float32)
        chs=["ch000","ch002","ch003","ch004","ch005","ch006"]
        drg=c.c/1e6/2/1e3
        w=n.repeat(1/4,4)
        w2=n.repeat(1/16,16)
        for ipp in range(n_ipp):
            rg=int(range_km[ipp]/drg)
            for ch in chs:
                z=n.convolve(d.read_vector_c81d(tx_idx[ipp]+rg-64,256,ch),w,mode="same")
                RTI[ipp,:]+=n.convolve(n.abs(z)**2.0,w2,mode="same")

        t0=n.min(tx_idx)/1e6

        dur=(n.max(tx_idx)-n.min(tx_idx))/1e6
        n_pad=int((dur*0.1)*1000000)

        m_txpa,m_txidxa,m_rnga,m_dopa,m_snra,m_beam=cmf.read_mf_output(dmf,n.min(tx_idx)-n_pad,n.max(tx_idx)+n_pad)

        if n.max(snr)>100:
            plt.subplot(221)
            plt.plot(tx_idx/1e6-t0,range_km,".")
            plt.plot(m_txidxa/1e6-t0,m_rnga,"x")
            plt.ylabel("Range (km)")
            plt.xlabel("Time (s)")        
            plt.subplot(222)
            plt.plot(tx_idx/1e6-t0,doppler_ms/1e3,".")
            plt.plot(m_txidxa/1e6-t0,m_dopa/1e3,"x")
            plt.ylabel("Doppler (km/s)")
            plt.xlabel("Time (s)")
            plt.subplot(223)
            plt.scatter(tx_idx/1e6-t0,10.0*n.log10(snr),s=1,c=n.mod(beam_idx,5),cmap="berlin",vmin=0,vmax=4)
            plt.scatter(m_txidxa/1e6-t0,10.0*n.log10(m_snra),s=1,c=n.mod(m_beam,5),cmap="berlin",vmin=0,vmax=4)

            plt.ylabel("SNR (dB)")
            plt.xlabel("Time (s)")        
            
            cb=plt.colorbar(location="top")
            cb.set_label("Beam index")
            plt.subplot(224)
            plt.pcolormesh(tx_idx/1e6-t0,n.arange(256),10.0*n.log10(RTI.T),cmap="plasma")
            plt.colorbar(location="top")
            plt.xlabel("Time (s)")                
            plt.tight_layout()
            plt.show()
        


