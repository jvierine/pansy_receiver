
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

def cluster(tx_idx,
            rg,
            dop,
            snr,
            min_det=5
            ):
    idx=n.where(n.abs(dop) > 3e3)[0]
    #tx_idx=tx_idx[gidx]
    #rg=rg[gidx]
    #dop=dop[gidx]
    tv=tx_idx/1e6
    #idx=n.arange(len(dop),dtype=n.int64)
    meteor_idxs=[]
    while len(idx) > min_det:
        # try to add measurements to peak snr obs
        i0=n.argmax(snr[idx])
        t0=tv[idx[i0]]
        # doppler migration
        rg_resid=rg[idx] - (rg[idx[i0]]+(0.5*dop[idx]+0.5*dop[idx[i0]])*(tv[idx]-t0))



        gidx=idx[n.where(n.abs(rg_resid)<3e3 )[0]]
        if len(gidx)>min_det:
            plt.plot(tv[idx]-t0,rg_resid,".")
            plt.show()

            meteor_idxs.append(gidx)
        idx=n.setdiff1d(idx,gidx)
    return(meteor_idxs)

def read_mf_output(dm_mf,i0,i1,snr_threshold=7,tx_pwr_threshold=1e9):
    txpa=n.array([],dtype=n.float32)
    txidxa=n.array([],dtype=n.uint64)
    rnga=n.array([],dtype=n.float32)
    dopa=n.array([],dtype=n.float32)
    snra=n.array([],dtype=n.float32)
    beam=n.array([],dtype=n.int32)         
    # read 20 pulse sequences   
    data_dict = dm_mf.read(i0, i1, ("tx_pwr","max_range","tx_idxs","max_dopvel","max_snr","beam_pos_idx"))

    if len(data_dict.keys()) == 0:
        print("no data")
    else:
        for k in data_dict.keys():
            data=data_dict[k]
            
            txp=data["tx_pwr"]
            # use this threshold for tx power. check that it is okay!!!
            if n.min(txp) > tx_pwr_threshold:
                txpa=n.concatenate((txpa,txp))
                txidxa=n.concatenate((txidxa,data["tx_idxs"]))
                rnga=n.concatenate((rnga,data["max_range"]))
                dopa=n.concatenate((dopa,data["max_dopvel"]))
                snra=n.concatenate((snra,data["max_snr"]))
                beam=n.concatenate((beam,data["beam_pos_idx"]))            
            else:
                print("low txpower. skipping")

        gidx=n.where(snra>snr_threshold)[0]
        txpa=txpa[gidx]
        txidxa=txidxa[gidx]
        rnga=rnga[gidx]
        dopa=dopa[gidx]
        snra=snra[gidx]
        beam=beam[gidx]
    return(txpa,txidxa,rnga,dopa,snra,beam)

mf_metadata_dir = "/media/archive/metadata/mf"
dm_mf = drf.DigitalMetadataReader(mf_metadata_dir)
db_mf = dm_mf.get_bounds()

dt=10000000
#n_min=int(n.floor((db_mf[1]-db_mf[0])/dt))
start_idx=db_mf[1]-2*60*60*1000000
#start_idx=db_mf[0]
n_min=int(n.floor((db_mf[1]-start_idx)/dt))
for i in range(n_min):
    i0=start_idx+i*dt
    i1=start_idx+i*dt+dt 
    txpa,txidxa,rnga,dopa,snra,beam=read_mf_output(dm_mf,i0,i1)

    cluster_idx=cluster(txidxa,rnga,dopa,snra)

    plt.subplot(311)
    plt.scatter((txidxa-i0)/1e6,rnga,s=1,c=10.0*n.log10(snra),vmin=13,vmax=30)

    for ci in range(len(cluster_idx)):
        plt.plot((txidxa[cluster_idx[ci]]-i0)/1e6,rnga[cluster_idx[ci]])                

    plt.xlim([0,dt/1e6])
    cb=plt.colorbar()
    plt.title("%s"%(stuffr.unix2datestr(i0/1e6)))
    cb.set_label("SNR (dB)")
    plt.ylabel("Range (km)")
    plt.subplot(312)
    plt.scatter((txidxa-i0)/1e6,dopa/1e3,s=1,c=10.0*n.log10(snra),vmin=13,vmax=30)
    plt.xlim([0,dt/1e6])
    cb=plt.colorbar()
    cb.set_label("SNR (dB)")
    plt.ylim([-100,10])
    plt.ylabel("Doppler velocity (km/s)")
    plt.subplot(313)
    plt.scatter((txidxa-i0)/1e6,10.0*n.log10(snra),s=1,c=beam%5,vmin=0,vmax=4,cmap="turbo")
    plt.xlim([0,dt/1e6])
    cb=plt.colorbar()
    cb.set_label("beam position")
    plt.ylabel("SNR (dB)")
    plt.xlabel("Time (s)")
    plt.ylim([7,50])
    plt.tight_layout()
    plt.show()
        
#        print(data.keys())
        
    
