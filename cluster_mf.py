
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
    idx=n.argsort(snr)[::-1]
    pairs=[]
    # look for measurement pair that best fits this measurement
    # has to be less than 10 ms apart and fit better than 500 meters
    # also, doppler can't change more than 15 km/s
    while len(idx)>1:
        i = idx[0]
        dt = (tx_idx[idx]-tx_idx[i])/1e6
        dr = rg[i]-(rg[idx]-dop[idx]*dt)
        ddop = dop[i]-dop[idx]

        fit_idx=n.where( (n.abs(dt) < 10e-3) & (n.abs(dr)<500) & (n.abs(ddop)<15e3))[0]
        pair_idx=idx[fit_idx]
        if len(pair_idx) > 1:
            print(pair_idx)
            pairs.append(pair_idx)
            if False:
                plt.subplot(131)
                plt.plot(dt[fit_idx],dr[fit_idx],".")
                plt.subplot(132)
                plt.plot(tx_idx[pair_idx],rg[pair_idx],".")
                plt.subplot(133)
                plt.plot(tx_idx[pair_idx],dop[pair_idx],".")
                plt.tight_layout()
                plt.show()
        else:
            print("no fit. removing")
            print(pair_idx)
            
        idx=n.setdiff1d(idx,pair_idx)
    pairs=n.array(pairs,dtype=n.int64)
    plt.subplot(121)
    for p in pairs:
        plt.plot(tx_idx[p],rg[p],".")
    plt.subplot(122)
    for p in pairs:
        plt.plot(tx_idx[p],dop[p],".")
    plt.show()
    

        
        
        

def read_mf_output(dm_mf,i0,i1,snr_threshold=7,tx_pwr_threshold=1e9):
    txpa=n.array([],dtype=n.float32)
    txidxa=n.array([],dtype=n.uint64)
    rnga=n.array([],dtype=n.float32)
    dopa=n.array([],dtype=n.float32)
    snra=n.array([],dtype=n.float32)
    beam=n.array([],dtype=n.int32)         
    # read 20 pulse sequences   
    try:
        data_dict = dm_mf.read(i0, i1, ("tx_pwr","max_range","tx_idxs","max_dopvel","max_snr","beam_pos_idx"))
    except:
        import traceback
        traceback.print_exc()
    

    if len(data_dict.keys()) == 0:
        print("no data")
    else:
        for k in data_dict.keys():
            data=data_dict[k]
            
            txp=data["tx_pwr"]
#            print(txp)
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
d=drf.DigitalRFReader("/media/archive/")
# tx channel bounds
b=d.get_bounds("ch007")
#start_idx=b[0]#db_mf[1]-2*60*60*1000000
start_idx=dt*int(n.floor(db_mf[0]/dt))#-2*60*60*1000000
#start_idx=db_mf[0]
n_min=int(n.floor((db_mf[1]-start_idx)/dt))
for i in range(n_min):
    i0=start_idx+i*dt
    i1=start_idx+i*dt+dt 
    txpa,txidxa,rnga,dopa,snra,beam=read_mf_output(dm_mf,i0,i1)

    gidx=n.where(n.abs(dopa)>3e3)[0]
    txpa=txpa[gidx]
    txidxa=txidxa[gidx]
    rnga=rnga[gidx]
    dopa=dopa[gidx]
    snra=snra[gidx]
    beam=beam[gidx]
    if len(txpa)<5:
        print("not enough data")
        continue
    
    cluster_idx=cluster(txidxa,rnga,dopa,snra)
    
