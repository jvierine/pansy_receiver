
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
import pansy_config as pc
import glob
import h5py 

def read_mf_output(dm_mf,i0,i1,snr_threshold=7,tx_pwr_threshold=1e9):

    txpa=n.array([],dtype=n.float32)
    txidxa=n.array([],dtype=n.uint64)
    rnga=n.array([],dtype=n.float32)
    dopa=n.array([],dtype=n.float32)
    snra=n.array([],dtype=n.float32)
    try:
        data_dict = dm_mf.read(i0, i1, ("tx_pwr","max_range","tx_idxs","max_dopvel","max_snr"))
    except:
        import traceback
        traceback.print_exc()

    if len(data_dict.keys()) == 0:
        pass
    else:
        for k in data_dict.keys():
            data=data_dict[k]
            
            txp=data["tx_pwr"]
            print(txp)
            if n.min(txp) > tx_pwr_threshold:
                txpa=n.concatenate((txpa,[txp]))
                txidxa=n.concatenate((txidxa,[data["tx_idxs"]]))
                rnga=n.concatenate((rnga,[data["max_range"]]))
                dopa=n.concatenate((dopa,[data["max_dopvel"]]))
                snra=n.concatenate((snra,[data["max_snr"]]))
            else:
                pass
#                print("low txpower. skipping")

        gidx=n.where(snra>snr_threshold)[0]
        txpa=txpa[gidx]
        txidxa=txidxa[gidx]
        rnga=rnga[gidx]
        dopa=dopa[gidx]
        snra=snra[gidx]
    return(txpa,txidxa,rnga,dopa,snra)


def analyze_until_now():
    dm_mf = drf.DigitalMetadataReader(pc.mf_isr_metadata_dir)
    b_mf = dm_mf.get_bounds()

    dt=10000000
    d=drf.DigitalRFReader("/media/archive/")
    b=d.get_bounds("ch007")

    start_idx=b_mf[0]
    n_min=int(n.floor((b_mf[1]-start_idx)/dt))

    for i in range(n_min):
        i0=start_idx+i*dt
        i1=start_idx+i*dt+dt

        txpa,txidxa,rnga,dopa,snra=read_mf_output(dm_mf,i0,i1)
        plt.scatter((txidxa-txidxa[0])/1e6,rnga,c=dopa,s=1,cmap="turbo")        
        plt.colorbar()
        plt.show()
if __name__ == "__main__":
    while True:
        try:
            analyze_until_now()
        except:
            import traceback
            traceback.print_exc()
        time.sleep(3600)
