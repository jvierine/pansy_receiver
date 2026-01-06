import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
import digital_rf as drf
import os
import stuffr
import time
import pansy_config as pc
import traceback
import h5py
import plot_simple_fits as psf
import process_cut_meteor as pcm
import healpix_radiant as hpr


def plot_pprof(t0,t1):
    """
    plot latest pmse
    """
    dm = drf.DigitalMetadataReader("/media/archive/metadata/xc_rvd")
    pprofs=[]
    dprofs=[]
    wprofs=[]
    data_dict = dm.read(t0, t1, ("i0","r0","r1","rvec","fvec"))
    n_b=0
    fmax=20.0
    keys=[]
    rvec=None
    for k in data_dict.keys():
        keys.append(k)
        print(k)
        rvec=data_dict[k]["rvec"]
        print(n.min(rvec))
        print(n.max(rvec))
        n_b+=1
    if n_b==0:
        print("no data")
        return(0)
    n_r=len(rvec)
    S=n.zeros([n_b,n_r],dtype=n.float32)
    M=n.zeros([n_b,n_r],dtype=n.float32)
    W=n.zeros([n_b,n_r],dtype=n.float32)
    
    for bi in range(n_b):
        print(bi)
        data_dict = dm.read(keys[bi]-10, keys[bi]+10, ("xc_arr","i0","i1","rvec","fvec"))
        for k in data_dict.keys():
            rvec=data_dict[k]["rvec"]
            xc=data_dict[k]["xc_arr"]
            fvec=data_dict[k]["fvec"]        
            fidx=n.where(n.abs(fvec)<fmax)[0]
            mean_pwr=n.sum(n.abs(data_dict[k]["xc_arr"][0:7,0,:,:]),axis=0)
            noise_floor=n.median(mean_pwr)
            snr=(mean_pwr-noise_floor)/noise_floor
            plt.pcolormesh(10.0*n.log10(snr[fidx,:].T))
            plt.show()
            snr_prof=n.zeros(len(rvec))
            mean_prof=n.zeros(len(rvec))
            width_prof=n.zeros(len(rvec))
            for i in range(len(rvec)):
                snr_prof[i]=n.trapz(snr[fidx,ri],fvec[fidx])
                mean_prof[i]=n.trapz(snr[fidx,ri]*fvec[fidx],fvec[fidx])/snr_prof[i]
                width_prof[i]=n.sqrt(n.trapz(snr[fidx,ri]*(fvec[fidx]-mean_prof[i])**2,fvec[fidx])/snr_prof[i])

            S[bi,:]=snr_prof
            M[bi,:]=mean_prof
            W[bi,:]=width_prof            
    tvec=n.array(keys)
    plt.pcolormesh(311)
    plt.pcolormesh(tvec/1e6,rvec,10.0*n.log10(S.T),cmap="plasma")
    plt.colorbar()
    plt.pcolormesh(312)
    plt.pcolormesh(tvec/1e6,rvec,M.T,cmap="turbo")
    plt.colorbar()
    plt.pcolormesh(313)
    plt.pcolormesh(tvec/1e6,rvec,W.T,cmap="plasma")
    plt.colorbar()
    
    plt.show()


tnow=time.time()
plot_pprof(int(tnow*1e6-5*3600*1e6),int(tnow*1e6))
