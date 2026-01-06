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
    data_dict = dm.read(t0, t1, ("i0"))
    n_b=0
    for k in data_dict.keys():
        print(k)
        n_b+=1
        
    exit(0)
    for bi in range(n_b):
        data_dict = dm.read(start_idx+bi*dt, start_idx+bi*dt+dt, ("xc_arr","i0","i1","r0","r1","f0","f1","n_fft","rdec"))
        for k in data_dict.keys():
            r0=data_dict[k]["r0"]
            r1=data_dict[k]["r1"]
            f0=data_dict[k]["f0"]
            f1=data_dict[k]["f1"]
            rdec=data_dict[k]["rdec"]        
            n_fft=data_dict[k]["n_fft"]
            fvec=n.fft.fftshift(n.fft.fftfreq(n_fft,d=5*1600/1e6))[f0:f1]
            tvs.append(stuffr.unix2date(data_dict[k]["i0"]/1e6))
            for i in range(5):
                mean_pwr=n.sum(n.abs(data_dict[k]["xc_arr"][0:7,i,:,:]),axis=0)
                dop_idx=n.argmax(mean_pwr,axis=0)
                pprofs[i].append(n.max(mean_pwr,axis=0))
                dprofs[i].append(fvec[dop_idx])
    pprofs=n.array(pprofs)
    dprofs=n.array(dprofs)
    rvec=n.arange(1600)*0.15
    rvec=rvec[r0:r1:rdec]
#    fig, axs = plt.subplots(nrows=5,ncols=1)
#    for i in range(5):
 #       ax=axs[i]    
  #      if i == 0:
    i=0
    dB=10.0*n.log10(pprofs[i,:,:].T)
    nfloor=n.nanmedian(dB)
    dB=dB-nfloor
    m=ax.pcolormesh(tvs,rvec,dB,vmin=0)
    ax.set_xlabel("Date (UTC)")
    ax.set_ylabel("Range (km)")
    fig.autofmt_xdate()
#    cb=fig.colorbar(m,ax=ax)
 #   cb.set_label("SNR (dB)")


tnow=time.time()
plot_pprof(int(tnow*1e6-24*3600*1e6),int(tnow*1e6))
