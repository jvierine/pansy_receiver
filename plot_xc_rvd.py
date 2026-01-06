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

def noise_floor_median(x):
    med = n.median(x)
    mad = n.median(n.abs(x - med))
    sigma = 1.4826 * mad  # robust std estimate
    return(med, sigma)

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
    tvs=[]
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
    os.system("rm tmp/spec*.png")
    for bi in range(n_b):
        print(bi)
        data_dict = dm.read(keys[bi]-10, keys[bi]+10, ("xc_arr","i0","i1","rvec","fvec"))
        tvs.append(stuffr.unix2date(keys[bi]/1e6))
        for k in data_dict.keys():
            rvec=data_dict[k]["rvec"]
            xc=data_dict[k]["xc_arr"]
            fvec=data_dict[k]["fvec"]        
            fidx=n.where(n.abs(fvec)<fmax)[0]
            mean_pwr=n.sum(n.abs(data_dict[k]["xc_arr"][0:7,0,:,:]),axis=0)
            plt.hist(10.0*n.log10(mean_pwr.flatten()),bins=100)
            plt.show()
    #            noise_floor=
            noise_floor, nf_sigma=noise_floor_median(mean_pwr)
            snr=(mean_pwr-noise_floor)/noise_floor
            psnr=n.copy(snr)
            plt.pcolormesh(fvec,rvec,psnr.T,cmap="plasma")
#            psnr[snr<0]=1e-3
            cb=plt.colorbar()
            cb.set_label("SNR (dB)")
            plt.title(stuffr.unix2datestr(keys[bi]/1e6))
            plt.xlim([-20,20])
            plt.ylim([75,100])
            plt.xlabel("Doppler (Hz)")
            plt.ylabel("Range (km)")
            plt.savefig("tmp/spec-%06d.png"%(bi))
            plt.close()
#            plt.show()
            # calculate moments
            snr_prof=n.zeros(len(rvec))
            mean_prof=n.zeros(len(rvec))
            width_prof=n.zeros(len(rvec))
            for i in range(len(rvec)):
                snr_prof[i]=n.trapz(snr[fidx,i],fvec[fidx])
                mean_prof[i]=n.trapz(snr[fidx,i]*fvec[fidx],fvec[fidx])/snr_prof[i]
                width_prof[i]=n.sqrt(n.trapz(snr[fidx,i]*(fvec[fidx]-mean_prof[i])**2,fvec[fidx])/snr_prof[i])

            S[bi,:]=snr_prof
            M[bi,:]=mean_prof
            W[bi,:]=width_prof            
    tvec=n.array(keys)
    print(S.shape)
    print(M.shape)
    print(W.shape)
    S[S<0]=1e-3
    min_snr=10

    import matplotlib.dates as mdates

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 8))

    # --- Power (dB) ---
    pcm0 = axes[0].pcolormesh(
        tvs, rvec, 10.0 * n.log10(S.T),
        cmap="plasma", vmin=0, vmax=50
    )
    fig.colorbar(pcm0, ax=axes[0], label="SNR (dB)")
    axes[0].set_ylim([75, 100])
    axes[0].set_ylabel("Range (km)")
  #  axes[0].set_title("SNR (dB)")
    axes[0].tick_params(labelbottom=False)
    
    # --- M ---
    M_masked = M.copy()
    M_masked[S < min_snr] = n.nan
    pcm1 = axes[1].pcolormesh(
        tvs, rvec, M_masked.T*3e8/47.5e6/2,
        cmap="seismic", vmin=-10, vmax=10
    )
    fig.colorbar(pcm1, ax=axes[1], label="Doppler mean (Hz)")
    axes[1].set_ylim([75, 100])
    axes[1].set_ylabel("Range (km)")
    #axes[1].set_title("Doppler mean (Hz)")
    axes[1].tick_params(labelbottom=False)
    
    # --- W ---
    W_masked = W.copy()
    W_masked[S < min_snr] = n.nan
    pcm2 = axes[2].pcolormesh(
        tvs, rvec, W_masked.T*3e8/47.5e6/2,
        cmap="plasma", vmin=0, vmax=6
    )
    fig.colorbar(pcm2, ax=axes[2], label="Doppler width (Hz)")
    axes[2].set_ylim([75, 100])
    axes[2].set_ylabel("Range (km)")
#    axes[2].set_title("Doppler width (Hz)")
#    axes[2].set_xlabel("Time")
    
    # --- Datetime formatting (bottom panel only) ---
    axes[2].xaxis.set_major_locator(mdates.AutoDateLocator())
    axes[2].xaxis.set_major_formatter(mdates.ConciseDateFormatter(
        axes[2].xaxis.get_major_locator()
    ))
    
    fig.tight_layout()
    fig.autofmt_xdate()
    plt.savefig("rvd-%06d.png"%(int(n.floor(t0/24/3600/1e6))))
    print("saving")
    plt.close()


    cmd="ffmpeg -framerate 10 -pattern_type glob -i \"tmp/spec*.png\" -c:v libx264 -pix_fmt yuv420p -profile:v high -level 4.1 -crf 23 -preset slow -movflags +faststart rvd-%06d.mp4"%(int(t0/24/3600/1e6))
    os.system(cmd)
    
    

dm = drf.DigitalMetadataReader("/media/archive/metadata/xc_rvd")
b=dm.get_bounds()
day0=int(n.floor(b[0]/1e6/24/3600))
day1=int(n.floor(b[1]/1e6/24/3600))

for di in range(day0,day1):
    print(di)
    plot_pprof(int(di*24*3600*1e6),int((di+1)*24*3600*1e6))
