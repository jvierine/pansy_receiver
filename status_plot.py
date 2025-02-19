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

def get_meteors(fig,ax,dt=24*3600*1000000):
    """
    plot latest meteors
    """
    dm = drf.DigitalMetadataReader(pc.detections_metadata_dir)
    b = dm.get_bounds()
    data_dict = dm.read(b[1]-dt, b[1], ("xhat","tx_idx","snr"))
    tv=[]
    v0s=[]
    r0s=[]
    durs=[]
    snrs=[]
    for k in data_dict.keys():
        data=data_dict[k]
        xhat=data["xhat"]
        snr=data["snr"]
        snrs.append(n.max(snr))
        dur=(n.max(data["tx_idx"])-n.min(data["tx_idx"]))/1e6
        r0=xhat[0]
        v0=xhat[1]
        t0=k
        tv.append(n.datetime64(stuffr.unix2date(t0/1e6)))
        v0s.append(v0)
        r0s.append(r0)
        durs.append(dur)
    v0s=-1*n.array(v0s)
    m=ax.scatter(tv,r0s,c=v0s,cmap="turbo",s=1,vmin=0,vmax=73)
    fig.autofmt_xdate()
    n_days=dt/3600/1000000/24
    ax.set_title("%1.1f meteors per day"%(len(r0s)/n_days))
    ax.set_xlabel("Date (UTC)")
    ax.set_ylabel("Range (km)")
#    cb=fig.colorbar(m,ax=ax)
 #   cb.set_label("Doppler (km/s)")
    return(t0,v0s,r0s)

def plot_pmse_modes(fig,ax,dt=24*3600*1000000):
    dm = drf.DigitalMetadataReader(pc.mesomode_metadata_dir)
    b = dm.get_bounds()
    start_idx=b[1]-dt
    dt2=1800*1000000
    n_t=int(n.floor((b[1]-start_idx)/dt2))
    for ti in range(n_t):
        dd=dmm.read(start_idx+ti*dt2,start_idx+ti*dt2+dt2)
        for k in kl:
            on.append(dd[k]["start"])
            off.append(dd[k]["end"])

def get_xc(fig,ax,dt=24*3600*1000000):
    """
    plot latest pmse
    """
    dm = drf.DigitalMetadataReader(pc.xc_metadata_dir)
    b = dm.get_bounds()
    start_idx=b[1]-dt
    n_b=int((b[1]-start_idx)/dt)
    pprofs=[]
    dprofs=[]
    for i in range(5):
        pprofs.append([])
        dprofs.append([])
    r0=0
    r1=1
    rdec=1
    tvs=[]
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
    
# how many receiver restarts in the log 
# sandra-rx events

# plot latest 24 hours of pmse

# plot latest 24 hours of meteors

d0=drf.DigitalMetadataReader(pc.mf_metadata_dir)
d1=drf.DigitalMetadataReader(pc.tx_metadata_dir)
d2=drf.DigitalMetadataReader(pc.detections_metadata_dir)
d3=drf.DigitalMetadataReader(pc.cut_metadata_dir)
d4=drf.DigitalMetadataReader(pc.mesomode_metadata_dir)
d5=drf.DigitalMetadataReader(pc.xc_metadata_dir)#meteor_cal_metadata_dir="/media/analysis/metadata/meteor_cal"
d6=drf.DigitalMetadataReader(pc.simple_fit_metadata_dir)
dr=drf.DigitalRFReader(pc.raw_voltage_dir)

tnow=time.time()*1e6
mfb=d0.get_bounds()
latest_mf=stuffr.unix2datestr(mfb[1]/1e6)

txb=d1.get_bounds()
latest_tx=stuffr.unix2datestr(txb[1]/1e6)

detb=d2.get_bounds()
latest_det=stuffr.unix2datestr(detb[1]/1e6)

cutb=d3.get_bounds()
latest_cut=stuffr.unix2datestr(cutb[1]/1e6)

modeb=d4.get_bounds()
latest_mode=stuffr.unix2datestr(modeb[1]/1e6)

xcb=d5.get_bounds()
latest_xc=stuffr.unix2datestr(xcb[1]/1e6)

fitb=d6.get_bounds()
latest_fit=stuffr.unix2datestr(fitb[1]/1e6)

b=dr.get_bounds("ch000")
latest_raw=stuffr.unix2datestr(b[1]/1e6)

try:
    h=h5py.File("/tmp/last_rem.h5","r")
    last_del=h["last_del"][()]
    print("removed non-meso mode up to %s (%1.0f s behind)"%(stuffr.unix2datestr(last_del/1e6),(tnow-last_del)/1e6))
    h.close()
except:
    pass


def plot_status():
    raw_delay=(tnow-b[1])/1e6
    print("raw voltage extent %s-%s (%1.0f s behind)"%(stuffr.unix2datestr(b[0]/1e6),latest_raw,raw_delay))
    tx_delay=(tnow-txb[1])/1e6
    print("latest tx %s (%1.0f s behind)"%(latest_tx,tx_delay))
    mf_delay=(tnow-mfb[1])/1e6
    print("latest mf %s (%1.0f s behind)"%(latest_mf,mf_delay))
    det_delay=(tnow-detb[1])/1e6
    print("latest det %s (%1.0f s behind)"%(latest_det,det_delay))
    cut_delay=(tnow-cutb[1])/1e6
    print("latest cut %s (%1.0f s behind)"%(latest_cut,cut_delay))
    mode_delay=(tnow-modeb[1])/1e6
    print("latest mode %s (%1.0f s behind)"%(latest_mode,mode_delay))
    xc_delay=(tnow-xcb[1])/1e6
    print("latest xc %s (%1.0f s behind)"%(latest_xc,xc_delay))
    fit_delay=(tnow-fitb[1])/1e6
    print("latest fit %s (%1.0f s behind)"%(latest_fit,fit_delay))

    labels=["Raw voltage","Transmit pulse detect","Match function","Clustering","Cutting","Mode boundaries","Cross-spectra"]
    delays=n.array([raw_delay,tx_delay,mf_delay,det_delay,cut_delay,mode_delay,xc_delay])/3600.0

    fig,(ax0,ax1)=plt.subplots(2,1,sharex=True,figsize=(8,8))
    get_xc(fig,ax0)
    get_meteors(fig,ax1)
    fig.tight_layout()
    plt.savefig("status.png")
    plt.close()
    #plt.show()

    fig,ax=plt.subplots(1,1,figsize=(8,4))
    ax.bar(labels,delays,0.6)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_ylabel("Processing delay (hours)")
    fig.tight_layout()
    plt.savefig("processing.png")
    plt.close()
    os.system("scp processing.png status.png j@4.235.86.214:/var/www/html/pansy/")


while True:
    try:
        plot_status()
    except:
        import traceback
        traceback.print_exc()
    time.sleep(3600)
