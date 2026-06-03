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
import pansy_plot

def metadata_bounds(path, label):
    try:
        reader = drf.DigitalMetadataReader(path)
        return reader, reader.get_bounds()
    except Exception as exc:
        print("status_plot: %s metadata unavailable at %s: %s"%(label, path, exc))
        return None, [-1, -1]

def raw_bounds(path, channel):
    try:
        reader = drf.DigitalRFReader(path)
        return reader, reader.get_bounds(channel)
    except Exception as exc:
        print("status_plot: raw voltage unavailable at %s/%s: %s"%(path, channel, exc))
        return None, [-1, -1]

def latest_label(bounds):
    if bounds[1] == -1:
        return "missing"
    return stuffr.unix2datestr(bounds[1]/1e6)

def extent_label(bounds):
    if bounds[0] == -1 or bounds[1] == -1:
        return "missing"
    return "%s-%s"%(stuffr.unix2datestr(bounds[0]/1e6), stuffr.unix2datestr(bounds[1]/1e6))

def delay_hours(tnow, bounds, scale=1.0):
    if bounds[1] == -1:
        return n.nan
    return scale*(tnow-bounds[1])/3600.0/1e6

def plot_unavailable(ax, title, exc):
    ax.text(0.5, 0.5, "%s unavailable\n%s"%(title, exc), ha="center", va="center", transform=ax.transAxes)
    ax.set_title(title)

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

def write_clock_offset():
    os.system("mkdir -p %s"%(pc.clock_metadata_dir))
    subdirectory_cadence_seconds = 3600
    file_cadence_seconds = 60
    samples_per_second_numerator = 1000000
    samples_per_second_denominator = 1
    file_name = "clock"
    dmw = drf.DigitalMetadataWriter(
        pc.clock_metadata_dir,
        subdirectory_cadence_seconds,
        file_cadence_seconds,
        samples_per_second_numerator,
        samples_per_second_denominator,
        file_name,
    )
    tnow=time.time()
    dr=drf.DigitalRFReader(pc.raw_voltage_dir)
    b=dr.get_bounds("ch000")
    data={"pc_tnow":tnow,"raw_b1":b[1]/1e6}
    try:
        dmw.write(int(tnow*1000000),data)
    except:
        import traceback
        traceback.print_exc()

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


def plot_status():
    #hpr.plot_last_48()
    try:
        psf.plot_latest_fits(save_png=True)
    except Exception as exc:
        print("status_plot: latest fit plot unavailable: %s"%(exc))
    try:
        pcm.plot_last()
    except Exception as exc:
        print("status_plot: latest meteor plot unavailable: %s"%(exc))

    d0,mfb=metadata_bounds(pc.mf_metadata_dir, "mf")
    d0_isr,mfb_isr=metadata_bounds(pc.mf_isr_metadata_dir, "mf_isr")
    d1,txb=metadata_bounds(pc.tx_metadata_dir, "tx")
    d2,detb=metadata_bounds(pc.detections_metadata_dir, "detections")
    d3,cutb=metadata_bounds(pc.cut_metadata_dir, "cut")
    d4,modeb=metadata_bounds(pc.mesomode_metadata_dir, "mesomode")
    d5,xcb=metadata_bounds(pc.xc_metadata_dir, "xc")
    d6,fitb=metadata_bounds(pc.simple_fit_metadata_dir, "simple_fit")
    dr,b=raw_bounds(pc.raw_voltage_dir, "ch000")
    #dgps=drf.DigitalMetadataReader(pc.gpslock_metadata_dir)
    #dclock=drf.DigitalMetadataReader(pc.clock_metadata_dir)

    if b[1] != -1:
        write_clock_offset()
    else:
        print("status_plot: skipping clock offset write because raw bounds are unavailable")


    tnow=time.time()*1e6
    b_mf = max(mfb[1], mfb_isr[1])
    mf_bounds = [-1, b_mf]
    latest_mf=latest_label(mf_bounds)
    latest_tx=latest_label(txb)
    latest_det=latest_label(detb)
    latest_cut=latest_label(cutb)
    latest_mode=latest_label(modeb)
    latest_xc=latest_label(xcb)
    latest_fit=latest_label(fitb)
    latest_raw=latest_label(b)

    #bg=dgps.get_bounds()
    #gdata=dgps.read(bg[1]-24*3600*1000000,bg[1])
    #holdover=0
    #holdovert=-1
    #holdovers=[]
    #for k in gdata:
    #    holdover=gdata[k]["holdover"]
    #    holdovert=k
    #    holdovers.append(holdover)
    #mean_holdover=n.mean(holdovers)

    try:
        h=h5py.File("/tmp/last_rem.h5","r")
        last_del=h["last_del"][()]
        print("removed non-meso mode up to %s (%1.0f s behind)"%(stuffr.unix2datestr(last_del/1e6),(tnow-last_del)/1e6))
        h.close()
    except:
        pass


    raw_delay=(tnow-b[1])/1e6 if b[1] != -1 else n.nan
    print("raw voltage extent %s (%1.0f s behind)"%(extent_label(b),raw_delay))
    tx_delay=(tnow-txb[1])/1e6 if txb[1] != -1 else n.nan
    print("latest tx %s (%1.0f s behind)"%(latest_tx,tx_delay))
    mf_delay=(tnow-b_mf)/1e6 if b_mf != -1 else n.nan
    print("latest mf %s (%1.0f s behind)"%(latest_mf,mf_delay))
    det_delay=(tnow-detb[1])/1e6 if detb[1] != -1 else n.nan
    print("latest det %s (%1.0f s behind)"%(latest_det,det_delay))
    cut_delay=(tnow-cutb[1])/1e6 if cutb[1] != -1 else n.nan
    print("latest cut %s (%1.0f s behind)"%(latest_cut,cut_delay))
    mode_delay=(tnow-modeb[1])/1e6 if modeb[1] != -1 else n.nan
    print("latest mode %s (%1.0f s behind)"%(latest_mode,mode_delay))
    xc_delay=(tnow-xcb[1])/1e6 if xcb[1] != -1 else n.nan
    print("latest xc %s (%1.0f s behind)"%(latest_xc,xc_delay))
    fit_delay=(tnow-fitb[1])/1e6 if fitb[1] != -1 else n.nan
    print("latest fit %s (%1.0f s behind)"%(latest_fit,fit_delay))
#    print("gpslock %s (%1.0f s holdover)"%(stuffr.unix2datestr(holdovert/1e6),holdover))

    labels=["Raw (s)","TX det","MF","Events","Cutting","Modes","FXC","fit"]
    delays=n.array([
        delay_hours(tnow,b,scale=3600.0),
        delay_hours(tnow,txb),
        delay_hours(tnow,mf_bounds),
        delay_hours(tnow,detb),
        delay_hours(tnow,cutb),
        delay_hours(tnow,modeb),
        delay_hours(tnow,xcb),
        delay_hours(tnow,fitb),
    ])

    fig,(ax0,ax1)=plt.subplots(2,1,sharex=True,figsize=(8,8))
    try:
        get_xc(fig,ax0)
    except Exception as exc:
        plot_unavailable(ax0, "Latest PMSE", exc)
    ax0.set_title(stuffr.unix2datestr(time.time()))
    try:
        get_meteors(fig,ax1)
    except Exception as exc:
        plot_unavailable(ax1, "Latest Meteors", exc)
    fig.tight_layout()
    plt.savefig("status.png")
    plt.close()
    #plt.show()

    fig,ax=plt.subplots(1,1,figsize=(8,4))
    ax.bar(labels,delays,0.6)
    ax.set_title(stuffr.unix2datestr(time.time()))
    ax.set_xticklabels(labels, rotation=90)
    ax.set_ylabel("Processing delay (hours)")
    fig.tight_layout()
    plt.savefig("processing.png")
    plt.close()
    try:
        pansy_plot.plot_raw(show_plot=False,fname="/tmp/raw.png")
    except Exception as exc:
        print("status_plot: raw plot unavailable: %s"%(exc))
    os.system("rsync -avz --bwlimit 1 /tmp/fit_data.h5 /tmp/raw.png /tmp/latest_meteor.png /tmp/latest_radiants.png /tmp/latest_hist.png status.png processing.png j@4.235.86.214:/var/www/html/pansy/")


while True:
    try:
        plot_status()
    except:
        import traceback
        traceback.print_exc()
    time.sleep(300)
