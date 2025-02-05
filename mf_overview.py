import numpy as n
import matplotlib.pyplot as plt
import digital_rf as drf
import os
import stuffr
import time
import pansy_config as pc

os.system("mkdir -p /media/analysis/overview/")

def create_overview(i0,md,dt=15*60*1000000):
    fname="/media/analysis/overview/overview-%d.png"%(i0/1e6)
    if os.path.exists(fname):
        print("%s already exists"%(fname))
        return
    L=24*3600*1000000
    n_windows=int(n.floor(L/dt))
    rvec=n.arange(0,150)
    n_r=len(rvec)
    rti=n.zeros([n_windows,n_r],dtype=n.float32)
    n_avg=n.zeros([n_windows,n_r],dtype=n.float32)
    tv_hr=n.arange(n_windows)*dt/1e6/3600
    for i in range(n_windows):
        mi0=i0+i*dt
        mi1=i0+i*dt+dt
        print(stuffr.unix2datestr(mi0/1e6))
        data_dict = md.read(mi0, mi1, ("tx_pwr","max_range","tx_idxs","max_dopvel","max_snr","beam_pos_idx"))
        for k in data_dict.keys():
            ridx=int(n.round(data_dict[k]["max_range"][0]))
        #        print(ridx)
            rti[i,ridx]+=data_dict[k]["max_snr"][0]
            n_avg[i,ridx]+=1.0
    rti=rti/(n_avg+1.0)
    plt.title(stuffr.unix2datestr(i0/1e6))
    plt.pcolormesh(tv_hr,rvec,rti.T,vmin=4,vmax=10)
    cb=plt.colorbar()
    cb.set_label("SNR")
    plt.xlabel("Time (UTC hour)")
    plt.ylabel("Range (km)")
    plt.tight_layout()
    plt.ylim([50,150])
    fname="/media/analysis/overview/overview-%d.png"%(i0/1e6)
    print("saving %s"%(fname))
    plt.savefig("/media/analysis/overview/overview-%d.png"%(i0/1e6))
    plt.close()


while True:
    mf_metadata_dir=pc.mf_metadata_dir
    md = drf.DigitalMetadataReader(mf_metadata_dir)
    b = md.get_bounds()
    day0=int(n.floor(b[0]/24/3600/1e6))
    day1=int(n.floor(b[1]/24/3600/1e6))
    for day in range(day0,day1):
        i0=day0*24*3600*1000000
        create_overview(i0,md)

    time.sleep(1800)