import numpy as n
import matplotlib.pyplot as plt

import digital_rf as drf
import sky_noise_cal as snc
import pansy_modes as pm
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import scipy.signal as ss
import stuffr
from datetime import datetime

mmode=pm.get_m_mode()
beam_pos=mmode["beam_pos_az_za"]

nlu=snc.noise_lookup()


#dm = drf.DigitalMetadataReader("/mnt/data/juha/pansy/metadata/xc2")
dm = drf.DigitalMetadataReader("/media/analysis/metadata/xc2")
b = dm.get_bounds()
print(b)

day0=n.ceil(b[0]/24/3600/1000000)+100
b0=day0*24*3600*1000000

n_days=int(n.floor((b[1]-b[0])/24/3600/1000000))
print(n_days)
beam_names=["Zenith","North","East","South","West"]
r0=None
r1=None
rdec=None
for di in range(n_days):
    print(di)

    nfloors=[]
    tms=[]
    max_noise=[]
    gdm_t=[]
    for i in range(5):
        gdm_t.append([])
    
    dd = dm.read(b0+di*24*3600*1000000, b0+(di+1)*24*3600*1000000)
    uts0=24*3600*n.floor((b0+di*24*3600*1000000)/1e6/24/3600)
    date_str = datetime.utcfromtimestamp(uts0).strftime("%Y-%m-%d")
#    print(len(dd.keys()))
    for k in dd.keys():
#        print(k)
#        print(dd[k].keys())
        chp=dd[k]["ch_pairs"]
        if r0 == None:
            r0 = dd[k]["r0"]
        if r1 == None:
            r1 = dd[k]["r1"]            
        if rdec == None:
            rdec = dd[k]["rdec"]    
            
        for i in range(5):
            gdm_t[i].append(nlu.get_noise_lookup(i,k/1e6))
 #           gdm_t[i].append(snc.get_noise(beam_pos[i][0],90-beam_pos[i][1],k/1e6))
        tms.append(k/1e6)
 #       print(chp[0:7])
        nfloor=n.zeros(5)
        max_noise.append(n.max(n.real(dd[k]["xc_arr"][0,0,:,:]),axis=0))
        for i in range(5):
            nfloor[i]=n.median(n.real(dd[k]["xc_arr"][0,i,:,:]))
        nfloors.append(nfloor)

    nfloors=n.array(nfloors)
    tms=n.array(tms)
    rti=n.array(max_noise)
    tms=tms.astype("datetime64[s]")
   # print(nfloors.shape)
    gdm_t=n.array(gdm_t)
    #   print(gdm_t.shape)
    #    fig,axes=plt.subplots(nrows=6,sharex=True,figsize=(10,6))
    fig=plt.figure(figsize=(6.5,5))
    gs=gridspec.GridSpec(nrows=4,ncols=1,height_ratios=[1,1,0.5,0.5],hspace=0)
    axes=[]
    for i in range(4):
        if i>0 and i<3:
            ax0=fig.add_subplot(gs[i],sharex=axes[0])
            axes.append(ax0)            
        else:
            ax0=fig.add_subplot(gs[i])            
            axes.append(ax0)

 #       if i<5:
#            ax0.set_
            
    scale=1
    axes[2].axis("off")
    axes[3].axis("off")    
    i=0
    nfloors0=n.max(gdm_t[i,:])*nfloors[:,i]/n.max(ss.medfilt(nfloors[:,i],5))
    if i==0:
        scale=n.max(gdm_t[i,:])/n.max(ss.medfilt(nfloors[:,i],5))
    axes[1].plot(tms,nfloors0,color="black",alpha=0.5)
    axes[1].plot(tms,gdm_t[i,:],color="black")
    axes[1].text(tms[int(len(tms)/4)],15e3,"$T_{\\mathrm{sys}}=%1.0f, %1.0f, %1.0f$ K (min, median, max)"%(n.min(nfloors0),n.median(nfloors0),n.max(nfloors0)))
    axes[1].set_ylabel(r"$T_{\mathrm{sys}}$ (K)")
    axr=axes[i].twinx()
    axes[1].set_ylim([3e3,20e3])
    axr.set_ylabel(beam_names[i])#Beam %d"%(i+1))
    axr.set_yticks([])
    axes[0].tick_params(labelbottom=False)
    axes[1].set_xlabel("Time (UTC)")
    axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

    rvec=n.arange(1600)*0.15
    #           ipp=1600
    #          rvec=n.arange(ipp)*drg
    #         ridx=n.where( (rvec >
            
    c=axes[0].pcolormesh(tms,rvec[r0:r1:rdec],
                         n.log10(scale*rti.T),vmin=3,vmax=5)
    cb=fig.colorbar(c,ax=axes[3],orientation="horizontal")#(axes[i].colorbar()
    axes[0].set_ylabel("Range (km)")
    cb.set_label("Noise temperature ($\log_{10}$ K)")
    axes[1].set_yticks([5e3,10e3,15e3])
    plt.title(date_str)
    fig.subplots_adjust(hspace=0)
    plt.savefig("figs/noise_temp.png",dpi=300,bbox_inches="tight")
    plt.show()
