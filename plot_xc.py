import numpy as n
import digital_rf as drf
import matplotlib.pyplot as plt
import stuffr

import pansy_config as pc



dm = drf.DigitalMetadataReader(pc.xc_metadata_dir)
b = dm.get_bounds()

dt=60*1000000
n_b=int((b[1]-b[0])/dt)
pprofs=[]
dprofs=[]
for i in range(5):
    pprofs.append([])
    dprofs.append([])
r0=0
r1=1
tvs=[]
for bi in range(n_b):
    data_dict = dm.read(b[0]+bi*dt, b[0]+bi*dt+dt, ("xc_arr","i0","i1","r0","r1","f0","f1","n_fft"))
    for k in data_dict.keys():
        r0=data_dict[k]["r0"]
        r1=data_dict[k]["r1"]
        f0=data_dict[k]["f0"]
        f1=data_dict[k]["f1"]
        n_fft=data_dict[k]["n_fft"]
        fvec=n.fft.fftshift(n.fft.fftfreq(n_fft,d=5*1600/1e6))[f0:f1]
        tvs.append(stuffr.unix2date(data_dict[k]["i0"]/1e6))
        for i in range(5):

            mean_pwr=n.sum(n.abs(data_dict[k]["xc_arr"][0:7,i,:,:]),axis=0)
            plt.pcolormesh(10.0*n.log10(mean_pwr))
            plt.colorbar()
            plt.show()
            dop_idx=n.argmax(mean_pwr,axis=0)
            pprofs[i].append(n.max(mean_pwr,axis=0))
            dprofs[i].append(fvec[dop_idx])#n.max(n.sum(n.abs(data_dict[k]["xc_arr"][0:7,i,:,:]),axis=0),axis=0))

pprofs=n.array(pprofs)
dprofs=n.array(dprofs)

print(pprofs.shape)
#print(pprofs[i,:,:])
rvec=n.arange(1600)*0.15
rvec=rvec[r0:r1]
fig, axs = plt.subplots(nrows=5,ncols=1)
for i in range(5):
    ax=axs[i]    
    if i == 0:
        m=ax.pcolormesh(tvs,rvec,dprofs[i,:,:].T,cmap="seismic",vmin=-5,vmax=5)
    else:
        m=ax.pcolormesh(tvs,rvec,dprofs[i,:,:].T,cmap="seismic")

    fig.autofmt_xdate()
    cb=fig.colorbar(m,ax=ax)
plt.show()
fig, axs = plt.subplots(nrows=5,ncols=1)
for i in range(5):
    ax=axs[i]    
    m=ax.pcolormesh(tvs,rvec,10.0*n.log10(pprofs[i,:,:].T))
    fig.autofmt_xdate()
    cb=fig.colorbar(m,ax=ax)
plt.show()
#    ax=axs[1]    
 #   m=ax.pcolormesh(tvs,rvec,10.0*n.log10(pprofs[i,:,:].T))
  ## cb=fig.colorbar(m,ax=ax)

#cb.set_label("Doppler (km/s)")
#ax.set_xlabel("Time (unix)")
#plt.show()

#plt.pcolormesh(tvs,rvec,10.0*n.log10(pprofs[i,:,:].T))
#plt.colorbar()
#plt.show()
    


