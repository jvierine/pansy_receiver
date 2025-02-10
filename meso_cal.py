import numpy as n
import digital_rf as drf
import matplotlib.pyplot as plt
import stuffr
import h5py
import pansy_config as pc

# assume beam pointing, and use echoes to determine what the phase calibration is

dm = drf.DigitalMetadataReader("/media/analysis/metadata/xc2/")
b = dm.get_bounds()

dt=60*1000000
pprofs=[]
dprofs=[]
for i in range(5):
    pprofs.append([])
    dprofs.append([])
r0=0
r1=1
tvs=[]
cals=[]
start_idx=stuffr.date2unix(2025,1,30,0,0,0)*1000000
end_idx=b[1]#stuffr.date2unix(2025,1,30,12,0,0)*1000000

n_b=int((end_idx-start_idx)/dt)
ch_pairs=None
noise_powers=[]
for bi in range(n_b):
    data_dict = dm.read(start_idx+bi*dt, start_idx+bi*dt+dt, ("xc_arr","i0","i1","r0","r1","f0","f1","n_fft","ch_pairs"))
    for k in data_dict.keys():
        print(stuffr.unix2datestr(k/1e6))
        r0=data_dict[k]["r0"]
        r1=data_dict[k]["r1"]
        f0=data_dict[k]["f0"]
        f1=data_dict[k]["f1"]
        ch_pairs=data_dict[k]["ch_pairs"]
        n_fft=data_dict[k]["n_fft"]
        fvec=n.fft.fftshift(n.fft.fftfreq(n_fft,d=5*1600/1e6))[f0:f1]
        zidx=n.argmin(n.abs(fvec))
        tvs.append(stuffr.unix2date(data_dict[k]["i0"]/1e6))
        # zenith beam
        xc=data_dict[k]["xc_arr"][:,0,:,:]
        mean_pwr=n.sum(n.abs(data_dict[k]["xc_arr"][0:7,0,:,:]),axis=0)
        
        ri=n.argmax(mean_pwr[zidx])

        #mi,mj=n.unravel_index(n.argmax(mean_pwr),shape=mean_pwr.shape)
        noise_floor=n.median(mean_pwr)
        dcpwr=(mean_pwr[zidx,:]-noise_floor)/noise_floor
        gidx=n.where(dcpwr > 30)[0]
        
#        if len(gidx)>0:
 #           power_ratios=n.zeros(7,dtype=n.float32)
  #          for i in range(7):
   #             power_ratios[i]=
    #        noise_powers.append(power_ratios)
        
        for gi in gidx:
            cals.append(xc[:,zidx,gi])
            noise_powers.append(n.abs(xc[:,zidx,gi])/n.abs(xc[0,zidx,gi]))
noise_powers=n.array(noise_powers)
#for i in range(7):
  #  plt.plot(noise_powers[:,i],".")
 #   plt.show()
cals=n.array(cals)
ho=h5py.File("data/mesocal.h5","w")
ho["cals"]=n.mean(cals,axis=0)
ho["pwr"]=n.mean(noise_powers,axis=0)
ho.close()
exit(0)
# stop
print(cals.shape)
for i in range(7,cals.shape[1]):
    rho=n.mean(cals[:,i])
    plt.subplot(211)
    plt.plot(n.angle(cals[:,i]),".")
    plt.title(ch_pairs[i])
    plt.axhline(n.angle(rho),color="red")
    plt.subplot(212)
    plt.hist(n.angle(cals[:,i]),bins=50)
    plt.axvline(n.angle(rho),color="red")
    plt.show()
#for i in range(cals.shape)

