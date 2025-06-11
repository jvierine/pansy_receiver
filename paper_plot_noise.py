import numpy as n
import matplotlib.pyplot as plt

import digital_rf as drf
import sky_noise_cal as snc
import pansy_modes as pm
mmode=pm.get_m_mode()
beam_pos=mmode["beam_pos_az_za"]

nlu=snc.noise_lookup()


dm = drf.DigitalMetadataReader("/mnt/data/juha/pansy/metadata/xc2")
b = dm.get_bounds()
print(b)


n_days=int(n.floor((b[1]-b[0])/24/3600/1000000))
print(n_days)
for di in range(n_days):
    print(di)

    nfloors=[]
    tms=[]
    gdm_t=[]
    for i in range(5):
        gdm_t.append([])
    
    dd = dm.read(b[0]+di*24*3600*1000000, b[0]+(di+1)*24*3600*1000000)
    print(len(dd.keys()))
    for k in dd.keys():
        print(k)
#        print(dd[k].keys())
        chp=dd[k]["ch_pairs"]
        for i in range(5):
            gdm_t[i].append(nlu.get_noise_lookup(i,k/1e6))
 #           gdm_t[i].append(snc.get_noise(beam_pos[i][0],90-beam_pos[i][1],k/1e6))
        tms.append(k/1e6)
 #       print(chp[0:7])
        nfloor=n.zeros(5)
        for i in range(5):
            nfloor[i]=n.median(n.real(dd[k]["xc_arr"][0,i,:,:]))
        nfloors.append(nfloor)

    nfloors=n.array(nfloors)
   # print(nfloors.shape)
    gdm_t=n.array(gdm_t)
    print(gdm_t.shape)

    for i in range(5):
        nfloors0=n.max(gdm_t[i,:])*nfloors[:,i]/n.max(nfloors[:,i])
        plt.plot(tms,nfloors0)
        plt.plot(tms,gdm_t[i,:])
        plt.show()
