import numpy as n
import matplotlib.pyplot as plt

import digital_rf as drf
import sky_noise_cal as snc
import pansy_modes as pm
mmode=pm.get_m_mode()
beam_pos=mmode["beam_pos_az_za"]

dm = drf.DigitalMetadataReader("/mnt/data/juha/pansy/metadata/xc2")
b = dm.get_bounds()
print(b)

nfloors=[]
n_days=int(n.floor((b[1]-b[0])/24/3600/1000000))
print(n_days)
tms=[]
gdm_t=[]
for i in range(5):
    gdm_t.append([])
for di in range(n_days):
    dd = dm.read(b[0]+di*24*3600*1000000, b[0]+(di+1)*24*3600*1000000)
    for k in dd.keys():
        print(dd[k].keys())
        chp=dd[k]["ch_pairs"]
        for i in range(5):
            gdm_t[i].append(get_noise(beam_pos[i][0],90-beam_pos[i][1],k/1e6))
        tms.append(k/1e6)
        print(chp[0])
        nfloor=n.zeros(5)
        for i in range(5):
            nfloor[i]=n.median(n.real(dd[k]["xc_arr"][0,i,:,:]))
        nfloors.append(nfloor)

nfloors=n.array(nfloors)
print(nfloors.shape)
gdm_t=n.array(gdm_t)
print(gdm_t.shape)
for i in range(5):
#    nfloors=n.max(gdm_t)*nfloors/n.max(nfloors)
    plt.plot(tms,nfloors[:,i])

plt.show()
