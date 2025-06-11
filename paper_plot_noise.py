import numpy as n
import matplotlib.pyplot as plt

import digital_rf as drf

dm = drf.DigitalMetadataReader("/mnt/data/juha/pansy/metadata/xc2")
b = dm.get_bounds()
print(b)

nfloors=[]
n_days=int(n.floor((b[1]-b[0])/24/3600/1000000))
print(n_days)
tms=[]
for di in range(n_days):
    dd = dm.read(b[0]+di*24*3600*1000000, b[0]+(di+1)*24*3600*1000000)
    for k in dd.keys():
        print(dd[k].keys())
        chp=dd[k]["ch_pairs"]
        tms.append(k/1e6)
        print(chp[0])
        nfloor=n.zeros(5)
        for i in range(5):
            nfloor[i]=n.median(n.real(dd[k]["xc_arr"][0,i,:,:]))
        nfloors.append(nfloor)

nfloors=n.array(nfloors)
print(nfloors.shape)
for i in range(5):
    plt.plot(tms,nfloors[:,i])
plt.show()
