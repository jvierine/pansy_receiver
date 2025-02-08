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
for i in range(5):
    pprofs.append([])
r0=0
r1=1
tvs=[]
for bi in range(n_b):
    data_dict = dm.read(b[0]+bi*dt, b[0]+bi*dt+dt, ("xc_arr","i0","i1","r0","r1"))
    for k in data_dict.keys():
        r0=data_dict[k]["r0"]
        r1=data_dict[k]["r1"]
        tvs.append([data_dict[k]["i0"],data_dict[k]["i1"]])
        for i in range(5):
            pprofs[i].append(n.max(n.sum(n.abs(data_dict[k]["xc_arr"][i,0:7,:,:]),axis=0),axis=0))
pprofs=n.array(pprofs)
print(pprofs.shape)
#print(pprofs[i,:,:])
rvec=n.arange(1600)*0.15
rvec=rvec[r0:(r1+1)]
for i in range(5):
    plt.pcolormesh(tv,rvec,10.0*n.log10(pprofs[i,:,:].T))
    plt.colorbar()
    plt.show()
    


