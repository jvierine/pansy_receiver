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

for bi in range(n_b):
    data_dict = dm.read(b[0]+bi*dt, b[0]+bi*dt+dt, ("xc_arr","i0","i1"))
    for k in data_dict.keys():
        for i in range(5):
            pprofs[i].append(n.max(data_dict[k]["xc_arr"][i,0,:,:],axis=0))
pprofs=n.array(pprofs)
for i in range(5):
    plt.pcolormesh(pprofs[i])
    plt.colorbar()
    plt.show()
    


