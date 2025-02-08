import numpy as n
import digital_rf as drf
import matplotlib.pyplot as plt
import stuffr

import pansy_config as pc



dm = drf.DigitalMetadataReader(pc.xc_metadata_dir)
b = dm.get_bounds()

dt=60*1000000
n_b=int((b[1]-b[0])/dt)
pprof0=[]
pprof1=[]
pprof2=[]
pprof3=[]
pprof4=[]

for bi in range(n_b):
    data_dict = dm.read(b[0]+bi*dt, b[0]+bi*dt+dt, ("xc_arr","i0","i1"))
    for k in data_dict.keys():
        plt.plot(n.max(data_dict[k]["xc_arr"][0,0,:,:],axis=0))
        plt.show()
