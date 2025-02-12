import numpy as n
import pansy_config as pc
import matplotlib.pyplot as plt
import digital_rf as drf




dm = drf.DigitalMetadataReader(pc.simple_fit_metadata_dir)
b = dm.get_bounds()
dt=100000000
n_block=int((b[1]-b[0])/dt)
hgts=[]
t0s=[]
vgs=[]
for i in range(n_block):
    data=dm.read(b[0]+i*dt,b[0]+i*dt + dt)
    for k in data.keys():
#                    dout["r0"]=r0
 #           dout["v0"]=v0            
  #          dout["std"]=[eres,nres,ures]

        r0=data[k]["r0"]
        v0=data[k]["v0"]
        std=data[k]["std"]
        if n.max(std*1e3) < 500.0:
            hgts.append(r0[2])
            vgs.append(n.linalg.norm(v0))
            t0s.append(k/1e6)

plt.scatter(t0s,hgts,c=vgs,s=2)
plt.colorbar()
plt.show()

plt.plot(vgs,hgts,".")
plt.show()


