import numpy as n
import digital_rf as drf
import matplotlib.pyplot as plt

det_md_dir = "/media/archive/metadata/detections"

dm = drf.DigitalMetadataReader(det_md_dir)
b = dm.get_bounds()

data_dict = dm.read(b[0], b[0]+3600*1000000, ("xhat"))

tv=[]
v0s=[]
r0s=[]
for k in data_dict.keys():
    print(k)
    data=data_dict[k]
    print(data)
    r0=data[0]
    v0=data[1]
    t0=k
    tv.append(t0/1e6)
    v0s.append(v0)
    r0s.append(r0)

plt.scatter(tv,r0s,c=v0s,vmin=0,vmax=73)
cb=plt.colorbar()
cb.set_label("Doppler (km/s)")
plt.xlabel("Time (unix)")
plt.show()