import numpy as n
import digital_rf as drf
import matplotlib.pyplot as plt
import stuffr

det_md_dir = "/media/archive/metadata/detections"

dm = drf.DigitalMetadataReader(det_md_dir)
b = dm.get_bounds()

data_dict = dm.read(b[0], b[1], ("xhat","tx_idx"))

tv=[]
v0s=[]
r0s=[]
durs=[]
for k in data_dict.keys():
    print(k)
    data=data_dict[k]
    xhat=data["xhat"]
    dur=(n.max(data["tx_idx"])-n.min(data["tx_idx"]))/1e6
    r0=xhat[0]
    v0=xhat[1]
    t0=k
    tv.append(stuffr.unix2date(t0/1e6))
    v0s.append(v0)
    r0s.append(r0)
    durs.append(dur)

plt.scatter(tv,r0s,c=durs,s=1,cmap="turbo")
plt.title("%d meteors"%(len(tv)))
cb=plt.colorbar()
cb.set_label("Duration (s)")
plt.xlabel("Time (unix)")
plt.show()


plt.scatter(tv,r0s,c=v0s,vmin=-73,vmax=0,s=1,cmap="turbo")
plt.title("%d meteors"%(len(tv)))
cb=plt.colorbar()
cb.set_label("Doppler (km/s)")
plt.xlabel("Time (unix)")
plt.show()

plt.scatter(tv,v0s,c=r0s,vmin=80,vmax=140,s=1,cmap="turbo")
plt.title("%d meteors"%(len(tv)))
cb=plt.colorbar()
cb.set_label("Height (km)")
plt.xlabel("Time (unix)")
plt.show()