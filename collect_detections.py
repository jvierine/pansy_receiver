import numpy as n
import digital_rf as drf
import matplotlib.pyplot as plt
import stuffr

import pansy_config as pc

det_md_dir = pc.detections_metadata_dir#"/media/archive/metadata/detections"

dm = drf.DigitalMetadataReader(det_md_dir)
b = dm.get_bounds()

#data_dict = dm.read(b[1]-24*3600*1000000, b[1], ("xhat","tx_idx","snr"))
data_dict = dm.read(b[0], b[1], ("xhat","tx_idx","snr"))

tv=[]
v0s=[]
r0s=[]
durs=[]
snrs=[]
for k in data_dict.keys():
    print(k)
    data=data_dict[k]
    xhat=data["xhat"]
    snr=data["snr"]
    snrs.append(n.max(snr))
    dur=(n.max(data["tx_idx"])-n.min(data["tx_idx"]))/1e6
    r0=xhat[0]
    v0=xhat[1]
    t0=k
    tv.append(n.datetime64(stuffr.unix2date(t0/1e6)))
    v0s.append(v0)
    r0s.append(r0)
    durs.append(dur)
v0s=-1*n.array(v0s)
plt.subplot(131)
plt.hist(10.0*n.log10(snrs),bins=50)
plt.xlabel("Signal-to-noise ratio (dB)")
plt.subplot(132)

plt.hist(durs,bins=50)
plt.xlabel("Event duration (s)")
plt.subplot(133)

plt.hist(r0s,bins=50)
plt.xlabel("Height (km)")
plt.tight_layout()
plt.show()

tv=n.array(tv)
#gidx=n.where(v0s>10)[0]
gidx=n.arange(len(tv),dtype=n.int64)

r0s=n.array(r0s)
fig, ax = plt.subplots()
m=ax.scatter(tv[gidx],r0s[gidx],c=v0s[gidx],s=1,cmap="turbo",vmax=72,vmin=-5)
fig.autofmt_xdate()
ax.set_title("%d meteors"%(len(tv)))
cb=fig.colorbar(m,ax=ax)
cb.set_label("Doppler (km/s)")
ax.set_xlabel("Time (unix)")
ax.set_ylabel("Range (km)")
plt.show()

fig, ax = plt.subplots()
m=ax.scatter(v0s,r0s,c=10.0*n.log10(snrs),s=1,cmap="turbo",vmin=10,vmax=42)
#fig.autofmt_xdate()
ax.set_title("%d meteors"%(len(tv)))
cb=fig.colorbar(m,ax=ax)
cb.set_label("SNR (dB)")
ax.set_xlabel("Radial velocity (km/s)")
plt.show()


if False:
    plt.scatter(tv,r0s,c=v0s,vmin=0,vmax=72,s=1,cmap="turbo")
    plt.title("%d meteors"%(len(tv)))
    cb=plt.colorbar()
    cb.set_label("Doppler (km/s)")
    plt.xlabel("Time (unix)")
    plt.show()

    plt.scatter(tv,r0s,c=10.0*n.log10(snrs),s=1,cmap="turbo")
    plt.title("%d meteors"%(len(tv)))
    cb=plt.colorbar()
    cb.set_label("SNR (dB)")
    plt.xlabel("Time (unix)")
    plt.show()

    plt.scatter(tv,v0s,c=r0s,vmin=80,vmax=140,s=1,cmap="turbo")
    plt.title("%d meteors"%(len(tv)))
    cb=plt.colorbar()
    cb.set_label("Height (km)")
    plt.xlabel("Time (unix)")
    plt.show()
