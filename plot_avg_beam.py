import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
import digital_rf as drf
import os
import stuffr


dm = drf.DigitalMetadataReader("/Users/j/src/pansy_test_data/metadata/simple_meteor_fit")
b = dm.get_bounds()
d=dm.read(b[0],b[1])

ew=[]
ns=[]
snr=[]
mfs=[]

for k in d.keys():
    print(k)
    vg=n.linalg.norm(d[k]["v0"])
    print(vg)
    if n.max(d[k]["std"]) < 400.0 and vg > 10: 
        gidx=n.where(d[k]["mfs"]/21>0.9)[0]
        ew=n.concatenate((ew,d[k]["ew"][gidx]))
        ns=n.concatenate((ns,d[k]["ns"][gidx]))
        snr=n.concatenate((snr,d[k]["snr"][gidx]))
        mfs=n.concatenate((mfs,d[k]["mfs"][gidx]))

plt.subplot(121)
#plt.scatter(ew,ns,c=10.0*n.log10(snr),s=1,vmin=10,vmax=20)
plt.scatter(ew,ns,c=10.0*n.log10(snr),vmin=10,vmax=20,s=1,cmap="gist_yarg",alpha=0.1)#,vmin=10,vmax=20)

plt.xlabel("East-West (km)")
plt.ylabel("North-South (km)")
#cb=plt.colorbar()
#cb.set_label("SNR (dB)")
#plt.show()
plt.subplot(122)

plt.hist2d(ew,ns,bins=[n.linspace(-30,30,num=250),n.linspace(-30,30,num=250)],weights=mfs)
plt.title("Histogram of detections")
plt.xlabel("East-West (km)")
plt.ylabel("North-South (km)")
plt.tight_layout()
plt.show()
