import pansy_config as pc 
import numpy as n
import matplotlib.pyplot as plt
import digital_rf as drf
import stuffr

dm = drf.DigitalMetadataReader(pc.gpslock_metadata_dir)
b = dm.get_bounds()

dd=dm.read(b[0],b[1])

tv=[]
holdover=[]

for k in dd.keys():
    holdover.append(dd[k]["holdover"])
    tv.append(stuffr.unix2date(k/1e6))

    #m=ax.pcolormesh(tvs,rvec,dB,vmin=0)
    #ax.set_xlabel("Date (UTC)")
    #ax.set_ylabel("Range (km)")
    #fig.autofmt_xdate()
fig,ax=plt.subplots(1,1)
ax.plot(tv,holdover)
ax.set_title("GPS log")
ax.set_xlabel("Time (UT)")
ax.set_ylabel("Holdover time (s)")
ax.autofmt_xdate()
plt.show()