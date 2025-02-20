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


plt.plot(tv,holdover)
plt.title("GPS log")
plt.xlabel("Time (UT)")
plt.ylabel("Holdover time (s)")
plt.show()