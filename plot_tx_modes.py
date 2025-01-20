import numpy as n
import pansy_detect as pd
import matplotlib.pyplot as plt
import pansy_modes as pm
import scipy.signal.windows as sw
import scipy.constants as c
import digital_rf as drf
import os
import stuffr
import time


metadata_dir = "/media/archive/metadata/tx"
if not os.path.exists(metadata_dir):
    print("metadata directory doesn't exist. exiting")
    exit(0)

dmr = drf.DigitalMetadataReader(metadata_dir)
db = dmr.get_bounds()
print(db)
i0=db[0]
n_hours=int(n.floor((db[1]-db[0])/1e6/3600))
for i in range(n_hours):
    data_dict = dmr.read(i0,i0+3600*1000000, "id")
    times=[]
    i0+=3600*1000000

    for key in data_dict.keys():
        times.append(key)
    times=n.array(times)
    plt.subplot(131)
    plt.plot(times,".")
    plt.subplot(132)
    plt.plot(n.diff(times)/1e6,".")
    plt.subplot(133)
    plt.plot(n.diff(times),".")
    plt.ylim([0,100000])
    plt.title(stuffr.unix2datestr(i0/1e6))
    plt.savefig("mode-%06d.png"%(i))
    print("mode-%06d.png"%(i))
    plt.close()


        

    
    
    
if __name__ == "__main__":
    meteor_search()
    
