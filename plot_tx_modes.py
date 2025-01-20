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
data_dict = dmr.read(db[0], db[0]+24*3600*1000000, "id")
times=[]
for key in data_dict.keys():
    print((key, data_dict[key]))
    times.append(key)
times=n.array(times)
plt.plot(times,".")
plt.show()


plt.plot(n.diff(times),".")
plt.show()

        

    
    
    
if __name__ == "__main__":
    meteor_search()
    
