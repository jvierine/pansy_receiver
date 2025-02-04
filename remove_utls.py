import numpy as n
import digital_rf as drf
import os
import stuffr
import time


 # this is where the data is
d=drf.DigitalRFReader("/media/archive/")
# tx channel bounds
b=d.get_bounds("ch007")
print(d._get_file_list(b[0],b[0]+100,1000000,3600,1000))

# transmit metadata
metadata_dir = "/media/archive/metadata/tx"
dmr = drf.DigitalMetadataReader(metadata_dir)
db = dmr.get_bounds()

dt=1000000

s0=int(n.max((n.floor(b[0]/dt),n.floor(db[0]/dt))))

print(s0)
s1=int(n.floor(db[1]/dt))

for i in range(s0,s1):
    #    print(i)
    data_dict = dmr.read(i*dt-1600*20, i*dt+dt, "id")
    print("%d n_ipp %d"%(i,len(data_dict.keys())))
    os.system("ls /media/archive/ch000/")
    #for k in data_dict.keys():
     #   print(k)
