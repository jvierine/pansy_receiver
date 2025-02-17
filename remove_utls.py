import numpy as n
import digital_rf as drf
import os
import stuffr
import time
import h5py

import pansy_config as pc

def clean_utls():
    # this is where the data is
    d=drf.DigitalRFReader("/media/archive/")
    # tx channel bounds
    b=d.get_bounds("ch007")

    ch=["ch000","ch001","ch002","ch003","ch004","ch005","ch006","ch007"]

    # transmit metadata
    metadata_dir=pc.tx_metadata_dir
    #metadata_dir = "/media/archive/metadata/tx"
    dmr = drf.DigitalMetadataReader(metadata_dir)
    db = dmr.get_bounds()

    dt=1000000

    start_idx=b[0]
    if os.path.exists("/tmp/last_rem.h5"):
        ho=h5py.File("/tmp/last_rem.h5","r")
        start_idx=ho["last_rem"][()]
        ho.close()

    s0=int(n.max((n.floor(start_idx/dt),n.floor(db[0]/dt))))

    #print(s0)
    s1=int(n.floor(db[1]/dt))

    for i in range(s0,s1):
        #    print(i)
        data_dict = dmr.read(i*dt-1600*20, i*dt+dt, "id")
    #    print("%d n_ipp %d"%(i,len(data_dict.keys())))
        n_ipp=len(data_dict.keys())
        if n_ipp == 0:
            print("no ipps %s. deleting."%(stuffr.unix2datestr(i)))
            fl=d._get_file_list(i*dt+10,i*dt+1000,1000000,3600,1000)
            for f in fl:
                for c in ch:
                    fname="%s/%s/%s"%("/media/archive/",c,f)
    #                print(fname)
                    os.system("rm -f %s"%(fname))
    ho=h5py.File("/tmp/last_rem.h5","w")
    ho["last_del"]=s1*dt+dt
    ho.close()
        #os.system("ls /media/archive/ch000/")
        #for k in data_dict.keys():
        #   print(k)

if __name__ == "__main__":
    while True:
        clean_utls()
        # sleep for a long time!
        time.sleep(12*3600)