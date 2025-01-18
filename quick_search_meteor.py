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

def meteor_search():

    d=drf.DigitalRFReader("/media/archive")
    print(d.get_bounds("ch000"))

    stm=pm.get_m_mode()
    ipp=stm["ipp_us"]
    beam_pos=stm["beam_pos_az_za"]
    n_beam=len(beam_pos)
    n_codes=len(stm["codes"])
    # get code vector
    codes=pm.get_vector(stm,ncodes=n_codes)

    metadata_dir = "/media/archive/metadata/tx"
    if not os.path.exists(metadata_dir):
        print("metadata directory doesn't exist. exiting")
        exit(0)

    dmr = drf.DigitalMetadataReader(metadata_dir)
    db = dmr.get_bounds()
    print(db)
    b=d.get_bounds("ch000")
    n_blocks=int(n.floor((db[1]-b[0])/(ipp*n_codes)))
    
    RTI = n.zeros([n_beam,n_codes,ipp],dtype=n.float32)

    N=20*1600
    for bi in range(n_blocks):
        #b=d.get_bounds("ch000")
        
        i0=bi*ipp*n_codes + b[0]
        i1=bi*ipp*n_codes + b[0] + ipp*n_codes + ipp

        if (i0 > b[0]) & (i1 < b[1]):
            data_dict = dmr.read(i0, i1, "id")
            for key in data_dict.keys():
                print((key, data_dict[key]))
        

    
    
    
if __name__ == "__main__":
    meteor_search()
    
