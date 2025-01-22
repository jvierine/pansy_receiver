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
import scipy.fftpack as fp

class range_doppler_search:
    def __init__(self,txlen=130,
                 rg=n.arange(400,950,2,dtype=n.int64)):
        self.idx=n.arange(130,dtype=n.int64)
        self.n_rg=len(rg)
        self.txlen=txlen
        self.idx_mat=n.zeros([self.n_rg,self.txlen],dtype=n.int64)
        self.idx=n.arange(self.txlen,dtype=n.int64)
        for ri in range(self.n_rg):
            self.idx_mat[ri,:]=self.idx+rg[ri]

    def mf(self,z,z_tx):
        """
        z = echo
        z_tx = transmit code
        """
        plt.subplot(121)
        plt.plot(z.real)
        plt.plot(z.imag)
        plt.subplot(122)
        plt.plot(z_tx.real)
        plt.plot(z_tx.imag)
        plt.show()
        z_tx=n.conj(z_tx)
        # decode each range gate
        Z=z[self.idx_mat]*z_tx[None,:]

        plt.pcolormesh(n.real(Z))
        plt.show()
        ZF=fp.fft(Z,axis=0)
        plt.pcolormesh(n.abs(ZF)**2.0)
        plt.show()
        
        

    

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

    start_idx=db[1]-1000000
    n_blocks=int(n.floor((db[1]-start_idx)/(ipp*n_codes)))
    
    RTI = n.zeros([n_beam,n_codes,ipp],dtype=n.float32)

    rds=range_doppler_search()
    N=20*1600
    for bi in range(n_blocks):

        
        i0=bi*ipp*n_codes + start_idx
        i1=bi*ipp*n_codes + start_idx + ipp*n_codes + ipp

        if (i0 > b[0]) & (i1 < b[1]):
            data_dict = dmr.read(i0, i1, "id")
            for key in data_dict.keys():
                print((key, data_dict[key]))
                z=d.read_vector_c81d(key,1600*20,"ch000")
                z_tx=d.read_vector_c81d(key,1600*20,"ch007")                
                for ti in range(20):
                    rds.mf(z[(0+ti*1600):(1600+ti*1600)],z_tx[(0+ti*1600):(rds.txlen+ti*1600)])
        
                plt.plot(z_tx.real)
                plt.plot(z_tx.imag)
                plt.show()
        

    
    
    
if __name__ == "__main__":
    meteor_search()
    
