import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
import digital_rf as drf
import os
import stuffr
import time
import pansy_config as pc
import cluster_mf as cmf
import traceback
import scipy.fftpack as fp
fft=fp.fft
class range_doppler_search:
    def __init__(self,
                 txlen,
                 echolen,
                 n_channels,
                 interp=20):
        self.txlen=txlen*interp
        
        self.fdec=8*interp
        self.fftlen=int(self.txlen*2/self.fdec)
        self.echolen=echolen*interp
        
        self.n_channels=n_channels
        self.interp=interp
        self.n_rg=self.echolen-self.txlen
        self.rg=n.arange(self.n_rg,dtype=n.int64)
        self.idx_mat=n.zeros([self.n_rg,self.txlen],dtype=n.int64)
        self.idx=n.arange(self.txlen,dtype=n.int64)
        for ri in range(self.n_rg):
            self.idx_mat[ri,:]=self.idx+self.rg[ri]
            
        self.z_rx=n.zeros([self.n_channels,self.echolen],dtype=n.complex64)

        drg=c.c/2/1e6/1e3
        self.rangev=self.rg*drg
        self.frad=47.5e6
        self.fvec=n.fft.fftshift(n.fft.fftfreq(self.fftlen,d=self.fdec/(interp*1e6)))
        self.dopv=self.fvec*c.c/2.0/self.frad

        
        
    def decim(self,Z):
        new_width=int(self.txlen/self.fdec)
        Z2=n.zeros([self.n_rg,new_width],dtype=n.complex64)
        idx=n.arange(new_width,dtype=n.int64)
        idx2=n.arange(new_width,dtype=n.int64)*self.fdec

        for i in range(self.fdec):
            Z2[:,idx]+=Z[:,idx2+i]
        return(Z2)

    def mf(self,tx,echo,debug=False):
            
        txi=n.repeat(tx,self.interp)
#        txi=tx#,self.interp)        
        for chi in range(self.n_channels):
            self.z_rx[chi,:]=n.repeat(echo[chi,:],self.interp)
#            self.z_rx[chi,:]=echo[chi,:]
        
        z_tx=n.conj(txi)
        
        MF=n.zeros([self.n_rg,self.fftlen],dtype=n.float32)
        for chi in range(self.n_channels):
            Z=self.z_rx[chi,self.idx_mat]*txi[None,:]
#            plt.pcolormesh(Z.real)
#            plt.show()
            
            # decimate
            ZD=self.decim(Z)
            ZF=n.fft.fftshift(fft(ZD,self.fftlen,axis=1),axes=1)
            MF+=ZF.real**2.0 + ZF.imag**2.0
        noise_floor=n.median(MF)
        pprof=n.max(MF,axis=1)
        peak_dopv=self.dopv[n.argmax(MF,axis=1)]
        
        if False:
            plt.pcolormesh(self.dopv,self.rangev,n.abs(ZF)**2.0)
            plt.show()
        return(MF,pprof,peak_dopv,noise_floor)





mddir="/home/j/src/pansy_receiver/test_data/metadata/cut"
dm = drf.DigitalMetadataReader(mddir)
b = dm.get_bounds()
dt=10000000
n_block=int((b[1]-b[0])/dt)
for bi in range(n_block):
    data=dm.read(b[0]+bi*dt,b[0]+bi*dt+dt)
    for k in data.keys():
        z_rx=n.array(data[k]["zrx_echoes_re"],dtype=n.complex64)+n.array(data[k]["zrx_echoes_im"],dtype=n.complex64)*1j
        z_tx=n.array(data[k]["ztx_pulses_re"],dtype=n.complex64)+n.array(data[k]["ztx_pulses_im"],dtype=n.complex64)*1j
        delays=data[k]["delays"]
        beam_id=data[k]["beam_id"]
        tx_idx=data[k]["tx_idx"]
        

#@        print(z_rx.shape)
  #      plt.pcolormesh(n.angle(z_rx[:,0,:]*n.conj(z_rx[:,1,:])).T,cmap="hsv")
   #     plt.colorbar()
    #    plt.show()

        interp=5
        txlen=z_tx.shape[1]
        echolen=z_rx.shape[2]        

        rds=range_doppler_search(txlen=txlen,echolen=echolen,interp=interp,n_channels=z_rx.shape[1])
        RTI=n.zeros([z_rx.shape[0],rds.n_rg],dtype=n.float32)

   #     for i in range(7):
  #          plt.pcolormesh(z_rx[:,i,:].real)
 #           plt.colorbar()
#            plt.show()
        print(z_rx.shape)
        for ti in range(z_rx.shape[0]):
            MF,pprof,peak_dopv,noise_floor=rds.mf(z_tx[ti,:],z_rx[ti,:,:])
            RTI[ti,:]=pprof
        plt.pcolormesh(RTI.T)
        plt.colorbar()
        plt.show()
            

        
        
        print(data[k])


