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
import itertools
import pansy_interferometry as pint
import h5py

# each antenna needs to be multiplied with these phases (n.exp(-1j*phasecal))
h=h5py.File("data/phases.h5","r")
phasecal=h["zenith"][()]
h.close()
# each antenna needs to be multiplied with these amplitudes
h=h5py.File("data/mesocal.h5","r")
amps=h["pwr"][()]
amps=1/n.sqrt(amps[0:7])

u,v,w=pint.uv_coverage(N=500,max_zenith_angle=20.0)
antpos=pint.get_antpos()
ch_pairs=n.array(list(itertools.combinations(n.arange(7),2)))
dmat=pint.pair_mat(ch_pairs,antpos)

fft=fp.fft

class range_doppler_search:
    def __init__(self,
                 txlen,
                 echolen,
                 n_channels,
                 interp=2):
        self.txlen=txlen*interp
        
        self.fdec=8*interp
        self.fftlen=512#int(self.txlen*32/self.fdec)
        if int(self.txlen*interp/self.fdec)>self.fftlen:
            self.fftlen=int(self.txlen*interp/self.fdec)
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

        # interferometry
        self.ch_pairs=list(itertools.combinations(n.arange(n_channels),2))
        self.n_pairs=len(self.ch_pairs)
        
        
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
        #print("zrx",self.z_rx.shape)
        z_tx=n.conj(txi)
        
        MF=n.zeros([self.n_rg,self.fftlen],dtype=n.float32)
        XC=n.zeros([self.n_pairs,self.n_rg,self.fftlen],dtype=n.complex64)
        XC2=n.zeros([self.n_pairs,self.n_rg],dtype=n.complex64)
        for chi in range(self.n_channels):
            Z=self.z_rx[chi,self.idx_mat]*txi[None,:]            
            # decimate
            ZD=self.decim(Z)
            ZF=n.fft.fftshift(fft(ZD,self.fftlen,axis=1),axes=1)
            MF+=ZF.real**2.0 + ZF.imag**2.0
     #   plt.pcolormesh(MF)
      #  plt.show()

        for chip in range(self.n_pairs):
            Z1=self.z_rx[self.ch_pairs[chip][0],self.idx_mat]*txi[None,:]            
            Z2=self.z_rx[self.ch_pairs[chip][1],self.idx_mat]*txi[None,:]            
            # decimate
            ZD1=self.decim(Z1)
            ZD2=self.decim(Z2)
            ZF1=n.fft.fftshift(fft(ZD1,self.fftlen,axis=1),axes=1)
            ZF2=n.fft.fftshift(fft(ZD2,self.fftlen,axis=1),axes=1)
            XC[chip,:,:]=ZF1*n.conj(ZF2)

        noise_floor=n.median(MF)
        pprof=n.max(MF,axis=1)
        dop_idx=n.argmax(MF,axis=1)

        # collect xc for peak doppler bin
        for i in range(self.n_rg):
            XC2[:,i]=XC[:,i,dop_idx[i]]
        
        peak_dopv=self.dopv[n.argmax(MF,axis=1)]
        
        return(MF,pprof,peak_dopv,noise_floor,XC2)

#mddir="../pansy_test_data/metadata/cut"
dm = drf.DigitalMetadataReader(pc.cut_metadata_dir)
b = dm.get_bounds()
dt=10000000
n_block=int((b[1]-b[0])/dt)
os.system("mkdir -p caldata")
for bi in range(n_block):
    data=dm.read(b[0]+bi*dt,b[0]+bi*dt+dt)
    for k in data.keys():
        z_rx=n.array(data[k]["zrx_echoes_re"],dtype=n.complex64)+n.array(data[k]["zrx_echoes_im"],dtype=n.complex64)*1j
        z_tx=n.array(data[k]["ztx_pulses_re"],dtype=n.complex64)+n.array(data[k]["ztx_pulses_im"],dtype=n.complex64)*1j
        delays=data[k]["delays"]
        beam_id=data[k]["beam_id"]
        tx_idx=data[k]["tx_idx"]
        c_snr=data[k]["c_snr"]
        #print(n.max(c_snr))
        if n.max(c_snr) < 100:# or (0 not in beam_id):
#            print(n.max(c_snr))
            continue
#        print("mf")
        interp=1
        txlen=z_tx.shape[1]
        echolen=z_rx.shape[2]        
        rds=range_doppler_search(txlen=txlen,echolen=echolen,interp=interp,n_channels=z_rx.shape[1])

        RTI=n.zeros([z_rx.shape[0],rds.n_rg],dtype=n.float32)
        XCT=n.zeros([rds.n_pairs,z_rx.shape[0],rds.n_rg],dtype=n.complex64)
        xct=n.zeros([rds.n_pairs,z_rx.shape[0]],dtype=n.complex64)

        peak_rg=[]
        peak_dop=[]
        drg=c.c/1e6/2/1e3
        snr=[]
        for ti in range(z_rx.shape[0]):
            MF,pprof,peak_dopv,noise_floor,XC=rds.mf(z_tx[ti,:],z_rx[ti,:,:])
            max_rg=n.argmax(pprof)
            peak_rg.append((delays[ti]+max_rg/interp)*drg)
            peak_dop.append(peak_dopv[max_rg])
            RTI[ti,:]=(pprof-noise_floor)/noise_floor
            XCT[:,ti,:]=XC
            xct[:,ti]=XC[:,max_rg]
            snr.append((pprof[max_rg]-noise_floor)/noise_floor)
        snr=n.array(snr)
        peak_rg=n.array(peak_rg)
        peak_dop=n.array(peak_dop)

        for bi in range(5):
            idx=n.where( (snr>7) & (beam_id == bi) & (n.abs(peak_dop) > 10e3) )[0]        
            if len(idx)>5:
                fidx0=tx_idx[idx[0]]
                ofname="caldata/meteor-%d-%d.h5"%(bi,fidx0)
                print(ofname)
                ho=h5py.File(ofname,"w")
                ho["xc"]=xct[:,idx]
                ho["bi"]=bi
                ho.close()
        #xct[:,idx]
        #pc.meteor_cal_metadata_dir

        if False:
            #for i in range(rds.n_pairs):
            plt.subplot(311)
            plt.plot(tx_idx[idx],peak_rg[idx],".")
            plt.title(n.unique(beam_id[idx]))
            plt.xlim([n.min(tx_idx),n.max(tx_idx)])
            plt.subplot(312)
            plt.plot(tx_idx[idx],peak_dop[idx],".")
            plt.xlim([n.min(tx_idx),n.max(tx_idx)])
            plt.subplot(313)
            plt.scatter(tx_idx[idx],n.angle(xct[0,idx]),c=beam_id[idx],s=1,cmap="turbo")
            plt.xlim([n.min(tx_idx),n.max(tx_idx)])
            plt.show()

            mu=[]
            mv=[]

            for ii in idx:
                print(ii)
                xc=xct[:,ii]
                xc=n.exp(1j*n.angle(xc))
                M=pint.mf(xc,dmat,u,v,w)
                mi,mj=n.unravel_index(n.argmax(M),shape=M.shape)
                mu.append(u[mi,mj])
                mv.append(v[mi,mj])

            mu=n.array(mu)
            mv=n.array(mv)
            mw=n.sqrt(1-mu**2.0-mv**2.0)
    #        peak_rg
            plt.subplot(131)
            plt.plot((tx_idx[idx]-tx_idx[0])/1e6,peak_rg[idx]*mv,".")
            plt.xlabel("Time (s)")
            plt.ylabel("North-South (km)")
            plt.subplot(132)
            plt.plot((tx_idx[idx]-tx_idx[0])/1e6,peak_rg[idx]*mu,".")
            plt.xlabel("Time (s)")
            plt.ylabel("East-West (km)")
            plt.subplot(133)
            plt.plot((tx_idx[idx]-tx_idx[0])/1e6,peak_rg[idx]*mw,".")
            plt.xlabel("Time (s)")
            plt.ylabel("Up (km)")
            plt.show()


