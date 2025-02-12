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

import meteor_fit as metfit

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


# each antenna needs to be multiplied with these phases (n.exp(-1j*phasecal))
h=h5py.File("data/phases.h5","r")
phasecal=h["phasecal"][()]
h.close()

u,v,w=pint.uv_coverage(N=500,max_zenith_angle=15.0)
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
        # interpolate transmit pulse
        txi=n.repeat(tx,self.interp)
        # interpolate echo
        for chi in range(self.n_channels):
            self.z_rx[chi,:]=n.repeat(echo[chi,:],self.interp)
        # conjugate transmit pulse sample
        z_tx=n.conj(txi)

        # match function output
        MF=n.zeros([self.n_rg,self.fftlen],dtype=n.float32)
        # cross-spectrum
        XC=n.zeros([self.n_pairs,self.n_rg,self.fftlen],dtype=n.complex64)
        
        # cross-spectrum at peak doppler and range-gate
        #XC2=n.zeros([self.n_pairs,self.n_rg],dtype=n.complex64)
        for chi in range(self.n_channels):
            Z=self.z_rx[chi,self.idx_mat]*txi[None,:]            
            # decimate
            ZD=self.decim(Z)
            ZF=n.fft.fftshift(fft(ZD,self.fftlen,axis=1),axes=1)
            MF+=ZF.real**2.0 + ZF.imag**2.0
        
            
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
        rgmax=n.argmax(pprof)
        xc=XC[:,rgmax,dop_idx[rgmax]]
        
        peak_dopv=self.dopv[n.argmax(MF,axis=1)]
        
        return(MF,pprof,peak_dopv,noise_floor,xc,rgmax)

def process_cut(data,
                dmw,
                interp=4,
                plot=True):


    # read phase calibration
    ho=h5py.File("data/phases.h5","r")
    phcal=ho["phasecal"][()]
    ho.close()

    # create a uniform grid of u,v,w values
    u,v,w=pint.uv_coverage(N=500,max_zenith_angle=10.0)
    # antenna positions
    antpos=pint.get_antpos()
    # channel pairs used for interferometry
    ch_pairs=n.array(list(itertools.combinations(n.arange(7),2)))
    # matrix of antenna direction vector pairs 
    # one row is a vector for each interferometry channel pair
    dmat=pint.pair_mat(ch_pairs[:,:],antpos)
    
    
    z_rx=n.array(data["zrx_echoes_re"],dtype=n.complex64)+n.array(data["zrx_echoes_im"],dtype=n.complex64)*1j
    z_tx=n.array(data["ztx_pulses_re"],dtype=n.complex64)+n.array(data["ztx_pulses_im"],dtype=n.complex64)*1j
    delays=data["delays"]
    beam_id=data["beam_id"]
    tx_idx=data["tx_idx"]
    c_snr=data["c_snr"]
#    print(n.max(c_snr))
    if n.max(c_snr) < 20 or (0 not in beam_id):
        
#        print("low snr and no zenith beam")
        return

#    interp=1
    if True:#try:
        txlen=z_tx.shape[1]
        echolen=z_rx.shape[2]        
        rds=range_doppler_search(txlen=txlen,echolen=echolen,interp=interp,n_channels=z_rx.shape[1])

        n_ipp=z_rx.shape[0]
        RTI=n.zeros([n_ipp, rds.n_rg],dtype=n.float32)
        xct=n.zeros([n_ipp,rds.n_pairs],dtype=n.complex64)

        peak_rg=[]
        peak_dop=[]
        drg=c.c/1e6/2/1e3
        snr=[]
        for ti in range(z_rx.shape[0]):
            MF,pprof,peak_dopv,noise_floor,xc,rgmax=rds.mf(z_tx[ti,:],z_rx[ti,:,:])
            max_rg=n.argmax(pprof)
            peak_rg.append((delays[ti]+max_rg/interp)*drg)
            peak_dop.append(peak_dopv[max_rg])
            RTI[ti,:]=(pprof-noise_floor)/noise_floor
            xct[ti,:]=xc
            snr.append((pprof[max_rg]-noise_floor)/noise_floor)
            
        snr=n.array(snr)
        peak_rg=n.array(peak_rg)
        peak_dop=n.array(peak_dop)

        gidx=n.where( (beam_id == 0) & (snr > 10))[0]
        if len(gidx)>2:
            mu,mv,mfs,mis,mjs=pint.image_points(phcal,xct[gidx,:],ch_pairs,u,v,w,dmat)
            mw=n.sqrt(1-mu**2-mv**2)
            
            ew=peak_rg[gidx]*mu
            ns=peak_rg[gidx]*mv
            up=peak_rg[gidx]*mw


            r0,v0,model,eres,nres,ures=metfit.simple_fit_meteor(tx_idx[gidx]/1e6,ew,ns,up)
            dout={}
            dout["r0"]=r0
            dout["v0"]=v0            
            dout["std"]=[eres,nres,ures]
            try:
                dmw.write(tx_idx[gidx[0]],dout)
            except:
                import traceback
                traceback.print_exc()

            if plot:
                nm=len(ew)
                plt.figure(figsize=(12,8))
                plt.subplot(231)
                plt.scatter(tx_idx[gidx]/1e6,ew,c=mfs/21)
                plt.plot(tx_idx[gidx]/1e6,model[0:nm],color="blue")
                plt.title(r"$|v_g|$=%1.1f km/s"%(n.linalg.norm(v0)))
     #           cb=plt.colorbar()
    #            cb.set_label("Coherence")

    #            plt.plot(tx_idx[gidx]/1e6,ew,".")
                plt.xlabel("Time (s)")
                plt.ylabel("EW (km)")            

                plt.subplot(232)
                plt.scatter(tx_idx[gidx]/1e6,ns,c=mfs/21)
                plt.title(r"$\sigma_e=%1.1f$ $\sigma_n=%1.1f$, $\sigma_u=%1.1f$ m"%(1e3*eres,1e3*nres,1e3*ures))
                plt.plot(tx_idx[gidx]/1e6,model[nm:(2*nm)],color="blue")            
    #            cb=plt.colorbar()
     #           cb.set_label("Coherence")

    #            plt.plot(tx_idx[gidx]/1e6,ns,".")
                plt.xlabel("Time (s)")
                plt.ylabel("NS (km)")                        
                plt.subplot(233)
                plt.scatter(tx_idx[gidx]/1e6,up,c=mfs/21)
                plt.title(stuffr.unix2datestr(tx_idx[gidx[0]]/1e6))
                plt.plot(tx_idx[gidx]/1e6,model[(2*nm):(3*nm)],color="blue")                        
      #          cb=plt.colorbar()
       #         cb.set_label("Coherence")
    #            plt.plot(tx_idx[gidx]/1e6,up,".")
                plt.xlabel("Time (s)")
                plt.ylabel("Up (km)")            
                plt.subplot(234)            
                plt.scatter(ew,ns,c=mfs/21)
                cb=plt.colorbar()
                cb.set_label("Coherence")
                plt.xlabel("EW (km)")
                plt.ylabel("NW (km)")            
                plt.subplot(235)
                plt.scatter(tx_idx[gidx]/1e6,10.0*n.log10(snr[gidx]),c=mfs/21)
                plt.xlabel("Time (s)")
                plt.ylabel("SNR (dB)")                        
                plt.subplot(236)
                plt.scatter(tx_idx[gidx]/1e6,peak_dop[gidx]/1e3,c=mfs/21)
                plt.xlabel("Time (s)")
                plt.ylabel("Doppler (km/s)")
                plt.tight_layout()
                plt.savefig("/media/analysis/fits/meteor-%1.1f.png"%(float(tx_idx[gidx[0]]/1e6)))
                plt.close()



mddir=pc.cut_metadata_dir
#mddir="../pansy_test_data/metadata/cut"
dm = drf.DigitalMetadataReader(mddir)
b = dm.get_bounds()
dt=10000000
n_block=int((b[1]-b[0])/dt)
os.system("mkdir -p caldata")
start_idx=b[0]
#start_idx=1737912526407585


subdirectory_cadence_seconds = 3600
file_cadence_seconds = 60
samples_per_second_numerator = 1000000
samples_per_second_denominator = 1
file_name = "fit"
os.system("mkdir -p %s"%(pc.simple_fit_metadata_dir))
dmw = drf.DigitalMetadataWriter(
    pc.simple_fit_metadata_dir,
    subdirectory_cadence_seconds,
    file_cadence_seconds,
    samples_per_second_numerator,
    samples_per_second_denominator,
    file_name,
)


for bi in range(n_block):
    data=dm.read(start_idx+bi*dt,start_idx+bi*dt+dt)
    kl=list(data.keys())
    for ki in range(rank,len(kl),size):
        k=kl[ki]
        print("%d %s"%(rank,stuffr.unix2datestr(k/1e6)))
        process_cut(data[k],dmw)



