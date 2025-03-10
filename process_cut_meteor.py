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
import pansy_modes as pmm
import meteor_fit as metfit
import simple_radiant as sr

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

phasecal = n.zeros([5,7])
imaging=[]
mmdict=pmm.get_m_mode()

def beam_pixmap():
    u0,v0,w0=pint.uv_coverage(N=400,max_zenith_angle=20.0)
    ew=u0*100
    ns=v0*100
    bitmap=n.zeros(u0.shape)
    bitmap[:,:]=0
    for i in range(5):
        az0=mmdict["beam_pos_az_za"][i][0]
        el0=90-mmdict["beam_pos_az_za"][i][1]
        p_h=n.cos(n.pi*el0/180.0)
        pw0=-n.sin(n.pi*el0/180.0)
        pv0=p_h*n.cos(-n.pi*az0/180)
        pu0=-p_h*n.sin(-n.pi*az0/180)
        angle=180*n.arccos(u0*pu0 + v0*pv0 + pw0*w0)/n.pi
        bitmap[angle<10]+=1
    return(ew,ns,bitmap)
# backdrop for plot, to show the five beam pointing directions
ewbm,nsbm,bitmap=beam_pixmap()

for i in range(5):
    # each antenna needs to be multiplied with these phases (n.exp(-1j*phasecal))
    h=h5py.File("data/phases%d.h5"%(i),"r")
    phasecal[i,:]=h["phasecal"][()]
    h.close()
    u,v,w=pint.uv_coverage2(N=500,
                            az0=mmdict["beam_pos_az_za"][i][0],
                            el0=90.0-mmdict["beam_pos_az_za"][i][1],
                            max_angle=10.0)
    if False:
        plt.subplot(121)
        plt.pcolormesh(u)
        plt.title("beam %d"%(i))
        plt.colorbar()
        plt.subplot(122)
        plt.pcolormesh(v)
        plt.colorbar()
        plt.tight_layout()
        plt.show()

    imaging.append({"u":u,"v":v,"w":w})

antpos=pint.get_antpos()
ch_pairs=n.array(list(itertools.combinations(n.arange(7),2)))
dmat=pint.pair_mat(ch_pairs,antpos)

fft=fp.fft

class range_doppler_search:
    """
    Range-Doppler matched filterbank
    """
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
                interp=1,
                write_dm=True,
                plot=False):


    # create a uniform grid of u,v,w values
#    u,v,w=pint.uv_coverage(N=500,max_zenith_angle=10.0)
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
    if n.max(c_snr) < 20:# or (0 not in beam_id):        
#        print("low snr and no zenith beam")
        return

    if True:
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
        print("range-doppler matched filter %d"%(z_rx.shape[0]))
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

        beam_ids=n.unique(beam_id)

        ews=[]
        nss=[]
        ups=[]
        txidxs=[]
        mfss=[]
        peak_dops=[]
        snrs=[]
        beam_idss=[]

        #        for bi in beam_ids:
        gidx=n.where( snr > 10 )[0]
        if len(gidx)>2:
            mu,mv,mfs=pint.image_points1(phasecal,xct[gidx,:],ch_pairs,dmat,snr[gidx],beam_id[gidx],tx_idx[gidx])
            mw=n.sqrt(1-mu**2-mv**2)
            ew=peak_rg[gidx]*mu
            ns=peak_rg[gidx]*mv
            up=peak_rg[gidx]*mw
            ews=n.concatenate((ews,ew))
            nss=n.concatenate((nss,ns))
            ups=n.concatenate((ups,up))
            txidxs=n.concatenate((txidxs,tx_idx[gidx]))
            mfss=n.concatenate((mfss,mfs))
            snrs=n.concatenate((snrs,snr[gidx]))
            peak_dops=n.concatenate((peak_dops,peak_dop[gidx]))
            beam_idss=n.concatenate((beam_idss,beam_id[gidx]))

            r0,v0,model,eres,nres,ures=metfit.simple_fit_meteor(txidxs/1e6,ews,nss,ups)

            rres=sr.get_radiant(r0*1e3,txidxs[0]/1e6,v0)

            dout={}
            dout["r0"]=r0
            dout["v0"]=v0            
            dout["std"]=[eres,nres,ures]
            dout["ew"]=ews
            dout["ns"]=nss
            dout["up"]=ups
            dout["snr"]=snrs
            dout["mfs"]=mfss
            dout["txidx"]=txidxs
            dout["eclat"]=rres["lat"]
            dout["eclon"]=rres["lon"]

            print("%d beam %d %s %1.2f km %1.2f km/s %1.2f %1.2f %1.2f"%(rank,bi,stuffr.unix2datestr(txidxs[0]/1e6),
                                                                         n.linalg.norm(r0),n.linalg.norm(v0),eres*1e3,nres*1e3,ures*1e3))
            if write_dm:
                try:
                    dmw.write(txidxs[0],dout)
                except:
                    import traceback
                    traceback.print_exc()

            if plot:
                nm=len(ews)
                plt.figure(figsize=(12,8))
                plt.subplot(231)
                plt.scatter(txidxs/1e6,ews,c=mfss/21)
                plt.plot(txidxs/1e6,model[0:nm],color="blue")
                plt.title(r"$|v_g|$=%1.1f km/s"%(n.linalg.norm(v0)))
                plt.xlabel("Time (s)")
                plt.ylabel("EW (km)")            
                plt.subplot(232)
                plt.scatter(txidxs/1e6,nss,c=mfss/21)
                plt.title(r"$\sigma_e=%1.1f$ $\sigma_n=%1.1f$, $\sigma_u=%1.1f$ m"%(1e3*eres,1e3*nres,1e3*ures))
                plt.plot(txidxs/1e6,model[nm:(2*nm)],color="blue")            
                plt.xlabel("Time (s)")
                plt.ylabel("NS (km)")                        
                plt.subplot(233)
                plt.scatter(txidxs/1e6,ups,c=mfss/21)
                plt.title(stuffr.unix2datestr(txidxs[0]/1e6))
                plt.plot(txidxs/1e6,model[(2*nm):(3*nm)],color="blue")
                plt.xlabel("Time (s)")
                plt.ylabel("Up (km)")            

                plt.subplot(234)

                plt.pcolormesh(ewbm,nsbm,bitmap,cmap="gist_yarg",vmax=5)
                plt.scatter(ews,nss,c=beam_idss,s=1,cmap="rainbow",vmin=0,vmax=4)
                sidx=n.argmin(txidxs)
                plt.text(ews[sidx],nss[sidx],r"$t_0$")

                plt.xlabel("EW (km)")
                plt.ylabel("NS (km)")     
                plt.ylim([-35,35])
                plt.xlim([-35,35])

                plt.subplot(235)
                plt.scatter(txidxs/1e6,10.0*n.log10(snrs),c=beam_idss,cmap="rainbow",vmin=0,vmax=4)
                cb=plt.colorbar()
                cb.set_label("beam")

                plt.xlabel("Time (s)")
                plt.ylabel("SNR (dB)")                        
                plt.subplot(236)
                plt.scatter(txidxs/1e6,peak_dops/1e3,c=mfss/21)
                cb=plt.colorbar()
                cb.set_label("Coherence")
                plt.xlabel("Time (s)")
                plt.ylabel("Doppler (km/s)")
                plt.tight_layout()
                #plt.show()
                plt.savefig("/tmp/pansy/meteor-%1.1f.png"%(float(tx_idx[0]/1e6)))
                plt.close()


                #fig=plt.figure(figsize=(16,9))
                if False:
                    fig,ax=plt.subplots(figsize=(16,9))
                    ax.pcolormesh(ewbm,nsbm,bitmap,cmap="gist_yarg",vmax=6)
                    ax.set_title("PANSY Meteor Head Echoes %s"%(stuffr.unix2datestr(tx_idx[0]/1e6)))
                    m=ax.scatter(ews,nss,c=10.0*n.log10(snrs),cmap="plasma",vmin=10,vmax=30)
                    cb=fig.colorbar(m,ax=ax)
                    cb.set_label("SNR (dB)")
                    ax.set_aspect("equal")#,s=1
                    #plt.axes().set_aspect('equal')
                    sidx=n.argmin(txidxs)
                    ax.text(ews[sidx],nss[sidx],r"$t_0$")
                    ax.set_xlabel("EW (km)")
                    ax.set_ylabel("NS (km)")     
                    ax.set_ylim([-35,35])
                    ax.set_xlim([-35,35])
                    plt.tight_layout()
#                    plt.show()
                    plt.savefig("/tmp/pansy/hor-%1.1f.png"%(float(tx_idx[0]/1e6)))
                    plt.close()


def process_latest():
    mddir=pc.cut_metadata_dir
 #   mddir="../pansy_test_data/metadata/cut"
    dm = drf.DigitalMetadataReader(mddir)
    b = dm.get_bounds()
    dt=120000000
#    os.system("mkdir -p caldata")

    start_idx=b[0]

    try:
        dmf=drf.DigitalMetadataReader(pc.simple_fit_metadata_dir)
        fitb=dmf.get_bounds()
        start_idx=fitb[1]
    except:
        print("no fit boundary found")
    
    n_block=int(n.ceil((b[1]-start_idx)/dt))
    print(stuffr.unix2datestr(b[1]/1e6))
    print(stuffr.unix2datestr(start_idx/1e6))

    subdirectory_cadence_seconds = 3600
    file_cadence_seconds = 60
    samples_per_second_numerator = 1000000
    samples_per_second_denominator = 1
    file_name = "fit"
    omddir=pc.simple_fit_metadata_dir
#    omddir="/tmp/simple_fit"

    os.system("mkdir -p %s"%(omddir))#pc.simple_fit_metadata_dir))

    dmw = drf.DigitalMetadataWriter(
        omddir,
        subdirectory_cadence_seconds,
        file_cadence_seconds,
        samples_per_second_numerator,
        samples_per_second_denominator,
        file_name,
    )

    for bi in range(rank,n_block,size):
        data=dm.read(start_idx+bi*dt,start_idx+bi*dt+dt)
        kl=list(data.keys())
        for ki in range(len(kl)):
            k=kl[ki]
            try:
                print(stuffr.unix2datestr(k/1e6))
                process_cut(data[k],dmw)
            except:
                import traceback
                traceback.print_exc()
            


if __name__ == "__main__":
    while True:
        process_latest()
        comm.Barrier()
        time.sleep(60)