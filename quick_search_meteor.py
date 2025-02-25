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
import pyfftw 
import h5py
import pansy_config as pc
import traceback

fft = pyfftw.interfaces.scipy_fftpack.fft
#fft=fp.fft

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

class range_doppler_search:
    def __init__(self,
                 txlen=132,
                 rg=n.arange(400,1000,2,dtype=n.int64),
                 fdec=8, # how much do we decimate before fft. can save a lot of compute
                 fftlen=256
                 ):
        
        self.idx=n.arange(txlen,dtype=n.int64)
        self.n_rg=len(rg)
        self.txlen=txlen
        if fftlen < txlen/fdec:
            print("too short fftlen. increasing")
            fftlen=n.ceil(txlen/fdec)
        self.fftlen=fftlen
        self.idx_mat=n.zeros([self.n_rg,self.txlen],dtype=n.int64)
        self.idx=n.arange(self.txlen,dtype=n.int64)
        self.rg=rg
        drg=c.c/2/1e6/1e3
        self.rangev=self.rg*drg
        self.frad=47.5e6
        self.fvec=n.fft.fftshift(n.fft.fftfreq(self.fftlen,d=fdec/1e6))
        self.dopv=self.fvec*c.c/2.0/self.frad
        self.fdec=fdec

        for ri in range(self.n_rg):
            self.idx_mat[ri,:]=self.idx+rg[ri]

    def decim(self,Z,fdec=4):
        new_width=int(self.txlen/self.fdec)
        Z2=n.zeros([self.n_rg,new_width],dtype=n.complex64)
        idx=n.arange(new_width,dtype=n.int64)
        idx2=n.arange(new_width,dtype=n.int64)*fdec

        for i in range(fdec):
            Z2[:,idx]+=Z[:,idx2+i]
        return(Z2)

    def mf(self,z,z_tx,debug=False):
        """
        z = echo
        z_tx = transmit code
        """
        if False:
            plt.subplot(121)
            plt.plot(z.real)
            plt.plot(z.imag)
            plt.subplot(122)
            plt.plot(z_tx.real)
            plt.plot(z_tx.imag)
            plt.show()

        z_tx=n.conj(z_tx)

        n_rx=z.shape[0]
        
        MF=n.zeros([self.n_rg,self.fftlen],dtype=n.float32)
        for rxi in range(n_rx):
            # decode each range gate
            #Z2=n.zeros([self.n_rg,self.txlen])
            ##   Z2[i,:]=z[(self.rg[i]):(self.rg[i]+self.txlen)]*z_tx

            Z=z[rxi,self.idx_mat]*z_tx[None,:]
            # decimate
            ZD=Z
            if self.fdec>1:            
                ZD=self.decim(Z,fdec=self.fdec)
            ZF=n.fft.fftshift(fft(ZD,self.fftlen,axis=1),axes=1)
            MF+=ZF.real**2.0 + ZF.imag**2.0#n.real(ZF*n.conj(ZF))
        noise_floor=n.median(MF)
        pprof=n.max(MF,axis=1)
        peak_dopv=self.dopv[n.argmax(MF,axis=1)]
        
        if debug:
            plt.pcolormesh(self.dopv,self.rangev,n.abs(ZF)**2.0)
            plt.show()
        return(MF,pprof,peak_dopv,noise_floor)


def process_m_mode(key,d,rds,dmw,dm_mf2,chs=["ch000","ch001","ch002","ch003","ch004","ch005","ch006"],debug=False):
    """
    process 20 ipps of meso mode starting at "key"
    """
    n_ch=len(chs)
    z=n.zeros([n_ch,1600*20],dtype=n.complex64)
    for chi in range(n_ch):
        z[chi,:]=d.read_vector_c81d(key,1600*20,chs[chi])
    z_tx=d.read_vector_c81d(key,1600*20,"ch007")
    RTI=n.zeros([20,rds.n_rg],dtype=n.float32)
    RTIV=n.zeros([20,rds.n_rg],dtype=n.float32)
    noise_floors=[]
    tx_pwrs=[]
    max_dops=[]
    max_snrs=[]
    max_ranges=[]
    tx_idxs=[]
    for ti in range(20):
        tx_pwr=n.sum(n.abs(z_tx[(0+ti*1600):(rds.txlen+ti*1600)])**2.0)
        tx_idxs.append(key+ti*1600)
        tx_pwrs.append(tx_pwr)
        MF,pprof,dop_prof,nf=rds.mf(z[:,(0+ti*1600):(1600+ti*1600)],z_tx[(0+ti*1600):(rds.txlen+ti*1600)])
        noise_floors.append(nf)
        RTI[ti,:]=pprof
        RTIV[ti,:]=dop_prof
    tv=n.arange(20)*1.6e-3
    noise_floor=n.median(noise_floors)
    snr= (RTI-noise_floor)/noise_floor
    snr[snr<=0]=0.0001
    for ti in range(20):
        max_rg=n.argmax(snr[ti,:])
        max_dop=RTIV[ti,max_rg]
        max_snr=snr[ti,max_rg]
        if debug:
            print("%s snr=%1.0f range=%1.1f km doppler=%1.1f km/s txp=%1.1f"%(stuffr.unix2datestr((int(key)+ti*1600)/1e6),max_snr,rds.rangev[max_rg],max_dop,tx_pwrs[ti]))
        max_dops.append(max_dop)
        max_ranges.append(rds.rangev[max_rg])
        max_snrs.append(max_snr)
    odata_dict={}
    
    tx_idxs=n.array(tx_idxs)
    #if (key >= i0) & (key<i1):                
    odata_dict["beam_pos_idx"]=[n.arange(20,dtype=n.int8)]
    odata_dict["tx_std"]=[n.repeat(n.std(tx_pwrs),20)]
    odata_dict["tx_pwr"]=[tx_pwrs]
    odata_dict["max_snr"]=[max_snrs]
    odata_dict["max_range"]=[max_ranges]
    odata_dict["max_dopvel"]=[max_dops]
    odata_dict["noise_floor"]=[noise_floors]
    odata_dict["tx_idxs"]=[tx_idxs]

    # write timestamp of last detection in tmp file, so that we know how far the analysis has reached.
    last_fname="/tmp/meteor_mf_%d.h5"%(rank)
    ho=h5py.File(last_fname,"w")
    ho["latest"]=key
    ho.close()
    # check if we have done this already?
    mf2out=dm_mf2.read(key-100,key+100,["beam_pos_idx"])
    if len(mf2out.keys())==0:
        dmw.write([key],odata_dict)
    else:
        print("%d %s looks like this is already processed."%(rank,stuffr.unix2datestr(key/1e6)))


def process_isr_mode(key,d,rds,dmw,chs=["ch000","ch001","ch002","ch003","ch004","ch005","ch006"]):
    """
    process 1 ipp of isr mode starting at "key"
    """
    n_ch=len(chs)
    ipp=12500
    gc=600
    z=n.zeros([n_ch,ipp],dtype=n.complex64)
    for chi in range(n_ch):
        z[chi,:]=d.read_vector_c81d(key,ipp,chs[chi])
    # gc cancel
    z[:,0:gc]=0.0
    z_tx=d.read_vector_c81d(key,ipp,"ch007")
    
    tx_pwr=n.sum(n.abs(z_tx[0:540])**2.0)
    tx_idx=key
    MF,pprof,dop_prof,nf=rds.mf(z,z_tx,debug=True)
    snr= (pprof-nf)/nf
    max_rg=n.argmax(pprof)
    max_range=rds.rangev[max_rg]
    max_dop=dop_prof[max_rg]
    snr[snr<=0]=0.0001
    max_snr=n.max(snr)
    
    if max_snr>10:
        print("%s snr=%1.0f range=%1.1f km doppler=%1.1f km/s txp=%1.1f"%(stuffr.unix2datestr((int(key)/1e6)),max_snr,rds.rangev[max_rg],max_dop,tx_pwr))
        odata_dict={}
        tx_idxs=n.array(tx_idxs)
        odata_dict["tx_pwr"]=tx_pwr
        odata_dict["max_snr"]=max_snr
        odata_dict["max_range"]=max_range
        odata_dict["max_dopvel"]=max_dop
        odata_dict["noise_floor"]=nf
        odata_dict["tx_idxs"]=tx_idx
        dmw.write([key],odata_dict)


def meteor_search(debug=False):

    # this is where the existing metadata lives
    mf_metadata_dir=pc.mf_metadata_dir
    #mf_metadata_dir = "/media/archive/metadata/mf"
    db_mf=[-1,-1]
    dm_mf=None
    if os.path.exists(mf_metadata_dir):
        print("metadata directory exists. searching for last timestamp")
        try:
            dm_mf = drf.DigitalMetadataReader(mf_metadata_dir)
            db_mf = dm_mf.get_bounds()
            print(db_mf)
        except:
            print("couldn't read mf metadata")
    else:
        os.system("mkdir -p %s"%(mf_metadata_dir))

    # setup the directory and file cadence.
    # use 1 MHz, as this is the sample-rate and thus a
    # natural resolution for timing.
    subdirectory_cadence_seconds = 3600
    file_cadence_seconds = 60
    samples_per_second_numerator = 1000000
    samples_per_second_denominator = 1
    file_name = "mf"

    dmw = drf.DigitalMetadataWriter(
        mf_metadata_dir,
        subdirectory_cadence_seconds,
        file_cadence_seconds,
        samples_per_second_numerator,
        samples_per_second_denominator,
        file_name,
    )


    # setup the directory and file cadence.
    # use 1 MHz, as this is the sample-rate and thus a
    # natural resolution for timing.
    subdirectory_cadence_seconds = 3600
    file_cadence_seconds = 60
    samples_per_second_numerator = 1000000
    samples_per_second_denominator = 1
    file_name = "mf"
    os.system("mkdir -p %s"%(pc.mf_isr_metadata_dir))

    dmw_isr = drf.DigitalMetadataWriter(
        pc.mf_isr_metadata_dir,
        subdirectory_cadence_seconds,
        file_cadence_seconds,
        samples_per_second_numerator,
        samples_per_second_denominator,
        file_name,
    )


    # see if results already exist
    dm_mf2 = drf.DigitalMetadataReader(mf_metadata_dir)

    # see if results already exist
    dm_mf_isr = drf.DigitalMetadataReader(pc.mf_isr_metadata_dir)

    # this is where the data is
    d=drf.DigitalRFReader("/media/archive/")
    # tx channel bounds
    b=d.get_bounds("ch007")

    # this is the first sample in data
    i0=b[0]
    if db_mf[1] != -1:
        # start where we left off, instead of the start
        i0=db_mf[1]+10
    
    print("starting at %s"%(stuffr.unix2datestr(i0/1e6)))

    # 100 seconds per analysis window
    #dt=60000000
    #n_windows = int(n.floor((b[1]-i0)/dt))
    
#    d=drf.DigitalRFReader("/media/archive")
 #   print(d.get_bounds("ch000"))

    stm=pm.get_m_mode()
    ipp=stm["ipp_us"]
    beam_pos=stm["beam_pos_az_za"]
    n_beam=len(beam_pos)
    n_codes=len(stm["codes"])
    # get code vector
    codes=pm.get_vector(stm,ncodes=n_codes)

    # transmit metadata
    metadata_dir=pc.tx_metadata_dir
    #metadata_dir = "/media/archive/metadata/tx"
    if not os.path.exists(metadata_dir):
        print("metadata directory doesn't exist. exiting")
        exit(0)

    dmr = drf.DigitalMetadataReader(metadata_dir)
    db = dmr.get_bounds()
    print(db)
    b=d.get_bounds("ch000")

    d_analysis=file_cadence_seconds*1000000
    # start analysis where the previous one left off
    start_idx=d_analysis*int(n.ceil((db_mf[1]-6*3600*1000000)/d_analysis))
    # stay 6 minutes behind realtime to avoid underfull metadata files
    end_idx=d_analysis*int(n.ceil(db[1]/d_analysis))-6*d_analysis

    if end_idx < start_idx:
        print("end before start! waiting for more data.")
        return
#        exit(0)

    # minutes since 1970
    end_minute=int(n.floor(end_idx/d_analysis))
    start_minute=int(n.floor(start_idx/d_analysis))
    n_minutes=end_minute-start_minute

    rds=range_doppler_search()
    N=20*1600
    beam_pos_idx=n.arange(20,dtype=n.int8)%5

    rds_isr=range_doppler_search(txlen=540)

    # analyze in parallel. one minute for each thread
    for bi in range(start_minute,end_minute):
        if bi%size != rank:
#            print("rank %d skipping minute for rank %d"%(rank,bi%size))
            continue
        
        cput0=time.time()                      
#        i0=start_minute*60*1000000 + bi*60*1000000
 #       i1=start_minute*60*1000000 + bi*60*1000000 + 60*1000000
        i0=bi*60*1000000 
        i1=bi*60*1000000 + 60*1000000
            
        b=d.get_bounds("ch000")

        # only process if we have raw voltage data in ringbuffer
        if (i0 > b[0]) & (i1 < b[1]):
            print("%d processing %s"%(rank,stuffr.unix2datestr(i0/1e6)))
            data_dict = dmr.read(i0, i1, "id")

            mf2out=dm_mf2.read(i0,i1,["beam_pos_idx"])
            if len(mf2out.keys())>0:
                print("rank %d already processed %d-%d %d results"%(rank,i0,i1,len(mf2out.keys())))
                continue

            mf2out=dm_mf_isr.read(i0,i1,["tx_pwr"])
            if len(mf2out.keys())>0:
                print("rank %d already processed %d-%d %d results"%(rank,i0,i1,len(mf2out.keys())))
                continue

            keys=data_dict.keys()
            n_keys=len(keys)
            #print("%d processing %d pulses"%(rank,20*n_keys))
            for key in keys:
                keyi=int(key)
                try:
                    if  data_dict[key] == 1:
                        print("meso mode %s"%(stuffr.unix2datestr(key/1e6)))
                        process_m_mode(key,d,rds,dmw,dm_mf2,chs=["ch000","ch001","ch002","ch003","ch004","ch005","ch006"])
                    elif data_dict[key] == 2:
                        print("isr mode %s"%(stuffr.unix2datestr(key/1e6)))
                        process_isr_mode(key,d,rds_isr,dmw_isr,chs=["ch000","ch001","ch002","ch003","ch004","ch005","ch006"])
                    else:
                        print("unknown mode %d"%(data_dict[key]))

                except:
                    traceback.print_exc()


        cput1=time.time()
        if n_keys > 0:
            print("rank %d %d pulses %s cputime/realtime %1.2f"% (rank,n_keys*20,stuffr.unix2datestr(i0/1e6), (cput1-cput0)/(size*n_keys*20*1.6e-3)))

    
    
    
if __name__ == "__main__":
    while True:
        meteor_search()
        print("wait for all threads to finnish")
        comm.Barrier()
        time.sleep(60)
    
