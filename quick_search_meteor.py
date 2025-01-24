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

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


class range_doppler_search:
    def __init__(self,txlen=130,
                 rg=n.arange(400,950,2,dtype=n.int64)):
        
        self.idx=n.arange(130,dtype=n.int64)
        self.n_rg=len(rg)
        self.txlen=txlen
        self.idx_mat=n.zeros([self.n_rg,self.txlen],dtype=n.int64)
        self.idx=n.arange(self.txlen,dtype=n.int64)
        self.rg=rg
        self.rangev=self.rg*0.15
        self.frad=47.5e6
        self.fvec=n.fft.fftshift(n.fft.fftfreq(self.txlen,d=1/1e6))
        self.dopv=self.fvec*c.c/2.0/self.frad

        for ri in range(self.n_rg):
            self.idx_mat[ri,:]=self.idx+rg[ri]

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

        # decode each range gate
        #Z2=n.zeros([self.n_rg,self.txlen])
        ##   Z2[i,:]=z[(self.rg[i]):(self.rg[i]+self.txlen)]*z_tx
        Z=z[self.idx_mat]*z_tx[None,:]

        # df = 2*f*v/c
        # df*c/2/f

        ZF=n.fft.fftshift(fp.fft(Z,axis=1),axes=1)
        pwr=n.real(ZF*n.conj(ZF))
        noise_floor=n.median(pwr)
        pprof=n.max(pwr,axis=1)
        peak_dopv=self.dopv[n.argmax(pwr,axis=1)]

        if debug:
            plt.pcolormesh(self.dopv,self.rangev,n.abs(ZF)**2.0)
            plt.show()
        return(pwr,pprof,peak_dopv,noise_floor)

def meteor_search(debug=False):

    mf_metadata_dir = "/media/archive/metadata/mf"
    db_mf=[-1,-1]
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
    metadata_dir = "/media/archive/metadata/tx"
    if not os.path.exists(metadata_dir):
        print("metadata directory doesn't exist. exiting")
        exit(0)

    dmr = drf.DigitalMetadataReader(metadata_dir)
    db = dmr.get_bounds()
    print(db)
    b=d.get_bounds("ch000")

    start_idx=i0#db[1]-200000000

    d_analysis=file_cadence_seconds*1000000
    end_minute=int(n.floor(db[1]/d_analysis))
    start_minute=int(n.floor(start_idx/d_analysis))

    n_minutes=end_minute-start_minute
    #    n_blocks=int(n.floor((db[1]-start_idx)/(ipp*n_codes)))
    #   RTI = n.zeros([n_beam,n_codes,ipp],dtype=n.float32)

    rds=range_doppler_search()
    N=20*1600
    beam_pos_idx=n.arange(20,dtype=n.int8)%5
    # analyze in parallel. one minute for each thread
    for bi in range(rank,n_minutes,size):
        cput0=time.time()                        
        i0=start_minute*60*1000000 + bi*60*1000000
        i1=start_minute*60*1000000 + bi*60*1000000 + 60*1000000
        print("rank %d processing minute %d/%d"%(rank,bi,n_minutes))
        db_mf = dm_mf.get_bounds()
        #       if db_mf[1] >= i0:
        #          print("skipping block %d, because it is already processed"%(bi))
        #         continue
        #    print("processing %d"%(bi))
            
        b=d.get_bounds("ch000")
        # if we have raw voltage

        # only process if we have raw voltage data in ringbuffer
        if (i0 > b[0]) & (i1 < b[1]):
            data_dict = dmr.read(i0, i1, "id")
            keys=data_dict.keys()
            print("processing %d pulses"%(20*len(keys)))
            for key in keys:
                if key < db_mf[1]:
                    print("already processed %d. skipping"%(key))
                
                print((key, data_dict[key]))
                z=d.read_vector_c81d(key,1600*20,"ch000")
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
#                    print(tx_pwr)
                    MF,pprof,dop_prof,nf=rds.mf(z[(0+ti*1600):(1600+ti*1600)],z_tx[(0+ti*1600):(rds.txlen+ti*1600)])
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
                data_dict={}
                
                tx_idxs=n.array(tx_idxs)
                data_dict["beam_pos_idx"]=n.arange(20,dtype=n.int8)
                data_dict["tx_std"]=n.repeat(n.std(tx_pwrs),20)
                data_dict["tx_pwr"]=tx_pwrs
                data_dict["max_snr"]=max_snrs
                data_dict["max_range"]=max_ranges
                data_dict["max_dopvel"]=max_dops
                data_dict["noise_floor"]=noise_floors
                dmw.write(tx_idxs,data_dict)

        cput1=time.time()
        print("%s cputime/realtime %1.2f"% (stuffr.unix2datestr(i0/1e6), (cput1-cput0)/(size*60.0)))

    
    
    
if __name__ == "__main__":
    while True:
        meteor_search()
        time.sleep(1)
    
