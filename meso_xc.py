import numpy as n
import digital_rf as drf
import matplotlib.pyplot as plt
import pansy_modes as pm
import scipy.signal.windows as sw
import scipy.constants as c
import stuffr
import time
import pansy_config as pc
import itertools
import os
import traceback
import h5py

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def analyze_block(i0,
                  i1,
                  rx_ch=["ch000","ch001","ch002","ch003","ch004","ch005","ch006"],
                  tx_ch="ch007",
                  r0=75,r1=110, # range interval to store
                  max_dop=27.0, # largest Doppler shift (Hz) stored
                  txlen=132,
                  plot=True,
                  n_cycles=312): # how many 20 ipp cycles are stored in one spectrum
    """
    find contiguous blocks of mesosphere mode
    """

    subdirectory_cadence_seconds = 3600
    file_cadence_seconds = 60
    samples_per_second_numerator = 1000000
    samples_per_second_denominator = 1
    file_name = "xc"
    os.system("mkdir -p %s"%(pc.xc_metadata_dir))

    dmw = drf.DigitalMetadataWriter(
        pc.xc_metadata_dir,
        subdirectory_cadence_seconds,
        file_cadence_seconds,
        samples_per_second_numerator,
        samples_per_second_denominator,
        file_name,
    )

    dmt = drf.DigitalMetadataReader(pc.tx_metadata_dir)
    dd=dmt.read(i0,i1)
    
    n_ch=len(rx_ch)
    # how many ipps 
    n_ipp=4*n_cycles
    # how many tx beams 
    n_beams=5
    beams=n.arange(5)
    # xc and self correlations
    
    # what pairs do we store
    ch_pairs=[]
    for i in range(n_ch):
        ch_pairs.append((i,i))
    ch_pairs2=list(itertools.combinations(n.arange(7),2))
    for i in range(len(ch_pairs2)):
        ch_pairs.append(ch_pairs2[i])

    n_xc=len(ch_pairs)
        
    drg=c.c/1e6/2/1e3
    ipp=1600
    rvec=n.arange(ipp)*drg

    ridx=n.where( (rvec > r0) & (rvec < r1))[0]
    ri0=n.min(ridx)
    ri1=n.max(ridx)
    wfun=sw.hann(n_ipp)
    fvec=n.fft.fftshift(n.fft.fftfreq(n_ipp,d=(5*ipp/1e6)))
    fidx=n.where( n.abs(fvec)<max_dop)[0]
    fi0=n.min(fidx)
    fi1=n.max(fidx)
    Z=n.zeros([n_ch,n_beams,n_ipp,ipp],dtype=n.complex64)
    S=n.zeros([n_ch,n_beams,n_ipp,ipp],dtype=n.complex64)
    XC=n.zeros([n_xc,n_beams,n_ipp,ipp],dtype=n.complex64)

    W=n.zeros([n_xc,n_beams,ipp])
    WS=n.zeros([n_xc,n_beams,ipp])    

    txpulse=n.zeros(1600,dtype=n.complex64)
    z_echo=n.zeros(1600,dtype=n.complex64)
    ipp_idx0=0
    print("%d %s"%(rank,stuffr.unix2datestr(i0/1e6)))
    for k in dd.keys():

        ztx=d.read_vector_c81d(k,1600*20,tx_ch)
        
        for chi in range(n_ch):
            z=d.read_vector_c81d(k,1600*20,rx_ch[chi])
            for bi in range(n_beams):
                for ti in range(4):
                    si0=ti*5*1600 + bi*1600
                    txpulse[0:txlen]=ztx[si0:(si0+txlen)]
                    z_echo[:]=z[si0:(si0+1600)]
                    # gc remove
                    z_echo[0:(txlen+100)]=0.0
                    Z[chi,bi,ti+ipp_idx0,:]=n.fft.ifft(n.conj(n.fft.fft(txpulse))*n.fft.fft(z_echo))
        # we get for ipps for each 20 pulse cycle
        ipp_idx0+=4
        if ipp_idx0 >= n_ipp:
            for chi in range(n_ch):
                for bi in range(n_beams):
                    for ri in range(1600):
                        S[chi,bi,:,ri]=n.fft.fftshift(n.fft.fft(wfun*Z[chi,bi,:,ri]))
                        
            for pi in range(n_xc):
                for bi in range(n_beams):
                    ch0=ch_pairs[pi][0]
                    ch1=ch_pairs[pi][1]
                    W[pi,bi,:]=n.sum(n.abs(S[ch0,bi,:,:]*n.conj(S[ch1,bi,:,:])),axis=0)
                    XC[pi,bi,:,:]+=S[ch0,bi,:,:]*n.conj(S[ch1,bi,:,:])/W[pi,bi,None,:]
                    WS[pi,bi,:]+=1.0/W[pi,bi,:]
            ipp_idx0=0
            if False:
                plt.subplot(121)
                plt.pcolormesh(fvec,rvec,10.0*n.log10(n.abs(XC[0,1,:,:].T)),vmin=noise_floor,vmax=noise_floor+20)
                plt.title("%s"%(stuffr.unix2datestr(i0/1e6)))
                plt.ylim([r0,r1])
                plt.xlim([-max_dop,max_dop])
                plt.xlabel("Doppler (Hz)")
                plt.ylabel("Range (km)")
                plt.colorbar()
                plt.subplot(122)            
                plt.pcolormesh(fvec,rvec,n.angle(XC[7,1,:,:].T),cmap="hsv")
                plt.title("beam %d (%d,%d)"%(1,ch_pairs[7][0],ch_pairs[7][1]))
                plt.ylim([r0,r1])
                plt.xlim([-max_dop,max_dop])
                plt.xlabel("Doppler (Hz)")
                plt.ylabel("Range (km)")
                plt.colorbar()
                plt.tight_layout()
                plt.show()
            

    for pi in range(n_xc):
        for bi in range(n_beams):
            XC[pi,bi,:,:]=XC[pi,bi,:,:]/WS[pi,bi,None,:]
    data_out={}
    data_out["xc"]=n.array(XC[:,:,fi0:fi1,ri0:ri1],dtype=n.complex64)
    data_out["rvec"]=rvec[ri0:ri1]
    data_out["fvec"]=rvec[fi0:fi1]
    data_out["ch_pairs"]=n.array(ch_pairs,dtype=n.int8)
    data_out["rx_ch"]=n.array(rx_ch)
    data_out["i0"]=i0
    data_out["i1"]=i1    
    data_out["beams"]=n.array(beams,dtype=n.int8)
    ho=h5py.File("/media/analysis/pmse/xc-%d.h5"%(i0/1e6),"w")
    for okey in data_out.keys():
        ho[okey]=data_out[okey]
    ho.close()
    

    if plot:
        pi=7
        bi=1
        dB=10.0*n.log10(n.abs(XC[pi,bi,:,:].T))
        if bi < 7:
            noise_floor=n.nanmedian(dB)
        else:
            noise_floor=-3
        #10.0*n.log10(n.abs(XC[pi,bi,:,:].T))
        plt.subplot(121)
        plt.pcolormesh(fvec,rvec,dB,vmin=noise_floor,vmax=noise_floor+20)
        plt.title("%s"%(stuffr.unix2datestr(i0/1e6)))
        plt.ylim([r0,r1])
        plt.xlim([-max_dop,max_dop])
        plt.xlabel("Doppler (Hz)")
        plt.ylabel("Range (km)")
        plt.colorbar()
        plt.subplot(122)            
        plt.pcolormesh(fvec,rvec,n.angle(XC[pi,bi,:,:].T),cmap="hsv")
        plt.title("beam %d (%d,%d)"%(bi,ch_pairs[pi][0],ch_pairs[pi][1]))
        plt.ylim([r0,r1])
        plt.xlim([-max_dop,max_dop])
        plt.xlabel("Doppler (Hz)")
        plt.ylabel("Range (km)")
        plt.colorbar()
        plt.tight_layout()
        plt.savefig("xc-%03d-%03d-%d.png"%(pi,bi,i0/1e6))
        plt.close()
                #                plt.show()
        if False:
            for pi in range(n_xc):
                for bi in range(n_beams):
                    dB=10.0*n.log10(n.abs(XC[pi,bi,:,:].T))
                    if bi < 8:
                        noise_floor=n.nanmedian(dB)
                    else:
                        noise_floor=-3
                    #10.0*n.log10(n.abs(XC[pi,bi,:,:].T))
                    plt.subplot(121)
                    plt.pcolormesh(fvec,rvec,dB,vmin=noise_floor,vmax=noise_floor+20)
                    plt.title("%s"%(stuffr.unix2datestr(i0/1e6)))
                    plt.ylim([r0,r1])
                    plt.xlim([-max_dop,max_dop])
                    plt.xlabel("Doppler (Hz)")
                    plt.ylabel("Range (km)")
                    plt.colorbar()
                    plt.subplot(122)            
                    plt.pcolormesh(fvec,rvec,n.angle(XC[pi,bi,:,:].T),cmap="hsv")
                    plt.title("beam %d (%d,%d)"%(bi,ch_pairs[pi][0],ch_pairs[pi][1]))
                    plt.ylim([r0,r1])
                    plt.xlim([-max_dop,max_dop])
                    plt.xlabel("Doppler (Hz)")
                    plt.ylabel("Range (km)")
                    plt.colorbar()
                    plt.tight_layout()
                    plt.savefig("xc-%03d-%03d-%d.png"%(pi,bi,i0))
                    plt.close()
    #                plt.show()


                    
                
            


    
    #mesomode_metadata_dir="/media/analysis/metadata/mesomode"

if __name__ == "__main__":
    dmm = drf.DigitalMetadataReader(pc.mesomode_metadata_dir)
    d = drf.DigitalRFReader(pc.raw_voltage_dir)
    b=dmm.get_bounds()
    t0=stuffr.date2unix(2025,1,30,0,0,00)*1000000
    dd=dmm.read(t0,t0+24*3600*1000000)
    kl=list(dd.keys())
    for ki in range(rank,len(kl),size):#dd.keys():
        k=kl[ki]
        i0=dd[k]["start"]
        i1=dd[k]["end"]
        rb=d.get_bounds("ch000")
        if (i1-i0)>60*1000000 and (rb[0]< i0):
            print("long enough")
#            print("rank %d %s - %s"%(rank,stuffr.unix2datestr(i0/1e6),stuffr.unix2datestr(i1/1e6)))
            analyze_block(i0,i1)

    

