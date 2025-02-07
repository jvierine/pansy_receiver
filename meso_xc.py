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

def analyze_block(i0,i1,
                  rx_ch=["ch000","ch001","ch002","ch003","ch004","ch005","ch006"],
                  tx_ch="ch007",
                  r0=60,r1=130, # range interval to store
                  max_dop=27.0, # largest Doppler shift (Hz) stored
                  txlen=132,
                  n_cycles=64): # how many 20 ipp cycles are stored in one spectrum
    """
    find contiguous blocks of mesosphere mode
    """
    dmt = drf.DigitalMetadataReader(pc.tx_metadata_dir)
    dd=dmt.read(i0,i1)
    
    n_ch=len(rx_ch)
    # how many ipps 
    n_ipp=4*n_cycles
    # how many tx beams 
    n_beams=5
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
    wfun=sw.hann(n_ipp)
    fvec=n.fft.fftshift(n.fft.fftfreq(n_ipp,d=(5*ipp/1e6)))
    Z=n.zeros([n_ch,n_beams,n_ipp,ipp],dtype=n.complex64)
    S=n.zeros([n_ch,n_beams,n_ipp,ipp],dtype=n.complex64)
    XC=n.zeros([n_xc,n_beams,n_ipp,ipp],dtype=n.complex64)

    txpulse=n.zeros(1600,dtype=n.complex64)
    z_echo=n.zeros(1600,dtype=n.complex64)
    ipp_idx0=0
    for k in dd.keys():
        print("%d %s"%(ipp_idx0,stuffr.unix2datestr(k/1e6)))
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
                    Z[chi,bi,ti+ipp_idx0,:]=n.fft.ifft(n.fft.fft(n.conj(txpulse))*n.fft.fft(z_echo))
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
                    XC[pi,bi,:,:]+=S[ch0,bi,:,:]*n.conj(S[ch1,bi,:,:])
            ipp_idx0=0

    for pi in range(n_xc):
        for bi in range(n_beams):
            plt.subplot(121)
            plt.pcolormesh(fvec,rvec,n.abs(XC[pi,0,:,:].T))
            plt.title("%s"%(stuffr.unix2datestr(i0/1e6)))
            plt.ylim([r0,r1])
            plt.xlim([-max_dop,max_dop])
            plt.xlabel("Doppler (Hz)")
            plt.ylabel("Range (km)")
            plt.colorbar()
            plt.subplot(122)            
            plt.pcolormesh(fvec,rvec,n.angle(XC[pi,0,:,:].T),cmap="hsv")
            plt.title("beam %d (%d,%d)"%(bi,ch_pairs[pi][0],ch_pairs[pi][1]))
            plt.ylim([r0,r1])
            plt.xlim([-max_dop,max_dop])
            plt.xlabel("Doppler (Hz)")
            plt.ylabel("Range (km)")
            plt.colorbar()
            plt.tight_layout()
            plt.show()


                    
                
            


    
    #mesomode_metadata_dir="/media/analysis/metadata/mesomode"

if __name__ == "__main__":
    dmm = drf.DigitalMetadataReader(pc.mesomode_metadata_dir)
    d = drf.DigitalRFReader(pc.raw_voltage_dir)
    b=dmm.get_bounds()
    t0=stuffr.date2unix(2025,1,31,11,00,00)*1000000
    dd=dmm.read(t0,t0+3600*1000000)
    for k in dd.keys():
        i0=dd[k]["start"]
        i1=dd[k]["end"]
        rb=d.get_bounds("ch000")
        if (i1-i0)>60*1000000 and (rb[0]< i0):
            print("long enough")
            analyze_block(i0,i1)

    
    #mesomode_metadata_dir="/media/analysis/metadata/mesomode"
    

def analyze_m_mode(d,
                   dt=10,
                   debug=False,
                   b0=None,
                   b1=None,
                   channels=["ch000","ch001","ch002","ch003","ch004","ch005","ch006"],
                   alias=False):
    # The mesospheric mode is a special sequence of
    # 20 16-bit complementary codes
    # two codes are sent in each direction
    # 
    
    # simple range-Doppler power spectrum for the M-mode
    b=d.get_bounds("ch000")
    if b0 == None:
        b0=b[0]
    if b1==None:
        b1=b[1]

    stm=pm.get_m_mode()
    ncodes=20
    codes=pm.get_vector(stm,ncodes=len(stm["codes"]))
    codes.shape=(ncodes,1600) 
    C=n.copy(codes)
    
    for i in range(codes.shape[0]):
        if debug:
            plt.plot(codes[i,:].real)
            plt.show()

        C[i,:]=n.conj(n.fft.fft(codes[i,:]))

    # how many integration periods
    n_ints=int(n.floor((b1-b0)/1e6/dt))
    for i in range(rank,n_ints,size):
        i0=b0+i*dt*1000000
        i1=b0+(i+1)*dt*1000000
        m_start_idxs=find_m_mode_start(d,i0,i1)
        n_beam=len(stm["beam_pos_az_za"])
        n_rg=1600
        n_pulse=4
        n_reps=32
        wf=sw.hann(n_reps
                   *n_pulse)
        
        n_channels=len(channels)
        
        Z=n.zeros([n_channels,n_beam,n_rg,n_pulse*n_reps],dtype=n.complex64)
        S=n.zeros([n_channels,n_beam,n_rg,n_pulse*n_reps],dtype=n.float32)
        
        for chi in range(len(channels)):
            ch=channels[chi]
            rvec=n.arange(n_rg)*0.15
            if alias:
                rvec+=2*1600*0.15
            # 2*f*v/c = df
            # df*c/f/2 = v
            fvec=n.fft.fftshift(n.fft.fftfreq(n_reps*n_pulse,d=n_beam*1600/1e6))*c.c/47.5e6/2.0
            if len(m_start_idxs) > 1:
                print("%1.1f found %d sequences"%(i*dt,len(m_start_idxs)))
                n_cycles=int(n.floor(len(m_start_idxs)/n_reps))
                print(n_cycles)
                for ci in range(n_cycles):
                    print("cycle %d %s"%(ci,ch))
                    for ri in range(n_reps):
                        i0=m_start_idxs[ci*n_reps + ri]
                        for pi in range(n_pulse):
                            for bi in range(n_beam):
                                iread = i0 + pi*1600*n_beam + bi*1600 #+ 20
                                z=d.read_vector_c81d(iread,1600,ch)
                                if alias:
                                    z2=d.read_vector_c81d(iread+2*1600,1600,ch)
                                    z[200:1600]=z2[200:1600]
                                codei = (pi*n_beam + bi)%20
                                Z[chi,bi,:,ri*n_pulse + pi]=n.fft.ifft(n.fft.fft(z)*C[codei,:])
                                #plt.plot(Z[bi,:,ri*n_pulse + pi].real)
                                #plt.plot(Z[bi,:,ri*n_pulse + pi].imag)
                                #plt.show()
                                # phase shift based on transmit pulse phase 
                                phase = n.angle(Z[chi,bi,0,ri*n_pulse + pi])
                                Z[chi,bi,:,ri*n_pulse + pi]=Z[chi,bi,:,ri*n_pulse + pi]*n.exp(-1j*phase)
        #                        if True:
        #                           plt.plot(Z[bi,:,pi].real)
        #                          plt.plot(Z[bi,:,pi].imag)
        #                         plt.show()

                    for bi in range(n_beam):
                        if False:
                            plt.pcolormesh(Z[chi,bi,:,:].real)
                            plt.title("beam %d"%(bi))
                            plt.colorbar()
                            plt.show()
                        for ri in range(n_rg):
                            S[chi,bi,ri,:]+=n.abs(n.fft.fftshift(n.fft.fft(wf*Z[chi,bi,ri,:])))**2.0
    #            for bi in range(n_beam):
     #               for ri in range(n_rg):
      #                  S[bi,ri,:]=S[bi,ri,:]/n.median(n.abs(S[bi,ri,:]-n.median(S[bi,ri,:])))

                plt.figure(figsize=(16,9))
                plt.subplot(231)
                dB=10.0*n.log10(S[chi,0,:,:])
                plt.title("Beam 1")
                dB=dB-n.nanmedian(dB)
                plt.pcolormesh(fvec,rvec,dB,vmin=-3,vmax=20)
                plt.ylabel("Range (km)")
                plt.xlabel("Doppler (Hz)")
                plt.xlim([-200,200])

                plt.subplot(232)
                dB=10.0*n.log10(S[chi,1,:,:])
                dB=dB-n.nanmedian(dB)
             #   plt.ylim([0,40])
                plt.pcolormesh(fvec,rvec,dB,vmin=-3,vmax=20)
                plt.ylabel("Range (km)")
                plt.xlabel("Doppler (Hz)")

                plt.title("Beam 2")
                plt.xlim([-200,200])

                plt.subplot(233)
                dB=10.0*n.log10(S[chi,2,:,:])
                dB=dB-n.nanmedian(dB)
              #  plt.ylim([0,40])
                plt.title("Beam 3")
                plt.xlim([-200,200])

                plt.pcolormesh(fvec,rvec,dB,vmin=-3,vmax=20)
                plt.ylabel("Range (km)")
                plt.xlabel("Doppler (Hz)")

                plt.subplot(234)
                dB=10.0*n.log10(S[chi,3,:,:])
                dB=dB-n.nanmedian(dB)
            #    plt.ylim([0,40])
                plt.title("Beam 4")
                plt.xlim([-200,200])


                plt.pcolormesh(fvec,rvec,dB,vmin=-3,vmax=20)
                plt.ylabel("Range (km)")
                plt.xlabel("Doppler (Hz)")

                plt.subplot(235)
                dB=10.0*n.log10(S[chi,4,:,:])
                dB=dB-n.nanmedian(dB)
                plt.pcolormesh(fvec,rvec,dB,vmin=-3,vmax=20)
                plt.ylabel("Range (km)")
                plt.xlabel("Doppler (Hz)")

                plt.xlim([-200,200])

    #            plt.ylim([0,40])
                plt.title("Beam 5")

                plt.subplot(236)
                plt.title(stuffr.unix2datestr(m_start_idxs[0]/1e6))
                A=n.zeros([Z.shape[2],Z.shape[3]],dtype=n.float32)
                for bi in range(5):
                    A+=n.abs(Z[chi,bi,:,:])**2.0
                # range average a bit to get all details
                for ti in range(A.shape[1]):
                    A[:,ti] = n.convolve(n.repeat(1,20),A[:,ti],mode="same")
                nf=n.nanmedian(A[200:1400,:])
                A=(A-nf)/nf
                plt.pcolormesh(A,vmin=0,vmax=5)
    #            plt.ylabel("Range (km)")
    #            plt.xlabel("Doppler (Hz)")

    #            plt.pcolormesh(Z[0,:,:].real)
                plt.tight_layout()
                if alias:
                    plt.savefig("amrd-%d-%s.png"%(int(m_start_idxs[0]/1e6),ch))
                else:
                    plt.savefig("mrd-%d-%s.png"%(int(m_start_idxs[0]/1e6),ch))

                #plt.show()    
                plt.close()
                plt.clf()
