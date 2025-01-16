import numpy as n
import digital_rf as drf
import matplotlib.pyplot as plt
import pansy_modes as pm
import scipy.signal.windows as sw
import scipy.constants as c
import stuffr
import time

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def plot_overview(dirname="test_data/test_data/mixmode",m_mode_start=n.array([],dtype=int),st_mode_start=n.array([],dtype=int)):
    d=drf.DigitalRFReader(dirname)

    b=d.get_bounds("ch000")

    L=b[1]-b[0]
    n_ipps=int(n.floor(L/1600))
    N_max=4000
    step=int(n.floor(n_ipps/N_max))
    print(step)
    S=n.zeros([N_max,1600],dtype=n.float32)
    for i in range(N_max):
        print(b[0]+i*1600*step)
        z=d.read_vector_c81d(b[0]+i*1600*step,1600,"ch000")
   
        S[i,:]=n.abs(z)**2.0

    tvec=n.arange(N_max)*step*1600/1e6
    ftime=n.arange(1600)
    dB=10.0*n.log10(S.T)
    dB=dB-n.nanmedian(dB)
    plt.figure(figsize=(16,9))
    plt.pcolormesh(tvec,ftime,dB,vmin=-3)
    cb=plt.colorbar()
    cb.set_label("Power (dB)")
    # plot detected sequence starts
    plt.plot( (m_mode_start-b[0])/1e6, n.mod(m_mode_start,1600), "." ,color="red")
    # plot detected sequence starts
    plt.plot( (st_mode_start-b[0])/1e6, n.mod(st_mode_start,1600), "." ,color="white")
    plt.title(dirname)
    plt.xlabel("Time (s)")
    plt.ylabel("Time (us)")
    plt.tight_layout()
    plt.savefig("%s/overview.png"%(dirname))
    plt.show()
    plt.close()

def find_st_mode_start(d,i0=None,i1=None,ch="ch000", rxp_max=3e8):
    """ 
    use 10 first pulses to find the start of sequence 
    it appears that the 10 first pulses does not repeat otherwise.
    also detect the minima before the transmit pulse, which is due to the 
    receiver protection. this seems to avoid strong airplanes echoes being detected
    (not guaranteed)
    """
    stm=pm.get_st_mode()
    ncodes=10
    codes=pm.get_vector(stm,ncodes=ncodes)

    # this one looks at nulls
    rxp_window = n.tile(n.concatenate((n.repeat(1,16),n.repeat(0,304))),ncodes)
    # add 20 samples of space in front of tx pulse to align with rxp protect window
    codes=n.roll(codes,20)
    tx_window = n.copy(codes)
    tx_window[n.abs(tx_window)>0.1]=1.0  
    N=len(rxp_window)

    b=d.get_bounds(ch)
    if i0==None:
        i0=b[0]
    if i1==None:
        i1=(b[1]-(160+1)*320)

    if False:
        plt.subplot(211)
        plt.plot(rxp_window.real)
        z=d.read_vector_c81d(i0+80,ncodes*320,ch)
        plt.plot(codes.real)
        plt.subplot(212)
        plt.plot(z.real)
        plt.plot(z.imag)
        plt.plot(1/n.abs(z))
        plt.show()

    R=n.conj(n.fft.fft(rxp_window,(160+1)*320))
    C=n.conj(n.fft.fft(codes,(160+1)*320))
    T=n.conj(n.fft.fft(tx_window,(160+1)*320))

    st_idxs=[]
    prev_idx=0
    while i0 < i1:
        z=d.read_vector_c81d(i0,(160+1)*320,ch)
        pwr=z.real**2.0 + z.imag**2.0
        pwr_cc=n.abs(n.fft.ifft(n.fft.fft(pwr)*R))
        tx_pwr=n.abs(n.fft.ifft(n.fft.fft(pwr)*T))
        code_ccc=n.fft.ifft(n.fft.fft(z)*C)
        code_cc=(code_ccc.real**2.0+code_ccc.imag**2.0)/tx_pwr

        # are we in a region where the rx protect is low enough
        mask=(pwr_cc < rxp_max)
        #peak_cc=n.max(code_cc)
        #plt.plot(pwr_cc)
        #plt.show()    
        #rxp_normalized=peak_cc*pwr_cc/n.max(pwr_cc)
        mf=code_cc*mask# - rxp_normalized
        mi=n.argmax(mf)
           
        if n.max(mf) > 100:
            if False:
                print("%f %f %d"%( (i0-b[0])/1e6, mf[mi], mi+i0 - prev_idx))
            st_idxs.append(mi+i0)
            prev_idx=mi+i0
            #z=d.read_vector_c81d(i0+mi,1600,"ch000")
            #plt.plot(z[0:64].real)
            #plt.plot(z[0:64].imag)
            #plt.plot(n.abs(z[0:64]))
            #plt.show()

            if False:
                plt.plot(peak_cc*pwr_cc/n.max(pwr_cc))
                plt.plot(code_cc)
                plt.plot(mf)
                plt.show()
        i0+=160*320
    return(n.array(st_idxs,dtype=int))

def find_m_mode_start(d,
                      i0=None,
                      i1=None,
                      ch="ch007", # channel 007 is the transmit sample
                      debug=False):
    # look for AA, which happens at the end of the cycle!
    # use the fact that the cross-correlation AB is zero
    # where AA maximizes
    # todo: estimate pulse-to-pulse phase change
    stm=pm.get_m_mode()
    codes=pm.get_vector(stm,ncodes=2)
    codeaa=n.concatenate((codes[0:1600],codes[0:1600]))
    codebb=n.concatenate((codes[1600:3200],codes[1600:3200]))
    # look for A A. This only happens at end of cycle 
    N=len(codeaa)

    # select power
    power=n.copy(codeaa)
    power[n.abs(codeaa)>0.1]=1.0

    b=d.get_bounds(ch)
    if i0==None:
        i0=b[0]
    if i1==None:
        i1=(b[1]-12*N)
        #        i1=(b[1]-(160+1)*320)
    CAA=n.conj(n.fft.fft(codeaa,11*N))
    CBB=n.conj(n.fft.fft(codebb,11*N))
    P=n.conj(n.fft.fft(power,11*N))

    start_idxs=[]
    prev_idx=0
    plot=False
    while i0 < i1:
        z=d.read_vector_c81d(i0,11*N,ch)
        Z=n.fft.fft(z)
        ccaaa=n.fft.ifft(Z*CAA)
        ccaa=ccaaa.real**2.0+ccaaa.imag**2.0
        ccbbb=n.fft.ifft(Z*CBB)
        ccbb=ccbbb.real**2.0+ccbbb.imag**2.0

        # regularize with +1 to avoid division by zero
        zp=z.real**2.0 + z.imag**2.0#n.abs(z)**2.0
        # correlate power pattern
        pwr=n.abs(n.fft.ifft(n.fft.fft(zp)*P))+1
        
        # the maxima of AA and minima of AB cross-correlations coincide
        # make use of this!
        pfr=(ccaa/pwr - ccbb/pwr)
        mi=n.argmax(pfr)
        #        plt.plot(pfr)
        #        plt.show()
#        if debug:
#            print(pfr[mi])        
        if pfr[mi]>150 and (i0+mi - prev_idx) != 0:
            # found new start
            if debug:
                print("%f %d %f"%( (i0-b[0])/1e6, i0+mi - prev_idx, pfr[mi]))
            # 
            # the idx points to the start of last ipp in a sequence of 20 ipps. we need to adjust to the start of the first ipp.
            start_idxs.append(i0+mi - 19*1600)
            if plot:
                plt.plot(pfr,label=r"$m_t$")
                plt.plot(ccaa/pwr,label=r"$c^A_{-t}*z_t/P_t$")
                plt.plot(ccbb/pwr,label=r"$c^B_{-t}*z_t/P_t$")
                plt.xlabel("Time (samples)")
                plt.legend()
                plt.show()
                z=d.read_vector_c81d(i0+mi,N,ch)
                plt.subplot(211)
                plt.plot(z.real*codeaa.real)
                plt.plot(z.imag*codeaa.real)
                plt.subplot(212)
                plt.plot(z.real*codebb.real)
                plt.plot(z.imag*codebb.real)
                plt.show()
            prev_idx=i0+mi
        # go forward by one code sequence
        i0+=10*N
    return(n.array(start_idxs,dtype=int))

def analyze_st_mode(d,
                    dt=10,
                    b0=None,
                    b1=None,
                    debug_phase=False):
    b=d.get_bounds("ch000")
    if b0==None:
        b0=b[0]
    if b1==None:
        b1=b[1]        

    stm=pm.get_st_mode()
    ncodes=10
    codes=pm.get_vector(stm,ncodes=len(stm["codes"]))
    codes.shape=(160,320) 
    C=n.copy(codes)
    
    for i in range(codes.shape[0]):
        C[i,:]=n.conj(n.fft.fft(codes[i,:]))

    # how many integration periods
    n_ints=int(n.floor((b1-b0)/1e6/dt))
    for i in range(rank,n_ints,size):
        i0=b0+i*dt*1000000
        i1=b0+(i+1)*dt*1000000
        st_start_idxs=find_st_mode_start(d,i0,i1)
        n_beam=len(stm["beam_pos_az_za"])
        n_rg=320
        n_pulse=32
        n_reps=32
        wf=sw.hann(n_reps*n_pulse)
        Z=n.zeros([n_beam,n_rg,n_pulse*n_reps],dtype=n.complex64)
        S=n.zeros([n_beam,n_rg,n_pulse*n_reps],dtype=n.float32)

        rvec=n.arange(n_rg)*0.15
        # 2*f*v/c = df
        # df*c/f/2 = v
        fvec=n.fft.fftshift(n.fft.fftfreq(n_reps*n_pulse,d=n_beam*320/1e6))*c.c/47.5e6/2.0
        if len(st_start_idxs) > 0:
            print("%1.1f found %d sequences"%(i*dt,len(st_start_idxs)))
            n_cycles=int(n.floor(len(st_start_idxs)/n_reps))
            print(n_cycles)
            for ci in range(n_cycles):
                print("cycle %d"%(ci))
                for ri in range(n_reps):
                    i0=st_start_idxs[ci*n_reps + ri]
                    for pi in range(n_pulse):
                        for bi in range(n_beam):
                            iread = i0 + pi*320*n_beam + bi*320 + 20
                            z=d.read_vector_c81d(iread,320,"ch000")
                            codei = (pi*n_beam + bi)%160
                            Z[bi,:,ri*n_pulse + pi]=n.fft.ifft(n.fft.fft(z)*C[codei,:])
                            if debug_phase:
                                plt.subplot(211)
                                plt.plot(Z[bi,:,ri*n_pulse + pi].real)
                                plt.plot(Z[bi,:,ri*n_pulse + pi].imag)
                                plt.xlabel(r"Time ($\mu s$)")
                                plt.title("Uncorrected phase")
                                plt.ylabel("Complex voltage")
                            # phase shift based on transmit pulse phase 
                            phase = n.angle(Z[bi,0,ri*n_pulse + pi])
                            # rx blank
#                            z[0:40]=0.0
 #                           Z[bi,:,ri*n_pulse + pi]=n.fft.ifft(n.fft.fft(z)*C[codei,:])

                            Z[bi,:,ri*n_pulse + pi]=Z[bi,:,ri*n_pulse + pi]*n.exp(-1j*phase)
                            if debug_phase:
                                plt.subplot(212)
                                plt.plot(Z[bi,:,ri*n_pulse + pi].real)
                                plt.plot(Z[bi,:,ri*n_pulse + pi].imag)
                                plt.xlabel(r"Time ($\mu s$)")
                                plt.title("Corrected phase")
                                plt.ylabel("Complex voltage")
                                plt.tight_layout()
                                plt.show()

    #                        if True:
    #                           plt.plot(Z[bi,:,pi].real)
    #                          plt.plot(Z[bi,:,pi].imag)
    #                         plt.show()

                for bi in range(n_beam):
                    if False:
                        plt.pcolormesh(Z[bi,:,:].real)
                        plt.title("beam %d"%(bi))
                        plt.colorbar()
                        plt.show()
                    for ri in range(n_rg):
                        S[bi,ri,:]+=n.abs(n.fft.fftshift(n.fft.fft(wf*Z[bi,ri,:])))**2.0
#            for bi in range(n_beam):
 #               for ri in range(n_rg):
  #                  S[bi,ri,:]=S[bi,ri,:]/n.median(n.abs(S[bi,ri,:]-n.median(S[bi,ri,:])))

            plt.figure(figsize=(16,9))
            plt.subplot(231)
            dB=10.0*n.log10(S[0,:,:])
            plt.title("Beam 1")
            dB=dB-n.nanmedian(dB)
            plt.pcolormesh(fvec,rvec,dB,vmin=-3,vmax=20)
            plt.xlim([-200,200])
            plt.ylim([0,40])

            plt.subplot(232)
            dB=10.0*n.log10(S[1,:,:])
            dB=dB-n.nanmedian(dB)
            plt.xlim([-200,200])
            plt.ylim([0,40])
            plt.pcolormesh(fvec,rvec,dB,vmin=-3,vmax=20)
            plt.title("Beam 2")
            
            plt.subplot(233)
            dB=10.0*n.log10(S[2,:,:])
            dB=dB-n.nanmedian(dB)
            plt.xlim([-200,200])
            plt.ylim([0,40])
            plt.title("Beam 3")

            plt.pcolormesh(fvec,rvec,dB,vmin=-3,vmax=20)
            plt.subplot(234)
            dB=10.0*n.log10(S[3,:,:])
            dB=dB-n.nanmedian(dB)
            plt.xlim([-200,200])
            plt.ylim([0,40])
            plt.title("Beam 4")


            plt.pcolormesh(fvec,rvec,dB,vmin=-3,vmax=20)
            plt.subplot(235)
            dB=10.0*n.log10(S[4,:,:])
            dB=dB-n.nanmedian(dB)
            plt.pcolormesh(fvec,rvec,dB,vmin=-3,vmax=20)
            plt.xlim([-200,200])
            plt.ylim([0,40])
            plt.title("Beam 5")
            plt.subplot(236)
            plt.title(stuffr.unix2datestr(st_start_idxs[0]/1e6))

            plt.pcolormesh(Z[1,:,:].real,vmin=-50e3,vmax=50e3)

            plt.tight_layout()
            plt.savefig("strd-%d.png"%(int(st_start_idxs[0]/1e6)))
            #plt.show()    
            plt.close()
            plt.clf()


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
            if len(m_start_idxs) > 2:
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

def test_mode_detection():
    # test mode finding with the test data
    d=drf.DigitalRFReader("test_data/mixmode")
    st_start_idxs=find_st_mode_start(d)
    m_start_idxs=find_m_mode_start(d)
    plot_overview(dirname="test_data/mixmode", m_mode_start=m_start_idxs, st_mode_start=st_start_idxs)

    d=drf.DigitalRFReader("test_data/stmode")
    start_idxs=find_st_mode_start(d)
    plot_overview(dirname="test_data/stmode",st_mode_start=start_idxs)

    d=drf.DigitalRFReader("test_data/mmode")
    start_idxs=find_m_mode_start(d)
    plot_overview(dirname="test_data/mmode",m_mode_start=start_idxs)


if __name__ == "__main__":
    if False:
        d=drf.DigitalRFReader("/media/archive/")
        b=d.get_bounds("ch007")
        # 10 seconds later then start of buffer
        idx_end=b[1]-3600*1000000
        print(b)
        while True:
            b=d.get_bounds("ch007")
            if idx_end < b[0]:
                print("cannot keep up. adjust start")
                idx_end=b[0]+10000000
            i1=idx_end+60000000
            if i1>b[1]:
                i1=b[1]
            start_idx=find_m_mode_start(d,
                                        i0=idx_end,
                                        i1=i1,
                                        ch="ch007", # channel 007 is the transmit sample
                                        debug=True)
            print(len(start_idx))
            if len(start_idx) == 0:
                idx_end=i1
            else:
                # we found 
                idx_end=start_idx[-1]+1600*10

            print(idx_end)
            time.sleep(1)
    
    d=drf.DigitalRFReader("/media/archive/")
    b=d.get_bounds("ch000")
    analyze_m_mode(d,b0=b[1]-3600*1000000,b1=b[1])
    
    b=d.get_bounds("ch000")
    print(b)
    analyze_st_mode(d,b0=b[1]-3600*1000000,b1=b[1])

    
    #
    #d=drf.DigitalRFReader("test_data/mmode")
    #analyze_m_mode(d)
    #
    #d=drf.DigitalRFReader("test_data/mixmode")
    #analyze_m_mode(d)
    #analyze_st_mode(d)





        
