import numpy as n
import digital_rf as drf
import matplotlib.pyplot as plt
import pansy_modes as pm

def plot_overview(dirname="test_data/test_data/mixmode",mmode_start=n.array([],dtype=int),stmode_start=n.array([],dtype=int)):
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
    plt.plot( (mmode_start-b[0])/1e6, n.mod(mmode_start,1600), "." ,color="red")
    # plot detected sequence starts
    plt.plot( (stmode_start-b[0])/1e6, n.mod(stmode_start,1600), "." ,color="green")
    plt.title(dirname)
    plt.xlabel("Time (s)")
    plt.ylabel("Time (us)")
    plt.tight_layout()
    plt.savefig("%s/overview.png"%(dirname))
    plt.show()
    plt.close()

def find_stmode_start(d,i0=None,ch="ch000"):
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
    while (i0 < (b[1]-(160+1)*320)):
        z=d.read_vector_c81d(i0,(160+1)*320,ch)
        pwr=z.real**2.0 + z.imag**2.0
        pwr_cc=n.abs(n.fft.ifft(n.fft.fft(pwr)*R))
        tx_pwr=n.abs(n.fft.ifft(n.fft.fft(pwr)*T))
        code_ccc=n.fft.ifft(n.fft.fft(z)*C)
        code_cc=(code_ccc.real**2.0+code_ccc.imag**2.0)/tx_pwr

        peak_cc=n.max(code_cc)
        rxp_normalized=peak_cc*pwr_cc/n.max(pwr_cc)
        mf=code_cc - rxp_normalized
        mi=n.argmax(mf)
           
        if n.max(mf) > 90:
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

def find_mmode_start(d,i0=None,ch="ch000"):
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

    CAA=n.conj(n.fft.fft(codeaa,11*N))
    CBB=n.conj(n.fft.fft(codebb,11*N))
    P=n.conj(n.fft.fft(power,11*N))

    start_idxs=[]
    prev_idx=0
    plot=False
    while (i0 < (b[1]-12*N)):
        z=d.read_vector_c81d(i0,11*N,ch)
        Z=n.fft.fft(z)
        ccaaa=n.fft.ifft(Z*CAA)
        ccaa=ccaaa.real**2.0+ccaaa.imag**2.0
        ccbbb=n.fft.ifft(Z*CBB)
        ccbb=ccbbb.real**2.0+ccbbb.imag**2.0

        # regularize with +1 to avoid division by zero
        zp=z.real**2.0 + z.imag**2.0#n.abs(z)**2.0
        pwr=n.abs(n.fft.ifft(n.fft.fft(zp)*P))+1
        # the maxima of AA and minima of AB cross-correlations coincide
        # make use of this!
        pfr=(ccaa/pwr - ccbb/pwr)
        mi=n.argmax(pfr)
        #print(pfr[mi])        
        if pfr[mi]>150 and (i0+mi - prev_idx) != 0:
            # found new start
            print("%f %d %f"%( (i0-b[0])/1e6, i0+mi - prev_idx, pfr[mi]))
            # the idx points to the start of last ipp in a sequence of 20 ipps
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
        i0+=10*N
    return(n.array(start_idxs,dtype=int))


plot_overview(dirname="test_data/stmode")
d=drf.DigitalRFReader("test_data/stmode")
start_idxs=find_stmode_start(d)
plot_overview(dirname="test_data/stmode",stmode_start=start_idxs)

if False:
    d=drf.DigitalRFReader("test_data/stmode")
    b=d.get_bounds("ch000")
    z=d.read_vector_c81d(b[0],1600,"ch000")
    plt.plot(z.real)
    plt.plot(z.imag)
    plt.show()

d=drf.DigitalRFReader("test_data/mmode")
start_idxs=find_mmode_start(d)
plot_overview(dirname="test_data/mmode",mmode_start=start_idxs)

d=drf.DigitalRFReader("test_data/mixmode")
start_idxs=find_mmode_start(d)
plot_overview(dirname="test_data/mixmode",mmode_start=start_idxs)


#b=d.get_bounds("ch000")
#plt.plot((start_idxs-b[0])/1e6,n.mod(start_idxs,1600),".")
#plt.show()
#plt.plot(n.diff(start_idxs),".")
#plt.show()



if False:
    plot_overview(dirname="test_data/mixmode")
    plot_overview(dirname="test_data/mmode")
    plot_overview(dirname="test_data/stmode")