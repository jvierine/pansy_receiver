import numpy as n
import digital_rf as drf
import matplotlib.pyplot as plt
import pansy_modes as pm

def plot_overview(dirname="test_data/test_data/mixmode",mmode_start=n.array([],dtype=int)):
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
    plt.title(dirname)
    plt.xlabel("Time (s)")
    plt.ylabel("Time (us)")
    plt.tight_layout()
    plt.savefig("%s/overview.png"%(dirname))
    plt.show()
    plt.close()

def find_stmode_start(d,i0=None,ch="ch000"):
    stm=pm.get_st_mode()
    codes=pm.get_vector(stm,ncodes=2)

    # this one looks at nulls
    rxp_window = n.concatenate((n.repeat(1,20),n.repeat(0,320),n.repeat(1,20),n.repeat(0,320)))
    codes=n.roll(codes,-20)
    plt.plot(codes.real)
    plt.show()
    N=len(rxp_window)

    b=d.get_bounds(ch)
    if i0==None:
        i0=b[0]

    R=n.conj(n.fft.fft(rxp_window,2*N))
    while (i0 < (b[1]-2*N)) and not_found:
        z=d.read_vector_c81d(i0,2*N,ch)
        pwr=z.real**2.0 + z.imag**2.0
        cc=n.abs(n.fft.ifft(n.fft.fft(pwr)*R))
        
        i0+=N

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



d=drf.DigitalRFReader("test_data/stmode")
start_idxs=find_stmode_start(d)
#plot_overview(dirname="test_data/stmode",mmode_start=start_idxs)

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