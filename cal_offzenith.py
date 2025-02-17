import numpy as n
import h5py
import matplotlib.pyplot as plt
import itertools
import scipy.optimize as so

def fit_beam_cal(xc0,xcn,ch_pairs):
    """
    xc0 = xc on zenith beam
    xcn = xc on a non-zenith beam
    phasecal = zenith beam calibration
    weight = goodness of xc
    ch_pairs = pairs of antenna channels cross-correlated

    Assume that z_a^1  = z_a^0 e^{i\phi^1_a}.
    """
    # each antenna needs to be multiplied with these phases (n.exp(1j*phasecal))
    h=h5py.File("data/phases.h5","r")
    phasecal=h["phasecal"][()]
    print(phasecal)
    h.close()
    def ss(phases):
        #phasecalz=n.exp(1j*(phasecal[ch_pairs[:,0]]-phasecal[ch_pairs[:,1]]))
        zenith_model = n.exp(1j*n.angle(xc0))*n.exp(1j*(phasecal[ch_pairs[:,0]]-phasecal[ch_pairs[:,1]]))
        #xprint(zenith_model)
        offzenith_model=zenith_model*n.exp(1j*(phases[ch_pairs[:,0]]-phases[ch_pairs[:,1]]))
        #print(offzenith_model)
        res=n.array(offzenith_model-n.exp(1j*n.angle(xcn)),dtype=n.complex128)
        #print(res)
        s=n.sum(res.real**2.0 + res.imag**2.0)
        #print(s)
        print(s)
        return(s)
    
    phases=so.fmin(ss,n.zeros(7,n.float32))
    phases=so.fmin(ss,phases)#n.zeros(7,n.float32))
    phases=so.fmin(ss,phases)#n.zeros(7,n.float32))

    zenith_model = n.exp(1j*n.angle(xc0))*n.exp(1j*(phasecal[ch_pairs[:,0]]-phasecal[ch_pairs[:,1]]))
    offzenith_model=zenith_model*n.exp(1j*(phases[ch_pairs[:,0]]-phases[ch_pairs[:,1]]))
    zd=n.angle(n.exp(1j*n.angle(xcn))*n.conj(offzenith_model))
    for i in range(21):
        plt.plot(zd[:,i],".")
        plt.show()
#    plt.plot(xc0.real,".")
 #   plt.plot(xc0.imag,".")
  #  plt.show()

h=h5py.File("caldata/cal_all.h5","r")
xc0=h["beam1/xc0"][()]
xcn=h["beam1/xcn"][()]
h.close()
#print(xcn.shape)
#xcn2=n.median(n.exp(1j*n.angle(xcn)),axis=0)
if False:
    for i in range(xcn.shape[1]):
        plt.plot(n.angle(xc0[:,i]*n.conj(xcn[:,i])),".")
        plt.show()

print(xcn.shape)
ch_pairs=n.array(list(itertools.combinations(n.arange(7,dtype=n.int64),2)))
fit_beam_cal(xc0,xcn,ch_pairs)
