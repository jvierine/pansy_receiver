import numpy as n
import matplotlib.pyplot as plt
import scipy.signal.windows as sw
import scipy.constants as c
import digital_rf as drf
import os
import stuffr
import time
import scipy.fftpack as fp
#import pyfftw 

rg=n.arange(400,6500,10)
n_rg=len(rg)
fftlen=1024
R=n.zeros([fftlen,n_rg])



# this is where the data is
d=drf.DigitalRFReader("/media/archive/")
# tx channel bounds
b=d.get_bounds("ch007")
n_samp=100*1000000
i0=b[1]-n_samp + 10736
ipp=12500
txlen=527
n_ipp=int(n.floor(n_samp/ipp))
for i in range(n_ipp):
    print(i)
    z=d.read_vector_c81d(i0+i*ipp,ipp,"ch000")
    #plt.plot(z.real)
    #plt.plot(z.imag)
    #plt.show()
    ztx=d.read_vector_c81d(i0+i*ipp,ipp,"ch007")
    z[0:(txlen+800)]=0
    for ri in range(n_rg):
        S=n.fft.fftshift(fp.fft(z[(rg[ri]):(rg[ri]+txlen)]*n.conj(ztx[0:txlen]),fftlen))
        R[:,ri]+=S.real**2.0 + S.imag**2.0
dop=n.fft.fftshift(n.fft.fftfreq(fftlen,d=1/1e6))
plt.pcolormesh(dop,rg*0.15,10.0*n.log10(R.T))
plt.colorbar()
plt.show()

