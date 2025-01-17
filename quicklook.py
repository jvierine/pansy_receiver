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

d=drf.DigitalRFReader("/media/archive")
b=d.get_bounds("ch007")

N=37500
dec=10
ipp=1600
N_dec=int(N/dec)

RTI=n.zeros([N_dec,ipp],dtype=n.float32)

n_blocks=int(n.floor((b[1]-b[0])/(N*ipp)))

b=d.get_bounds("ch007")
for i in range(n_blocks):
    print(i)
    z=d.read_vector_c81d(b[0]+i*N*ipp,N*ipp,"ch007")
    for j in range(N_dec):
        RTI[j,:]=n.abs(z[(j*ipp*dec):(j*ipp*dec+ipp)])
    plt.pcolormesh(RTI)
    plt.title(stuffr.unix2datestr( (b[0]+i*N*ipp)/1e6 ))
    plt.tight_layout()
    plt.savefig("ql-%d.png"%((b[0]+i*N*ipp)/1e6))
    plt.close()
#    plt.show()
                
