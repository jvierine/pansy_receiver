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

d=drf.DigitalRFReader("/media/archive")
b=d.get_bounds("ch000")
i0=b[0]+10000000
n_windows=2000
dt=int(n.floor((b[1]-i0)/n_windows))

window_len=1600

S=n.zeros([n_windows,window_len],dtype=n.float64)

tv=[]
for i in range(n_windows):
    print(i)
    read_idx=i0+dt
    z=d.read_vector_c81d(read_idx,window_len,"ch000")
    print(z[0])
    tv.append(n.datetime64(int(read_idx/1e6), 's'))
    S[i,:]=n.abs(z)

plt.pcolormesh(tv,n.arange(1600),S.T)
plt.colorbar()
plt.show()
plt.pcolormesh(S)
plt.colorbar()
plt.show()
