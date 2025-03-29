import numpy as n
import matplotlib.pyplot as plt
import digital_rf as drf
import h5py
h=h5py.File("data/mesocal.h5","r")
pwr=h["pwr"][()]
scale=n.real(n.sqrt(pwr[0:7])/n.sqrt(pwr[0:7]))
h.close()
print(scale)
d=drf.DigitalRFReader("/media/archive")
import sys
channels=["ch000","ch001","ch002","ch003","ch004","ch005","ch006","ch007"]
b=d.get_bounds("ch000")
print(b)
fig,((ax0,ax1),(ax2,ax3),(ax4,ax5),(ax6,ax7))=plt.subplots(4,2,sharex=True)
axs=[ax0,ax1,ax2,ax3,ax4,ax5,ax6,ax7]
for i in range(8):
    ch=channels[i]
    z=d.read_vector_c81d(b[1]-1000000,100000,ch)
    if i<7:
        z=z*scale[i]
    #plt.subplot(4,2,i+1)
    axs[i].set_ylabel(i+1)
    axs[i].plot(z[0:20000].real)
    axs[i].plot(z[0:20000].imag)
    axs[i].set_title("noise std=%1.1f"%(n.median(n.abs(z))))
plt.tight_layout()
plt.show()
