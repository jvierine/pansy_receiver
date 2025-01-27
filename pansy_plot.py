import numpy as n
import matplotlib.pyplot as plt
import digital_rf as drf

d=drf.DigitalRFReader("/media/buffer")
import sys
channels=["ch000","ch001","ch002","ch003","ch004","ch005","ch006","ch007"]
b=d.get_bounds("ch000")
print(b)
for i in range(8):
    ch=channels[i]
    z=d.read_vector_c81d(b[1]-1000000,100000,ch)
    plt.subplot(4,2,i+1)
    plt.ylabel(i+1)
    plt.plot(z[0:2000].real)
    plt.plot(z[0:2000].imag)
    plt.title("noise std=%1.1f"%(n.median(n.abs(z))))
plt.tight_layout()
plt.show()
