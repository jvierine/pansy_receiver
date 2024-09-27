import numpy as n
import matplotlib.pyplot as plt
import digital_rf as drf

d=drf.DigitalRFReader("/media/buffer")
import sys

b=d.get_bounds("ch000")
print(b)
z=d.read_vector_c81d(b[1]-1000000,100000,sys.argv[1])
plt.plot(z[0:10000].real)
plt.plot(z[0:10000].imag)
plt.show()
