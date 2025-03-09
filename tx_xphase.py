import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
import digital_rf as drf
import os
import stuffr
import time
import pansy_config as pc
import cluster_mf as cmf
import traceback
import scipy.fftpack as fp
import itertools
import pansy_modes as pmm


if __name__ == "__main__":
    mddir=pc.tx_metadata_dir
    dm = drf.DigitalMetadataReader(mddir)
    rdir=pc.raw_voltage_dir
    d = drf.DigitalRFReader(rdir)

    b = d.get_bounds("ch000")
    dt=10000000
    start_idx=b[0]
    n_block=int(n.ceil((b[1]-start_idx)/dt))

    subdirectory_cadence_seconds = 3600
    file_cadence_seconds = 60
    samples_per_second_numerator = 1000000
    samples_per_second_denominator = 1
    file_name = "tx_xphase"

    for bi in range(n_block):
        data=dm.read(start_idx+bi*dt,start_idx+bi*dt+dt,"id")
        kl=list(data.keys())
        for ki in range(len(kl)):
            k=kl[ki]
            if data[k] == 1:
                try:
                    z0=d.read_vector_c81d(k,120,"ch000")
                    z1=d.read_vector_c81d(k,120,"ch001")
                    print(n.angle(n.mean(z0*n.conj(z1))))
#                    plt.subplot(121)
 #                   plt.plot(z0.real)
  #                  plt.plot(z0.imag)
   #                 plt.subplot(122)
    #                plt.plot(z1.real)
     #               plt.plot(z1.imag)
      #              plt.show()
                except:
                    import traceback
                    traceback.print_exc()
            



