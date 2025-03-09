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

    # setup the directory and file cadence.
    # use 1 MHz, as this is the sample-rate and thus a
    # natural resolution for timing.
    subdirectory_cadence_seconds = 3600
    file_cadence_seconds = 60
    samples_per_second_numerator = 1000000
    samples_per_second_denominator = 1
    file_name = "txphase"
    os.system("mkdir -p %s"%(pc.phase_metadata_dir))

    dmw = drf.DigitalMetadataWriter(
        pc.phase_metadata_dir,
        subdirectory_cadence_seconds,
        file_cadence_seconds,
        samples_per_second_numerator,
        samples_per_second_denominator,
        file_name,
    )


    b = d.get_bounds("ch000")
    dt=60000000
    start_idx=b[0]
    n_block=int(n.ceil((b[1]-start_idx)/dt))

    channels=["ch000","ch001","ch002","ch003","ch004","ch005","ch006","ch007"]
    xphase=n.zeros(8,dtype=n.complex64)
    for bi in range(n_block):
        data=dm.read(start_idx+bi*dt,start_idx+bi*dt+dt,"id")
        kl=list(data.keys())
        if len(kl) > 0:
            ki=0
#        for ki in range(len(kl)):
            k=kl[ki]
            if data[k] == 1:
                try:
                    z0=d.read_vector_c81d(k,120,"ch000")
                    for j in range(8):                    
                        z1=d.read_vector_c81d(k,120,channels[j])
                        xphase[j]=n.mean(z0*n.conj(z1))
                    print(n.angle(xphase))
                    print(n.abs(xphase))
                    dout={"xphase":xphase}
                    dmw.write(k,dout)
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
            



