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


def tx_xphase():
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

    try:
        dmp = drf.DigitalMetadataReader(pc.phase_metadata_dir)
        dmp_b=dmp.get_bounds()
        start_idx=dmp_b[1]
    except:
        import traceback
        trackback.print_exc()

    n_block=int(n.ceil((b[1]-start_idx)/dt))

    channels=["ch000","ch001","ch002","ch003","ch004","ch005","ch006","ch007"]
    xphase=n.zeros(8,dtype=n.complex64)
    for bi in range(n_block):
        print(stuffr.unix2datestr((start_idx+bi*dt)/1e6))

        data=dm.read(start_idx+bi*dt,start_idx+bi*dt+1000000,"id")
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
                        if False:
                            #plt.subplot(121)
                            plt.plot((z0*n.conj(z1)).real)
                            plt.plot((z0*n.conj(z1)).imag)

                            #plt.plot(z0.imag)
                            #plt.subplot(122)
                            #plt.plot(z1.real)
                            #plt.plot(z1.imag)
                            plt.show()
                    if n.mean(n.abs(xphase))>1e7:
                        print("%s %1.1f %1.1f %1.1f %1.1f %1.1f %1.1f %1.1f %1.1f"%(stuffr.unix2datestr(k/1e6),n.angle(xphase[0]),n.angle(xphase[1]),n.angle(xphase[2]),n.angle(xphase[3]),n.angle(xphase[4]),n.angle(xphase[5]),n.angle(xphase[6]),n.angle(xphase[7])))
                        dout={"xphase":xphase}
                        dmw.write(k,dout)
                    else:
                        print("%s %1.1g"%(stuffr.unix2datestr(k/1e6),n.mean(n.abs(xphase))))

                except:
                    pass
#                    print()
#                    import traceback
#                    traceback.print_exc()
            


if __name__ == "__main__":
    while True:
        tx_xphase()
        time.sleep(3600)
