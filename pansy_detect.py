import numpy as n
import digital_rf as drf
import matplotlib.pyplot as plt

def plot_overview(dirname="test_data/test_data/mixmode"):
    d=drf.DigitalRFReader(dirname)

    b=d.get_bounds("ch000")

    L=b[1]-b[0]
    n_ipps=int(n.floor(L/1600))
    N_max=4000
    step=int(n.floor(n_ipps/N_max))
    print(step)
    S=n.zeros([N_max,1600],dtype=n.float32)
    for i in range(N_max):
        print(b[0]+i*1600*step)
        z=d.read_vector_c81d(b[0]+i*1600*step,1600,"ch000")
        S[i,:]=n.abs(z)**2.0

    tvec=n.arange(N_max)*step*1600/1e6
    ftime=n.arange(1600)
    dB=10.0*n.log10(S.T)
    dB=dB-n.nanmedian(dB)
    plt.figure(figsize=(16,9))
    plt.pcolormesh(tvec,ftime,dB,vmin=-3)
    cb=plt.colorbar()
    cb.set_label("Power (dB)")

    plt.title(dirname)
    plt.xlabel("Time (s)")
    plt.ylabel("Time (us)")
    plt.tight_layout()
    plt.savefig("%s/overview.png"%(dirname))
    plt.show()
    plt.close()

plot_overview(dirname="test_data/mixmode")
plot_overview(dirname="test_data/mmode")
plot_overview(dirname="test_data/stmode")