import numpy as n
import matplotlib.pyplot as plt
import digital_rf as drf
import pansy_config as pc # type: ignore
import h5py

# each antenna needs to be multiplied with these phases (n.exp(-1j*phasecal))
h=h5py.File("data/phases.h5","r")
phasecal=h["zenith"][()]
h.close()
# each antenna needs to be multiplied with these amplitudes
h=h5py.File("data/mesocal.h5","r")
amps=h["pwr"][()]
amps=1/n.sqrt(amps[0:7])
# ch_pairs[0] = (0,1) => z_0 * z_1^* = xc[0]*amps[0]*amps[1]*n.exp(-1j*(phasecal[0]-phasecal[1]))
print(phasecal)

dmd = drf.DigitalMetadataReader("../pansy_test_data/metadata/xc2/")
b = dmd.get_bounds()
dt=60*1000000

n_blocks=int((b[1]-b[0])/dt)

for bi in range(n_blocks):
    dd=dmd.read(b[0]+bi*dt,b[0]+bi*dt+dt,["xc_arr","r0","r1","rdec","f0","f1","n_fft","ipp","ch_pairs"])
    for k in dd:
        xca=dd[k]["xc_arr"]
        chp=dd[k]["ch_pairs"]
        xc=xca[7:xca.shape[0],0,:,:]
        xcc=n.copy(xc)
        chp=chp[7:chp.shape[0],:]
        #print(xcc.shape)
        print(chp.shape)
        for i in range(chp.shape[0]):
            xcc[i,:,:]=xc[i,:,:]*amps[chp[i][0]]*amps[chp[i][1]]*n.exp(-1j*(phasecal[chp[i][0]]-phasecal[chp[i][1]]))

        noise_floor=n.median(n.abs(xca[0,0,:,:]))
        snr=(n.abs(xca[0,0,:,:])-noise_floor)/noise_floor
        gidx=n.where(snr>100)
        for i in range(len(gidx[0])):
            fi=gidx[0][i]
            ri=gidx[1][i]
            plt.pcolormesh(n.angle(xcc[i,:,:]),cmap="hsv")
            plt.show()
            print(snr[fi,ri])
            print(xc.shape)
#        print(len(gidx))