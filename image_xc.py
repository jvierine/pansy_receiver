import numpy as n
import matplotlib.pyplot as plt
import digital_rf as drf
import pansy_config as pc # type: ignore
import h5py
import pansy_interferometry as pint
import itertools
import rgbim

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

u,v,w=pint.uv_coverage(N=500,max_zenith_angle=10.0)
antpos=pint.get_antpos()
ch_pairs=n.array(list(itertools.combinations(n.arange(7),2)))
dmat=pint.pair_mat(ch_pairs,antpos)


for bi in range(n_blocks):
    dd=dmd.read(b[0]+bi*dt,b[0]+bi*dt+dt,["xc_arr","r0","r1","rdec","f0","f1","n_fft","ipp","ch_pairs"])
    for k in dd:
        f0=dd[k]["f0"]
        f1=dd[k]["f1"]
        fv=n.arange(f0,f1)
        n_f=len(fv)
        I=n.zeros([u.shape[0],u.shape[1]])#,n_f])
        #I[:,:]=0.0
        xca=dd[k]["xc_arr"]
        chp=dd[k]["ch_pairs"]
        xc=xca[7:xca.shape[0],0,:,:]
        xcc=n.copy(xc)
        print(chp)

        chp=chp[7:chp.shape[0],:]
        print(chp)
        #print(xcc.shape)
        print(chp.shape)
        for i in range(chp.shape[0]):
            xcc[i,:,:]=xc[i,:,:]*amps[chp[i][0]]*amps[chp[i][1]]*n.exp(-1j*(phasecal[chp[i][0]]-phasecal[chp[i][1]]))

        noise_floor=n.median(n.abs(xca[0,0,:,:]))
        snr=(n.abs(xca[0,0,:,:])-noise_floor)/noise_floor
#        gidx=n.where(snr>100)
        plt.pcolormesh(snr.T,vmin=0,vmax=100)
        plt.show()
        plt.pcolormesh(n.angle(xcc[0,:,:].T),cmap="hsv")
        plt.show()
        for fi in range(xc.shape[1]):
            if snr[fi,15]>30:
                M=pint.mf(xc[:,fi,15],dmat,u,v,w)
                plt.pcolormesh(n.abs(M))
                plt.show()
        #rgbim.rgb_image(I,peak_fraction=1.0)
        #plt.show()
#        plt.pcolormesh(n.abs(I))
 #       plt.colorbar()
  #      plt.savefig("img-%d.png"%(k/1e6))
   #     plt.close()
#        plt.show()
#        print(snr[fi,ri])
 #       print(xc.shape)
#        print(len(gidx))