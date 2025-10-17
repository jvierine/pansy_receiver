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
import pansy_interferometry as pint
import h5py
import pansy_modes as pmm
import meteor_fit as metfit
#import simple_radiant as sr
import process_cut_meteor as pcm
import tensorflow as tf
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mc

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def process_latest():
#    mddir=pc.cut_metadata_dir
    mddir="/data1/pansy/metadata/cut"
    dm = drf.DigitalMetadataReader(mddir)
    b = dm.get_bounds()
    dt=120000000
#    os.system("mkdir -p caldata")

    start_idx=b[0]

    n_block=int(n.ceil((b[1]-start_idx)/dt))
    print(stuffr.unix2datestr(b[1]/1e6))
    print(stuffr.unix2datestr(start_idx/1e6))
    I=n.zeros([256,256,3],dtype=n.float32)

    cmap=matplotlib.colormaps["rainbow"]
    #    B=cmap(hv[sort_idx])[:,0:3]*vvv[:,None]
    for bi in range(rank,n_block,size):
        data=dm.read(start_idx+bi*dt,start_idx+bi*dt+dt)
        kl=list(data.keys())
        for ki in range(len(kl)):
            k=kl[ki]
            try:
                img,dimg=pcm.process_cut(data[k],None,do_cnn_image=True,plot=False,write_dm=False)
                dimg=-1*dimg
                dimg[dimg<0]=0
                dimg[dimg>72e3]=72e3
                dimg=(dimg)/72e3
#                plt.pcolormesh(dimg,cmap="turbo")
 #               plt.colorbar()
  #              plt.show()
                if True:
                    # dB intensity scale
                    print(img.shape)
                    dB=10*n.log10(img)
                    dB[n.isnan(dB)]=-3
                    dB[dB<-3]=-3
                    #                    dB.shape=(dB.shape[0],dB.shape[1],1)
                    v=(dB-n.min(dB))/(n.max(dB)-n.min(dB))

                   # linear intensity scale                    
#                    v=(img-n.min(img))/(n.max(img)-n.min(img))
                    
                    print(dB.shape)
                    print(v.shape)
#                    print(rgbim.shape)
                    rgbim=cmap(dimg)[:,:,0:3]*v[:,:,None]
                    
                    dB2=tf.image.resize(rgbim,(256,256),antialias=True)
#                    plt.imshow(dB2[::-1,:,:])
 #                   plt.show()
                    
                    if False:
                        plt.imshow(dB2[::-1,:,0])
                        plt.colorbar()
                        plt.show()

                    tf.keras.utils.save_img("/data1/pansy/cnn_images/cnn-%d.png"%(k),dB2[::-1,:,:],scale=True)
            except:
                import traceback
                traceback.print_exc()

process_latest()
