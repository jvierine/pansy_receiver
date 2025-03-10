import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
import digital_rf as drf
import os
import stuffr


def plot_blocks(block=3600*1000000):
    dm = drf.DigitalMetadataReader("/Users/j/src/pansy_test_data/metadata/simple_meteor_fit")
    b = dm.get_bounds()

    n_blocks=int(n.floor((b[1]-b[0])/block))

    for bi in range(n_blocks):
        ew=[]
        ns=[]
        up=[]
        snr=[]
        mfs=[]
        d=dm.read(b[0]+bi*block,b[0]+bi*block+block)

        for k in d.keys():
            print(k)
  #          vg=n.linalg.norm(d[k]["v0"])
 #           print(vg)
#            if n.max(d[k]["std"]) < 400.0 and vg > 10: 
            gidx=n.where(d[k]["mfs"]/21>0.1)[0]
            ew=n.concatenate((ew,d[k]["ew"][gidx]))
            ns=n.concatenate((ns,d[k]["ns"][gidx]))
            up=n.concatenate((up,d[k]["up"][gidx]))
            snr=n.concatenate((snr,d[k]["snr"][gidx]))
            mfs=n.concatenate((mfs,d[k]["mfs"][gidx]))
        if len(ew)>0:
            ew=n.array(ew)
            ns=n.array(ns)
            up=n.array(up)
            rng=n.sqrt(ew**2.0+ns**2.0+up**2.0)
            u=ew/rng
            v=ns/rng
     #       plt.hist2d(ns,up,bins=[n.linspace(-30,30,num=100),n.linspace(75,135,num=100)],weights=mfs)
    #        plt.title("Histogram of detections")
   #         plt.xlabel("North-South (km)")
  #          plt.ylabel("Up (km)")
 #           plt.tight_layout()
#            plt.show()


    #        plt.subplot(121)
            #plt.scatter(ew,ns,c=10.0*n.log10(snr),s=1,vmin=10,vmax=20)
            plt.figure(figsize=(16,8))
            plt.subplot(121)
            plt.scatter(u,v,c=mfs/21,vmin=0,vmax=1,s=1,cmap="viridis",alpha=1.0)#,vmin=10,vmax=20)
            plt.colorbar()
            plt.xlim([-0.3,0.3])
            plt.ylim([-0.3,0.3])

#            plt.xlim([-30,30])
 #           plt.ylim([-30,30])
            plt.title(stuffr.unix2datestr((b[0]+bi*block)/1e6))
            plt.xlabel("East-West (km)")
            plt.ylabel("North-South (km)")
            plt.tight_layout()
            plt.subplot(122)
#            plt.gca().set_aspect('equal') 
            plt.scatter(ns,up,c=mfs/21,vmin=0,vmax=1,s=1,cmap="viridis",alpha=1.0)#,vmin=10,vmax=20)
            plt.colorbar()
            plt.xlim([-30,30])
            plt.ylim([75,135])
#            plt.title(stuffr.unix2datestr((b[0]+bi*block)/1e6))
            plt.xlabel("North-South (km)")
            plt.ylabel("Up (km)")
            plt.tight_layout()
            plt.gca().set_aspect('equal') 
            plt.savefig("/tmp/uanim-%d.png"%(b[0]+bi*block))
            plt.close()
            plt.show()
 #           plt.show()

        #cb=plt.colorbar()
        #cb.set_label("SNR (dB)")
        #plt.show()
 #       plt.subplot(122)

plot_blocks()