import numpy as n
import pansy_config as pc
import matplotlib.pyplot as plt
import digital_rf as drf
import stuffr
import numpy as n
import h5py

def plot_latest_fits(save_png=False):
    dm = drf.DigitalMetadataReader(pc.simple_fit_metadata_dir)
  #  dm = drf.DigitalMetadataReader("/tmp/simple_fit")

    b = dm.get_bounds()

    dd=dm.read(b[1]-24*3600*1000000,b[1])
    vgs=[]
    slats=[]
    slons=[]
    ea=[]
    no=[]
    sn=[]
    tv=[]

    for k in dd.keys():
        vgs.append(n.linalg.norm(dd[k]["v0"]))
        snr=dd[k]["snr"]
        ew=dd[k]["ew"]
        ns=dd[k]["ns"]
        mi=n.argmax(snr)
        ea.append(ew[mi])
        no.append(ns[mi])
        sn.append(snr[mi])
       # print(snr[mi])
        slats.append(dd[k]["eclat"])
        slons.append(dd[k]["eclon"])
        tv.append(k/1e6)

    slons=n.array(slons)
    slats=n.array(slats)
    vgs=n.array(vgs)
    tv=n.array(tv)
    ea=n.array(ea)
    no=n.array(no)
    sn=n.array(sn)

    ho=h5py.File("/tmp/fit_data.h5","w")
    ho["north"]=no
    ho["east"]=ea
    ho["snr"]=sn
    ho["vg"]=vgs
    ho["slon"]=slons
    ho["slat"]=slats
    ho["time"]=tv
    ho.close()
    print(len(sn))
    print(len(ew))
    print(len(ns))

    plt.figure(figsize=(8,4))
    if True:
        ax = plt.subplot(121, projection="lambert")
        sp=ax.scatter(n.angle(n.exp(1j*n.pi*slons/180.0)*n.exp(1j*n.pi/2)),n.pi*slats/180.0,c=vgs,vmin=10,vmax=72,s=0.5,cmap="turbo")
        ax.set_title("%s\n%s"%(stuffr.unix2datestr(n.min(tv)),stuffr.unix2datestr(n.max(tv))))
        #frame1 = plt.gca()
        ax.xaxis.set_ticklabels([])
        cb=plt.colorbar(sp)
        cb.set_label("Geocentric velocity (km/s)")
        ax.grid(True)
        ax.set_xlabel("Apex-centered ecliptic longitude (deg)")
        ax.set_ylabel("Ecliptic latitude (deg)")

    ax = plt.subplot(122)
    sp=ax.scatter(ea,no,c=10.0*n.log10(sn),s=0.5,cmap="gist_yarg",vmin=10,vmax=25)
    ax.set_xlim([-30,30])
    ax.set_ylim([-30,30])
    #cb=plt.colorbar(sp)
    #cb.set_label("SNR (dB)")
    ax.set_xlabel("East-West (km)")
    ax.set_ylabel("North-South (km)")

    plt.tight_layout()
    if save_png:
        plt.savefig("/tmp/latest_radiants.png")
        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    plot_latest_fits()
