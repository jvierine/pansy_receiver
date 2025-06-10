import numpy as np
import healpy as hp
import matplotlib.pyplot as plt 
import h5py 
import numpy as n
import stuffr
import pansy_config as pc

import digital_rf as drf

from astropy.time import Time
from astropy.coordinates import get_sun, HeliocentricTrueEcliptic
import astropy.units as u

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()



def read_block(i0,i1,dm):
    d=dm.read(i0,i1)
    eclats=[]
    eclons=[]
    vgs=[]
    tv=[]
    h0=[]
    h1=[]    
    ew0=[]
    ew1=[]
    ns0=[]
    ns1=[]
    mf=[]
    snr=[]
    std=[]
    ewp=[]
    nsp=[]
    upp=[]
    dur=[]    
    

    for k in d.keys():
#        print(d[k].keys())
#        print(k)
        eclat=d[k]["eclat"]
        eclon=d[k]["eclon"]
        tunix=n.min(d[k]["txidx"])/1e6
        
        vg=n.linalg.norm(d[k]["v0"])
        h0.append(n.min(d[k]["up"]))
        h1.append(n.max(d[k]["up"]))
        
        ew0.append(n.min(d[k]["ew"]))
        ew1.append(n.max(d[k]["ew"]))
        
        ns0.append(n.min(d[k]["ns"]))
        ns1.append(n.max(d[k]["ns"]))
        
        eclats.append(eclat)
        eclons.append(eclon)

        mi=n.argmax(d[k]["snr"])
        snr.append(d[k]["snr"][mi])
        mf.append(d[k]["mfs"][mi])

        ewp.append(d[k]["ew"][mi])
        nsp.append(d[k]["ns"][mi])
        upp.append(d[k]["up"][mi])

        dur.append( (n.max(d[k]["txidx"])-n.min(d[k]["txidx"]))/1e6)
        
        vgs.append(vg)
        tv.append(tunix)
        
    return({"eclat":n.array(eclats),
            "eclon":n.array(eclons),
            "vg":n.array(vgs),
            "tunix":n.array(tv),
            "h0":n.array(h0),
            "h1":n.array(h1),
            "ew0":n.array(ew0),
            "ew1":n.array(ew1),
            "ns0":n.array(ns0),
            "ns1":n.array(ns1),
            "up0":n.array(h0),
            "up1":n.array(h1),
            "ewp":n.array(ewp),
            "nsp":n.array(nsp),            
            "upp":n.array(upp),
            "dur":n.array(dur),            
            "snr":n.array(snr),
            "mf":n.array(mf)}
           )

def radiant_dist(lats,lons,vgs,tv,title="Sun centered ecliptic",savefig=True,nside=16,fname="/tmp/latest_radiants.png"):
    # Convert latitude to colatitude (theta) in radians
    gidx=n.where(n.isnan(lats)!=True)[0]
    lats=lats[gidx]
    lons=lons[gidx]
    vgs=vgs[gidx]
    tv=tv[gidx]
    print(n.where(lats<=-90))
    print(n.where(lats>=90))
    theta = np.radians(90 - lats)  # Colatitude: 90° - latitude
    print(n.min(theta))
    print(n.max(theta))

    phi = np.radians(lons+90)  # Longitude in radians

    # Set HEALPix resolution (higher Nside means higher resolution)
    nside = nside  # Must be a power of 2

    # Convert lat/lon to HEALPix pixel indices
    pixels = hp.ang2pix(nside, theta, phi)

    # Create a HEALPix map and count occurrences in each pixel
    histogram = np.bincount(pixels, minlength=hp.nside2npix(nside))

    slon=solar_ecliptic_longitude(tv[0])
    mean_vel=n.zeros(len(histogram))

    for i in range(len(pixels)):
        mean_vel[pixels[i]]+=vg[i]
    mean_vel=mean_vel/histogram
    mean_vel[histogram<5]=n.nan
#    histogram[histogram<5]=1
    # Plot the histogram as a HEALPix map
    if True:
        hp.mollview(histogram, title=r"$\lambda_{\mathrm{sun}}=%1.1f^{\circ}$ %s"%(slon,title), unit="Counts",cmap="turbo",flip="geo",norm="linear",max=100,min=0)
        hp.projtext(-90., 0., '0°', lonlat=True, coord='geo',color="white")
        hp.projtext(0., 0., '270°', lonlat=True, coord='geo',color="white")
        hp.projtext(90., 0., '180°', lonlat=True, coord='geo',color="white")
        hp.graticule(color="white",alpha=0.2,dpar=10,verbose=True)

    # Set custom x-ticks (longitude)
    if False:
        plt.text(0,0,'270°',color="white")
        plt.text(-n.pi/2,0,'0°',color="white")
        plt.text(n.pi/2,0,'90°',color="white")

    if savefig:
        plt.savefig(fname)
        plt.close()
    else:
        plt.show()

#    titlestr=stuffr.unix2datestr((b[0])/1e6)
 #   lons=180*n.angle(n.exp(1j*n.pi*270/180)*n.exp(1j*-n.angle(n.exp(1j*n.pi*lons/180)*n.exp(-1j*n.pi*270/180))))/n.pi
  #  radiant_dist(lats,lons,vg,tv,title=titlestr,savefig=True,nside=32)


if __name__ == "__main__":
    mddir="/mnt/data/juha/pansy/metadata/simple_meteor_fit"
    dm = drf.DigitalMetadataReader(mddir)
    b = dm.get_bounds()
    bd=read_block(b[0],b[1],dm)
    ho=h5py.File("pansy_simple_pars.h5","w")
    for k in bd.keys():
        ho[k]=bd[k]
    ho.close()
