import numpy as n
import pansy_config as pc
import matplotlib.pyplot as plt
import digital_rf as drf
import stuffr
import numpy as n
import h5py
import healpy as hp

from astropy.time import Time
from astropy.coordinates import get_sun, HeliocentricTrueEcliptic
import astropy.units as u

def solar_ecliptic_longitude(unix_time):
    # Convert Unix time to Astropy Time object
    time = Time(unix_time, format='unix', scale='utc')
    
    # Get Sun's position in the sky
    sun = get_sun(time)

    sp=get_sun(time)
    sun_pos=sp.transform_to('geocentricmeanecliptic')
    
    # Return ecliptic longitude in degrees
    return(sun_pos.lon.deg)

def plot_latest_fits(save_png=False):
    dm = drf.DigitalMetadataReader(pc.simple_fit_metadata_dir)
    #dm = drf.DigitalMetadataReader("../pansy_test_data/metadata/simple_meteor_fit")

    b = dm.get_bounds()

    dd=dm.read(b[1]-2*24*3600*1000000,b[1])
    vgs=[]
    slats=[]
    slons=[]
    ea=[]
    no=[]
    up=[]
    sn=[]
    tv=[]
    mfs=[]

    for k in dd.keys():
        vgs.append(n.linalg.norm(dd[k]["v0"]))
        snr=dd[k]["snr"]
        ew=dd[k]["ew"]
        ns=dd[k]["ns"]
        upp=dd[k]["up"]

        mf=dd[k]["mfs"]
        mi=n.argmax(snr)
        ea.append(ew[mi])
        no.append(ns[mi])
        up.append(upp[mi])
        sn.append(snr[mi])
        mfs.append(mf[mi])
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
    up=n.array(up)
    sn=n.array(sn)
    mfs=n.array(mfs)

    ho=h5py.File("/tmp/fit_data.h5","w")
    ho["north"]=no
    ho["east"]=ea
    ho["snr"]=sn
    ho["vg"]=vgs
    ho["slon"]=slons
    ho["slat"]=slats
    ho["time"]=tv
    ho["mfs"]=mfs

    ho.close()
    print(len(sn))
    print(len(ew))
    print(len(ns))


        #ax = plt.subplot(121)

    gidx=n.where(n.isnan(slats)!=True)[0]
    slats=slats[gidx]
    # flip slon.
    slons=180*n.angle(n.exp(-1j*n.pi*slons[gidx]/180))/n.pi
    theta = n.radians(90 - slats)  # Colatitude: 90° - latitude
    title="%s - %s"%(stuffr.unix2datestr(n.min(tv)),stuffr.unix2datestr(n.max(tv)))
    phi = n.radians(slons+90)  # Longitude in radians
    nside = 32  
    pixels = hp.ang2pix(nside, theta, phi)
    histogram = n.bincount(pixels, minlength=hp.nside2npix(nside))
    slon=solar_ecliptic_longitude(tv[0])
    hp.mollview(histogram, title=r"$\lambda_{\mathrm{sun}}=%1.1f^{\circ}$ %s"%(slon,title), unit="Counts",cmap="turbo",flip="astro",norm="linear")

    hp.projtext(-90., 0., '0°', lonlat=True, coord='astro',color="white")
    hp.projtext(0., 0., '270°', lonlat=True, coord='astro',color="white")
    hp.projtext(90., 0., '180°', lonlat=True, coord='astro',color="white")
    hp.graticule(color="white",alpha=0.2,dpar=10,verbose=True)

    #hp.graticule(color="white",alpha=0.1)
    plt.savefig("/tmp/latest_radiants.png")
    plt.close()
    
    plt.figure(figsize=(8,4))
    if True:
        ax = plt.subplot(121, projection="lambert")
        sp=ax.scatter(n.angle(n.exp(-1j*n.pi*slons/180.0)*n.exp(1j*n.pi/2)),n.pi*slats/180.0,c=vgs,vmin=10,vmax=72,s=0.5,cmap="turbo")
        ax.set_title("%s\n%s"%(stuffr.unix2datestr(n.min(tv)),stuffr.unix2datestr(n.max(tv))))
        #frame1 = plt.gca()
        ax.xaxis.set_ticklabels([])
        cb=plt.colorbar(sp)
        cb.set_label("Geocentric velocity (km/s)")
        ax.grid(True)
        ax.set_xlabel("Apex-centered ecliptic longitude (deg)")
        ax.set_ylabel("Ecliptic latitude (deg)")

    ax = plt.subplot(122)
    L=n.sqrt(up**2.0+ea**2.0+no**2.0)
    u=ea/L
    v=no/L

    sp=ax.scatter(u,v,c=mfs/21,s=0.5,cmap="brg",vmin=0.5,vmax=1,alpha=0.2)

    ax.set_xlim([-0.3,0.3])
    ax.set_ylim([-0.3,0.3])
    cb=plt.colorbar(sp)
    cb.set_label("Match function")
    ax.set_aspect("equal")
    ax.set_xlabel("East-West (unit)")
    ax.set_ylabel("North-South (unit)")

    plt.tight_layout()
    if save_png:
        plt.savefig("/tmp/latest_hist.png")
        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    plot_latest_fits()
