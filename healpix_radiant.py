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

def solar_ecliptic_longitude(unix_time):
    # Convert Unix time to Astropy Time object
    time = Time(unix_time, format='unix', scale='utc')
    
    # Get Sun's position in the sky
    sun = get_sun(time)

    sp=get_sun(time)
    sun_pos=sp.transform_to('geocentricmeanecliptic')
    
    # Return ecliptic longitude in degrees
    return(sun_pos.lon.deg)


def read_block(i0,i1,dm):
    d=dm.read(i0,i1)
    eclats=[]
    eclons=[]
    vgs=[]
    tv=[]

    for k in d.keys():
        print(k)
        eclat=d[k]["eclat"]
        eclon=d[k]["eclon"]
        tunix=n.min(d[k]["txidx"])/1e6
        vg=n.linalg.norm(d[k]["v0"])
        eclats.append(eclat)
        eclons.append(eclon)
        vgs.append(vg)
        tv.append(tunix)
    return(n.array(eclats),n.array(eclons),n.array(vgs),n.array(tv))

def radiant_dist(lats,lons,vgs,tv,title="Sun centered ecliptic",savefig=True,nside=16):
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
        hp.mollview(histogram, title=r"$\lambda_{\mathrm{sun}}=%1.1f^{\circ}$ %s"%(slon,title), unit="Counts",cmap="turbo",flip="geo",norm="linear")
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
        plt.savefig("/tmp/latest_radiants.png")
        plt.close()
    else:
        plt.show()

def plot_last_48():
    #simple_fit_metadata_dir
    dm = drf.DigitalMetadataReader(pc.simple_fit_metadata_dir)
    b = dm.get_bounds()
    lats,lons,vg,tv=read_block(b[1]-48*3600*1000000,b[1],dm)
    titlestr=stuffr.unix2datestr((b[0])/1e6)
    lons=180*n.angle(n.exp(1j*n.pi*270/180)*n.exp(1j*-n.angle(n.exp(1j*n.pi*lons/180)*n.exp(-1j*n.pi*270/180))))/n.pi
    radiant_dist(lats,lons,vg,tv,title=titlestr,savefig=True,nside=32)


if __name__ == "__main__":
    dm = drf.DigitalMetadataReader("/Users/j/src/pansy_test_data/metadata/simple_meteor_fit")
    b = dm.get_bounds()

    lats,lons,vg,tv=read_block(b[0],b[1],dm)
    #lats,lons,vg,tv=read_block(b[0],b[0]+24*3600*1000000,dm)
    #plt.hist(vg)
    #plt.show()
    print(len(lats))
    titlestr=stuffr.unix2datestr((b[0])/1e6)
    #plt.hist(lats)
    #plt.show()
    print("plot")
    # looks like lon is flipped aroud 270
    lons=180*n.angle(n.exp(1j*n.pi*270/180)*n.exp(1j*-n.angle(n.exp(1j*n.pi*lons/180)*n.exp(-1j*n.pi*270/180))))/n.pi

    radiant_dist(lats,lons,vg,tv,title=titlestr,savefig=False,nside=64)

    w=2*24*3600*1000000
    dt=24*3600*1000000

    n_days=int(n.floor((b[1]-b[0])/dt))
    for di in range(n_days):
        lats,lons,vg,tv=read_block(b[0]+di*dt,b[0]+di*dt+w,dm)
        titlestr=stuffr.unix2datestr((b[0]+di*dt)/1e6)
        #plt.hist(lats)
        #plt.show()
        print("plot")
        radiant_dist(lats,lons,vg,tv,title=titlestr)
        #lats = np.random.uniform(-90, 90, num_points)  # Latitude in degrees
        #lons = np.random.uniform(-180, 180, num_points)  # Longitude in degrees


#mean_vel[mean_vel<0]=0
#mean_vel[mean_vel>72]=72

#hp.mollview(mean_vel, title="Histogram of radiants", unit="Counts",cmap="turbo",flip="geo",norm="linear")
#hp.graticule(color="white",alpha=0.1)
#plt.show()