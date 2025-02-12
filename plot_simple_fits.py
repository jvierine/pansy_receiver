import numpy as n
import pansy_config as pc
import matplotlib.pyplot as plt
import digital_rf as drf
import stuffr

import jcoord
from astropy import units as u
from astropy.coordinates import AltAz
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
import astropy.coordinates as ac
from astropy.time import Time
import numpy as n

def get_radiant(p0,t0,u0):
    """
    initial position (enu)
    time (unix)
    u0 unit vector indication velocity vector of meteor (enu)
    """
    p0=jcoord.enu2ecef(pc.lat,pc.lon,0,p0[0],p0[1],p0[2])
    u0=jcoord.enu2ecef(pc.lat,pc.lon,0,u0[0],u0[1],u0[2])
    print(p0)
    print(u0)
    
    # tbd: convert to astropy
    llh=jcoord.ecef2geodetic(p0[0], p0[1], p0[2])
    print(llh)
    llh2=jcoord.ecef2geodetic(p0[0]-0.1*u0[0], p0[1]-0.1*u0[1], p0[2]-0.1*u0[2])
    print(llh2)
    aar=jcoord.geodetic_to_az_el_r(llh[0],llh[1],llh[2], llh2[0], llh2[1], llh2[2])
#    print(llh[2]/1e3)
    
    loc = EarthLocation(lat=llh[0]*u.deg, lon=llh[1]*u.deg, height=llh[2]*u.m)
    obs_time=Time(t0,format="unix")
    aa = AltAz(location=loc, obstime=obs_time)
    sky = SkyCoord(alt = aar[1]*u.deg, az = aar[0]*u.deg, obstime = obs_time, frame = 'altaz', location = loc)

    c_ecl=sky.transform_to('geocentricmeanecliptic')

    sp=ac.get_sun(obs_time)
    sun_pos=sp.transform_to('geocentricmeanecliptic')

    sc_gc_lon = 180.0*n.angle(n.exp(-1j*n.pi*sun_pos.lon.deg/180.0)*n.exp(1j*n.pi*c_ecl.lon.deg/180.0))/n.pi
    sc_gc_lat = c_ecl.lat.deg

    return(sky,sc_gc_lat,sc_gc_lon, sun_pos.lon.deg)

dm = drf.DigitalMetadataReader(pc.simple_fit_metadata_dir)
b = dm.get_bounds()
dt=100000000
n_block=int((b[1]-b[0])/dt)
hgts=[]
t0s=[]
vgs=[]
els=[]
azs=[]
slats=[]
slons=[]

for i in range(n_block):
    data=dm.read(b[0]+i*dt,b[0]+i*dt + dt)
    for k in data.keys():
#                    dout["r0"]=r0
 #           dout["v0"]=v0            
  #          dout["std"]=[eres,nres,ures]

        r0=data[k]["r0"]
        v0=data[k]["v0"]

        # radiant direction
        v1=-v0
        v_h=n.sqrt(v1[0]**2.0+v1[1]**2.0)
        el=180*n.arctan(n.abs(v1[2])/v_h)/n.pi
        az=180*n.arccos(v1[1]/v_h)/n.pi
        
        r,sc_lat,sc_lon,sun_lon=get_radiant(r0*1e3,k/1e6,v0/n.linalg.norm(v0))
        slats.append(sc_lat)
        slons.append(sc_lon)
        
 #       print(sc_lat)
#        print(sc_lon)
                
        std=data[k]["std"]
        if n.max(std*1e3) < 500.0:
            hgts.append(r0[2])
            vgs.append(n.linalg.norm(v0))
            t0s.append(stuffr.unix2date(k/1e6))
            els.append(el)
            azs.append(az)

plt.scatter(slons,slats,c=vgs,vmin=0,vmax=73)
plt.colorbar()
plt.show()
plt.scatter(t0s,azs,c=vgs,s=2,vmin=0,vmax=73)
plt.xlabel("Time (unix)")
plt.ylabel("Azimuth (deg)")
cb=plt.colorbar()
cb.set_label("$v_g$ (km/s)")
plt.show()

plt.scatter(t0s,els,c=vgs,s=2,vmin=0,vmax=73)
plt.xlabel("Time (unix)")
plt.ylabel("El (deg)")
cb=plt.colorbar()
cb.set_label("$v_g$ (km/s)")
plt.show()

            
plt.scatter(t0s,hgts,c=vgs,s=2,vmin=0,vmax=73)
plt.xlabel("Time (unix)")
plt.ylabel("Height (km)")
cb=plt.colorbar()
cb.set_label("$v_g$ (km/s)")
plt.show()

plt.plot(vgs,hgts,".")
plt.xlabel("Geocentric velocity (km/s)")
plt.ylabel("Height (km)")
plt.show()


