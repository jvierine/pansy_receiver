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
    u0=u0/n.linalg.norm(u0)
    p0=jcoord.enu2ecef(pc.lat,pc.lon,0,p0[0],p0[1],p0[2])
    u0=jcoord.enu2ecef(pc.lat,pc.lon,0,u0[0],u0[1],u0[2])
  #  print(p0)
   # print(u0)
    
    # tbd: convert to astropy
    llh=jcoord.ecef2geodetic(p0[0], p0[1], p0[2])
#    print(llh)
    llh2=jcoord.ecef2geodetic(p0[0]-0.1*u0[0], p0[1]-0.1*u0[1], p0[2]-0.1*u0[2])
 #   print(llh2)
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

    return({"sky":sky,"lat":sc_gc_lat,"lon":sc_gc_lon, "sun_lon":sun_pos.lon.deg})