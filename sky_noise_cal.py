import numpy as n
import matplotlib.pyplot as plt
from pygdsm import GlobalSkyModel
import pansy_config as pc
from datetime import datetime

from pygdsm import GSMObserver
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.time import Time
import astropy.units as u
import h5py

import pansy_modes as pm
mmode=pm.get_m_mode()
beam_pos=mmode["beam_pos_az_za"]


# Inputs
azimuth = 120 * u.deg        # Azimuth (measured from North towards East)
elevation = 45 * u.deg       # Elevation (Altitude)
latitude = pc.lat           # Example: Tromsø, Norway
longitude = pc.lon
elevation_m = 100              # Optional: Observer's elevation in meters
unix_time = 1718112000       # Example UNIX timestamp (2024-06-11 00:00:00 UTC)

# Define observer location and time
location = EarthLocation(lat=latitude * u.deg, lon=longitude * u.deg, height=elevation_m * u.m)

# sidereal day in seconds
d_sidereal=86164.0905

gsm = GlobalSkyModel()
gsm.generate(47)

def az_el_to_unit_vector(az_deg, el_deg):
    # Convert degrees to radians
    az = n.deg2rad(az_deg)
    el = n.deg2rad(el_deg)

    # Spherical to Cartesian conversion
    x = n.cos(el) * n.sin(az)  # East component
    y = n.cos(el) * n.cos(az)  # North component
    z = n.sin(el)               # Up component

    return n.array([x, y, z])

def rotate_vector(v, angle1_deg, axis1, angle2_deg, axis2):
    # Convert angles to radians
    theta1 = n.deg2rad(angle1_deg)
    theta2 = n.deg2rad(angle2_deg)

    # Define rotation matrices
    def rotation_matrix(axis, theta):
        if axis == 'x':
            return n.array([
                [1, 0, 0],
                [0, n.cos(theta), -n.sin(theta)],
                [0, n.sin(theta),  n.cos(theta)]
            ])
        elif axis == 'y':
            return n.array([
                [n.cos(theta), 0, n.sin(theta)],
                [0, 1, 0],
                [-n.sin(theta), 0, n.cos(theta)]
            ])
        elif axis == 'z':
            return n.array([
                [n.cos(theta), -n.sin(theta), 0],
                [n.sin(theta),  n.cos(theta), 0],
                [0, 0, 1]
            ])
        else:
            raise ValueError("Axis must be 'x', 'y', or 'z'")

    # Apply first rotation
    R1 = rotation_matrix(axis1, theta1)
    v_rot1 = R1 @ v

    # Apply second rotation
    R2 = rotation_matrix(axis2, theta2)
    v_rot2 = R2 @ v_rot1

    return v_rot2

def unit_vector_to_az_el(v):
    x, y, z = v

    # Normalize (just in case)
    v_norm = v / n.linalg.norm(v)
    x, y, z = v_norm

    # Elevation (altitude) angle
    el = n.arcsin(z)

    # Azimuth (from North, increasing toward East)
    az = n.arctan2(x, y)  # Note: atan2(x, y) gives az from North

    # Convert to degrees
    az_deg = n.degrees(az) % 360
    el_deg = n.degrees(el)

    return az_deg, el_deg

def get_noise(az,el,unix_time,bw=5):
    time = Time(unix_time, format='unix')

    # AltAz frame at given time and location
    altaz_frame = AltAz(obstime=time, location=location)
    # Create SkyCoord in AltAz
    az_del=n.linspace(-bw,bw,num=10)
    el_del=n.linspace(-bw,bw,num=10)
    temps=[]
    for i in range(len(az_del)):
        for j in range(len(el_del)):
            v=az_el_to_unit_vector(az, el)
            v_rot = rotate_vector(v, az_del[i], 'x', el_del[j], 'y')
            azr,elr=unit_vector_to_az_el(v_rot)
#            print(azr,elr)
 #           print(az,el)
            altaz_coord = SkyCoord(az=azr*u.deg, alt=elr*u.deg, frame=altaz_frame)
            temps.append(gsm.get_sky_temperature(altaz_coord))
    return(n.mean(temps))


def calc_beam_lookup(n_times=48*2,t0=1718112000):
    temp=n.zeros([5,n_times])
    tsid=n.arange(n_times)*(d_sidereal/n_times)
    for i in range(5):
        for ti in range(n_times):
            tnow=t0+tsid[ti]
            temp[i,ti]=get_noise(beam_pos[i][0],90-beam_pos[i][1],tnow)
#        plt.plot(tsid,temp[i,:])
 #       plt.show()
    ho=h5py.File("data/noise_lookup.h5","w")
    ho["temp"]=temp
    ho["t0"]=t0
    ho["tsid"]=tsid
    ho.close()
#    print("done")

    #nlu=noise_lookup()
    #temp2=n.zeros([5,n_times])
    #for i in range(5):    
    #    for ti in range(n_times):
    #        tnow=t0+tsid[ti]+1
     #       temp2[i,ti]=nlu.get_noise_lookup(i,tnow)
#        plt.plot(tsid,temp2[i,:])
 #       plt.show()

class noise_lookup():
    def __init__(self):
        h=h5py.File("data/noise_lookup.h5","r")
        self.temp=h["temp"][()]
        self.n_times=self.temp.shape[1]
        self.t_epoch=h["t0"][()]
        self.t_sid=h["tsid"][()]
        h.close()
    def get_noise_lookup(self,beam_id,tnow):
        t_sid_now=(tnow-self.t_epoch)%d_sidereal
        ti=n.argmin(n.abs(t_sid_now - self.t_sid))
        return(self.temp[beam_id,ti])
    
    
            
if __name__ == "__main__":
    calc_beam_lookup()
#    times=n.arange(24*2*2)*3600/2/2 + 1718112000
 #   temps=[]
  #  for t in times:
   #     temps.append(get_noise(0,90,t))
   # plt.plot(times,temps)
    #plt.show()
