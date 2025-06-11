import numpy as n
import matplotlib.pyplot as plt
from pygdsm import GlobalSkyModel
import pansy_config as pc
from datetime import datetime

from pygdsm import GSMObserver
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.time import Time
import astropy.units as u

# Inputs
azimuth = 120 * u.deg        # Azimuth (measured from North towards East)
elevation = 45 * u.deg       # Elevation (Altitude)
latitude = pc.lat           # Example: Tromsø, Norway
longitude = pc.lon
elevation_m = 100              # Optional: Observer's elevation in meters
unix_time = 1718112000       # Example UNIX timestamp (2024-06-11 00:00:00 UTC)

# Define observer location and time
location = EarthLocation(lat=latitude * u.deg, lon=longitude * u.deg, height=elevation_m * u.m)


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
def get_noise(az,el,unix_time,bw=3):
    time = Time(unix_time, format='unix')

    # AltAz frame at given time and location
    altaz_frame = AltAz(obstime=time, location=location)
    # Create SkyCoord in AltAz
    az_del=n.linspace(-bw,bw,num=3)
    el_del=n.linspace(-bw,bw,num=3)
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

if __name__ == "__main__":
    times=n.arange(24*2*2)*3600/2/2 + 1718112000
    temps=[]
    for t in times:
        temps.append(get_noise(0,90,t))
    plt.plot(times,temps)
    plt.show()
