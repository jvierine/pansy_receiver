import numpy as np
import healpy as hp
import matplotlib.pyplot as plt 
import h5py 

# Generate some random latitude and longitude points (or use real data)
#num_points = 1000
h=h5py.File("slon.h5","r")
lons=h["slon"][()]
lats=h["slat"][()]
h.close()
print(len(lons))
#lats = np.random.uniform(-90, 90, num_points)  # Latitude in degrees
#lons = np.random.uniform(-180, 180, num_points)  # Longitude in degrees

# Convert latitude to colatitude (theta) in radians
theta = np.radians(90 - lats)  # Colatitude: 90° - latitude
phi = np.radians(lons+90)  # Longitude in radians

# Set HEALPix resolution (higher Nside means higher resolution)
nside = 32  # Must be a power of 2

# Convert lat/lon to HEALPix pixel indices
pixels = hp.ang2pix(nside, theta, phi)

# Create a HEALPix map and count occurrences in each pixel
histogram = np.bincount(pixels, minlength=hp.nside2npix(nside))

# Plot the histogram as a HEALPix map
hp.mollview(histogram, title="Histogram of radiants", unit="Counts",cmap="turbo")
hp.graticule(color="white",alpha=0.1)
plt.show()