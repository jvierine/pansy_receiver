"""Validate geometry, Doppler units, spectral moments and blending in a sidecar."""
import argparse
import h5py
import numpy as np
from render_five_beam_hsv import colors


def moments(spectrum,velocity):
    power=spectrum.sum(axis=-1)
    mean=np.divide(spectrum@velocity,power,out=np.zeros_like(power),where=power>0)
    variance=np.divide(spectrum@(velocity**2),power,out=np.zeros_like(power),where=power>0)-mean**2
    return np.array([power,mean,np.sqrt(np.maximum(variance,0))])


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('source')
    a=p.parse_args()
    with h5py.File(a.source) as h:
        assert h.attrs['geometry']=='constant_altitude_plane'
        east,north=np.meshgrid(h['east_km'][()],h['north_km'][()])
        r=h['slant_range_km'][()]; height=h['height_km'][()]
        np.testing.assert_allclose(r**2,east**2+north**2+height**2)
        np.testing.assert_allclose(-h['direction_uvw'][()][:,:,2]*r,height)
        velocity=h['velocity_mps'][()]
        np.testing.assert_allclose(h['doppler_frequency_hz'][()],
            2*h.attrs['radar_frequency_hz']*velocity/h.attrs['speed_of_light_mps'])
        weights=h['beam_weights'][()].reshape(5,-1)
        cover=h['coverage_mask'][()].ravel()
        combined_idx=h['spectra/combined/flat_pixel_index'][()]
        np.testing.assert_array_equal(combined_idx,np.flatnonzero(cover))
        for frame in range(len(h['time_unix_us'])):
            blend=np.zeros((cover.size,len(velocity)))
            for beam in range(5):
                bg=h[f'spectra/beam{beam}']
                idx=bg['flat_pixel_index'][()]
                spectra=bg['power'][frame].astype(float)
                assert np.isfinite(spectra).all() and (spectra>=0).all()
                recomputed=moments(spectra,velocity)
                saved=h['beam_moments'][frame,beam].reshape(3,-1)[:,idx]
                np.testing.assert_allclose(recomputed,saved,rtol=2e-5,atol=1e-5)
                blend[idx]+=spectra*weights[beam,idx,None]
            blend=blend[cover]/weights.sum(axis=0)[cover,None]
            saved_spectrum=h['spectra/combined/power'][frame].astype(float)
            np.testing.assert_allclose(blend,saved_spectrum,rtol=2e-5,atol=1e-7)
            recomputed=moments(saved_spectrum,velocity)
            saved=np.array([h[key][frame].ravel()[cover] for key in
                           ['power_density_proxy','centroid_mps','width_mps']])
            np.testing.assert_allclose(recomputed,saved,rtol=2e-5,atol=1e-5)
        rgb_args=[h[k][()] for k in ['centroid_mps','width_mps','relative_power_db','coverage_mask']]
        settings=dict(doppler_max=9.568,width_min=0,width_max=6.379,saturation_min=.2,db_min=-30,db_max=0)
        original=colors(*rgb_args,**settings)
        adjusted=colors(*rgb_args,**dict(settings,doppler_max=5,width_max=15))
        assert np.isfinite(original).all() and not np.allclose(original,adjusted)
        print('PASS: all frames, five beam spectra, blend, moments, altitude, Doppler units, recoloring')
        print('Centroid percentiles:',np.percentile(h['centroid_mps'][()][:,cover.reshape(east.shape)],[5,50,95]))
        print('Width percentiles:',np.percentile(h['width_mps'][()][:,cover.reshape(east.shape)],[5,50,95]))


if __name__=='__main__':
    main()
