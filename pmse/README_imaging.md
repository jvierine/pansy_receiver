# Five-beam altitude-plane HSV imaging

The new combined product is distinct from the existing first-image plots.
It uses the selected zenith gate's height, 85.290954301 km **above the radar**,
on a common flat horizontal plane. Off zenith it interpolates reconstructed
spectra at longer slant ranges, not at a fixed range gate.

## Compute on Revontuli

Run from `/home/j/src/pansy_receiver` (reference contains `cfg/` and `data/`):

```sh
OPENBLAS_NUM_THREADS=8 /usr/bin/python3 pmse/animate_five_beams.py \
  --root /mnt/data/juha/pansy/metadata/xc2 --reference . \
  --start 2025-01-30T04:00:00 --end 2025-01-30T04:30:00 \
  --height-km 85.290954301 --output /tmp/five_beam_hsv
```

For local Python use `conda run -n base python`. Dependencies: NumPy, h5py,
Matplotlib, Pillow. No CSV products are written; the geometry CSV is an
existing instrument input. All generated scientific data are HDF5.

## Recolor without source XC2 or SVD

Download the [imaging sidecar](https://juha.no/share/pansy_pmse_2025-01-30_85km_five_beam_hsv.h5)
as `five_beam_hsv.h5`, then run:

```sh
conda run -n base python pmse/render_five_beam_hsv.py five_beam_hsv.h5 \
  --output figures/recolored --doppler-max 7 --width-max 10 \
  --saturation-min 0.2 --db-min -25 --db-max 0
conda run -n base python pmse/check_imaging_sidecar.py five_beam_hsv.h5
```

PNG uses 166 mm width, 11 pt labels and 300 dpi. GIF uses the same layout at
screen resolution. `--frame-ms` changes playback speed, not the sample times.
Each render writes a `.style.h5` recording limits and the data-sidecar SHA256.
No numerical data in the source sidecar are modified during recoloring.

## Sidecar schema (v1)

- `spectra/beam0` through `spectra/beam4` and `spectra/combined` contain `power`
  in `(time, covered_pixel, Doppler)` order and `flat_pixel_index`, the C-order
  indices into `(north, east)`. Values are float32 with lossless gzip/shuffle.
- `velocity_mps`, `doppler_frequency_hz`, `time_unix_us` provide spectral/time axes.
- `power_density_proxy`, `relative_power_db`, `centroid_mps`, `width_mps` are
  `(time, north, east)` float64 arrays. `beam_moments` stores corresponding
  unblended values in `(time, beam, moment, north, east)` order.
- `east_km`, `north_km`, `height_km`, `slant_range_km`, `direction_uvw`,
  `coverage_mask`, `beam_weights` specify the common geometry and blend.
- Calibration, antenna positions, channel pairs, singular values, source paths
  and timestamp groups, and algorithm settings are stored alongside the data.

No clipped RGB values are used as the scientific source. The saved full
spectra permit recomputing moments using different Doppler windows without
an inversion. Changing altitude, imaging mask, or calibration requires a new
inversion; this sidecar does not store a full three-dimensional image volume.

## Scientific limitations

Uses positive meteor calibration/steering signs, with a distinct phasecal for
each TX beam. Calibration is not independently validated for this historical
observation. The 18-mode inverse and 1/(s+1) denominator match the earlier
regularization, now on beam-centered domains on a common Cartesian grid.

Squared inversion magnitudes are normalized over each beam's image area and
weighted by nonnegative excess spectral power/noise. This is a relative
spectral density proxy, **not absolute SNR or calibrated echo brightness**.
Cosine-squared overlap tapers are combination weights, not measured TX gains.
The noise-floor/positive-power clipping can bias weak-signal moments.
The full archived Doppler swath is included; color clipping does not change
the moment integration window. Black outside coverage is not interpolated.
An animation frame represents one recorded integration, not an instantaneous
simultaneous transmission: the transmitter cycles through the five beams.
