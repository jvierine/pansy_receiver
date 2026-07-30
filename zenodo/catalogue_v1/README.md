# PANSY meteor head-echo catalogue, Version 1

Join all products with `event_id`, a Unix sample index at 1 MHz. Float missing
values are `NaN`. Dataset attributes give units, frames, parameter order,
schema version, release version, and processing commit.

## Files

| Path | Contents |
|---|---|
| `pansy_level1_example_1757458600402300.h5` | Figure 3 raw-voltage cut |
| `level2/pansy_level2_YYYY-MM.h5` | Time-resolved measurements |
| `level3/pansy_level3_YYYY-MM.h5` | Radiants and orbits |
| `release_summary.json` | Coverage and row counts |
| `SHA256SUMS` | File checksums |
| `example_level2_radiant.py` | Straight-line Level 2 radiant example |
| `verify_level2_level3.py` | Independent Level 2 orbit check |

Level 1 voltage is stored as separate real and imaginary signed ADC counts.
It is an event cut, not continuous voltage, and is not calibrated to volts.

## Level 2

`/events/measurement_start` and `/events/measurement_count` index contiguous
rows in `/measurements`.

| Dataset | Unit |
|---|---|
| `event_id`, `time_sample_idx` | sample at 1 MHz |
| `time_offset_s` | s |
| `east_km`, `north_km`, `up_km`, `range_km` | km |
| `snr_db` | dB |
| `doppler_velocity_mps` | m s-1 |
| `beam_id` | index |
| `selection_keep` | boolean |

Position is local east-north-up at PANSY.

## Level 2 example

Requires Python, NumPy, h5py, and Astropy:

```sh
python example_level2_radiant.py level2/pansy_level2_2025-07.h5
```

The script reads one event, fits a constant-velocity straight line to retained
ENU positions, and prints its apparent GCRS radiant.

## Level 3

The file contains `v0`, `vg`, GCRS and ecliptic radiants, Sun-centered
ecliptic radiants, solar longitude, two GCRS Cartesian states, fit-quality
fields, and:

| Quantity | Order |
|---|---|
| Kepler elements, standard deviations, covariance | `a,e,i,Omega,omega,nu,q` |
| Local fit parameters and covariance | `east,north,up,ve,vn,vu,log10_beta` |
| Cartesian states | `x,y,z,vx,vy,vz` |

`initial_state_gcrs_m_mps` is the fitted state. 
`geocentric_state_gcrs_m_mps` replaces its velocity with the geocentric
radiant velocity after zenith-attraction correction.

The local ENU fit covariance is not a covariance of either GCRS state. A
zenith-attraction-corrected Cartesian-state covariance is unavailable in the
compact source catalogue. The Kepler covariance comes from the orbit
uncertainty ensemble.

## Independent check

Requires Python, NumPy, h5py, and Astropy:

```sh
python verify_level2_level3.py . --month 2025-07
```

The script fits Level 2 positions, transforms the state to GCRS, and calculates
heliocentric Keplerian elements without zenith-attraction correction.

## Deposit metadata

Add the article citation, DOI, creators, ORCIDs, affiliations, funding,
contact, and data license in Zenodo.
