# PANSY Meteor Extension



https://github.com/user-attachments/assets/2e66e59c-f924-4e65-8808-d4a50f88fefd



The <a href="https://pansy.eps.s.u-tokyo.ac.jp/en/">PANSY</a> radar is a powerful MST radar operating in Antarctica at Syowa station. It is primarily used to study the troposphere and the lower thermosphere using radar echoes of turbulence. The radar is extremely capable in terms of sensitivity also for meteors. For the purpose of extending the scientific outcome of PANSY, a receiver capable of detecting meteor head echoes was shipped to Syowa station in late 2024, and became operation in January 2025. This will allow creating a continuous record of meteor head echoes during the remainder of the operational lifetime of the PANSY radar. This is a unique opportunity to observe the southern hemisphere meteor flux via the head echo technique, which is not likely to occur anything soon.

Software and engineering documantation related to the PANSY meteor head echo receiver. This is a software defined radio receiver that operates independently of the PANSY receiver, using the transmit pulse leakthrough as a phase and radar experiment sequence reference. 

<img width="681" alt="Screenshot 2025-02-06 at 22 45 44" src="https://github.com/user-attachments/assets/96143e2d-e476-4a23-aa50-909337f93215" />

<img width="913" alt="Screenshot 2025-02-04 at 22 42 37" src="https://github.com/user-attachments/assets/23685b9a-d38c-402e-ae0a-d1c42ffa0ac4" />

The video below shows the troposphere-stratosphere mode analyzed using the receiver measurements connected to the 1 kW 19 antenna mini PANSY in Shigaraki. The radar sees mainly airplanes and ground clutter.

https://github.com/user-attachments/assets/849bf707-b011-4c54-bbf1-1a141d7f4ec6

## Instructions

Run project Python programs in the base conda environment:

```bash
conda run -n base python PROGRAM.py
```

Runtime analysis overrides belong in
`~/.config/pansy-receiver/pansy-analysis.env`, or in the file selected by
`PANSY_ANALYSIS_CONFIG`.

## Important programs

This repository contains experimental studies alongside the operational
pipeline. The following files are the main entry points and shared modules.

### Receiver and detection chain

- `find_mode_starts.py`: identifies PANSY experiment modes from transmit
  pulses and writes transmit-mode metadata.
- `mesomode_boundary.py`: converts transmit-mode records into mesosphere-mode
  observing intervals.
- `tx_xphase.py`: estimates and monitors transmit-pulse cross-phase alignment.
- `quick_search_meteor.py`: performs the first range-Doppler meteor search over
  mesosphere-mode data.
- `cluster_mf.py`: clusters matched-filter detections into candidate events.
- `cut_meteors.py`: writes complex-voltage event cuts around clustered
  detections.
- `process_cut_meteor.py`: performs the initial per-cut measurement extraction.
- `receiver/systemd/`: user services for the continuously running receiver
  processing chain.
- `receiver/README.md`: deployment, configuration, and service-operation notes.

### Trajectory and orbit analysis

- `run_meteor_analysis_pipeline.py`: standard science-quality post-cut analysis
  entry point; supports MPI event parallelism.
- `plot_interferometric_disambiguation.py`: interferometric alias generation,
  constant-velocity trajectory ranking, event diagnostics, and event plots.
- `pansy_orbit.py`: orbit determination and DASST integration for the winning
  trajectory.
- `orbit_metadata_table.py`: canonical compact HDF5 schemas and conversion
  helpers for event, path, and alias catalogue metadata.
- `pansy_interferometry.py`: shared interferometric response and direction
  calculations.
- `pansy_ballistic.py`: ballistic and atmospheric trajectory model helpers.

### Catalogue and paper figures

- `plot_orbit_catalogue_statistics.py`: catalogue count, height, velocity, and
  uncertainty statistics used by the paper.
- `plot_paper_radiant_results.py`: aggregate, corrected, snapshot, and
  shower-radiant paper products.
- `plot_fitted_radiant_distribution.py`: static and HEALPix radiant
  distributions from orbit metadata.
- `plot_meteor_position_histograms.py`: beam-centered meteor-position
  histograms using the interferometry search grid.
- `plot_height_band_spatial_frequency_decomposition.py`: native four-panel
  height-band and spherical-harmonic radiant figure.
- `healpix_hammer.py`: common true-boundary HEALPix renderer for Hammer
  projections.
- `radiant_spatial_frequency_filter.py`: spherical-harmonic filtering and
  iterative positive-residual extraction.

### Interactive review tools

- `annotate_event_plot_issues.py`: fast keyboard-driven event-plot quality
  annotation.
- `height_band_web_gui.py`: browser interface for selecting height-velocity
  populations.
- `select_height_bands.py`: Matplotlib height-band selection and comparison
  utilities.
- `interactive_radiant_frequency_filter.py`: Mac-first Qt GUI for tuning the
  spherical-harmonic decomposition; filtering runs only when **Recompute** is
  pressed.

### Monitoring and operations

- `status_plot.py`: builds receiver, processing, radiant, and operations status
  products for the PANSY web monitor.
- `receiver/install_user_service.sh`: installs and updates the user systemd
  services.
- `receiver/scripts/`: service launchers, wait helpers, watchdog, and raw-data
  retention scripts.

Scientific and intermediate data products are HDF5. Do not use CSV or NPZ for
project data.

## Meteor cut analysis

The standard post-cut meteor analysis pipeline is the interferometric
disambiguation pipeline in `run_meteor_analysis_pipeline.py`. It recomputes
range-Doppler observables for each cut, constructs all high-coherence
interferometer alias candidates, filters path hypotheses by visibility,
linearity, and descent, then ranks surviving solutions with the ballistic fit
and transmit-beam consistency diagnostics.

Use this pipeline for science-quality meteor cut processing and for generating
the accompanying candidate, ranking, orbit, HDF5 metadata, and LaTeX summary
products.
