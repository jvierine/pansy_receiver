# XC retirement audit

XC retirement is a producer change, not deletion of historical XC/xc2 data.
Do not change `xc_metadata_dir` to point to sparse data: its schema, normalization,
frequency resolution, and selection are different. The old readers must remain
available for historical reproducibility.

## Consumer inventory (repository audit, 2026-10-19)

| Consumer | Dependency and disposition |
| --- | --- |
| `status_plot.py` | Active monitor reads XC and reports FXC lag. Must be migrated before stopping XC. New PMSE panels render independently at `/rvd/`; meteor/raw/status monitoring must continue. |
| `paper_plot_noise.py:read_xc2_day` | Legacy dense XC-based noise/peak-range figure. Retain for historical inputs. Never give it a sparse channel or reuse its XC power scale. |
| `save_cut_noise_measurements.py`, `fit_cut_noise_pygdsm.py` | Adopted noise-fit workflow uses guarded raw meteor cutouts and HDF5 measurements, not XC. Unaffected. New `pmse.export_noise` produces compatible per-module measurement datasets from dense sparse-archive noise fields. |
| `meso_cal.py` | Offline zenith calibration reads XC; leave historical implementation intact. `pmse.sparse_readers.cross_products` supplies all 28 equivalent pair products at retained cells with an explicit calibration sign. A sparse calibration selection/fit still needs real-data validation. |
| `plot_xc.py`, `image_xc.py` | Historical/offline XC examples. Preserve legacy support; use explicit sparse adapter for new work rather than manufacturing dense arrays. |
| `pmse/image_first_block.py`, `image_five_beams.py`, `animate_five_beams.py` | Historical January 2025 imaging reproduction. Their original XC2 input remains available and is not replaced. |
| `pmse/plot_pansy_pmse_{day,spectra,falsecolor}.py` | Historical XC2 figure generation. Preserve original input and results. New viewer consumes sparse data explicitly. |
| `meso_xc_rvd.py`, `plot_xc_rvd.py` | Separate legacy `xc_rvd` path, not the active `xc2` producer. Must not be confused with the new sparse archive or stopped based on its name. |
| `cal_zenith.py` | Local variable called `xc2` is meteor correlation, not a reader of the XC2 product. Not a dependency. |
| Receiver, TX-mode finding, mesomode boundaries, meteor MF/search/cut/fit | Upstream/parallel chains do not require XC; do not restart them for this migration. |
| Metadata backup/mirror | Keep old XC backup/history. Add new sparse channel separately; no `--delete`, permission/owner/group preservation disabled. |

## Information and normalization limits

Dense `noise_power[7,5,range]` and its uncertainty are stored for **every** full
mode, even when no sparse signal cells are retained. Noise-model selection is
therefore independent of detected echoes. Unit-energy measured-TX decoding and
unit-L2 Hann weighting have unit white-noise gain. Export these noise estimates
with `python -m pmse.export_noise --metadata ROOT --day YYYY-MM-DD --output noise.h5`,
then pass that HDF5 file to the existing fitter with `--measurements noise.h5`.
It has mode-midpoint sampling, not cutout-pulse sampling. Validate colored-noise
filter response and compare with contemporaneous cut-padding measurements before
reusing an absolute Kelvin calibration. Historical XC uses TX-amplitude-dependent
and reciprocal-weighted powers; its numerical scale is not interchangeable.

At a retained cell, the seven spectra from every coherent look reconstruct
`C_cd = mean_l(S_l,c * conj(S_l,d))` exactly. This preserves interferometric phase
and permits later calibration. Current five-beam imaging uses the positive
calibration factor `exp(+i*(phi_c-phi_d))`; callers must request the sign explicitly.
An independent known-direction check remains necessary. The archive stores
uncalibrated spectra, geometry and phase-calibration snapshots, not overwritten
calibrated voltages.

Sparse selection is intentionally lossy: below-threshold spectra, full-band
spatial noise covariance and the old finer coherent frequency resolution cannot
be recovered from it. Dense per-receiver noise powers support the adopted
system-noise model, but are not a substitute for measured 7x7 sky covariance.
Do not describe absent bins as zero or threshold-censored moments as unbiased.

## Required cutover gates

1. Rust/Python same-raw-block parity: same retained cells, axes, actual support,
   noise and phase, with documented floating-point tolerances.
2. New live records advance, include empty/noise-only modes, respect resource and
   disk budgets, and survive interrupted writes without duplicate publication.
3. New native files replicate to Revontuli and are independently readable there.
4. `/rvd/` is verified with real data and actual browser rendering. The established
   monitor no longer reads XC for its live panel or freshness indicator.
5. Existing noise-model measurements/calibrations and historical imaging inputs
   remain unchanged. New noise estimates pass contemporaneous raw-noise checks.
6. Only then stop/disable **`pansy-meso-xc.service`**. Do not delete its data, stop
   mesomode/TX services, or alter receiver acquisition. Keep the old service unit
   as a rollback path. Check unaffected services and fresh output after cutover.

Until these gates pass, XC remains enabled. An audit or successful synthetic
test alone is not permission to claim a completed operational replacement.
