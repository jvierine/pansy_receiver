# PANSY PMSE / PMWE sparse spectral archive

Status: implementation and numerical tests complete; real-data commissioning
is separate from the existing cross-spectral service. The new archive must
pass measured noise, signal-injection, coverage, and throughput checks before
it replaces any existing product. No old data are removed.

## Scientific target

Preserve the information needed to search for extreme mesospheric vertical
drafts, including angular checks against sidelobes. Chau et al. (2021),
[doi:10.1029/2021GL094918](https://doi.org/10.1029/2021GL094918), report PMSE
vertical drafts near +/-50 m/s. Chau, Berera and Huyghebaert (2025),
[doi:10.3389/fspas.2025.1716224](https://doi.org/10.3389/fspas.2025.1716224),
discuss events above 12.5 m/s and five standard deviations of the observed
vertical-velocity distribution. The 12.5 m/s value is a screening flag here,
not a declaration that a spectral cell is a rogue vertical draft.

High radial velocities require arrival-angle and interference checks. An
off-zenith echo can project horizontal motion onto the line of sight. The
archive keeps all seven complex receiver channels for later imaging and
calibration. A final draft catalogue should use persistent spectral features,
local background variability, and spatial consistency, rather than isolated
high-speed bins. The physical upward/downward sign still needs an independent
known-target check; the stored algebraic convention is f_D=2*f_radar*v/c.

## Observing and processing choices

- One record per **closed measured mesomode**, start to stop. Current examples
  on 2026-10-19 last 85.9--86.5 seconds. The format also accommodates a 120 s
  mode; fixed two-minute windows must not cross observing-mode changes.
- Five beams: zenith, north, east, south, west; tilted beams have 80 deg
  elevation. All seven receiver channels share the same sparse cell mask.
- Select nominal beam-center altitudes 50--120 km above the radar. Tilted
  beams need slant ranges beyond 120 km. Store slant range and nominal height
  separately; a sidelobe return's true height depends on its recovered angle.
- Sampling is 1 MHz, IPP 1600 samples, and the same beam repeats every 8 ms.
  The unambiguous velocity interval is approximately +/-199.33 m/s at 47 MHz.
- Default coherent FFT length 256 pulses per beam: 2.048 s and 1.557 m/s
  Doppler spacing. Retain the nearest bins spanning +/-90 m/s (actual outer
  bin centers +/-90.321 m/s). Do not reduce the Nyquist limit by pulse
  presumming, notch the zero-Doppler region, or discard extreme Dopplers.
- Decode with the recorded transmit pulse for each IPP, not an ideal code.
  Normalize the measured waveform to unit pulse energy. This is equivalent
  to unit mean pulse power followed by correlator division by sqrt(132).
  White receiver-noise gain is then one, independent of reference amplitude.
  Store measured transmit power summaries and the normalization convention.
- Use linear cross-correlation, a Hann coherent window with unit L2 norm,
  and an unnormalized forward FFT. Preserve the complex channel phases.
  Range spacing defaults to four samples, approximately 600 m, matching the
  old quicklook spacing. This spacing is distinct from range resolution.
- Incoherently average **power**, preserving the complex spectrum from every
  coherent look. Never average complex spectra from independent looks into
  a single voltage value. Keep actual look bounds, pulse counts, valid
  coverage and missing-look reasons; never synthesize samples across gaps.

## Detection and uncertainty

Estimate decoded spectral noise separately for each receiver, beam and range.
Use the median of the incoherently averaged spectrum in guard bins outside
the archive's +/-90 m/s band, correcting the median for the Gamma distribution
of K averaged exponential powers. The full FFT is available for this noise
estimate, although complex voltages are archived only in selected cells.
RFI or atmospheric signal in guard bins can bias the estimate, generally
upwards; retain noise diagnostics and validate this on measured data.

The detector averages the seven noise-whitened receiver powers. Define
SNR = (total_power - noise_power)/noise_power. Thus 3 dB means the **excess
signal** is approximately twice the noise, or total power approximately
three times noise. It does not mean total power/noise = 3 dB.

Select cells with SNR >=3 dB and a conservative thermal-noise significance
gate. The latter uses K, rather than 7*K, as the Gamma shape, allowing fully
correlated receiver noise. A finite-noise-estimate F-tail approximation and
Bonferroni factor target 1e-3 false detections per complete record under the
model. This is an approximate thermal-noise model, **not an empirical bound
on false atmospheric events**: interference, clutter, meteors, and decoder
sidelobes require separate checks. At least 16 coherent looks are required.

For K=42, a conservative noise-only relative power standard deviation is
1/sqrt(K)=15.4% (84 degrees of freedom). Record each selected cell's measured
power variance and standard error across looks, plus the noise power and
its estimated standard error. These are linear-power errors; SNR is stored
and displayed in dB. Missing or insufficient data produce an explicit quality
record; they must not be interpreted as a confirmed absence of PMSE/PMWE.

Selection is performed once on the full ensemble. Retain all seven channels
and every valid coherent look at each selected cell, including weak channel
values. Thresholding each channel or look independently would bias the
cross-correlations and prevent faithful imaging.

## Digital Metadata v1

Standard Digital Metadata channel, sample rate 1,000,000/1 Hz, hourly
directories and 600 s files, basename `pmse`. Record key is the actual mode
start sample. Files remain readable by `digital_rf.DigitalMetadataReader`.
Scientific arrays use lossless HDF5 gzip level 1 and shuffle. A writer lock
serializes output. `complete=1` is committed only after the record and its
compressed arrays are flushed; readers reject incomplete records.
Re-running a completed mode is idempotent. Only the producer's interrupted,
explicitly incomplete groups may be repaired in place.

| Field | Shape / purpose |
| --- | --- |
| `sparse_index` | N x 3: beam, full-FFT bin, range-axis index |
| `spectra` | K x 7 x N complex64 measured-TX decoded coherent spectra |
| `mean_power` | 7 x N float32; average squared magnitude across looks |
| `power_variance`, `power_standard_error` | 7 x N float32; inter-look statistics |
| `noise_power`, `noise_power_standard_error` | 7 x 5 x R; retained even in empty records |
| `detector_snr_db` | N; excess signal/noise, with explicit definition |
| `detector_standard_error`, `null_probability` | N; detector diagnostics |
| `look_bounds_sample` | K x 2; start and exclusive stop |
| `pulses_per_beam_per_look` | K; distinguishes full and partial looks |
| `mode_start_sample`, `mode_stop_sample_exclusive` | actual observing boundaries |
| `mode_duration_seconds`, `coverage_seconds` | full interval and actually decoded support |
| `slant_range_km`, `nominal_altitude_km`, `altitude_valid` | range/height geometry |
| `doppler_frequency_hz`, `radial_velocity_mps` | complete FFT coordinate axes |
| `antenna_positions_m`, calibration snapshot and hashes | later interferometry inputs |
| `tx_mean_min_max_power` | K x 3; measured pulse-reference power diagnostics |
| `high_velocity_candidate` | N; |v_radial|>=12.5 m/s screening mask |
| configuration, normalization, quality and source metadata | reproducibility, limits, coverage |

Power units are decoded ADC-count squared per spectral bin under the stated
normalization, not an absolute radar reflectivity or temperature calibration.
The scalar noise estimate does not supply a measured spatial sky-noise
covariance. Channel covariance at selected cells can be recovered exactly as
mean_l(S_l,c * conj(S_l,d)); neither 21 separate cross-products nor the seven
autocorrelations need to be stored redundantly.

The absence of a sparse cell means it failed the archive's selection. It does
not mean zero power. Spectral moments formed only from selected cells are
censored estimates and must be labeled/treated accordingly. Recoloring and
interferometric reanalysis of retained cells are possible. Lowering the SNR
threshold or changing altitude/resolution after raw deletion cannot recover
discarded cells.

## Run and inspect

Use base conda locally. At Syowa use the existing processing virtual
environment `/home/radar/testenv/bin/python`; no base conda installation was
found there. Deploy source only through git into an isolated checkout.

```sh
conda run -n base python -m unittest discover -s tests -p test_pmse_sparse_record.py
conda run -n base python -m pmse.sparse_record --output /path/to/new/channel --latest 1
conda run -n base python -m pmse.inspect_sparse_record /path/to/new/channel --output detected_cells.png
```

For commissioning on Syowa, use one worker and low scheduling priority:

```sh
nice -n 15 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/radar/testenv/bin/python -m pmse.sparse_record \
  --output /media/analysis/metadata/pmse_sparse_trial --latest 1
```

The loop option processes newly closed modes and polls every 60 s. It is not
enabled automatically. In its initial form the discovery lookback defaults
to three hours; choose an explicit lookback/batch for a retention backfill.
Commission a quota, oldest-first catch-up policy, and a full recovery test
before using this as the only long-term archive.

## Resources and retained history

2026-10-19 live snapshot: Syowa has 16 logical CPUs, load about 9, 236 GiB
available memory, 1.7 TiB free on the analysis filesystem and 8.7 TiB on the
raw archive filesystem. Revontuli has 96 logical CPUs, load about 43, 279 GiB
available memory and 7.9 TiB free on `/mnt/data`. Raw voltages are at Syowa,
where the configured cleanup interval retains roughly ten days. Process
there first to avoid shipping full raw voltages over the limited link.

Uncompressed selected voltages cost 8*K*7*N bytes per mode. At K=42 this is
2352 bytes per selected cell, plus power/noise/index overhead. A 1% occupancy
of roughly 70,000 tested cells is about 1.6 MB per mode before compression;
actual occupancy and full-day rates must be measured in both summer and
winter before committing disk capacity. Dense geometry/noise metadata remains
even when N=0, making coverage and non-detections auditable.

Historical cross-spectra remain useful, but cannot recreate the seven complex
look spectra, unrecorded altitude gates, or unrecorded Doppler bins once the
original raw voltages have expired. This archive improves retention from its
start and can backfill only the raw data that still exist.
