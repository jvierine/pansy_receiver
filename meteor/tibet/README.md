# Tibet meteor and interferometry chain

Offline Rust DSP for the November 2026 15-channel campaign. It reads the external
raw files without modifying them. No acquisition or PANSY service is changed.
Python is used only for Digital Metadata packaging, independent verification,
and plotting. Range-Doppler filtering reuses the existing sibling PANSY Rust
crate through an in-repository Cargo dependency. This is the canonical Tibet
implementation, maintained together with PANSY on `main`.

## Commands

From the repository root:

```sh
cargo build --release --locked --manifest-path meteor/tibet/Cargo.toml
meteor/tibet/target/release/tibet-meteor inspect RAW_FOLDER
meteor/tibet/target/release/tibet-meteor process RAW_FOLDER NEW_RUN --write-legacy-cuts 0
meteor/tibet/target/release/tibet-meteor cut NEW_RUN NEW_CUTS
meteor/tibet/target/release/tibet-meteor verify-cuts NEW_CUTS
```

Use a fresh output directory. Optional processing limits are `--max-files N`,
`--max-pulses N`, and `--threshold-db D`. The default statistic threshold is
12 dB excess over median matched-filter power. This statistic is not calibrated
SNR or false-alarm probability. At most three peaks per IPP are retained, so
the chain does not claim complete recovery of arbitrarily crowded/weak meteors.

The final ten-channel cuts use interleaved int16 with shuffle + deflate level 4.
Sparse low/high int16 overflow words reconstruct the original int32 exactly.
The cutter verifies storage round trips before replaying the matched filters.
`repack-cuts OLD_CUTS NEW_CUTS` upgrades existing version-1 ten-channel cuts
without needing the raw source drive. The shared Python voltage reader accepts
both schemas. The standalone MATLAB implementation is in `matlab/tibet`.

Publish standard Digital RF metadata without expanding the compact voltages:

```sh
conda run -n base python -m meteor_analysis.cut_metadata NEW_CUTS --run NEW_RUN --geometry ARRAY_GEOMETRY.h5
```

The metadata retains interleaved int16 and sparse overflow tables with shuffle
and deflate. Use `meteor_analysis.voltage.read_iq` to reconstruct voltages from
either cut files or Digital Metadata records. Legacy split int32 metadata is
still supported when publishing older cuts.

## Data contract and metadata evidence

- Header: 8192 bytes, little endian, guarded by `AAAA5555` and `AA5555AA`.
- Pulse: 8-byte marker/counter, 640 gates x 15 channels x interleaved int32 I/Q,
  then matching 8-byte marker/counter. All declared records are validated.
- Sampling 500 kHz; IPP 1280 us; 6144 records per file. Only those records are
  observations. Remaining file capacity is not decoded as additional pulses.
- Header time is 1970 with invalid subsecond words. Filename timestamps are
  interpreted as Beijing time, explicitly marked unverified. The filename has
  only one-second precision; integer pulse indices have exact relative timing
  only if the header IPP describes the exported records.
- `INACC=4` is retained, not interpreted as another decimation factor. Its
  firmware accumulation semantics require confirmation. The 32-sample measured
  leakage is used rather than assuming `width=8 us` is the whole coded waveform.
- TX alpha/beta and all RX phase and weight-index words are retained per source.
  A changed signature starts a new pointing segment; no track crosses a file
  boundary or acquisition gap. Joint calibration is withheld on signature changes.
  An undocumented change within a single-header file cannot be recovered from
  nonexistent per-IPP steering metadata.

Zero-based channel interpretation from the supplied 15-panel image and weight
indices:

| Channels | Interpretation | Amplitude mask indices |
|---|---|---|
| 0-4 | Full array S, N, zenith, E, W | 15 for all five |
| 5-7 | S0, S1, S2 (19-element groups) | 23, 24, 25 |
| 8-14 | L0-L6 (large groups) | 26-32 |

The signed phase increments are consistent with the first five receivers being
15 degrees S/N/Z/E/W using a 2.1 m x / 3.637 m y phase grid. This is an inference
from the phase codes, lattice and image labels, not verified firmware semantics.
The ten subarray RX phase codes and TX codes are zero in all five input files.
No transmit pointing changes were found in these files. The mask indices are
identifiers, not antenna lists. Geometry comes from the earlier XLSX-derived
selection and the supplied screenshot; an instrument mask-table export is still
needed to confirm it independently. All 931 positions are retained in the map,
with the geometrical 799-element active mask kept explicitly provisional.

## Signal processing

1. Extract the common complex leakage from channel 2, samples `[7,39)`, and
   normalize its energy. Use the same reference for every channel. Never
   self-match each receiver against its own leakage: that removes receiver phase.
2. Estimate noise from the median raw powers at gates 250-639. Search five
   full-array outputs with noise-normalized, native-sample measured-code FFT
   filters, FFT length 128. Preserve original receiver amplitudes by evaluating
   all 15 channels at each common lag/frequency. Sub-bin Doppler interpolation
   is applied before the final complex-amplitude evaluation.
3. Search nominal slant ranges about 70-180 km. Lags are relative to the measured
   reference start; hardware group delay and refractive path are uncalibrated.
4. Reject secondary peaks only when their full 15-channel complex response is
   explained by the stronger peak convolved with the measured code ambiguity
   function. Retain vetoed observations separately. Residual-power fraction
   0.35 is provisional; close same-direction simultaneous meteors can remain
   unresolved. Filtering does not prove head/trail/satellite classification.
5. Group locally by range, Doppler and receiver signature. Require at least
   eight pulses, 30 ms span, gaps <=50 ms, and range-rate/Doppler agreement
   within 10 km/s. These gates are explicit prototype settings, not validated
   detection completeness thresholds.
6. Cut every channel with 80 ms temporal padding clipped at each source-file
   boundary, and 16 range-gate margins around the detected echo window. Store
   original int32 I/Q without 16-bit requantization, pulse indices, measured TX,
   raw headers, range offsets, noise and matched complex observables.

## Calibration and imaging

Initial phase calibration uses one strongest point per independent selected
track, not every pulse as an independent meteor. Require a zenith-dominated
full-array response, >18 dB detection statistic and >3 dB zenith/other look ratio.
The circular phase estimate is relative to channel 8 (L0), with correction
`exp(-i phasecal)`. It covers channels 5-14; the differently pointed full-array
channels are not silently treated as interchangeable interferometer elements.

The estimate assumes the ensemble is centred on boresight. Absolute pointing
and a linear receiver phase ramp remain degenerate with the ensemble mean
direction. Between-event phase dispersion contains genuine angular spread as
well as calibration error. The result is explicitly provisional and never
overwrites or automatically corrects the raw cuts.

Native imaging sums the elements in each subarray and uses their phase response,
not only centroid baselines. It uses the provisional active mask, retains wide
20-degree-box angular score maps, and reports conditional fits in a +/-4-degree
direction-component box around zenith. A phase-only capped-weight match score
is not a noise-covariance likelihood or formal angle uncertainty. Receive-phase
sign/axis conventions and the boresight calibration gauge remain unconfirmed.
These fits do not establish unique all-sky directions or final meteor orbits.

## Outputs and verification

`catalogue.h5`, `cuts/event_*.h5`, `phasecal_boresight_provisional.h5`, and
`interferometry_provisional.h5` retain scientific results. Digital Metadata
channels are `mf`, `detections`, `cut`, `geometry`, `phasecal_provisional`, and
`angles_provisional`, readable with the actual `digital_rf` library. Their sample
rate is 500 kHz. `tx_idx` denotes IPP start; `tx_start_sample` adds the seven-gate
reference offset in cut metadata. Unique metadata keys may have a small collision
offset, recorded separately; observation timestamps are never shifted.

```sh
cargo test --manifest-path rust/meteor/Cargo.toml
conda run -n base python scripts/verify_tibet_native.py
```

The independent injected-voltage check exercises a moving target, known receiver
phases, code-sidelobe suppression, numerical agreement with direct complex NumPy
matched sums, and exact int32 cut reconstruction. Unit tests exercise signed I/Q,
header bounds, integer timing, counters and track-gap handling. Real-data outputs
must also be inspected; these checks do not prove population completeness.

## Ten-channel raw-voltage cuts

```sh
rust/meteor/target/release/tibet-meteor cut outputs/new_run outputs/new_cuts10
conda run -n base python scripts/export_tibet_voltage_cuts.py outputs/new_cuts10 --run outputs/new_run
conda run -n base python scripts/verify_tibet_voltage_cuts.py outputs/new_cuts10 --baseline outputs/new_run --reprocessed outputs/independent_rerun
conda run -n base python scripts/verify_tibet_pointing_cuts.py outputs/new_run
```

The separate `cut` command does not change the detector or the original catalogue.
It stores original int32 echo I/Q for source channels 5-14 only, in S0-S2, L0-L6
order. Each event gets 80 ms padding on both ends (clipped at source boundaries)
and 64 sample range margins. A single, separately sampled channel-2 TX reference
preserves cross-channel receiver phases. The current and preceding two IPPs'
references, decoded baud signs, validity flags, originating pulse beam IDs and
directions are retained; missing previous-file context is not invented.

Header, marker counters, exact pulse indices, sample origin, frequency, IPP,
waveform ID, width field, accumulation field, steering codes and weight-table IDs
are preserved. No raw phase correction is applied. Ten-channel phase metadata is
explicitly provisional, with L0 at cut index 3 (source index 8). Geometry retains
all 931 antenna positions and the provisional active-mask flag.

New-cut Digital Metadata has `cut`, `detections`, `geometry`, and
`phasecal_provisional` channels at 500 kHz. Unlike the legacy 15-channel export,
`tx_idx` means the sampled reference start (IPP start + 7 samples), `ipp_idx`
means record start, and `delays` is echo-window start relative to that reference.
The schema is `tibet_meteor_cut_v1`. Sampled leakage is not calibrated TX power.

For an externally documented within-file sequence, use
`cut RUN OUTPUT --pointing-schedule schedule.h5`. Each source group (`0000`,
`0001`, ...) must contain full-length one-dimensional arrays `tx_alpha_code`
(uint16), `tx_beta_code` (uint16), `tx_beam_id` (int32), `tx_azimuth_deg` (float64,
degrees clockwise from geographic north), and `tx_zenith_angle_deg` (float64).
Azimuth may be NaN for zenith. Schedule pulse zero is the first source record.
Nonzero physical pointings are not guessed from an unconfirmed phase-code sign.
The documented multi-frame cycle uses each frame header's current beam ID and
TX steering codes automatically; it does not imply within-file switching.
Only an actual within-file sequence needs the external schedule. Header-only
pointings are explicitly marked as lacking within-file telemetry; RX metadata remains
per-file and an undocumented within-file RX change cannot be recovered.

Every saved detection is replayed from the cut at its original delay/Doppler.
Diagnostics compare current/previous/two-previous TX references both at that
delay/Doppler and across the entire padded cut with 128-bin range-Doppler
filters. A wide-search maximum may be a different target or code sidelobe;
these diagnostics do not silently select a physical range or alter detections.
The alternating A/B code can test odd-IPP aliases, but two IPPs repeat the code:
the 191.87 km range interval does not prove unaliased range. Other pulse coding,
staggered IPP or independent geometry is needed for an absolute alias order.

Exact detector regression uses all original channels, as before. The ten-channel
cuts reproduce subarray complex amplitudes, not the original five-full-array
detection statistic: that statistic cannot be independently recomputed from
the ten cut channels without the omitted full-array voltages or known combiner
weights. Metadata retains its original value rather than claiming otherwise.
# MST Spectra

`tibet-meteor mst RAW_DIRECTORY DETECTION_RUN NEW_OUTPUT` reduces all declared
pulses into complete complementary A/B pairs for the five full-array receive
looks. It uses the sampled transmit reference and one-baud range
spacing (1.2 km). These observations have zenith transmission only.

`tibet-meteor mst-fine OUTPUT` adds 1024-lag, pair-count-normalized ACF spectra
using the reusable `pansy_receiver/meteor/spectra` FFT estimator. It must be run
only once per output. The native frequency spacing is 0.095 Hz; the Bartlett
window's first-zero resolution is approximately 0.38 Hz. File gaps are never
treated as observed samples.

Masks retain independent bits for a bad reference/pulse, extreme transient
voltages, known meteor tracks (including their code footprint), and other
strong transient targets. Persistent ground echoes are not called impulsive
RFI. The retained-voltage complex mean is removed separately at each range.
The central 0.5 Hz on either side of zero is not interpreted as atmospheric
power. Persistent narrow lines are not blindly notched because they may
contain the desired MST echo.

Thresholds are conservative working choices, not a validated classifier.
Satellites and RFI are not guaranteed absent. Range aliases remain unresolved;
low-range features can include residual clutter and code sidelobes.
