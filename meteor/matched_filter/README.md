# Rust meteor range-Doppler filter

Optional native kernel; not enabled by default in receiver services. Receiver calibration is
applied by the caller. Uses measured complex TX samples, forward FFTs and an
incoherent sum of receiver powers. FFT plans and scratch storage are reused.

`Filter::new(gates, txlen, decimation, nfft)` has two explicitly separate modes:

- `decimation=8`: legacy boxcar accumulation, including truncation of incomplete
  groups. Used for numerical parity, not as the preferred new algorithm.
- `decimation=1`: every TX sample and the full native Doppler band. With 132 TX
  samples this retains the four samples lost by the old eight-sample sums.

`Output::peaks` finds multiple local maxima with configurable range/frequency
suppression instead of discarding everything except one maximum per pulse.
Every peak retains complex amplitudes from all receivers at the *same* range
and Doppler cell, for subsequent calibrated interferometry. Frequency sign
matches Python's forward FFT. Range is the supplied sample lag; cut-relative
lags must have the recorded cut delay restored before reporting physical range.

The median-power ratio is retained for comparison, not advertised as calibrated
SNR or false-alarm probability. Removing boxcar accumulation changes the noise
distribution across Doppler; thresholds cannot simply be copied between modes.
Code sidelobes and spectral leakage can also create candidate peaks. Multiple
candidates are not proof of multiple meteors.

## Verification

```sh
cargo build --release --locked --manifest-path meteor/matched_filter/Cargo.toml
conda run -n base python meteor/matched_filter/preflight.py
cargo test --manifest-path meteor/matched_filter/Cargo.toml
conda run --no-capture-output -n base python meteor/matched_filter/compare_reference.py --cut-root data/metadata/cut
cargo clippy --manifest-path meteor/matched_filter/Cargo.toml --all-targets -- -D warnings
```

The Python test extracts the actual legacy filter class without importing its
MPI/service startup. It compares all power cells and the median, plus independent
complex FFT values at extracted peaks. Fixtures are temporary HDF5. Tests cover
three recorded cut pulses and two simultaneous injected targets, in both modes;
the full-sample case must recover the injected range/frequency cells, including
a target beyond the legacy decimated Nyquist band. This does not validate event
recall, a production threshold, or interferometric alias selection.

Timings use 20 warm calls, exclude I/O and peak extraction, and compare RustFFT
with the existing Python class using SciPy FFT (not the production pyFFTW
backend). Recorded cut gates cover only a narrow window, not the whole online
search. Do not extrapolate these timings to station throughput.

Next gates: noise-only/weak-injection false-alarm tests, resolution-aware
suppression settings, calibrated angle candidates connected to consensus,
and bounded raw-window shadow replay. Only after these checks should a
prebuilt, reversible service update be considered. No live service is changed
by building or testing this crate.

## First live integration

`quick_search_meteor.py` selects the ctypes adapter with `PANSY_MF_BACKEND=rust`.
The default remains `python`. The live adapter deliberately keeps the original
decimation, single-peak selection, calibration, SNR, metadata and MPI checkpoint
logic unchanged. Full-band and multi-peak detection are not enabled by this
switch. Unknown backend names and unavailable libraries fail explicitly.

Deploy through Git, build the shared library while the old service runs, then
run `preflight.py` with the receiver's Python environment. It tests 20 full-range
seven-receiver pulses against production pyFFTW, including range profiles,
Doppler selections and median noise. It opens no metadata and starts no MPI job.
Only after it passes should the quick-search service receive the environment
override and be restarted. Keep a rollback to `PANSY_MF_BACKEND=python` ready;
verify Rust startup messages and advancing MF metadata, not merely active status.
Never restart acquisition or other processing services for this change.

## Diagnostic peak interface

`RustMatchedFilter.peaks` exposes up to 256 peaks (default eight), each with
complex64 amplitudes at the identical range/Doppler cell for every receiver.
It is a separate call; production `process` and its ABI are unchanged.
Amplitudes are not phase-calibrated directions. Test with
`preflight_peaks.py`, which independently checks the complex FFT and power
normalization for two injected targets at both decimation settings.

`shadow_raw_peaks.py RAW MF BATCH_SAMPLE OUTPUT.h5` reads exactly one existing
20-pulse MF batch, applies the same receiver amplitude calibration and TX gate,
and checks every exported complex amplitude against NumPy. It retains only
the compact peak diagnostics, not raw voltages or full spectra. Run as a
resource-limited one-shot job. Additional peaks may be sidelobes, not meteors;
calibrated interferometric consistency and false-alarm validation remain open.
