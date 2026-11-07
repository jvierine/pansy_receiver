# Offline meteor consensus prototype

This is an isolated first step, not connected to the realtime services or their
metadata writers. No deployment or service restart is required to test it.

Uses a pinned Git revision of `rdigital_rf` with its numeric Digital Metadata
reader. Native HDF5 is required; `--features static-hdf5` builds HDF5 into the
binary when the build host and receiver have different HDF5 ABIs. No private
observations belong in Git.

```sh
cargo test --manifest-path meteor/consensus/Cargo.toml
cargo run --release --manifest-path meteor/consensus/Cargo.toml -- \
  /path/to/metadata/mf START_SAMPLE END_EXCLUSIVE
```

Reads explicit half-open windows up to 10 seconds, at most 4096 input
rows (oversized windows are errors, not silently truncated). Use `--detections`
to exercise the same reader on already selected legacy detections; that is an
I/O smoke test, not an independent detection-recall measurement. The command
prints a diagnostic summary. `--output FILE.h5` optionally writes diagnostic
points and membership atomically; it does not write production metadata.
Legacy floating-point pulse indices are accepted only when finite, integral,
nonnegative and below 2^53; they are never silently rounded.

Local smoke replay: detections channel, samples
`1746541806935733..1746541816935733`, contains two records and 104 points.
The initial constant radial-speed prototype fragmented these into five groups.
The quadratic model returns two groups (78 and 24 points) in approximately
3 ms on the development Mac. This is not ground truth or a receiver benchmark.

## Method and limits

- Deterministic, bounded pair consensus with local least-squares refitting.
  Range and radial velocity must agree with the same quadratic range model
  with constant radial acceleration. This is not the final 3-D trajectory
  fit, and radial acceleration is not an estimate of physical drag alone.
- Strong seeds use the existing linear SNR 7 threshold (8.45 dB). Growth uses
  linear SNR 2 (3.01 dB), at least four strong seeds, six distinct pulses and
  60 ms duration. These are provisional values, not tuned on real data.
- Range tolerance is 0.5 km and radial-velocity tolerance is 3 km/s. No minimum
  radial speed is imposed. Head/trail discrimination is still required later.
- Groups cannot bridge 150 ms gaps. Selected rows are exclusive across tracks;
  near-duplicate peaks within a selected pulse/model are suppressed. At most
  one point per pulse is retained within a track. Five-beam cycling is handled
  through the actual sample timestamps, not a presumed adjacent-pulse interval.
- At most 512 hypotheses per extracted track and 16 tracks per input window.
  This controls work, but does not prove real-time throughput or completeness.
- No receiver-phase observables are present in the current `mf` channel.
  Do not substitute transmit beam IDs for interferometric direction. Tracks
  indistinguishable in range/Doppler cannot yet be separated. The raw producer
  currently saves only one peak per pulse, so it may already have lost weaker
  simultaneous echoes before this prototype sees them.
- Quadratic range is a short-arc approximation. Higher-order curvature,
  threshold calibration, window-boundary reconciliation, and alias-aware
  angular continuity remain to be validated before a production cutover.

## Next gates

1. Compare old/new grouping on bounded raw `mf` windows with manually identified
   overlapping meteors and weak edges; record purity, recall and runtime in HDF5.
2. Preserve additional local peaks and calibrated complex receiver amplitudes
   in a separate shadow metadata channel. Compare interferometric consistency
   without prematurely choosing one alias.
3. Shadow-run on the receiver with bounded CPU/I/O and independent checkpoints.
   Do not write to `detections`, `cut`, or advance their progress markers.
4. Only after inspection and throughput checks, switch one service with a
   rollback path; any interruption must remain below 60 seconds. Never rebuild
   binaries or replay historical data while the production service is stopped.

## Finite station shadow replay

`shadow_batch.py METADATA START_SAMPLE END_EXCLUSIVE OUTPUT.h5` processes at
most five minutes in two-second windows. Reading and fitting are native Rust;
Python only runs bounded jobs and collects HDF5 diagnostics. Each window must
finish within 30 seconds. A failed window aborts the batch, not silently skips
data. Use a resource-limited transient service (CPUQuota=20%, Nice=10,
MemoryMax=256M, RuntimeMaxSec=300), without modifying any production service.

The MF reader includes the preceding 32 ms batch, clips pulse times exactly,
and applies the legacy TX-power gate (all pulses > 1e9) before fitting. Only
finite points above the growth threshold enter consensus. Float products are
float32, sample IDs int64, and membership IDs int32 (-1 means unassigned).
These windows are not reconciled across boundaries: track counts are window
fragments, NOT a meteor catalogue or detection-recall estimate. No angular
disambiguation or head/trail decision is made by this prototype.

Integration test: `conda run -n base python meteor/consensus/tests/check_cli.py`.
