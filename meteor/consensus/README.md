# Meteor consensus and interferometric grouping

The standalone diagnostic does not change production metadata. An optional
realtime grouping backend is described below; offline tests need no restart.

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

## Hourly diagnostic deployment

After installing the native binary, run `bash receiver/install_consensus_shadow.sh`.
This installs ONLY the separate shadow timer, never restarts acquisition or any
production processor. It samples the latest complete minute once per hour and
skips stale (>1 hour), future or insufficient metadata. Each job is limited to
20% of one CPU, 256 MiB and five minutes, with production paths read-only in
the service. No historical catch-up or network transfer is performed.

Only `hourly_*.h5` in `/media/analysis/shadow/meteor_consensus` are rotated:
at most 48 files or 64 MiB, whichever is reached first. Manually requested
`replay_*.h5` files are preserved. Failed/interrupted partial files are cleaned
at the next invocation. Exclusive locking prevents concurrent manual jobs.
Stored model residuals and per-track pulse uniqueness are checked before
publication, and summary counts/runtime are HDF5 attributes. These checks
do NOT establish meteor purity or completeness; visual/labelled validation
and interferometric consistency remain prerequisites for live grouping.

Inspect `journalctl --user -u pansy-consensus-shadow.service` and stop sampling
with `systemctl --user disable --now pansy-consensus-shadow.timer` (stop the
service too to cancel an in-progress sample). Python test:
`conda run -n base python -m unittest discover -s meteor/consensus/tests -p 'test_*.py'`.

## Interferometric grouping backend

`PANSY_GROUPING_BACKEND=rust_phase` replaces the grouping inside `cluster_mf.py`.
For each existing MF peak, it reads the corresponding 20-pulse raw batch,
recomputes only that range gate using the native filter, verifies the Doppler
bin against MF metadata, and retains seven complex receiver amplitudes. The
per-beam phase corrections use the same sign as `image_points1`. Common TX/
meteor phase cancels in the normalized inter-receiver signature comparison.
No phase unwrapping, sky-position choice, or interferometer alias selection is
done here. Array aliases with identical responses remain indistinguishable.

Native grouping forms local phase-compatible components, then uses RANSAC
with quadratic range and corresponding Doppler. Links require squared
phase-signature coherence >=0.65 and separation <=50 ms (including missing
echoes during five-beam cycling). Strong links are merged first; overlapping
components need majority phase consistency to prevent a single outlier from
bridging them. Tests cover indistinguishable radial tracks, common pulse
phase, a bridging outlier and 32-ms measurement gaps.

Sequential fragments can be joined across gaps up to 150 ms only when their
joint range/Doppler fit and independently extrapolated receiver phases agree.
Each of the 21 baselines is fitted directly as wrapped phasors over the adjacent
0.12-second endpoints, searching rates from -300 to 300 rad/s in 2 rad/s steps.
Endpoint fits must have total coherence weight >=7 and weighted phase agreement
>=0.85 at the gap midpoint. Robust radial trimming retains at least 90% of the
combined points and 80% of each fragment. Tests include missing-cycle bridges
and a different direction that must not be joined. This is grouping consistency,
not pulse-to-pulse Doppler refinement or a sky-position solution.

The first production integration retains the existing SNR>7, |radial speed|>
3 km/s and PMSE-range screening. It does NOT yet claim recovery of weak edges
or low-radial-speed meteors. It removes the legacy 300-point truncation, but
fails explicitly above 16384 points. Multi-peak production detection is also
not enabled; it still groups the single peak already stored per pulse.

Detection fields and event keys are unchanged. An independent atomic HDF5
checkpoint advances only after a complete ten-second window; failed writes
do not advance it, and duplicate records on replay are tolerated. Batch tails
are included and clipped to half-open windows. As in the existing detector,
tracks crossing ten-second boundaries are not reconciled yet. Calibration/
raw-data failures stop progress rather than silently falling back to phase-free
grouping. No raw voltages or extra phase arrays are written to live metadata.

Build the shared library (`--features static-hdf5` for the receiver), run
`tests/check_phase_ffi.py`, then `replay_phase.py MF RAW START END OUTPUT.h5`
for explicit <=10-second windows. `plot_phase_replay.py` renders comparisons.
Only after real-data checks, `bash receiver/install_grouping_rust.sh enable`
installs a drop-in and restarts the clustering service alone; `rollback` removes
it. The native cut-filter library is required for its amplitude-return ABI.
