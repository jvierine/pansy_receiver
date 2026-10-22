# Legacy XC2 / xc_rvd Capon candidates

Native Rust processing of existing dense Digital Metadata HDF5. No Python,
raw-voltage decoding, acquisition changes, spectral three-width clipping,
receiver-noise subtraction, forced unit covariance diagonal, or eigenvalue
clipping. Output is a compact HDF5 candidate catalogue; source spectra remain
unchanged. This is an exploratory radial-velocity candidate screen, not a
calibrated vertical-draft detector.

Build and test (requires native HDF5 development libraries):

```sh
cargo test --manifest-path pmse/legacy_capon/Cargo.toml --locked
cargo build --release --locked --jobs 1 --manifest-path pmse/legacy_capon/Cargo.toml
```

One-block pilot, zenith beam, 80–92 km, |radial velocity| >= 15 m/s:

```sh
pmse/legacy_capon/target/release/pansy-legacy-capon \
  --root /media/archive/metadata/xc_rvd --day 2026-01-30 \
  --reference . --output /tmp/pansy-capon-pilot.h5 --max-blocks 1
```

For XC2 choose `/media/old_analysis/metadata/xc2` or
`/media/analysis/metadata/xc2`. Multiple `--root` arguments are accepted;
duplicate integration starts use the first encountered copy and are logged.
An invalid first copy is not silently replaced with a second copy.

Keep the station's live processing checkout unchanged. Create a separate Git
worktree at `/home/radar/src/git/pansy_legacy_capon` for this analysis and place
the built binaries under its `pmse/legacy_capon/target/release/` directory.
Install `pmse/systemd/pansy-legacy-capon@.service` in the station user's systemd
directory, reload, and start `pansy-legacy-capon@2026-01-30.service`.
It has a **100% CPU quota (one core total)**, no worker pool, low CPU/I/O
priority, and a 1 GiB memory cap. No cron or automatic retries. Override the
input/output directories with `~/.config/pansy-legacy-capon.env` if needed.
Outputs and stale partial outputs are never overwritten automatically.

Summarize a completed result with `pansy-capon-summary RESULT.h5`; an optional
second argument bounds the maximum event speed shown, e.g. `RESULT.h5 80`.
This is a descriptive event-level view, not a new cell-level rejection rule.

## Axis compatibility

- Read `xc_arr` shape first: `[28,5,nfrequency,nrange]`.
- A stored `fvec` already matching `nfrequency` is authoritative. Legacy
  `xc_rvd` stores all 1,248 bins but `f1=1247`; treating this as an exclusive
  crop would lose a bin. This implementation retains the full stored axis.
- Without `fvec`, derive the FFT-shifted axis from `n_fft`, the 1 MHz sample
  rate (overridable), and `5*ipp`. Full, exclusive-crop, and inclusive-crop
  conventions must reconcile with the stored dimension; otherwise fail.
- Use stored `rvec` when available, otherwise derive slant range with the
  physical speed of light. Range is not mislabeled altitude.
- Validate channel-pair coverage, axis lengths/spacing, powers, and matrix
  inputs. Errors are explicit; input-error counts produce a nonzero exit.

## Method and efficiency

One beam and its full stored range span are read at a time; the full day is never
loaded. Subarray geometry and sensor steering vectors are cached once. Only
bins above the configured power/noise threshold and absolute-velocity limit
are imaged. The threshold retains the supplied legacy definition:
`10 log10(total power / estimated noise)`, not excess signal/noise. Noise is a
per-range median of outer-band summed autocorrelation power; this legacy
estimate can be contaminated by signal and is not an absolute calibration.

Construct `R_ij = XC_ij * exp(+i(phi_i-phi_j))`, with conjugate lower triangle,
and sensor steering `a_i = exp(-ik r_i.s)`, using the historical negative
vertical direction cosine. Add `epsilon * trace(R)/7` to the measured
diagonal, retaining noise. Require real positive Hermitian Cholesky pivots.
For each cached steering vector, evaluate `1/||L^-1 a||^2`; no explicit inverse
or eigendecomposition is performed. Pair-dependent historical averaging can
make legacy matrices indefinite: these cells are flagged invalid rather than
clipped into an apparently valid detection. Default loading is 0.01.

The default odd 61x61 grid includes zenith, within a 9-degree disk. Its maximum
must lie within 2 degrees to pass; search-edge peaks are flagged separately.
Acceptance is a peak-angle criterion, not an angular confidence interval or
proof that grating-lobe alternatives do not exist. Grid resolution and
loading are user-configurable and retained in the output.

Every tested extreme bin is retained with decision 0 (retained), 1
(off-zenith), 2 (invalid covariance), 3 (search boundary), or 4
(excluded as parent-compatible code contamination). The spectral
array is not copied. Source file/group, frequency/range indices, axes, phase
and geometry snapshots, SNR proxy, loading, and angles permit inspection.
The full raw above-threshold occupancy and accepted **extreme** occupancy are
stored separately on each block's native frequency axis: the latter is not a
full-velocity filtered histogram.

Retained cells are grouped as connected components with default links <=0.7 km
and <=2 m/s within a mode or the immediately previous successfully read mode,
with <=600 s separation between mode intervals and matching velocity sign.
Links are transitively closed. These configurable event definitions are not
independence assumptions. Mode counts distinguish single-mode candidates;
no meteor-duration rejection is silently applied. File/group failures break
temporal linkage. The output is flushed per block and atomically published at
the end; partial files remain inspectable if interrupted.

## Legacy code-cycle warning

The full xc_rvd band includes potential code-cycle range-Doppler sidelobes
near +/-99.665 m/s for a 32 ms four-pulse cycle. Complementary autocorrelation
sidelobes in an A,B,B,A cycle have a zero DC term but nonzero quarter-PRF terms;
they can therefore have a near-zenith Capon peak inherited from a real PMSE
source. An angle-only filter does not distinguish them from real high-speed
scatter. The synthetic regression demonstrates the mechanism, not the actual
historical TX sequence. No automatic quarter-PRF exclusion is applied. Treat
these bands as suspect pending source-spectrum/TX checks, not extreme winds.

### What can be corrected from legacy files

An exact conversion to complementary-pair spectra cannot in general be made
from these files. They retain `E[X_i(k) conj(X_j(k))]`, not the cross-frequency
terms. Pair addition followed by factor-two decimation mixes `X(k)` and
`X(k+N/2)` coherently. Its covariance therefore needs terms such as
`E[X_i(k) conj(X_j(k+N/2))]`. Those phases are absent, including in full-band
xc_rvd. A regression constructs two records with identical stored covariance
but different pair-summed power. Spectral folding is not a replacement for
complementary decoding; the narrower XC2 band loses additional information.

For legacy analysis, use a forward model of the actual code-cycle ambiguity
response to test possible parent echoes at the appropriate range and Doppler
offsets. Any power-domain deconvolution additionally assumes sufficiently
uncorrelated scattering cells and a known response/averaging rule; it is not
lossless recovery. Keep candidate flags and original spectra, and do not
silently subtract a template or notch the quarter-PRF bands. Exact re-decoding
requires the original complex pulse records (and measured TX templates).

For raw-data decoding, a simple adjacent complementary sum also has a
Doppler-dependent residual. In the equal-energy stop-and-hop model, if a
nonzero-lag code response is `+c,-c`, the main response is `L,L`, and the pulse
separation is `dt`, the residual-to-main power ratio is
`(c/L)^2 tan^2(pi*f*dt)`. Thus 16 ms output spacing alone does not demonstrate
sidelobe-free decoding. Validate a Doppler-aware decoder against the measured
TX sequence across the full retained velocity band, including its aliases,
before interpreting high-speed candidates as atmospheric motion. The current
tests establish these limitations; they do not claim a historical-data
correction or change the live processor.

## Parent-compatible exclusion (v2, default enabled)

The catalogue now excludes plausible code ghosts from its event membership,
without subtracting power, reconstructing voltages, altering source spectra,
or notching a velocity band. `--no-code-ghost-filter` reproduces the angle-only
selection for comparison. This is an explicit conservative exclusion rule,
not a calibrated posterior probability that an echo is artificial.

For each angularly accepted candidate, search all stored ranges for stronger
parents at the configured code's nonzero range-sidelobe lags, and frequencies
displaced by +/- PRF/4. Wrap frequencies at the original 8 ms Nyquist, never
at the narrower XC2 crop edges. Compute the harmonic power gain from the
aperiodic autocorrelations of the actual configured 16-chip A/B arrays with
8 us rectangular chips and ABBA sequencing (BAAB has identical harmonic power).
This is an equal-energy stop-and-hop model; historical measured TX shapes and
pair-dependent averaging may differ. The filter fails explicitly outside the
configured 1 MHz / 1600 sample IPP mode.

A parent must exceed 10 dB total/noise and the candidate's excess power. Its
predicted ghost excess power must explain at least 50% of the observed
candidate excess, and its 21 complex off-diagonal covariances must have real
cosine similarity >=0.95 with the candidate. Covariances are not phase-rotated
to force agreement. The power estimate subtracts the per-range noise proxy;
the Capon covariance itself still retains noise. Configure these thresholds
with `--ghost-parent-db`, `--ghost-min-fraction`, `--ghost-min-similarity`.

Decision 4 retains the strongest supporting parent's frequency/range indices,
code power gain, predicted fraction and spatial similarity in HDF5. Only
decision 0 enters the new event catalogue; source data and rejected candidates
are preserved. `ghost_parent_search_complete=0` records missing possible
parents outside the stored range/frequency extent (or a cell not tested due
to another decision / diagnostic opt-out). Such retained cells are uncertain,
not demonstrated free of ghosts. There is no automatic cross-mode tracking
requirement: tracking would lose single-mode contamination. Event grouping
still allows inspecting persistence of the surviving candidates.

Voltage-domain synthetic tests verify the code's response at all 128 sample
lags and all four cycle harmonics. Injection tests require that model ghosts
are excluded, different-direction echoes and echoes without parent support
survive, stronger independent excess survives, and cropped parent coverage is
not labeled clean. A real echo coincident in range, Doppler and spatial
covariance with a predicted ghost can also be excluded: the data cannot
distinguish them. These tests are not an empirical false-rejection rate on
atmospheric data. Compare threshold settings before interpreting extremes.

## Archive sweep and small plot transfers

`pansy-capon-sweep` visits populated days across all three supplied metadata
roots, prioritizes January 30, and invokes the processor sequentially. Each
completed day is atomically published; `sweep-index.h5` records pending,
completed and failed days. Existing v2 results are reused, never overwritten.
The `pansy-legacy-capon-sweep.service` quota is 99% of one CPU; the companion
`pansy-legacy-capon-plots.service` quota is 1%, keeping total analysis <=one
core. The companion exits when the sweep has finished and creates compressed
HDF5 plotting projections, not copies of the source spectra. No cron.

Keep full catalogues on Syowa. Transfer only plotting projections below 5 MiB
and rendered plots; use rsync `--max-size=5m --no-perms --no-owner --no-group`.
The Rustplotlib renderer in `pmse/legacy_plots` makes before/after Doppler-UTC
overviews and individual event tracks with native axes and colourbars. It
uses 10 pt type at 190 mm figure width, and saves both PNG and PDF. Each point
is the maximum power/noise across selected ranges in one integration and a
0.5 m/s bin. Blank central velocities were not selected, not measured absent.

## Full Doppler context: default plots

The selected-cell plots above are diagnostic only and are not full spectra.
The default plotting executable is now `pansy-full-doppler`. First run
`pansy-capon-full-doppler CATALOGUE.h5 SMALL_FULL_DOPPLER.h5` **on Syowa or
Revontuli**, with a one-core quota on Syowa. It reads the original source
autocorrelations referenced by the catalogue and retains every native
frequency bin, including zero Doppler, without velocity, angular, candidate
or SNR cuts. The plotted quantity is the maximum seven-channel total/noise
over zenith ranges 80--92 km. A colour limit clips colour only, never removes
bins. The projection retains original mode start/stop times and source
provenance in HDF5. Only small projections and plots leave the processing host.

For xc_rvd this is all 1248 bins over approximately +/-199.33 m/s, not just
the +/-100 m/s sidebands. XC2 can only show its narrower **stored** band;
missing frequencies are not reconstructed. The Rustplotlib renderer makes
a full-day overview, shared-full-axis candidate panels and individual
full-size candidate plots. No central-Doppler hole or zoomed velocity axes.

All computation and compilation belongs on Revontuli or Syowa. **juha.no is
a thin web/transfer host, not a build or processing host.** Build portable
station binaries on Revontuli using the musl target and `--features portable`
(bundled HDF5 and zlib), `CC_x86_64_unknown_linux_musl=musl-gcc`,
`LIBZ_SYS_STATIC=1`, and native `CMAKE=/usr/bin/cmake`. This avoids transferring
a compiler to Syowa or depending on Revontuli's newer glibc/HDF5 ABI.
