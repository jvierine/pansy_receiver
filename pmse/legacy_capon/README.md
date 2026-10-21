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

Install `pmse/systemd/pansy-legacy-capon@.service` in the station user's systemd
directory, reload, and start `pansy-legacy-capon@2026-01-30.service`.
It has a **100% CPU quota (one core total)**, no worker pool, low CPU/I/O
priority, and a 1 GiB memory cap. No cron or automatic retries. Override the
input/output directories with `~/.config/pansy-legacy-capon.env` if needed.
Outputs and stale partial outputs are never overwritten automatically.

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

One beam and the selected range slab are read at a time; the full day is never
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

Every tested extreme bin is retained with decision 0 (accepted angle), 1
(off-zenith), 2 (invalid covariance), or 3 (search boundary). The spectral
array is not copied. Source file/group, frequency/range indices, axes, phase
and geometry snapshots, SNR proxy, loading, and angles permit inspection.
The full raw above-threshold occupancy and accepted **extreme** occupancy are
stored separately on each block's native frequency axis: the latter is not a
full-velocity filtered histogram.

Accepted cells are grouped as connected components with default links <=0.7 km
and <=2 m/s within a mode or the immediately previous successfully read mode,
with <=600 s separation between mode intervals and matching velocity sign.
Links are transitively closed. These configurable event definitions are not
independence assumptions. Mode counts distinguish single-mode candidates;
no meteor-duration rejection is silently applied. File/group failures break
temporal linkage. The output is flushed per block and atomically published at
the end; partial files remain inspectable if interrupted.
