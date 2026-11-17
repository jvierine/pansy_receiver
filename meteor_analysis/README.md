# Shared Meteor Analysis

Offline, site-independent analysis with a Tibet cut adapter. Nothing here starts
acquisition, changes PANSY services, or overwrites the detection catalogue.

- `stations.py`: explicit geodetic sites, including Tibet MST (30.1 N, 90.5 E,
  4287 m). Array orientation is deliberately not inferred from station location.
- `measurements.py` and the Rust library: common-TX fractional-delay and Doppler
  refinement, energy-centred complex phases, and measured off-echo covariance.
- `geometry.py`: native full-element array responses with a checked interpolation
  table, preserving the provisional powered/dummy mask.
- `trajectory.py`: shared joint trajectory fitting and shrinking-radius profiles.
- `pulse_phase.py`: nonoverlapping velocity pairs with branch candidates, actual
  scattering-time intervals, pointing-change guards, and same-code diagnostics.
- `atmosphere.py`: explicit MSIS 2.1 inputs for the observation date and location;
  network tables are parsed in memory and retained as HDF5, not CSV products.
- `coordinates.py`: generic ENU-to-inertial conversion and the incoming
  Earth-hyperbola asymptote, including zenith attraction and Earth rotation.
- `meteor/analysis/src/plotting.rs`: reusable Rustplotlib event, channel-phase,
  and radiant figures at their final physical size and body font.

## Build And Run

From this repository, with HDF5 installed (Rustplotlib is pinned in Cargo):

```sh
cargo build --release --locked --manifest-path meteor/analysis/Cargo.toml --features plotting
conda run -n base python -m meteor_analysis.tibet CUTS --run DETECTOR_RUN --geometry GEOMETRY.h5 --output NEW_ANALYSIS
meteor/analysis/target/release/meteor-plot NEW_ANALYSIS NEW_ANALYSIS/figures
conda run -n base python -m meteor_analysis.report NEW_ANALYSIS
conda run -n base python -m meteor_analysis.metadata NEW_ANALYSIS
conda run -n base python -m meteor_analysis.verify NEW_ANALYSIS
```

Python requires numpy, scipy, h5py, astropy, pymsis and pypdf. The report requires
LaTeX. DSP and dynamics use the native Rust library through a small ctypes API;
SciPy performs optimization and Astropy handles reference frames. All scientific
data products are HDF5. Numerical inputs, source hashes and model assumptions are
stored with the results, not rendered as source-code credits.
Digital Metadata export uses `digital_rf`; PDF validation uses PyMuPDF. Run the
focused checks with `conda run -n base python -m pytest tests/test_meteor_analysis.py`
after building the native library, and use `cargo test` and `cargo clippy` with
the same manifest and `--all-features` for the Rust components.

## Interpretation

The Tibet antenna axes and receive phase convention are not confirmed. Its sky
map is explicitly illustrative under x-east/y-north assumptions; array-frame
radiants are also supplied. The nominal range branch and boresight-ensemble phase
calibration are conditional, not independent absolute calibrations.

The single-body model follows the existing PANSY Ceplecha implementation:
`dv/dt = -(3/4) rho_air |v| v / (rho_material r)` and
`dr/dt = -rho_air sigma |v|^3 / (4 rho_material)`, with optional local gravity.
The default material density is 3000 kg/m3 and sigma is 1e-8 kg/J. Winds, Coriolis,
fragmentation and a changing ablation coefficient are not modeled. Radius is at
the first observed point. A fixed-radius profile still permits radius to shrink
after that initial point. Upper-open profiles remain open.

Range and within-pulse Doppler retain their local 2x2 covariance. Angular fitting
uses the measured ten-channel covariance, shared-reference phase covariance, and
event-constant phase-offset nuisance parameters. Pulse-phase velocity branches
are selected jointly across each event using a robust within-pulse Doppler trend
and continuity, rather than independent rounding against a constant-velocity
trajectory. Ambiguous alternatives are retained as diagnostics; branches are
not globally marginalized. Adjacent-pulse pairs are
nonoverlapping; same-code two-IPP pairs are a separate, correlated diagnostic and
are not counted as additional independent data. Working scatter inflation is
estimated once and held fixed across each radius profile.

Profile thresholds and local covariance are conditional diagnostics, not a
validated total uncertainty or a marginal radius posterior. The explicit
density/radius scaling degeneracy means atmosphere and material assumptions can
dominate the physical radius uncertainty. Run simulated coverage studies and
independent calibration before interpreting narrow profile intervals literally.

## Tests

```sh
cargo test --locked --manifest-path meteor/analysis/Cargo.toml --all-features
cargo clippy --locked --manifest-path meteor/analysis/Cargo.toml --all-targets --all-features -- -D warnings
conda run -n base python -m pytest tests/test_meteor_analysis.py -q
```

Checks include an independent complex direct sum, an independent SciPy ODE and
time-step halving, synthetic phase-sign/branch recovery, nonoverlap and pointing
guards, Earth-hyperbola energy/direction, the apex-centred map convention, and
joint shrinking-radius recovery. These are correctness checks, not interval
coverage or a validated complete meteor catalogue.
