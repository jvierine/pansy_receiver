# Bin-by-bin XC2 screening

`pansy-koki-xc2` ports the supplied extreme-filter numerical selection into
single-threaded Rust. It is separate from the older dominant-main-peak screen.
Inputs are XC2 source references in the existing full-day HDF5 overviews.
No acquisition or source-data files are changed.

- Beam 0; nominal slant range 80–92 km using the supplied 0.15 km/sample rule.
- Per-range outer-band median noise; candidates at total/noise >=5 dB.
- Above a 30 dB peak, retain the primary-cluster centroid +/-3 sigma,
  using the same monotonic walk, baseline subtraction and half-bin width floor.
- Per-channel outer-band noise, `max(auto-noise, noise)` normalization, unit
  covariance diagonal, Hermitian PSD projection and 1% trace diagonal loading.
- Same 60x60 grid within 9 degrees, acceptance <=2 degrees. Conjugating both
  covariance and steering gives the shared kernel's equivalent sign convention.
- Histogram every surviving bin at 1 m/s spacing. Retain all accepted cells
  with |velocity| >=15 m/s, including wings underneath a near-zero dominant peak.
  These are spectral-cell counts, not independent-event counts.

## Storage and transfer contract

The supplied full-block export is intentionally replaced by compact event
tables to obey the explicit aggregate storage budget. Dense XC arrays are
**never copied**. Event tables retain native range/frequency indices, source
file/group, arrival direction, power/noise and noise estimates. Existing source
data can be reread by reference. HDF5 fields are compressed; no CSV is emitted.

The station's **entire sweep** has a default 400,000,000-byte budget, including
partial files, event tables and daily summaries. A conservative uncompressed
next-block bound and summary headroom are checked before writing. A budget stop
leaves partial evidence and exits; it does not discard events or report a
completed sweep. CLI limits above 450,000,000 bytes are rejected.

The collector permits at most 200,000,000 bytes of new Koki summaries and
figures. Its static web mirror is a subset of those files. Combined with the
station limit, this leaves substantial headroom below 1,000 MB, including
mirrored results. PNGs use 300 dpi at article dimensions; labels remain 10 pt.
Only small summaries and plots travel, sequentially with `rsync --bwlimit=7`,
leaving headroom for existing 1 KiB/s transfers under the 10 kB/s link budget.
Do not start parallel result collectors.

Use the user service `pmse/systemd/pansy-koki-xc2.service` on Syowa and
`pansy-koki-overview-plots.service` on Revontuli. Both have one-core quotas;
only the latter publishes static artifacts to juha.no. No cron or tunnels.

For the PMSE summer archive on Revontuli, use `pansy-koki-pmse.service`
and `pansy-koki-pmse-plots.service`. These share a one-core user slice. Native
XC2 roots are scanned directly for November–February dates; 2025-01-30 is
processed first. Dedicated XC2-only context overviews avoid substituting an
xc_rvd copy of the same day. Context HDF5 has a separate 100 MB ceiling;
400 MB results plus 100 MB context plus 200 MB collected figures/summaries
and a <=200 MB static mirror remain below the aggregate 1,000 MB budget.

Each detected integration also gets a compact PNG with +/-30 minutes of
full-native-Doppler context and the integration's range-maximum spectrum.
All accepted extreme bins are marked; nearby detections share a context but
are not declared independent physical events. A daily event index links all
integration plots. The spectrum is explicitly the maximum over 80–92 km,
not the spectrum at a single candidate range gate.

Tests cover complex PSD projection, preservation of rank-one PSD covariance,
leakage threshold behavior, and a complete HDF5 binwise fixture where a fast
wing passes despite a stronger near-zero main peak. The synthetic fixture
also checks that no XC array is exported.
