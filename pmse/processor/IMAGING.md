# Native sparse-product imaging

`src/imaging.rs` implements calibrated seven-channel covariances, a cached
complex TSVD operator, altitude-plane geometry, integration-window splitting,
and full-band circular Doppler moments in Rust. No Python numerical imaging
or Matplotlib rendering is used by these kernels.

The plane is 85.290954301 km above the radar, matching the historical imaging
plane. Every pixel requires its own slant range, sqrt(h*h+east*east+north*north),
not a common range gate for all beams. Five beam-specific operators and phase
calibrations must be used before combining their spatial supports.

The intended integration policy is one start-to-stop mesomode, split into
consecutive half-open windows of at most 120 seconds for long/continuous modes.
The final short window is retained. Live continuous-mode scheduling still needs
to be connected: merely splitting closed-mode records is not sufficient.

Covariance is mean(S_a conj(S_b)), with autocorrelations first, then increasing
unordered channel pairs. Historical calibration uses exp(+i(phi_a-phi_b)) and
positive steering with baseline position_b-position_a and negative vertical
direction cosine. The phase sign is mandatory. This convention is not an
independent validation of contemporary instrument calibration.

The default historical reconstruction uses 18 singular modes and 1/(sigma+1).
The operator exposes both parameters explicitly. Its squared magnitude is a
relative image proxy, not calibrated SNR or absolute echo brightness. Sparse
non-detections must remain censored, not be labeled measured zero power.

Status: unit-tested kernels only. Native Digital Metadata ingestion, sparse
range interpolation with coverage masks, five-beam blending, recolorable HDF5
sidecars, Rustplotlib display, and service deployment are not yet wired to this
module. The historical Python imaging scripts are unchanged.

Test with `cargo test --manifest-path pmse/processor/Cargo.toml`.
