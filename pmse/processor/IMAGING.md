# Native sparse-product imaging

`src/imaging.rs` implements calibrated seven-channel covariances, a cached
complex TSVD operator, altitude-plane geometry, integration-window splitting,
and full-band circular Doppler moments in Rust. No Python numerical imaging
or Matplotlib rendering is used by these kernels.

The plane is 85.290954301 km above the radar, matching the historical imaging
plane. Every pixel requires its own slant range, sqrt(h*h+east*east+north*north),
not a common range gate for all beams. Five beam-specific operators and phase
calibrations must be used before combining their spatial supports.

The integration policy is one start-to-stop mesomode, split into
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

Status: unit-tested kernels only. Existing-record transport, WebAssembly
bindings, sparse range interpolation with coverage masks, five-beam blending,
and Rustplotlib browser display are not yet wired to this module. The historical
Python imaging scripts are unchanged.

## Browser-only deployment contract

Imaging runs in the web browser through Rust compiled to WebAssembly, not on
Syowa or in a server-side imaging service. Rustplotlib owns browser rendering.
Only the current combined five-beam image is displayed. Existing sparse
complex spectra are the scientific source; do not introduce an image archive,
image HDF5 sidecars, server-rendered images, or another native data product.
Transfer the existing spectra plus the required geometry/calibration and dense
noise metadata. The browser computes covariances, cached TSVD operators,
altitude-plane projection, blending and Doppler color moments locally.

An ephemeral response selecting existing records for transport is not a new
persisted scientific product. Preserve complex phases, all incoherent looks,
missing-bin masks, integration bounds, and calibration provenance in transport.
Keep only bounded current-frame buffers and reuse geometry operators. Changing
colors reuses browser-resident numerical results without another data transfer.

Test with `cargo test --manifest-path pmse/processor/Cargo.toml`.
