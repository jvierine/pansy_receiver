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

The legacy historical reconstruction uses 18 singular modes and 1/(sigma+1).
The operator exposes both parameters explicitly. Its squared magnitude is a
relative image proxy, not calibrated SNR or absolute echo brightness. Sparse
non-detections must remain censored, not be labeled measured zero power.

The browser crate in `pmse/browser` calls these kernels through WebAssembly.
`pmse.serve_browser` selects original seven-channel spectra from an existing
record and serves them in an ephemeral binary response; it does not calculate
covariances or images. The historical Python imaging scripts are unchanged.

## Experimental positive sparse reconstruction (2026-10-19)

The Radar image box now defaults to **PMSE wavelet sparsity**. A checkbox enables
sparse reconstruction; its method selector chooses PMSE wavelets or compact
pixel L1. Unchecking it selects the legacy TSVD comparison. Changing method or
regularization reuses one browser-resident source packet, not another download.

`src/sparse_imaging.rs` fits calibrated complex covariance directly, with
nonnegative **power** pixels. It does not deconvolve a rasterized TSVD image and
does not square a recovered power estimate. Seven real autocorrelations and
sqrt(2)-weighted real/imaginary parts of the 21 cross-correlations pack the
Hermitian Frobenius norm into 49 real rows. These are not 49 statistically
independent spatial constraints.

For compact echoes the objective is
`0.5 ||W (A x - y)||² + sum(lambda_p x_p)`, with `x >= 0`.
Positive FISTA uses a cached forward matrix and a weighted row-Gram spectral
step bound. Optional reweighting adds one stage, with weights
`1.05*peak/(x + 0.05*peak)`. A nonnegative unpenalized power refit is restricted
to the selected support (>=2% of peak, at most 64 pixels). This reduces shrinkage
without activating rejected pixels. It does not perform off-grid refinement.

For PMSE the objective replaces pixel L1 by `lambda_w ||H E x||_1`, excluding
the single global-mean coefficient. `E` embeds each beam support in a zero-padded
power-of-two square; `H` is an orthonormal separable 2-D Haar pyramid. A
primal-dual forward-backward solver uses the exact adjoint `E^T H^T`; it uses
neither TV nor a point-source penalty. Zero extension can create support-edge
wavelet effects; this is a first experimental basis, not SARA/SASIR equivalence.
No compact-support refit or pixel reweighting is applied to wavelet images.

The browser subtracts dense channel noise powers from covariance diagonals,
then uses empirical gain scales `s_i=max(C_ii-N_i,N_i)`: its fitted data are
`(C_ij-delta_ij*N_i)/sqrt(s_i*s_j)`, with unit-modulus steering. Approximate
row standard deviations are `sqrt(C_ii*C_jj/(L*s_i*s_j))`, floored with the
dense noise powers. The compact threshold is
`strength*sqrt(2 ln(Npixels))*||W A_column||`; the wavelet threshold is the mean
of these pixel thresholds. This wavelet threshold is a heuristic, not an
independently calibrated coefficient-wise detection test.

Important assumptions: approximately independent proper-Gaussian looks,
mutually incoherent source components, direction-independent channel gain
normalization, and negligible unmodeled spatial noise correlation. Full noise
covariance is not archived, shared-baseline errors are correlated, sparse
selection is conditional on detection, and historical phases are not validated
for contemporary data. Thus residual RMS is an approximate diagnostic, not a
chi-square statistic or confidence level. Images remain relative power, not
calibrated SNR or unbiased source flux. Spatial ambiguities and off-grid sources
can remain even when an image looks sharp.

Each solver stage has a 400-iteration cap. Diagnostics below the image report
mean relative covariance residual, weighted RMS, iterations and cells that did
not reach the iterate-change tolerance. No tiny nonzero solution is normalized
to unit flux. The five-beam support weighting and circular Doppler moments are
retained. A superseded active worker is terminated, computation has a 60-second
timeout, and the host retains at most one 64 MiB source packet. Computation is
browser-only; there are no receiver/service changes or image archive products.

Build the checked-in WASM artifacts with `bash scripts/build_pmse_imaging.sh`
(Rust wasm32 target and wasm-bindgen-cli 0.2.129). `PANSY_CARGO` and
`PANSY_WASM_BINDGEN` can select installed executables. Run native kernel tests
with the command below and browser/real-data checks with
`npm run test:pmse-imaging`. The latter overlays local assets on read-only live
APIs by default; `PMSE_TEST_DEPLOYED=1` checks deployed assets instead.

References are also linked immediately below the radar image:

- Beck & Teboulle (2009), [FISTA](https://www.ceremade.dauphine.fr/~carlier/FISTA).
- Candès, Wakin & Boyd (2008), [reweighted L1](https://arxiv.org/abs/0711.1612).
- Girard et al. (2015), [sparse interferometric imaging](https://arxiv.org/abs/1504.03896).

## Browser-only deployment contract

Imaging runs in the web browser through Rust compiled to WebAssembly, not on
Syowa or in a server-side imaging service. Rustplotlib owns browser rendering.
Only the selected combined five-beam image is displayed. Existing sparse
complex spectra are the scientific source; do not introduce an image archive,
image HDF5 sidecars, server-rendered images, or another native data product.
Transfer the existing spectra plus the required geometry/calibration and dense
noise metadata. The browser computes covariances, cached sparse/legacy operators,
altitude-plane projection, blending and Doppler color moments locally.

An ephemeral response selecting existing records for transport is not a new
persisted scientific product. Preserve complex phases, all incoherent looks,
missing-bin masks, integration bounds, and calibration provenance in transport.
Keep only bounded current-frame buffers and reuse geometry operators. Changing
colors reuses browser-resident numerical results without another data transfer.

Clicking an RTI selects its mesomode and slant-range gate (converting the
displayed nominal height using the clicked beam elevation). All five beams
use that same gate. This is a spherical range shell, not a constant-altitude
plane: its east/north projection is shown with equal physical axis scales.
The reusable Rust geometry also supports the earlier constant-altitude plane,
but that is not the click-selection default. A worker keeps TSVD computation
off the browser UI thread. The 61 by 61 Cartesian grid and nine-degree supports
use the historical positive phase convention, 18 modes and ridge one.

Test with `cargo test --manifest-path pmse/processor/Cargo.toml`.
