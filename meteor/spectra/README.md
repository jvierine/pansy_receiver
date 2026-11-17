# Radar Spectra With Missing Data

Native complex-voltage ACF estimation with valid-pair counts at every lag.
Direct and FFT implementations use the same estimator. Missing observations
are not interpolated, and absent lags remain missing. The Bartlett-windowed
spectrum is signed: negative estimates are retained rather than truncated.

Run the tests with `cargo test --manifest-path meteor/spectra/Cargo.toml`.
They check missing-data amplitude, FFT/direct agreement, spectral frequency
sign, integrated power, and fully missing observations.

Pair-count normalization avoids attenuation from independent missing samples;
it does not restore lost information or guarantee an unbiased result when the
mask depends on the echo. The approximate noise penalty relative to a complete
integration is the square root of the ratio of complete to retained lag pairs.
Spectral covariance also depends on the lag window and the mask distribution.

The Tibet adapter retains masks, decoded voltages, removed coherent means,
signed before/after spectra, and lag counts. Combining complete complementary
A/B pairs avoids unequal code weights after masking. Stationary clutter
subtraction also removes stationary atmospheric signal. Strong targets can
retain range sidelobes and Doppler aliases; a cleaned spectrum alone does not
establish an atmospheric echo.
