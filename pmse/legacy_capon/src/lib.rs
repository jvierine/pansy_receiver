//! Legacy dense XC2/xc_rvd adapter and diagonally loaded Capon kernels.
//! Cross-products are C_ij = <z_i conj(z_j)>, positive baseline phase;
//! the matching sensor steering is exp(-ik r_i.s). No eigenvalue clipping.
// Fixed seven-sensor triangular kernels deliberately use explicit indices.
#![allow(clippy::needless_range_loop)]
use anyhow::{bail, ensure, Result};
use num_complex::Complex64;
use std::f64::consts::PI;
pub type C = Complex64;
pub type Matrix = [[C; 7]; 7];
pub const LIGHT: f64 = 299_792_458.;

#[derive(Clone, Debug)]
pub struct Axes {
    pub frequency_hz: Vec<f64>,
    pub range_km: Vec<f64>,
    pub frequency_source: &'static str,
}

/// The stored dimension wins over inconsistent historical f1 conventions.
/// In xc_rvd, fvec and xc_arr are full length but f1 is the LAST bin, not
/// an exclusive slice end. Never apply f0:f1 to an already matching fvec.
#[allow(clippy::too_many_arguments)]
pub fn axes(
    nf: usize,
    nr: usize,
    nfft: usize,
    ipp: usize,
    rate: f64,
    f0: usize,
    f1: usize,
    r0: usize,
    r1: usize,
    rdec: usize,
    fvec: Option<Vec<f64>>,
    rvec: Option<Vec<f64>>,
) -> Result<Axes> {
    ensure!(
        nf > 1
            && nr > 0
            && nfft >= nf
            && ipp > 0
            && rate.is_finite()
            && rate > 0.
            && rdec > 0
            && r1 > r0,
        "invalid axis metadata"
    );
    let df = rate / (5. * ipp as f64 * nfft as f64);
    let full = || {
        (0..nfft)
            .map(|j| (j as f64 - (nfft / 2) as f64) * df)
            .collect::<Vec<_>>()
    };
    let (frequency_hz, frequency_source) = match fvec {
        Some(v) if v.len() == nf => (v, "stored_fvec_matching_shape"),
        Some(v) if v.len() == nfft && f1 >= f0 && f1 - f0 == nf && f1 <= nfft => {
            (v[f0..f1].to_vec(), "stored_full_fvec_exclusive_crop")
        }
        Some(_) => bail!("fvec length is inconsistent with xc_arr"),
        None if nf == nfft && f0 == 0 && (f1 == nfft || f1 + 1 == nfft) => {
            (full(), "full_fft_shape_legacy_f1_safe")
        }
        None if f1 >= f0 && f1 - f0 == nf && f1 <= nfft => {
            (full()[f0..f1].to_vec(), "reconstructed_exclusive_crop")
        }
        None if f1 >= f0 && f1 - f0 + 1 == nf && f1 < nfft => {
            (full()[f0..=f1].to_vec(), "reconstructed_inclusive_crop")
        }
        None => bail!("cannot reconcile f0/f1/n_fft with stored spectral dimension"),
    };
    ensure!(
        frequency_hz.iter().all(|f| f.is_finite())
            && frequency_hz
                .windows(2)
                .all(|w| (w[1] - w[0] - df).abs() <= df * 1e-5)
            && frequency_hz
                .iter()
                .all(|f| f.abs() <= rate / (10. * ipp as f64) + df * 1e-5),
        "invalid/nonuniform Doppler axis or unexpected effective IPP"
    );
    let range_km = match rvec {
        Some(v) => {
            ensure!(v.len() == nr, "rvec length mismatch");
            v
        }
        None => {
            ensure!((r1 - r0).div_ceil(rdec) == nr, "range shape mismatch");
            (0..nr)
                .map(|j| (r0 + j * rdec) as f64 * LIGHT / (2. * rate * 1000.))
                .collect()
        }
    };
    ensure!(
        range_km.iter().all(|r| r.is_finite() && *r > 0.)
            && range_km.windows(2).all(|w| w[1] > w[0]),
        "invalid range axis"
    );
    Ok(Axes {
        frequency_hz,
        range_km,
        frequency_source,
    })
}

pub fn validate_pairs(pairs: &[[usize; 2]]) -> Result<[usize; 7]> {
    ensure!(pairs.len() == 28, "need all 28 channel pairs");
    let mut seen = [[false; 7]; 7];
    let mut autos = [0; 7];
    for (p, &[i, j]) in pairs.iter().enumerate() {
        ensure!(i < 7 && j < 7, "channel index outside 0..7");
        let (a, b) = (i.min(j), i.max(j));
        ensure!(!seen[a][b], "duplicate channel pair");
        seen[a][b] = true;
        if i == j {
            autos[i] = p;
        }
    }
    ensure!(
        (0..7).all(|i| (i..7).all(|j| seen[i][j])),
        "missing channel pair"
    );
    Ok(autos)
}

pub fn covariance(packed: &[C], pairs: &[[usize; 2]], phases: &[f64; 7]) -> Result<Matrix> {
    ensure!(
        packed.len() == 28 && pairs.len() == 28 && phases.iter().all(|p| p.is_finite()),
        "invalid covariance input"
    );
    let mut r = [[C::new(0., 0.); 7]; 7];
    for (&z, &[i, j]) in packed.iter().zip(pairs) {
        ensure!(
            i < 7 && j < 7 && z.re.is_finite() && z.im.is_finite(),
            "nonfinite covariance"
        );
        if i == j {
            ensure!(
                z.re > 0. && z.im.abs() <= 1e-4 * z.re,
                "invalid autocorrelation"
            );
            r[i][i] = C::new(z.re, 0.);
        } else {
            let calibrated = z * C::from_polar(1., phases[i] - phases[j]);
            r[i][j] = calibrated;
            r[j][i] = calibrated.conj();
        }
    }
    Ok(r)
}

pub struct Capon {
    pub directions: Vec<[f64; 3]>,
    steering: Vec<[C; 7]>,
    pub grid_step_sine: f64,
    radius_sine: f64,
}
pub struct Peak {
    pub theta_deg: f64,
    pub u: f64,
    pub v: f64,
    pub power: f64,
    pub boundary: bool,
    pub loading: f64,
}
impl Capon {
    pub fn new(
        positions: &[[f64; 3]; 7],
        frequency: f64,
        radius: f64,
        side: usize,
    ) -> Result<Self> {
        ensure!(
            frequency.is_finite()
                && frequency > 0.
                && radius.is_finite()
                && radius > 0.
                && radius < 90.
                && (3..=129).contains(&side)
                && side % 2 == 1
                && positions.iter().flatten().all(|x| x.is_finite()),
            "invalid Capon geometry"
        );
        let edge = radius.to_radians().sin();
        let step = 2. * edge / (side - 1) as f64;
        let k = 2. * PI * frequency / LIGHT;
        let mut directions = Vec::new();
        let mut steering = Vec::new();
        for j in 0..side {
            for i in 0..side {
                let u = -edge + i as f64 * step;
                let v = -edge + j as f64 * step;
                if u * u + v * v > edge * edge + 1e-15 {
                    continue;
                }
                let d = [u, v, -(1. - u * u - v * v).sqrt()];
                directions.push(d);
                steering.push(std::array::from_fn(|c| {
                    C::from_polar(1., -k * (0..3).map(|a| positions[c][a] * d[a]).sum::<f64>())
                }));
            }
        }
        Ok(Self {
            directions,
            steering,
            grid_step_sine: step,
            radius_sine: edge,
        })
    }
    pub fn peak(&self, r: Matrix, epsilon: f64) -> Result<Peak> {
        ensure!(epsilon.is_finite() && epsilon > 0., "invalid loading");
        let scale = (0..7).map(|i| r[i][i].re).sum::<f64>() / 7.;
        ensure!(scale.is_finite() && scale > 0., "invalid covariance trace");
        // Scale for numerical stability. Retain all measured noise and powers.
        let mut loaded = r.map(|row| row.map(|v| v / scale));
        for c in 0..7 {
            loaded[c][c] += C::new(epsilon, 0.);
        }
        // Complex Cholesky APIs may take a complex square root of a negative
        // pivot. Require a real positive Hermitian pivot explicitly instead.
        let mut lower = [[C::new(0., 0.); 7]; 7];
        for i in 0..7 {
            for j in 0..=i {
                ensure!(
                    (loaded[i][j] - loaded[j][i].conj()).norm() < 1e-8
                        && loaded[i][j].re.is_finite()
                        && loaded[i][j].im.is_finite(),
                    "non-Hermitian covariance"
                );
                let mut value = loaded[i][j];
                for k in 0..j {
                    value -= lower[i][k] * lower[j][k].conj();
                }
                if i == j {
                    ensure!(
                        value.re > 0. && value.re.is_finite() && value.im.abs() < 1e-8,
                        "loaded covariance is not positive definite"
                    );
                    lower[i][i] = C::new(value.re.sqrt(), 0.);
                } else {
                    lower[i][j] = value / lower[j][j].re;
                }
            }
        }
        let mut best = (f64::INFINITY, 0usize);
        // a^H R^-1 a = ||L^-1 a||^2: only one triangular solve per pixel.
        for (p, a) in self.steering.iter().enumerate() {
            let mut y = [C::new(0., 0.); 7];
            for i in 0..7 {
                let mut s = a[i];
                for j in 0..i {
                    s -= lower[i][j] * y[j];
                }
                y[i] = s / lower[i][i];
            }
            let den = y.iter().map(|z| z.norm_sqr()).sum::<f64>();
            if den < best.0 {
                best = (den, p);
            }
        }
        ensure!(best.0.is_finite() && best.0 > 0., "invalid Capon power");
        let d = self.directions[best.1];
        let sine = d[0].hypot(d[1]);
        Ok(Peak {
            theta_deg: sine.min(1.).asin().to_degrees(),
            u: d[0],
            v: d[1],
            power: scale / best.0,
            boundary: sine >= self.radius_sine - 1.5 * self.grid_step_sine,
            loading: epsilon * scale,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rvd_full_axis_does_not_drop_last_bin() {
        let f: Vec<_> = (0..1248).map(|j| (j as f64 - 624.) / 9.984).collect();
        let a = axes(
            1248,
            100,
            1248,
            1600,
            1e6,
            0,
            1247,
            334,
            733,
            4,
            Some(f.clone()),
            None,
        )
        .unwrap();
        assert_eq!(a.frequency_hz, f);
        assert_eq!(a.frequency_hz.len(), 1248);
        let b = axes(1248, 100, 1248, 1600, 1e6, 0, 1247, 334, 733, 4, None, None).unwrap();
        assert_eq!(b.frequency_hz.len(), 1248);
        assert!((b.frequency_hz[1247] - f[1247]).abs() < 1e-10);
    }
    #[test]
    fn xc2_crop_and_bad_axes() {
        assert_eq!(
            axes(538, 67, 1248, 1600, 1e6, 355, 893, 400, 667, 4, None, None)
                .unwrap()
                .frequency_hz
                .len(),
            538
        );
        assert!(axes(
            538,
            67,
            1248,
            1600,
            1e6,
            355,
            894,
            400,
            667,
            4,
            Some(vec![0.; 3]),
            None
        )
        .is_err());
    }
    #[test]
    fn calibrated_known_direction_and_noise() {
        let pos = std::array::from_fn(|i| [i as f64 * 9., (i * i) as f64 * 3., i as f64 * 0.1]);
        let c = Capon::new(&pos, 47e6, 9., 61).unwrap();
        let p = c
            .directions
            .iter()
            .position(|d| d[0] > 0.04 && d[1] < -0.025)
            .unwrap();
        let a = c.steering[p];
        let phases = [0., 0.7, -0.3, 1.2, -1., 0.1, 0.4];
        let raw: [C; 7] = std::array::from_fn(|i| a[i] * C::from_polar(1., -phases[i]));
        let pairs: Vec<_> = (0..7).flat_map(|i| (i..7).map(move |j| [i, j])).collect();
        validate_pairs(&pairs).unwrap();
        let packed: Vec<_> = pairs
            .iter()
            .map(|&[i, j]| {
                raw[i] * raw[j].conj()
                    + if i == j {
                        C::new(0.01 * (i + 1) as f64, 0.)
                    } else {
                        C::new(0., 0.)
                    }
            })
            .collect();
        let r = covariance(&packed, &pairs, &phases).unwrap();
        let peak = c.peak(r, 0.01).unwrap();
        assert!((peak.u - c.directions[p][0]).abs() < 1e-12);
        assert!((peak.v - c.directions[p][1]).abs() < 1e-12);
        // A rank-one covariance remains solvable with positive loading.
        assert!(c
            .peak(
                std::array::from_fn(|i| std::array::from_fn(|j| a[i] * a[j].conj())),
                0.01
            )
            .is_ok());
    }
    #[test]
    fn invalid_covariance_is_not_silently_repaired() {
        let c = Capon::new(&[[0.; 3]; 7], 47e6, 9., 11).unwrap();
        let mut r: Matrix = std::array::from_fn(|i| {
            std::array::from_fn(|j| C::new(if i == j { 1. } else { 0. }, 0.))
        });
        r[0][1] = C::new(2., 0.);
        r[1][0] = C::new(2., 0.);
        assert!(c.peak(r, 0.01).is_err());
        assert!(validate_pairs(&[[0, 0]; 28]).is_err());
    }
}
