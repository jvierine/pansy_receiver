//! Faithful numerical port of the supplied extreme-filter covariance and leakage rules.
//! Conjugating both covariance and steering preserves a^H R^-1 a: the shared
//! kernel uses negative steering, whereas the supplied implementation uses positive.
use crate::{Matrix, C};
use anyhow::{ensure, Result};

pub fn median(v: &mut [f64]) -> Result<f64> {
    ensure!(
        !v.is_empty() && v.iter().all(|x| x.is_finite()),
        "invalid median input"
    );
    v.sort_unstable_by(f64::total_cmp);
    let k = v.len() / 2;
    Ok(if v.len() % 2 == 0 {
        (v[k - 1] + v[k]) / 2.
    } else {
        v[k]
    })
}

/// Primary-cluster monotonic walk and baseline-subtracted second moment.
pub fn leakage(power: &[f64], frequency: &[f64], noise: f64) -> Result<Vec<bool>> {
    ensure!(
        power.len() == frequency.len()
            && power.len() > 1
            && noise > 0.
            && power.iter().all(|p| p.is_finite() && *p >= 0.),
        "invalid spectrum"
    );
    let pk = (0..power.len())
        .max_by(|&a, &b| power[a].total_cmp(&power[b]).then_with(|| b.cmp(&a)))
        .unwrap();
    if 10. * (power[pk] / noise.max(1e-12)).log10() <= 30. {
        return Ok(vec![true; power.len()]);
    }
    let mut l = pk;
    let mut r = pk;
    while l > 0 && power[l - 1] <= power[l] && power[l] >= 0.05 * power[pk] {
        l -= 1;
    }
    while r + 1 < power.len() && power[r + 1] <= power[r] && power[r] >= 0.05 * power[pk] {
        r += 1;
    }
    let base = power[l].min(power[r]);
    let mut w: Vec<_> = power[l..=r].iter().map(|p| (p - base).max(0.)).collect();
    if w.iter().sum::<f64>() <= 0. {
        w.copy_from_slice(&power[l..=r]);
    }
    let sum = w.iter().sum::<f64>();
    ensure!(sum > 0., "empty leakage cluster");
    let mean = w
        .iter()
        .zip(&frequency[l..=r])
        .map(|(w, f)| w * f)
        .sum::<f64>()
        / sum;
    let variance = w
        .iter()
        .zip(&frequency[l..=r])
        .map(|(w, f)| w * (f - mean).powi(2))
        .sum::<f64>()
        / sum;
    let sigma = variance
        .max(1e-12)
        .sqrt()
        .max(0.5 * (frequency[1] - frequency[0]).abs());
    Ok(frequency
        .iter()
        .map(|f| (f - mean).abs() <= 3. * sigma)
        .collect())
}

/// PSD projection via a symmetric 14x14 real representation of a 7x7 Hermitian
/// matrix. Cyclic Jacobi avoids external BLAS threads and explicit inversion.
pub fn psd(r: Matrix) -> Result<Matrix> {
    const N: usize = 14;
    let mut a = [[0.; N]; N];
    let mut v = [[0.; N]; N];
    for i in 0..N {
        v[i][i] = 1.;
    }
    for i in 0..7 {
        for j in 0..7 {
            ensure!(
                r[i][j].re.is_finite() && r[i][j].im.is_finite(),
                "nonfinite covariance"
            );
            let z = (r[i][j] + r[j][i].conj()) * 0.5;
            a[i][j] = z.re;
            a[i + 7][j + 7] = z.re;
            a[i + 7][j] = z.im;
            a[i][j + 7] = -z.im;
        }
    }
    let scale = a
        .iter()
        .flatten()
        .map(|x| x.abs())
        .fold(0., f64::max)
        .max(1.);
    let mut converged = false;
    for _ in 0..80 {
        let mut largest: f64 = 0.;
        for p in 0..N {
            for q in p + 1..N {
                let x = a[p][q];
                largest = largest.max(x.abs());
                if x.abs() < 1e-14 * scale {
                    continue;
                }
                let tau = (a[q][q] - a[p][p]) / (2. * x);
                let t = if tau >= 0. {
                    1. / (tau + (1. + tau * tau).sqrt())
                } else {
                    -1. / (-tau + (1. + tau * tau).sqrt())
                };
                let c = 1. / (1. + t * t).sqrt();
                let s = t * c;
                a[p][p] -= t * x;
                a[q][q] += t * x;
                a[p][q] = 0.;
                a[q][p] = 0.;
                for k in 0..N {
                    if k != p && k != q {
                        let kp = a[k][p];
                        let kq = a[k][q];
                        a[k][p] = c * kp - s * kq;
                        a[p][k] = a[k][p];
                        a[k][q] = s * kp + c * kq;
                        a[q][k] = a[k][q];
                    }
                    let kp = v[k][p];
                    let kq = v[k][q];
                    v[k][p] = c * kp - s * kq;
                    v[k][q] = s * kp + c * kq;
                }
            }
        }
        if largest < 1e-12 * scale {
            converged = true;
            break;
        }
    }
    ensure!(converged, "PSD eigensolver failed to converge");
    let mut out = [[C::new(0., 0.); 7]; 7];
    for i in 0..7 {
        for j in 0..7 {
            for k in 0..N {
                let d = a[k][k].max(0.);
                out[i][j] += C::new(
                    0.5 * (v[i][k] * v[j][k] + v[i + 7][k] * v[j + 7][k]),
                    0.5 * (v[i + 7][k] * v[j][k] - v[i][k] * v[j + 7][k]),
                ) * d;
            }
        }
    }
    Ok(out)
}

pub fn covariance(
    packed: &[C],
    pairs: &[[usize; 2]],
    phases: &[f64; 7],
    noise: &[f64; 7],
) -> Result<Matrix> {
    let autos = crate::validate_pairs(pairs)?;
    ensure!(packed.len() == 28, "invalid packed covariance");
    let signal: [f64; 7] =
        std::array::from_fn(|i| (packed[autos[i]].norm() - noise[i]).max(noise[i]));
    ensure!(
        signal.iter().all(|x| x.is_finite() && *x > 0.),
        "invalid channel noise"
    );
    let mut r = [[C::new(0., 0.); 7]; 7];
    for i in 0..7 {
        r[i][i] = C::new(1., 0.);
    }
    for (p, &[i, j]) in pairs.iter().enumerate() {
        if i != j {
            let z = packed[p] * C::from_polar(1., phases[i] - phases[j])
                / (signal[i] * signal[j]).sqrt().max(1e-12);
            r[i][j] = z;
            r[j][i] = z.conj();
        }
    }
    psd(r)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn clipping_indefinite_complex_matrix() {
        let mut r = [[C::new(0., 0.); 7]; 7];
        for i in 0..7 {
            r[i][i] = C::new(1., 0.);
        }
        r[0][1] = C::new(0., 2.);
        r[1][0] = C::new(0., -2.);
        let p = psd(r).unwrap();
        assert!((p[0][0].re - 1.5).abs() < 1e-10);
        assert!((p[0][1].im - 1.5).abs() < 1e-10);
        assert!((p[2][2].re - 1.).abs() < 1e-10);
    }
    #[test]
    fn psd_preserves_positive_rank_one() {
        let a: [C; 7] = std::array::from_fn(|i| C::from_polar(1., i as f64 * 0.4));
        let r = std::array::from_fn(|i| std::array::from_fn(|j| a[i] * a[j].conj()));
        let p = psd(r).unwrap();
        for i in 0..7 {
            for j in 0..7 {
                assert!((p[i][j] - r[i][j]).norm() < 1e-10);
            }
        }
    }
    #[test]
    fn leakage_masks_distant_wing_only_for_strong_peak() {
        let f: Vec<_> = (-5..=5).map(|x| x as f64).collect();
        let p = vec![1., 1., 1., 1., 100., 10000., 100., 1., 1., 1., 20.];
        let keep = leakage(&p, &f, 1.).unwrap();
        assert!(keep[5]);
        assert!(!keep[10]);
        assert!(leakage(&p, &f, 100.).unwrap().iter().all(|x| *x));
    }
}
