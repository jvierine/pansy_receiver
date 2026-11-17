//! Radar ACF estimates with explicit missing-sample pair counts.
use rustfft::{num_complex::Complex32 as C, Fft};

/// No mean is removed here. A mask must describe both ends of each valid product.
/// Data-dependent censoring can bias the ACF even with correct pair counts.
pub fn masked_acf(z: &[C], valid: &[bool], lags: usize) -> (Vec<C>, Vec<u32>) {
    assert_eq!(z.len(), valid.len());
    assert!(lags <= z.len());
    let mut a = vec![C::default(); lags];
    let mut count = vec![0; lags];
    for lag in 0..lags {
        for i in lag..z.len() {
            if valid[i] && valid[i - lag] {
                a[lag] += z[i] * z[i - lag].conj();
                count[lag] += 1;
            }
        }
        a[lag] = if count[lag] > 0 {
            a[lag] / count[lag] as f32
        } else {
            C::new(f32::NAN, f32::NAN)
        };
    }
    (a, count)
}

/// Signed, FFT-shifted spectrum. Negative estimates are not clipped to zero.
/// Frequency-bin spacing is 1/(FFT length * sample interval), not resolution.
pub fn bartlett_spectrum(a: &[C], fft: &dyn Fft<f32>) -> Vec<f32> {
    let nf = fft.len();
    assert!(nf >= 2 * a.len() && nf.is_multiple_of(2));
    let mut z = vec![C::default(); nf];
    z[0] = a[0];
    for l in 1..a.len() {
        let v = a[l] * (1. - l as f32 / a.len() as f32);
        z[l] = v;
        z[nf - l] = v.conj();
    }
    fft.process(&mut z);
    (0..nf).map(|i| z[(i + nf / 2) % nf].re).collect()
}

/// FFT implementation of the same valid-pair estimator, including mask correlation.
pub fn masked_acf_fft(
    z: &[C],
    valid: &[bool],
    lags: usize,
    forward: &dyn Fft<f32>,
    inverse: &dyn Fft<f32>,
) -> (Vec<C>, Vec<u32>) {
    assert_eq!(z.len(), valid.len());
    let n = forward.len();
    assert!(n >= 2 * z.len() && inverse.len() == n && lags <= z.len());
    let mut x = vec![C::default(); n];
    let mut mask = x.clone();
    for i in 0..z.len() {
        if valid[i] {
            x[i] = z[i];
            mask[i].re = 1.;
        }
    }
    forward.process(&mut x);
    forward.process(&mut mask);
    for v in &mut x {
        *v = C::new(v.norm_sqr(), 0.);
    }
    for v in &mut mask {
        *v = C::new(v.norm_sqr(), 0.);
    }
    inverse.process(&mut x);
    inverse.process(&mut mask);
    let count = (0..lags)
        .map(|i| (mask[i].re / n as f32).round().max(0.) as u32)
        .collect::<Vec<_>>();
    let a = (0..lags)
        .map(|i| {
            if count[i] > 0 {
                x[i] / (n as f32 * count[i] as f32)
            } else {
                C::new(f32::NAN, f32::NAN)
            }
        })
        .collect();
    (a, count)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fft_matches_direct_with_gaps() {
        let z = (0..500)
            .map(|i| C::new((i as f32 * 0.37).sin(), (i as f32 * 0.29).cos()))
            .collect::<Vec<_>>();
        let mask = (0..500)
            .map(|i| i % 7 != 0 && !(100..130).contains(&i))
            .collect::<Vec<_>>();
        let (a, n) = masked_acf(&z, &mask, 250);
        let mut planner = rustfft::FftPlanner::new();
        let (b, m) = masked_acf_fft(
            &z,
            &mask,
            250,
            planner.plan_fft_forward(1024).as_ref(),
            planner.plan_fft_inverse(1024).as_ref(),
        );
        assert_eq!(n, m);
        for (x, y) in a.iter().zip(&b) {
            assert!((*x - *y).norm() < 2e-6);
        }
    }
    #[test]
    fn random_gaps_do_not_attenuate_tone() {
        let z = (0..1000)
            .map(|i| C::from_polar(2., i as f32 * 0.13))
            .collect::<Vec<_>>();
        let good = (0..1000).map(|i| i % 4 != 0).collect::<Vec<_>>();
        let (a, count) = masked_acf(&z, &good, 50);
        for l in 0..50 {
            assert!(count[l] > 0);
            assert!((a[l] - C::from_polar(4., l as f32 * 0.13)).norm() < 0.001);
        }
    }
    #[test]
    fn no_valid_pairs_remain_missing() {
        let (a, count) = masked_acf(&[C::new(1., 0.); 8], &[false; 8], 4);
        assert!(a.iter().all(|z| z.re.is_nan()));
        assert_eq!(count, vec![0; 4]);
    }
    #[test]
    fn spectral_frequency_sign_and_power() {
        let z = (0..2048)
            .map(|i| C::from_polar(2., std::f32::consts::TAU * 0.125 * i as f32))
            .collect::<Vec<_>>();
        let (a, _) = masked_acf(&z, &vec![true; z.len()], 128);
        let fft = rustfft::FftPlanner::new().plan_fft_forward(512);
        let p = bartlett_spectrum(&a, fft.as_ref());
        let peak = (0..512).max_by(|&i, &j| p[i].total_cmp(&p[j])).unwrap();
        assert_eq!(peak, 320);
        assert!((p.iter().sum::<f32>() / 512. - 4.).abs() < 1e-4);
    }
}
