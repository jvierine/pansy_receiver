//! Conservative exclusion, not deconvolution. The configured, equal-energy
//! rectangular-chip ABBA response is a hypothesis for historical measured TX.
//! All decisions retain the supporting parent and model parameters.
use crate::C;
use std::f64::consts::PI;

const A: [f64; 16] = [
    -1., -1., -1., -1., -1., -1., 1., 1., 1., -1., -1., 1., 1., -1., 1., -1.,
];
const B: [f64; 16] = [
    -1., -1., -1., -1., 1., 1., -1., -1., 1., -1., -1., 1., -1., 1., -1., 1.,
];

/// Aperiodic matched-filter response of 8-sample rectangular chips.
fn ac(code: &[f64; 16], lag: usize) -> f64 {
    (0..128 - lag)
        .map(|n| code[n / 8] * code[(n + lag) / 8])
        .sum()
}

/// Power of one code-cycle harmonic relative to the parent's main response.
/// Stop-and-hop approximation within the 128 us pulse; no empirical scaling.
pub fn gain(delay_samples: usize, harmonic: usize) -> f64 {
    if delay_samples >= 128 || harmonic > 3 {
        return 0.;
    }
    let a = ac(&A, delay_samples);
    let b = ac(&B, delay_samples);
    let response: C = [a, b, b, a]
        .iter()
        .enumerate()
        .map(|(n, &x)| C::from_polar(x, -2. * PI * harmonic as f64 * n as f64 / 4.))
        .sum::<C>()
        / 4.;
    response.norm_sqr() / (128. * 128.)
}

/// Real cosine similarity of the off-diagonal covariance vectors. A common
/// calibration cancels in this comparison. No unit-diagonal transformation,
/// noise subtraction from covariance, fitted phase rotation or PSD clipping.
pub fn spatial_similarity(a: &[C], b: &[C]) -> f64 {
    if a.len() != b.len() || a.is_empty() {
        return f64::NAN;
    }
    let mut aa = 0.;
    let mut bb = 0.;
    let mut ab = 0.;
    for (&x, &y) in a.iter().zip(b) {
        aa += x.norm_sqr();
        bb += y.norm_sqr();
        ab += (x.conj() * y).re;
    }
    if !aa.is_finite() || !bb.is_finite() || !ab.is_finite() || aa <= 0. || bb <= 0. {
        return f64::NAN;
    }
    (ab / (aa.sqrt() * bb.sqrt())).clamp(-1., 1.)
}

#[derive(Clone, Copy, Debug)]
pub struct Evidence {
    pub parent_frequency: usize,
    pub parent_range: usize,
    pub gain: f64,
    pub predicted_fraction: f64,
    pub similarity: f64,
}

pub struct Search<'a> {
    pub frequencies: &'a [f64],
    pub ranges: &'a [f64],
    pub pulse_prf: f64,
    pub sample_rate: f64,
    pub power: &'a [f64], // frequency-major, full stored range span
    pub noise: &'a [f64],
    pub min_fraction: f64,
    pub min_similarity: f64,
    pub parent_threshold: f64, // total/noise, not excess SNR
    pub offsets: &'a [(isize, f64)],
}

/// Cache code response once per range grid, rather than per candidate.
pub fn offsets(range_step_km: f64, sample_rate: f64) -> Vec<(isize, f64)> {
    let samples = range_step_km * 2000. * sample_rate / crate::LIGHT;
    if !samples.is_finite() || samples <= 0. {
        return Vec::new();
    }
    let max = (127. / samples).floor() as isize;
    (-max..=max)
        .filter_map(|dr| {
            let delay = (dr as f64 * samples).abs();
            if (delay - delay.round()).abs() > 0.05 {
                return None;
            }
            let g = gain(delay.round() as usize, 1);
            (g > 1e-12).then_some((dr, g))
        })
        .collect()
}

impl Search<'_> {
    /// Match two quarter-PRF images to observed stronger parents. Spectral
    /// aliases wrap at the ORIGINAL 8 ms Nyquist, not at a cropped XC2 edge.
    /// Missing parent bins/ranges mark incomplete coverage, never a clean bill.
    pub fn find(
        &self,
        f: usize,
        r: usize,
        covariance: impl Fn(usize, usize) -> Vec<C>,
    ) -> (Option<Evidence>, bool) {
        let nr = self.ranges.len();
        let observed = (self.power[f * nr + r] - self.noise[r]).max(0.);
        if observed <= 0. {
            return (None, false);
        }
        let df = self.frequencies[1] - self.frequencies[0];
        let target = covariance(f, r);
        let mut complete = true;
        let mut best: Option<Evidence> = None;
        for &(dr, g) in self.offsets {
            let pr = r as isize + dr;
            if pr < 0 || pr >= nr as isize {
                complete = false;
                continue;
            }
            let pr = pr as usize;
            for harmonic in [1., 3.] {
                let pf = (self.frequencies[f] - harmonic * self.pulse_prf / 4.
                    + self.pulse_prf / 2.)
                    .rem_euclid(self.pulse_prf)
                    - self.pulse_prf / 2.;
                let index = ((pf - self.frequencies[0]) / df).round();
                if index < 0. || index >= self.frequencies.len() as f64 {
                    complete = false;
                    continue;
                }
                let index = index as usize;
                if (self.frequencies[index] - pf).abs() > 0.05 * df {
                    complete = false;
                    continue;
                }
                let parent = self.power[index * nr + pr];
                let excess = (parent - self.noise[pr]).max(0.);
                if parent < self.parent_threshold * self.noise[pr] || excess <= observed {
                    continue;
                }
                let fraction = g * excess / observed;
                if fraction < self.min_fraction {
                    continue;
                }
                let similarity = spatial_similarity(&target, &covariance(index, pr));
                if !similarity.is_finite() || similarity < self.min_similarity {
                    continue;
                }
                if best.is_none_or(|b| fraction > b.predicted_fraction) {
                    best = Some(Evidence {
                        parent_frequency: index,
                        parent_range: pr,
                        gain: g,
                        predicted_fraction: fraction,
                        similarity,
                    });
                }
            }
        }
        (best, complete)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn response_matches_voltage_domain_injection() {
        // Explicitly synthesize 32 decoded pulses of a single stationary
        // scatterer, then DFT; not a covariance-only self-consistency test.
        for delay in 0..128 {
            let voltage: Vec<f64> = (0..32)
                .map(|n| ac(if n % 4 == 0 || n % 4 == 3 { &A } else { &B }, delay))
                .collect();
            for h in 0..4 {
                let spectrum: C = voltage
                    .iter()
                    .enumerate()
                    .map(|(n, &v)| C::from_polar(v, -2. * PI * (h * 8) as f64 * n as f64 / 32.))
                    .sum();
                let measured = spectrum.norm_sqr() / (32. * 128_f64).powi(2);
                assert!((measured - gain(delay, h)).abs() < 1e-12);
            }
            assert!(gain(delay, 2) < 1e-25);
        }
    }
    #[test]
    fn exclusions_require_strength_spatial_match_and_visible_parent() {
        let frequencies: Vec<_> = (0..32).map(|i| (i as f64 - 16.) * 125. / 32.).collect();
        let ranges: Vec<_> = (0..100)
            .map(|i| 50. + i as f64 * crate::LIGHT * 4. / 2e9)
            .collect();
        let noise = vec![7.; 100];
        let mut power = vec![7.; 3200];
        let parent = 16 * 100 + 50;
        power[parent] += 70000.;
        let delay = (8..128).step_by(4).find(|&d| gain(d, 1) > 0.001).unwrap();
        let r = 50 + delay / 4;
        power[24 * 100 + r] += 70000. * gain(delay, 1);
        let same = vec![C::new(1., 0.); 21];
        let offsets = offsets(ranges[1] - ranges[0], 1e6);
        let search = Search {
            frequencies: &frequencies,
            ranges: &ranges,
            pulse_prf: 125.,
            sample_rate: 1e6,
            power: &power,
            noise: &noise,
            min_fraction: 0.5,
            min_similarity: 0.95,
            parent_threshold: 10.,
            offsets: &offsets,
        };
        let (hit, complete) = search.find(24, r, |_, _| same.clone());
        assert!(complete);
        let hit = hit.unwrap();
        assert_eq!((hit.parent_frequency, hit.parent_range), (16, 50));
        assert!((hit.predicted_fraction - 1.).abs() < 1e-10);
        // Genuine extra signal in the same cell, with a different direction.
        assert!(search
            .find(24, r, |f, _| if f == 16 {
                same.clone()
            } else {
                (0..7)
                    .flat_map(|i| {
                        ((i + 1)..7).map(move |j| C::from_polar(1., 0.7 * (i as f64 - j as f64)))
                    })
                    .collect()
            })
            .0
            .is_none());
        // No blanket high-velocity rejection, even in the ghost band.
        assert!(search.find(24, 90, |_, _| same.clone()).0.is_none());
        // An identical overlapping real echo is intrinsically ambiguous and
        // conservatively excluded; this is not a statistical certainty claim.
        let mut stronger = power.clone();
        stronger[24 * 100 + r] += 70000. * gain(delay, 1) * 4.;
        let search2 = Search {
            power: &stronger,
            ..search
        };
        assert!(search2.find(24, r, |_, _| same.clone()).0.is_none());
        let cropped = Search {
            frequencies: &frequencies[20..],
            power: &power[20 * 100..],
            ..search
        };
        let (hit, complete) = cropped.find(4, r, |_, _| same.clone());
        assert!(hit.is_none() && !complete);
        assert!(spatial_similarity(&[C::new(0., 0.)], &[C::new(1., 0.)]).is_nan());
    }
}
