//! Experimental positive covariance imaging. No image-domain deconvolution.
//! A real power distribution maps to a full Hermitian covariance. The packed
//! representation preserves its Frobenius norm (sqrt(2) on cross products).
use crate::imaging::channel_pairs;
use nalgebra::{Complex, DMatrix, DVector};
use std::f64::consts::PI;

type C = Complex<f64>;
const ROWS: usize = 49;

pub fn pack_covariance(c: &[C]) -> Result<Vec<f64>, &'static str> {
    if c.len() != 28 || c.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
        return Err("invalid covariance");
    }
    Ok(c[..7]
        .iter()
        .map(|z| z.re)
        .chain(c[7..].iter().map(|z| 2f64.sqrt() * z.re))
        .chain(c[7..].iter().map(|z| 2f64.sqrt() * z.im))
        .collect())
}

#[derive(Clone, Copy)]
pub struct SparseOptions {
    /// Multiplier of sqrt(2 ln(number of pixels)) in approximately whitened
    /// covariance units. This is a heuristic, not a calibrated false alarm rate.
    pub strength: f64,
    /// Analysis L1 in an orthonormal 2-D Haar basis, rather than pixel L1.
    pub wavelet: bool,
    pub reweight: bool,
    pub max_iterations: usize,
}
impl Default for SparseOptions {
    fn default() -> Self {
        Self {
            strength: 1.,
            wavelet: false,
            reweight: false,
            max_iterations: 200,
        }
    }
}

pub struct SparseFit {
    /// Source powers, NOT squared inverse amplitudes.
    pub power: Vec<f64>,
    pub residual_rms: f64,
    pub relative_residual: f64,
    pub iterations: usize,
    pub converged: bool,
}

pub struct SparseImager {
    forward: DMatrix<f64>,
    gram_rows: DMatrix<f64>,
    grid_coordinates: Vec<[usize; 2]>,
    wavelet_side: usize,
}

impl SparseImager {
    /// Grid coordinates embed the support in a zero-padded power-of-two grid.
    /// Wavelet analysis is H E and its exact adjoint is E^T H^T.
    pub fn new(
        positions: &[[f64; 3]; 7],
        directions: &[[f64; 3]],
        frequency: f64,
        grid_coordinates: &[[usize; 2]],
    ) -> Result<Self, &'static str> {
        if directions.is_empty()
            || directions.len() != grid_coordinates.len()
            || !frequency.is_finite()
            || frequency <= 0.
            || positions.iter().flatten().any(|v| !v.is_finite())
            || directions.iter().any(|d| {
                d.iter().any(|v| !v.is_finite())
                    || (d.iter().map(|v| v * v).sum::<f64>() - 1.).abs() > 1e-6
            })
        {
            return Err("invalid sparse imaging geometry");
        }
        let pairs = channel_pairs();
        let k = 2. * PI * frequency / 299792458.;
        let mut forward = DMatrix::zeros(ROWS, directions.len());
        for (p, d) in directions.iter().enumerate() {
            let c: Vec<C> = pairs
                .iter()
                .map(|&(a, b)| {
                    C::from_polar(
                        1.,
                        k * (0..3)
                            .map(|j| (positions[b][j] - positions[a][j]) * d[j])
                            .sum::<f64>(),
                    )
                })
                .collect();
            forward.set_column(p, &DVector::from_vec(pack_covariance(&c)?));
        }
        let lookup: std::collections::HashMap<_, _> = grid_coordinates
            .iter()
            .copied()
            .enumerate()
            .map(|(i, c)| (c, i))
            .collect();
        if lookup.len() != directions.len() {
            return Err("duplicate grid coordinates");
        }
        let edge = grid_coordinates.iter().flatten().copied().max().unwrap();
        if edge >= 128 {
            return Err("wavelet grid too large");
        }
        let wavelet_side = (edge + 1).next_power_of_two();
        let gram_rows = &forward * forward.transpose();
        Ok(Self {
            forward,
            gram_rows,
            grid_coordinates: grid_coordinates.to_vec(),
            wavelet_side,
        })
    }

    pub fn fit(
        &self,
        covariance: &[C],
        standard_deviation: &[f64],
        options: SparseOptions,
    ) -> Result<SparseFit, &'static str> {
        if standard_deviation.len() != ROWS
            || standard_deviation
                .iter()
                .any(|s| !s.is_finite() || *s <= 0.)
            || !options.strength.is_finite()
            || !(0.05..=5.).contains(&options.strength)
            || !(20..=1000).contains(&options.max_iterations)
            || (options.reweight && options.wavelet)
        {
            return Err("invalid sparse solver options or uncertainties");
        }
        let raw = pack_covariance(covariance)?;
        let weights: Vec<_> = standard_deviation.iter().map(|s| 1. / s).collect();
        let target: Vec<_> = raw.iter().zip(&weights).map(|(v, w)| v * w).collect();
        let n = self.forward.ncols();
        let gram = DMatrix::from_fn(ROWS, ROWS, |i, j| {
            self.gram_rows[(i, j)] * weights[i] * weights[j]
        });
        let lipschitz = gram.symmetric_eigen().eigenvalues.max() * 1.000001;
        if !lipschitz.is_finite() || lipschitz <= 0. {
            return Err("invalid weighted operator");
        }
        let norms: Vec<f64> = (0..n)
            .map(|p| {
                (0..ROWS)
                    .map(|r| (self.forward[(r, p)] * weights[r]).powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .collect();
        let threshold = options.strength * (2. * (n.max(2) as f64).ln()).sqrt();
        let base_penalty: Vec<f64> = norms.iter().map(|v| threshold * v).collect();
        let wavelet = if options.wavelet {
            base_penalty.iter().sum::<f64>() / n as f64
        } else {
            0.
        };
        let mut penalty = if options.wavelet {
            vec![0.; n]
        } else {
            base_penalty.clone()
        };
        let mut x = vec![0.; n];
        let mut iterations = 0;
        let mut converged = true;
        for stage in 0..if options.reweight { 2 } else { 1 } {
            if stage > 0 {
                let peak = x.iter().copied().fold(0., f64::max);
                if peak <= 0. {
                    break;
                }
                for p in 0..n {
                    penalty[p] = base_penalty[p] * (1.05 * peak) / (x[p] + 0.05 * peak);
                }
            }
            let (steps, ok) = self.solve(
                &target,
                &weights,
                &penalty,
                wavelet,
                lipschitz,
                options.max_iterations,
                &mut x,
            );
            iterations += steps;
            converged &= ok;
        }
        // Remove L1 amplitude shrinkage only on the selected compact support.
        // Do not erase the wavelet prior or activate pixels rejected by sparsity.
        if !options.wavelet {
            let peak = x.iter().copied().fold(0., f64::max);
            let active: Vec<_> = x
                .iter()
                .enumerate()
                .filter_map(|(p, &v)| (v > 0. && v >= peak * 0.02).then_some(p))
                .collect();
            if !active.is_empty() && active.len() <= 64 {
                let mut penalty = vec![f64::INFINITY; n];
                for p in &active {
                    penalty[*p] = 0.;
                }
                for p in 0..n {
                    if penalty[p].is_infinite() {
                        x[p] = 0.;
                    }
                }
                let (steps, ok) = self.solve(
                    &target,
                    &weights,
                    &penalty,
                    0.,
                    lipschitz,
                    options.max_iterations,
                    &mut x,
                );
                iterations += steps;
                converged &= ok;
            }
        }
        let prediction = &self.forward * DVector::from_column_slice(&x);
        let error = prediction
            .iter()
            .zip(&raw)
            .zip(&weights)
            .map(|((p, y), w)| ((p - y) * w).powi(2))
            .sum::<f64>();
        let target_norm = target.iter().map(|v| v * v).sum::<f64>();
        if x.iter().any(|v| !v.is_finite() || *v < 0.) || !error.is_finite() {
            return Err("nonfinite sparse fit");
        }
        Ok(SparseFit {
            power: x,
            residual_rms: (error / ROWS as f64).sqrt(),
            relative_residual: (error / target_norm.max(1e-30)).sqrt(),
            iterations,
            converged,
        })
    }

    /// FISTA for positive pixel L1; primal-dual forward-backward for wavelet L1.
    /// All large multiplies are matrix-free in image space (49 covariance rows).
    fn solve(
        &self,
        target: &[f64],
        weights: &[f64],
        penalty: &[f64],
        wavelet: f64,
        lipschitz: f64,
        max_iterations: usize,
        x: &mut Vec<f64>,
    ) -> (usize, bool) {
        let n = x.len();
        let mut z = x.clone();
        let mut t = 1.;
        let mut residual = vec![0.; ROWS];
        let mut gradient = vec![0.; n];
        let mut dual = vec![
            0.;
            if wavelet > 0. {
                self.wavelet_side * self.wavelet_side
            } else {
                0
            }
        ];
        let dual_step = 0.5 * lipschitz;
        let step = if wavelet > 0. {
            0.95 / lipschitz
        } else {
            1. / lipschitz
        };
        let mut stable = 0;
        for iteration in 0..max_iterations {
            for r in 0..ROWS {
                residual[r] = -target[r];
            }
            for p in 0..n {
                if z[p] != 0. {
                    for r in 0..ROWS {
                        residual[r] += self.forward[(r, p)] * weights[r] * z[p];
                    }
                }
            }
            for r in 0..ROWS {
                residual[r] *= weights[r];
            }
            let mut adjoint = dual.clone();
            if wavelet > 0. {
                haar_2d(&mut adjoint, self.wavelet_side, true);
            }
            for p in 0..n {
                gradient[p] = (0..ROWS)
                    .map(|r| self.forward[(r, p)] * residual[r])
                    .sum::<f64>()
                    + if wavelet > 0. {
                        let [a, b] = self.grid_coordinates[p];
                        adjoint[b * self.wavelet_side + a]
                    } else {
                        0.
                    };
            }
            let mut next = vec![0.; n];
            let mut delta = 0.;
            let mut size = 0.;
            for p in 0..n {
                next[p] = (z[p] - step * (gradient[p] + penalty[p])).max(0.);
                delta += (next[p] - x[p]).powi(2);
                size += next[p] * next[p];
            }
            let mut dual_change = 0.;
            let mut dual_size = 0.;
            if wavelet > 0. {
                let mut analysis = vec![0.; dual.len()];
                for (p, &[a, b]) in self.grid_coordinates.iter().enumerate() {
                    analysis[b * self.wavelet_side + a] = 2. * next[p] - x[p];
                }
                haar_2d(&mut analysis, self.wavelet_side, false);
                for (i, (q, v)) in dual.iter_mut().zip(analysis).enumerate() {
                    let next = if i == 0 {
                        0.
                    } else {
                        (*q + dual_step * v).clamp(-wavelet, wavelet)
                    };
                    dual_change += (next - *q).powi(2);
                    dual_size += next * next;
                    *q = next;
                }
                // Leave the global mean unpenalized so diffuse brightness is
                // not explicitly penalized as a point-source count.
                dual[0] = 0.;
                z.clone_from(&next);
            } else {
                let nt = (1. + (1f64 + 4. * t * t).sqrt()) / 2.;
                // Adaptive restart prevents oscillations on a settled support.
                if (0..n)
                    .map(|p| (z[p] - next[p]) * (next[p] - x[p]))
                    .sum::<f64>()
                    > 0.
                {
                    z.clone_from(&next);
                    t = 1.;
                } else {
                    for p in 0..n {
                        z[p] = next[p] + (t - 1.) / nt * (next[p] - x[p]);
                    }
                    t = nt;
                }
            }
            *x = next;
            if delta.sqrt() <= 1e-5 * size.sqrt().max(1e-8)
                && dual_change.sqrt() <= 1e-5 * dual_size.sqrt().max(1e-8)
            {
                stable += 1;
            } else {
                stable = 0;
            }
            if iteration >= 29 && stable >= 5 {
                return (iteration + 1, true);
            }
        }
        (max_iterations, false)
    }
}

/// Orthonormal separable Haar pyramid, in-place with a single row scratchpad.
fn haar_2d(data: &mut [f64], side: usize, inverse: bool) {
    let mut scratch = vec![0.; side];
    let c = 1. / 2f64.sqrt();
    let mut size = if inverse { 2 } else { side };
    while size >= 2 && size <= side {
        let half = size / 2;
        // Inversion reverses column/row ordering as well as scale ordering.
        for axis in 0..2 {
            for line in 0..size {
                let column = if inverse { axis == 0 } else { axis == 1 };
                let index = |p| {
                    if column {
                        p * side + line
                    } else {
                        line * side + p
                    }
                };
                if inverse {
                    for p in 0..half {
                        let a = data[index(p)];
                        let b = data[index(p + half)];
                        scratch[2 * p] = (a + b) * c;
                        scratch[2 * p + 1] = (a - b) * c;
                    }
                } else {
                    for p in 0..half {
                        let a = data[index(2 * p)];
                        let b = data[index(2 * p + 1)];
                        scratch[p] = (a + b) * c;
                        scratch[p + half] = (a - b) * c;
                    }
                }
                for p in 0..size {
                    data[index(p)] = scratch[p];
                }
            }
        }
        if inverse {
            size *= 2;
        } else {
            size /= 2;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn geometry() -> ([[f64; 3]; 7], Vec<[f64; 3]>, Vec<[usize; 2]>) {
        let positions = std::array::from_fn(|i| [i as f64 * 9., (i * i) as f64 * 3., 0.]);
        let directions = (-3..=3)
            .map(|i| {
                let u = i as f64 * 0.04;
                [u, 0.04, -(1. - u * u - 0.04f64.powi(2)).sqrt()]
            })
            .collect();
        (positions, directions, (0..7).map(|p| [p, 0]).collect())
    }
    fn observed(positions: &[[f64; 3]; 7], directions: &[[f64; 3]], power: &[f64]) -> Vec<C> {
        let k = 2. * PI * 47e6 / 299792458.;
        channel_pairs()
            .iter()
            .map(|&(a, b)| {
                directions
                    .iter()
                    .zip(power)
                    .map(|(d, p)| {
                        C::from_polar(
                            *p,
                            k * (0..3)
                                .map(|j| (positions[b][j] - positions[a][j]) * d[j])
                                .sum::<f64>(),
                        )
                    })
                    .sum()
            })
            .collect()
    }
    #[test]
    fn hermitian_packing_preserves_frobenius_norm() {
        let c: Vec<_> = (0..28)
            .map(|i| C::new(i as f64, if i < 7 { 0. } else { 0.5 }))
            .collect();
        let packed = pack_covariance(&c).unwrap();
        assert!(
            (packed.iter().map(|v| v * v).sum::<f64>()
                - c.iter()
                    .enumerate()
                    .map(|(i, v)| v.norm_sqr() * if i < 7 { 1. } else { 2. })
                    .sum::<f64>())
            .abs()
                < 1e-8
        );
    }
    #[test]
    fn point_source_returns_power_not_power_squared() {
        let (pos, dirs, coords) = geometry();
        let im = SparseImager::new(&pos, &dirs, 47e6, &coords).unwrap();
        let mut source = vec![0.; dirs.len()];
        source[4] = 2.;
        let c = observed(&pos, &dirs, &source);
        let fit = im
            .fit(
                &c,
                &[0.01; 49],
                SparseOptions {
                    reweight: true,
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(fit.relative_residual < 0.002, "{}", fit.relative_residual);
        assert!((fit.power[4] - 2.).abs() < 0.01, "{:?}", fit.power);
        assert!(fit
            .power
            .iter()
            .enumerate()
            .all(|(i, p)| i == 4 || *p < 0.01));
    }
    #[test]
    fn resolves_unequal_pair() {
        let (pos, dirs, coords) = geometry();
        let im = SparseImager::new(&pos, &dirs, 47e6, &coords).unwrap();
        let mut source = vec![0.; dirs.len()];
        source[1] = 1.;
        source[5] = 0.25;
        let fit = im
            .fit(
                &observed(&pos, &dirs, &source),
                &[0.001; 49],
                SparseOptions::default(),
            )
            .unwrap();
        assert!(fit.relative_residual < 0.005);
        assert!(
            (fit.power[1] - 1.).abs() < 0.02 && (fit.power[5] - 0.25).abs() < 0.02,
            "{:?}",
            fit.power
        );
    }
    #[test]
    fn zero_data_stays_zero_and_invalid_input_fails() {
        let (pos, dirs, coords) = geometry();
        let im = SparseImager::new(&pos, &dirs, 47e6, &coords).unwrap();
        let zero = vec![C::new(0., 0.); 28];
        assert!(im
            .fit(&zero, &[1.; 49], SparseOptions::default())
            .unwrap()
            .power
            .iter()
            .all(|v| *v == 0.));
        assert!(im.fit(&zero, &[0.; 49], SparseOptions::default()).is_err());
        assert!(im
            .fit(
                &zero,
                &[1.; 49],
                SparseOptions {
                    wavelet: true,
                    reweight: true,
                    ..Default::default()
                }
            )
            .is_err());
    }
    #[test]
    fn extended_wavelet_preserves_positive_patch() {
        let (pos, dirs, coords) = geometry();
        let im = SparseImager::new(&pos, &dirs, 47e6, &coords).unwrap();
        let source = vec![0., 0., 0.5, 1., 0.5, 0., 0.];
        let fit = im
            .fit(
                &observed(&pos, &dirs, &source),
                &[0.001; 49],
                SparseOptions {
                    wavelet: true,
                    max_iterations: 800,
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(fit.relative_residual < 0.02, "{}", fit.relative_residual);
        assert!(fit.power.iter().all(|p| p.is_finite() && *p >= 0.));
        assert!(
            fit.power[2] > 0.3 && fit.power[3] > 0.7 && fit.power[4] > 0.3,
            "{:?}",
            fit.power
        );
    }
    #[test]
    fn haar_is_orthonormal_and_inverse_is_adjoint() {
        let x: Vec<_> = (0..64).map(|i| (i as f64 * 0.7).sin()).collect();
        let y: Vec<_> = (0..64).map(|i| (i as f64 * 0.3).cos()).collect();
        let mut hx = x.clone();
        haar_2d(&mut hx, 8, false);
        let mut ity = y.clone();
        haar_2d(&mut ity, 8, true);
        let lhs = hx.iter().zip(&y).map(|(a, b)| a * b).sum::<f64>();
        let rhs = x.iter().zip(&ity).map(|(a, b)| a * b).sum::<f64>();
        assert!((lhs - rhs).abs() < 1e-10);
        assert!(
            (x.iter().map(|v| v * v).sum::<f64>() - hx.iter().map(|v| v * v).sum::<f64>()).abs()
                < 1e-10
        );
        haar_2d(&mut hx, 8, true);
        assert!(x.iter().zip(hx).all(|(a, b)| (a - b).abs() < 1e-12));
    }
    #[test]
    fn two_dimensional_wavelet_patch_and_negative_noise() {
        let (pos, _, _) = geometry();
        let mut directions = Vec::new();
        let mut coords = Vec::new();
        let mut source = Vec::new();
        for y in 0..7 {
            for x in 0..7 {
                let u = (x as f64 - 3.) * 0.04;
                let v = (y as f64 - 3.) * 0.04;
                directions.push([u, v, -(1. - u * u - v * v).sqrt()]);
                coords.push([x, y]);
                source.push((-0.5 * ((x as f64 - 3.).powi(2) + (y as f64 - 3.).powi(2))).exp());
            }
        }
        let im = SparseImager::new(&pos, &directions, 47e6, &coords).unwrap();
        let fit = im
            .fit(
                &observed(&pos, &directions, &source),
                &[0.01; 49],
                SparseOptions {
                    wavelet: true,
                    strength: 0.5,
                    max_iterations: 1000,
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(fit.relative_residual < 0.08, "{}", fit.relative_residual);
        let peak = fit.power.iter().copied().fold(0., f64::max);
        assert!(
            fit.power.iter().filter(|&&v| v > peak * 0.2).count() >= 8,
            "{:?}",
            fit.power
        );
        let mut negative = vec![C::new(0., 0.); 28];
        negative[..7].fill(C::new(-0.2, 0.));
        let noise = im
            .fit(
                &negative,
                &[1.; 49],
                SparseOptions {
                    wavelet: true,
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(noise.power.iter().all(|v| *v == 0.));
        assert!(noise.residual_rms > 0.);
    }
}
