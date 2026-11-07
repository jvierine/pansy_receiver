//! Offline prototype: joint range/Doppler consensus, never a production writer.
use anyhow::{ensure, Result};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug)]
pub struct Point {
    pub sample_idx: u64,
    pub range_km: f64,
    pub radial_velocity_km_s: f64,
    pub snr_linear: f64,
}
#[derive(Clone, Debug)]
pub struct Config {
    pub sample_rate_hz: f64,
    pub range_tolerance_km: f64,
    pub velocity_tolerance_km_s: f64,
    pub seed_snr_linear: f64,
    pub grow_snr_linear: f64,
    pub min_seed_points: usize,
    pub min_points: usize,
    pub min_duration_s: f64,
    pub max_gap_s: f64,
    pub trials: usize,
    pub max_tracks: usize,
    pub max_points: usize,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            sample_rate_hz: 1e6,
            range_tolerance_km: 0.5,
            velocity_tolerance_km_s: 3.0,
            seed_snr_linear: 7.0,
            grow_snr_linear: 2.0,
            min_seed_points: 4,
            min_points: 6,
            min_duration_s: 0.06,
            max_gap_s: 0.15,
            trials: 512,
            max_tracks: 16,
            max_points: 4096,
        }
    }
}
#[derive(Debug)]
pub struct Track {
    /// Original input indices; exclusive membership and at most one per pulse.
    pub members: Vec<usize>,
    pub epoch_sample: u64,
    pub range_at_epoch_km: f64,
    pub radial_velocity_km_s: f64,
    pub radial_acceleration_km_s2: f64,
}
type Model = (f64, f64, f64);
fn dt(p: &Point, epoch: u64, c: &Config) -> f64 {
    (i128::from(p.sample_idx) - i128::from(epoch)) as f64 / c.sample_rate_hz
}
fn fit(points: &[Point], ids: &[usize], epoch: u64, c: &Config) -> Option<Model> {
    use nalgebra::{Matrix3, Vector3};
    let n = ids.len() as f64;
    let tm = ids.iter().map(|&i| dt(&points[i], epoch, c)).sum::<f64>() / n;
    let scale = ids
        .iter()
        .map(|&i| (dt(&points[i], epoch, c) - tm).abs())
        .fold(0., f64::max);
    if scale <= 0. {
        return None;
    }
    let mut normal = Matrix3::zeros();
    let mut rhs = Vector3::zeros();
    for &i in ids {
        let t = (dt(&points[i], epoch, c) - tm) / scale;
        let r = Vector3::new(1., t, 0.5 * t * t) / c.range_tolerance_km;
        let v = Vector3::new(0., 1. / scale, t / scale) / c.velocity_tolerance_km_s;
        normal += r * r.transpose() + v * v.transpose();
        rhs += r * (points[i].range_km / c.range_tolerance_km)
            + v * (points[i].radial_velocity_km_s / c.velocity_tolerance_km_s);
    }
    let x = normal.cholesky()?.solve(&rhs);
    let v = x[1] / scale;
    let a = x[2] / (scale * scale);
    let model = (x[0] - v * tm + 0.5 * a * tm * tm, v - a * tm, a);
    if [model.0, model.1, model.2].iter().all(|x| x.is_finite()) {
        Some(model)
    } else {
        None
    }
}
fn residual(p: &Point, epoch: u64, model: Model, c: &Config) -> (f64, f64) {
    let t = dt(p, epoch, c);
    (
        (p.range_km - model.0 - model.1 * t - 0.5 * model.2 * t * t) / c.range_tolerance_km,
        (p.radial_velocity_km_s - model.1 - model.2 * t) / c.velocity_tolerance_km_s,
    )
}
fn consensus(
    points: &[Point],
    available: &[usize],
    epoch: u64,
    model: Model,
    c: &Config,
) -> Vec<Vec<usize>> {
    let mut pulse: BTreeMap<u64, (usize, f64)> = BTreeMap::new();
    for &i in available {
        let (r, v) = residual(&points[i], epoch, model, c);
        if r.abs() > 1. || v.abs() > 1. {
            continue;
        }
        let score = r * r + v * v;
        let entry = pulse.entry(points[i].sample_idx).or_insert((i, score));
        if score < entry.1 {
            *entry = (i, score);
        }
    }
    let mut groups: Vec<Vec<usize>> = vec![];
    for (_, (i, _)) in pulse {
        let split = groups.last().and_then(|g| g.last()).is_none_or(|&j| {
            (points[i].sample_idx - points[j].sample_idx) as f64 / c.sample_rate_hz > c.max_gap_s
        });
        if split {
            groups.push(vec![]);
        }
        groups.last_mut().unwrap().push(i);
    }
    groups
}
fn qualifies(points: &[Point], ids: &[usize], c: &Config) -> bool {
    ids.len() >= c.min_points
        && ids
            .iter()
            .filter(|&&i| points[i].snr_linear >= c.seed_snr_linear)
            .count()
            >= c.min_seed_points
        && (points[*ids.last().unwrap()].sample_idx - points[ids[0]].sample_idx) as f64
            / c.sample_rate_hz
            >= c.min_duration_s
}

pub fn split(points: &[Point], c: &Config) -> Result<Vec<Track>> {
    for v in [
        c.sample_rate_hz,
        c.range_tolerance_km,
        c.velocity_tolerance_km_s,
        c.seed_snr_linear,
        c.grow_snr_linear,
        c.min_duration_s,
        c.max_gap_s,
    ] {
        ensure!(
            v.is_finite() && v > 0.,
            "invalid positive configuration value"
        );
    }
    ensure!(
        c.seed_snr_linear >= c.grow_snr_linear
            && c.min_seed_points >= 2
            && c.min_points >= c.min_seed_points
            && c.trials > 0
            && c.max_tracks > 0,
        "invalid consensus configuration"
    );
    ensure!(
        points.len() <= c.max_points,
        "window too crowded: split the input window, never truncate by Doppler"
    );
    let mut available: Vec<_> = (0..points.len())
        .filter(|&i| {
            let p = points[i];
            p.range_km.is_finite()
                && p.range_km > 0.
                && p.radial_velocity_km_s.is_finite()
                && p.snr_linear.is_finite()
                && p.snr_linear >= c.grow_snr_linear
        })
        .collect();
    available.sort_by_key(|&i| (points[i].sample_idx, i));
    let mut tracks = vec![];
    let mut random = 0x52414e534143u64;
    for _ in 0..c.max_tracks {
        let seeds: Vec<_> = available
            .iter()
            .copied()
            .filter(|&i| points[i].snr_linear >= c.seed_snr_linear)
            .collect();
        if seeds.len() < c.min_seed_points {
            break;
        }
        let epoch = points[available[0]].sample_idx;
        let mut best: Option<(Vec<usize>, Model, f64)> = None;
        for trial in 0..c.trials {
            random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
            let i = seeds[trial % seeds.len()];
            let j = seeds[(random >> 32) as usize % seeds.len()];
            if points[i].sample_idx == points[j].sample_idx {
                continue;
            }
            let Some(model) = fit(points, &[i, j], epoch, c) else {
                continue;
            };
            for mut ids in consensus(points, &available, epoch, model, c) {
                if !qualifies(points, &ids, c) {
                    continue;
                }
                // Local optimization never recruits other temporal components.
                let Some(mut model) = fit(points, &ids, epoch, c) else {
                    continue;
                };
                let mut stable = false;
                for _ in 0..8 {
                    let before = ids.len();
                    ids.retain(|&i| {
                        let (r, v) = residual(&points[i], epoch, model, c);
                        r.abs() <= 1. && v.abs() <= 1.
                    });
                    if !qualifies(points, &ids, c) {
                        break;
                    }
                    if ids.len() == before {
                        stable = true;
                        break;
                    }
                    let Some(refit) = fit(points, &ids, epoch, c) else {
                        break;
                    };
                    model = refit;
                }
                if !stable
                    || ids.windows(2).any(|p| {
                        (points[p[1]].sample_idx - points[p[0]].sample_idx) as f64
                            / c.sample_rate_hz
                            > c.max_gap_s
                    })
                {
                    continue;
                }
                let score = ids
                    .iter()
                    .map(|&i| {
                        let (r, v) = residual(&points[i], epoch, model, c);
                        r * r + v * v
                    })
                    .sum();
                if best.as_ref().is_none_or(|(b, _, s)| {
                    ids.len() > b.len() || (ids.len() == b.len() && score < *s)
                }) {
                    best = Some((ids, model, score));
                }
            }
        }
        let Some((ids, model, _)) = best else {
            break;
        };
        let selected: std::collections::BTreeSet<_> = ids.iter().copied().collect();
        let pulses: std::collections::BTreeSet<_> =
            ids.iter().map(|&i| points[i].sample_idx).collect();
        available.retain(|i| {
            let (r, v) = residual(&points[*i], epoch, model, c);
            !selected.contains(i)
                && !(pulses.contains(&points[*i].sample_idx) && r.abs() <= 1. && v.abs() <= 1.)
        });
        let track_epoch = points[ids[0]].sample_idx;
        let shift = dt(&points[ids[0]], epoch, c);
        tracks.push(Track {
            members: ids,
            epoch_sample: track_epoch,
            range_at_epoch_km: model.0 + model.1 * shift + 0.5 * model.2 * shift * shift,
            radial_velocity_km_s: model.1 + model.2 * shift,
            radial_acceleration_km_s2: model.2,
        });
    }
    Ok(tracks)
}
