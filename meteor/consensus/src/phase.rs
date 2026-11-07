//! Local interferometric compatibility followed by robust range/Doppler fits.
//! Relative receiver phases are circular; no phase unwrapping or alias selection.
use crate::{split, Config, Point, Track};
use anyhow::{ensure, Result};
use std::collections::BTreeMap;

pub fn split_phase(
    points: &[Point],
    phase: &[[f32; 14]],
    coherence_min: f64,
    lag_s: f64,
    config: &Config,
) -> Result<Vec<Track>> {
    ensure!(
        points.len() == phase.len() && points.len() <= 16384,
        "phase shape/window limit"
    );
    ensure!(
        coherence_min > 0. && coherence_min <= 1. && lag_s > 0. && lag_s <= 0.15,
        "invalid phase settings"
    );
    let mut unit = Vec::with_capacity(phase.len());
    for row in phase {
        let mut normalized = [0f64; 14];
        for ch in 0..7 {
            let re = row[2 * ch] as f64;
            let im = row[2 * ch + 1] as f64;
            let amp = re.hypot(im);
            ensure!(amp.is_finite() && amp > 0., "missing receiver phase");
            normalized[2 * ch] = re / amp;
            normalized[2 * ch + 1] = im / amp;
        }
        unit.push(normalized);
    }
    let mut order: Vec<_> = (0..points.len()).collect();
    order.sort_by_key(|&i| points[i].sample_idx);
    let mut parent: Vec<_> = (0..points.len()).collect();
    let mut members: Vec<Vec<usize>> = (0..points.len()).map(|i| vec![i]).collect();
    let coherence = |i: usize, j: usize| {
        let mut re = 0.;
        let mut im = 0.;
        for ch in 0..7 {
            re += unit[i][2 * ch] * unit[j][2 * ch] + unit[i][2 * ch + 1] * unit[j][2 * ch + 1];
            im += unit[i][2 * ch + 1] * unit[j][2 * ch] - unit[i][2 * ch] * unit[j][2 * ch + 1];
        }
        (re * re + im * im) / 49.
    };
    let mut edges = Vec::new();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for (pos, &i) in order.iter().enumerate() {
        for &j in order[..pos].iter().rev() {
            let dt = (points[i].sample_idx - points[j].sample_idx) as f64 / config.sample_rate_hz;
            if dt > lag_s {
                break;
            }
            if dt == 0. {
                continue;
            }
            let dr = points[i].range_km
                - points[j].range_km
                - 0.5 * (points[i].radial_velocity_km_s + points[j].radial_velocity_km_s) * dt;
            if dr.abs() > config.range_tolerance_km * 2.
                || (points[i].radial_velocity_km_s - points[j].radial_velocity_km_s).abs()
                    > config.velocity_tolerance_km_s * 2.
            {
                continue;
            }
            // Common meteor/TX phase cancels; only inter-receiver structure remains.
            let score = coherence(i, j);
            if score >= coherence_min {
                edges.push((score, i, j));
            }
        }
    }
    edges.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    for (_, i, j) in edges {
        let mut a = root(&mut parent, i);
        let mut b = root(&mut parent, j);
        if a == b {
            continue;
        }
        if members[a].len() > members[b].len() {
            std::mem::swap(&mut a, &mut b);
        }
        // An isolated intermediate-phase point must not bridge two simultaneous
        // tracks. Compare nearby observations across the overlapping components.
        if members[a].len() >= 4 {
            let mut compared = 0;
            let mut good = 0;
            for &k in members[a].iter().step_by((members[a].len() / 64).max(1)) {
                let t = points[k].sample_idx;
                let pos = members[b].partition_point(|&l| points[l].sample_idx < t);
                let nearest = [pos.saturating_sub(1), pos.min(members[b].len() - 1)]
                    .into_iter()
                    .map(|p| members[b][p])
                    .min_by_key(|&l| points[l].sample_idx.abs_diff(t))
                    .unwrap();
                if points[nearest].sample_idx.abs_diff(t) as f64 / config.sample_rate_hz <= lag_s {
                    compared += 1;
                    if coherence(k, nearest) >= coherence_min {
                        good += 1;
                    }
                }
            }
            if compared >= 3 && (good as f64) / (compared as f64) < 0.6 {
                continue;
            }
        }
        parent[a] = b;
        let moved = std::mem::take(&mut members[a]);
        members[b].extend(moved);
        members[b].sort_by_key(|&k| points[k].sample_idx);
    }
    let mut components: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &i in &order {
        components.entry(root(&mut parent, i)).or_default().push(i);
    }
    let mut tracks = Vec::new();
    for ids in components.values() {
        if ids.len() < config.min_points {
            continue;
        }
        let subset: Vec<_> = ids.iter().map(|&i| points[i]).collect();
        for mut track in split(&subset, config)? {
            track.members = track.members.iter().map(|&i| ids[i]).collect();
            tracks.push(track);
        }
    }
    tracks.sort_by_key(|t| t.epoch_sample);
    Ok(tracks)
}

/// # Safety
/// Each input scalar buffer has n entries, phase has n*14 f32 entries and
/// labels has n writable i32 entries. No buffers overlap. No pointers are null.
#[no_mangle]
pub unsafe extern "C" fn pansy_consensus_phase(
    samples: *const u64,
    ranges: *const f32,
    velocities: *const f32,
    snr: *const f32,
    phase: *const f32,
    n: usize,
    coherence: f64,
    lag: f64,
    labels: *mut i32,
) -> i32 {
    if samples.is_null()
        || ranges.is_null()
        || velocities.is_null()
        || snr.is_null()
        || phase.is_null()
        || labels.is_null()
        || n > 16384
    {
        return -1;
    }
    std::panic::catch_unwind(|| {
        let s = std::slice::from_raw_parts(samples, n);
        let r = std::slice::from_raw_parts(ranges, n);
        let v = std::slice::from_raw_parts(velocities, n);
        let q = std::slice::from_raw_parts(snr, n);
        let p = std::slice::from_raw_parts(phase.cast::<[f32; 14]>(), n);
        let out = std::slice::from_raw_parts_mut(labels, n);
        out.fill(-1);
        let points: Vec<_> = (0..n)
            .map(|i| Point {
                sample_idx: s[i],
                range_km: r[i] as f64,
                radial_velocity_km_s: v[i] as f64 / 1000.,
                snr_linear: q[i] as f64,
            })
            .collect();
        let config = Config {
            max_points: 16384,
            ..Config::default()
        };
        match split_phase(&points, p, coherence, lag, &config) {
            Ok(tracks) => {
                for (id, t) in tracks.iter().enumerate() {
                    for &i in &t.members {
                        out[i] = id as i32;
                    }
                }
                tracks.len() as i32
            }
            Err(_) => -2,
        }
    })
    .unwrap_or(-3)
}
