//! Conditional phase-only imaging with full element sums, not centroid baselines.
use anyhow::{ensure, Result};
use hdf5::{File, Group, H5Type};
use rustfft::num_complex::Complex32;
use std::path::Path;

fn put<T: H5Type>(g: &Group, name: &str, shape: &[usize], v: &[T]) -> Result<()> {
    g.new_dataset::<T>()
        .shape(shape)
        .create(name)?
        .write_raw(v)?;
    Ok(())
}
struct Model {
    uv: [f32; 2],
    phases: [Complex32; 10],
}
fn models(xy: &[f64], mask: &[i64], active: &[u8], n: usize, width_deg: f32) -> Vec<Model> {
    let limit = width_deg.to_radians().sin();
    let mut result = Vec::with_capacity(n * n);
    for y in 0..n {
        for x in 0..n {
            let uv = [
                -limit + 2. * limit * x as f32 / (n - 1) as f32,
                -limit + 2. * limit * y as f32 / (n - 1) as f32,
            ];
            let mut response = [Complex32::default(); 10];
            for i in 0..active.len() {
                if active[i] == 0 {
                    continue;
                }
                let phase = 2. * std::f64::consts::PI / 5.99584916
                    * (xy[2 * i] * uv[0] as f64 + xy[2 * i + 1] * uv[1] as f64);
                let z = Complex32::from_polar(1., phase as f32);
                for j in 0..10 {
                    if mask[10 * i + j] != 0 {
                        response[j] += z;
                    }
                }
            }
            for z in &mut response {
                *z /= z.norm().max(1e-20);
            }
            result.push(Model {
                uv,
                phases: response,
            });
        }
    }
    result
}
fn score(model: &Model, voltage: &[Complex32; 10], weights: &[f32; 10]) -> f32 {
    let mut sum = Complex32::default();
    let total: f32 = weights.iter().sum();
    for j in 0..10 {
        sum += voltage[j] * model.phases[j].conj() * weights[j];
    }
    sum.norm_sqr() / (total * total).max(1e-20)
}
fn measurement(iq: &[f32], noise: &[f32], phase: &[f32]) -> ([Complex32; 10], [f32; 10]) {
    let mut z = [Complex32::default(); 10];
    let mut weights = [0f32; 10];
    for j in 0..10 {
        let channel = j + 5;
        let a = Complex32::new(iq[2 * channel], iq[2 * channel + 1]);
        z[j] = a / a.norm().max(1e-20) * Complex32::from_polar(1., -phase[channel]);
        weights[j] = (a.norm_sqr() / noise[channel]).sqrt().min(10.);
    }
    (z, weights)
}
pub fn run(root: &Path, geometry: &Path) -> Result<()> {
    let path = root.join("interferometry_provisional.h5");
    ensure!(!path.exists(), "imaging output exists");
    let cal = File::open(root.join("phasecal_boresight_provisional.h5"))?;
    ensure!(
        cal.dataset("selected_events")?.read_scalar::<u64>()? >= 3,
        "need at least three calibration events"
    );
    let phase = cal.dataset("phasecal_rad")?.read_raw::<f32>()?;
    ensure!(
        phase[5..].iter().all(|x| x.is_finite()),
        "nonfinite calibration"
    );
    let g = File::open(geometry)?;
    let xy = g.dataset("xy_m")?.read_raw::<f64>()?;
    let membership = g.dataset("membership")?.read_raw::<i64>()?;
    let active = g.dataset("provisional_active_mask")?.read_raw::<u8>()?;
    ensure!(
        xy.len() == 2 * active.len() && membership.len() == 10 * active.len(),
        "geometry shape mismatch"
    );
    let wide_n = 201;
    let fine_n = 141;
    let wide = models(&xy, &membership, &active, wide_n, 20.);
    let fine = models(&xy, &membership, &active, fine_n, 4.);
    let catalog = File::open(root.join("catalogue.h5"))?;
    let f = File::create(path)?;
    let status="PROVISIONAL: boresight-centred ensemble phase gauge; unconfirmed active mask and receive phase sign; phase-only score, not calibrated likelihood";
    put(&f, "status", &[status.len()], status.as_bytes())?;
    put(
        &f,
        "wide_uv",
        &[wide_n, wide_n, 2],
        &wide.iter().flat_map(|m| m.uv).collect::<Vec<_>>(),
    )?;
    put(
        &f,
        "fine_uv",
        &[fine_n, fine_n, 2],
        &fine.iter().flat_map(|m| m.uv).collect::<Vec<_>>(),
    )?;
    let events = f.create_group("events")?;
    for name in catalog.group("events")?.member_names()? {
        let event = catalog.group(&format!("events/{name}"))?;
        let iq = event
            .dataset("matched_amplitude_re_im")?
            .read_raw::<f32>()?;
        let noise = event.dataset("noise_power")?.read_raw::<f32>()?;
        let statistic = event
            .dataset("detector_excess_over_median_db")?
            .read_raw::<f32>()?;
        let peak = (0..statistic.len())
            .max_by(|&a, &b| statistic[a].total_cmp(&statistic[b]))
            .unwrap();
        let (z, w) = measurement(
            &iq[peak * 30..(peak + 1) * 30],
            &noise[peak * 15..(peak + 1) * 15],
            &phase,
        );
        let scores = wide.iter().map(|m| score(m, &z, &w)).collect::<Vec<_>>();
        let out = events.create_group(&name)?;
        put(&out, "peak_phase_score_map", &[wide_n, wide_n], &scores)?;
        let mut uv = Vec::new();
        let mut values = Vec::new();
        let mut broad_uv = Vec::new();
        for i in 0..statistic.len() {
            let (z, w) = measurement(
                &iq[i * 30..(i + 1) * 30],
                &noise[i * 15..(i + 1) * 15],
                &phase,
            );
            let best = fine
                .iter()
                .max_by(|a, b| score(a, &z, &w).total_cmp(&score(b, &z, &w)))
                .unwrap();
            uv.extend(best.uv);
            values.push(score(best, &z, &w));
        }
        let global = wide
            .iter()
            .max_by(|a, b| score(a, &z, &w).total_cmp(&score(b, &z, &w)))
            .unwrap();
        broad_uv.extend(global.uv);
        put(
            &out,
            "conditional_uv_within_4deg_box",
            &[statistic.len(), 2],
            &uv,
        )?;
        put(&out, "phase_score", &[statistic.len()], &values)?;
        put(&out, "global_peak_uv_20deg_box", &[2], &broad_uv)?;
        put(
            &out,
            "tx_time_us",
            &[statistic.len()],
            &event.dataset("tx_time_us")?.read_raw::<u64>()?,
        )?;
    }
    println!("Native provisional interferometry written for all tracks; 10 same-pointing subarray channels, full element-sum phase model.");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn recovers_injected_direction_after_gain_phase_correction() {
        let xy = (0..10)
            .flat_map(|j| [j as f64 * 7.3, (j * j % 7) as f64 * 8.1])
            .collect::<Vec<_>>();
        let mut membership = vec![0i64; 100];
        for j in 0..10 {
            membership[j * 10 + j] = 1;
        }
        let grid = models(&xy, &membership, &[1; 10], 5, 4.);
        let uv = [4f32.to_radians().sin() / 2., -4f32.to_radians().sin() / 2.];
        let mut iq = vec![0f32; 30];
        let mut cal = vec![0f32; 15];
        for j in 0..10 {
            cal[j + 5] = 0.1 * j as f32;
            let spatial = 2. * std::f64::consts::PI / 5.99584916
                * (xy[2 * j] * uv[0] as f64 + xy[2 * j + 1] * uv[1] as f64);
            let z = Complex32::from_polar(20., spatial as f32 + cal[j + 5]);
            iq[2 * (j + 5)] = z.re;
            iq[2 * (j + 5) + 1] = z.im;
        }
        let (z, w) = measurement(&iq, &[1.; 15], &cal);
        let best = grid
            .iter()
            .max_by(|a, b| score(a, &z, &w).total_cmp(&score(b, &z, &w)))
            .unwrap();
        assert!((best.uv[0] - uv[0]).abs() < 1e-6 && (best.uv[1] - uv[1]).abs() < 1e-6);
        assert!((score(best, &z, &w) - 1.).abs() < 1e-5);
    }
}
