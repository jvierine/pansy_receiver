//! Reusable phase-preserving meteor DSP and compact-sphere drag/ablation dynamics.
use rustfft::num_complex::Complex64;
use std::f64::consts::PI;

#[cfg(feature = "plotting")]
pub mod plotting;

fn sample(rx: &[Complex64], x: f64) -> Complex64 {
    let nearest = x.round();
    if (x - nearest).abs() < 1e-12 {
        return rx.get(nearest as usize).copied().unwrap_or_default();
    }
    let mut sum = Complex64::default();
    for i in x.floor() as isize - 16..=x.floor() as isize + 17 {
        if i < 0 || i >= rx.len() as isize {
            continue;
        }
        let d = x - i as f64;
        if d.abs() >= 17. {
            continue;
        }
        let sinc = (PI * d).sin() / (PI * d);
        let window = 0.42 + 0.5 * (PI * d / 17.).cos() + 0.08 * (2. * PI * d / 17.).cos();
        sum += rx[i as usize] * (sinc * window);
    }
    sum
}

/// Phases refer to the energy-weighted transmit fast-time centre, not gate zero.
pub fn matched(rx: &[Complex64], tx: &[Complex64], start: f64, fd: f64, fs: f64) -> Complex64 {
    let energy: f64 = tx.iter().map(|v| v.norm_sqr()).sum();
    let centre = tx
        .iter()
        .enumerate()
        .map(|(i, v)| i as f64 * v.norm_sqr())
        .sum::<f64>()
        / energy;
    tx.iter()
        .enumerate()
        .map(|(i, t)| {
            sample(rx, start + i as f64) * t.conj() / energy.sqrt()
                * Complex64::from_polar(1., -2. * PI * fd * (i as f64 - centre) / fs)
        })
        .sum()
}

/// # Safety
/// RX/TX interleaved arrays and output must span their declared dimensions.
#[no_mangle]
pub unsafe extern "C" fn meteor_matched(
    rx: *const f64,
    channels: usize,
    gates: usize,
    tx: *const f64,
    len: usize,
    start: f64,
    fd: f64,
    fs: f64,
    out: *mut f64,
) -> i32 {
    if rx.is_null()
        || tx.is_null()
        || out.is_null()
        || channels == 0
        || len == 0
        || fs <= 0.
        || !start.is_finite()
        || start < 17.
        || start + len as f64 + 17. >= gates as f64
    {
        return -1;
    }
    let raw = std::slice::from_raw_parts(rx, channels * gates * 2);
    let t = std::slice::from_raw_parts(tx, len * 2)
        .as_chunks::<2>()
        .0
        .iter()
        .map(|x| Complex64::new(x[0], x[1]))
        .collect::<Vec<_>>();
    if t.iter().map(|v| v.norm_sqr()).sum::<f64>() <= 0. {
        return -2;
    }
    let result = std::slice::from_raw_parts_mut(out, channels * 2);
    for ch in 0..channels {
        let z = raw[ch * gates * 2..(ch + 1) * gates * 2]
            .as_chunks::<2>()
            .0
            .iter()
            .map(|x| Complex64::new(x[0], x[1]))
            .collect::<Vec<_>>();
        let a = matched(&z, &t, start, fd, fs);
        result[2 * ch] = a.re;
        result[2 * ch + 1] = a.im;
    }
    0
}

fn density(h: f64, alt: &[f64], rho: &[f64]) -> f64 {
    if h < alt[0] || h > alt[alt.len() - 1] {
        return f64::NAN;
    }
    let j = alt.partition_point(|&a| a < h).max(1).min(alt.len() - 1);
    let f = (h - alt[j - 1]) / (alt[j] - alt[j - 1]);
    ((1. - f) * rho[j - 1].ln() + f * rho[j].ln()).exp()
}

fn rhs(
    y: &[f64; 7],
    alt: &[f64],
    rho: &[f64],
    station_alt: f64,
    material: f64,
    sigma: f64,
    gravity: bool,
) -> [f64; 7] {
    let h = station_alt + y[2] + (y[0] * y[0] + y[1] * y[1]) / (2. * 6_371_000.);
    let air = density(h, alt, rho);
    let speed = y[3..6].iter().map(|v| v * v).sum::<f64>().sqrt();
    let a = -0.75 * air / (material * y[6]) * speed;
    [
        y[3],
        y[4],
        y[5],
        a * y[3],
        a * y[4],
        a * y[5]
            - if gravity {
                9.80665 * (6_371_000. / (6_371_000. + h)).powi(2)
            } else {
                0.
            },
        -air * sigma * speed.powi(3) / (4. * material),
    ]
}

#[allow(clippy::too_many_arguments)]
pub fn propagate(
    initial: [f64; 7],
    times: &[f64],
    alt: &[f64],
    rho: &[f64],
    station_alt: f64,
    material: f64,
    sigma: f64,
    dt_max: f64,
    gravity: bool,
) -> Option<Vec<[f64; 7]>> {
    if times.iter().any(|v| !v.is_finite() || *v < 0.)
        || !times.windows(2).all(|v| v[1] >= v[0])
        || alt.len() < 2
        || alt.len() != rho.len()
        || dt_max <= 0.
        || initial[6] <= 0.
    {
        return None;
    }
    let mut y = initial;
    let mut t = 0.;
    let mut result = Vec::with_capacity(times.len());
    for &target in times {
        if target < t || !target.is_finite() {
            return None;
        }
        while target - t > 1e-12 {
            let dt = (target - t).min(dt_max);
            let k1 = rhs(&y, alt, rho, station_alt, material, sigma, gravity);
            let add = |k: &[f64; 7], f: f64| std::array::from_fn(|j| y[j] + dt * f * k[j]);
            let k2 = rhs(
                &add(&k1, 0.5),
                alt,
                rho,
                station_alt,
                material,
                sigma,
                gravity,
            );
            let k3 = rhs(
                &add(&k2, 0.5),
                alt,
                rho,
                station_alt,
                material,
                sigma,
                gravity,
            );
            let k4 = rhs(
                &add(&k3, 1.),
                alt,
                rho,
                station_alt,
                material,
                sigma,
                gravity,
            );
            for j in 0..7 {
                y[j] += dt * (k1[j] + 2. * k2[j] + 2. * k3[j] + k4[j]) / 6.;
            }
            if y[6] <= 1e-9 || !y.iter().all(|v| v.is_finite()) {
                return None;
            }
            t += dt;
        }
        result.push(y);
    }
    Some(result)
}

/// # Safety
/// Arrays must span their declared lengths. Output contains n rows of seven f64s.
#[no_mangle]
pub unsafe extern "C" fn meteor_propagate(
    initial: *const f64,
    times: *const f64,
    n: usize,
    alt: *const f64,
    rho: *const f64,
    nh: usize,
    station_alt: f64,
    material: f64,
    sigma: f64,
    dt: f64,
    gravity: i32,
    out: *mut f64,
) -> i32 {
    if initial.is_null()
        || times.is_null()
        || alt.is_null()
        || rho.is_null()
        || out.is_null()
        || nh < 2
        || dt <= 0.
        || material <= 0.
        || sigma < 0.
    {
        return -1;
    }
    let y: [f64; 7] = std::slice::from_raw_parts(initial, 7).try_into().unwrap();
    let a = std::slice::from_raw_parts(alt, nh);
    let r = std::slice::from_raw_parts(rho, nh);
    if !a.windows(2).all(|v| v[1] > v[0]) || !r.iter().all(|v| v.is_finite() && *v > 0.) {
        return -2;
    }
    let Some(result) = propagate(
        y,
        std::slice::from_raw_parts(times, n),
        a,
        r,
        station_alt,
        material,
        sigma,
        dt,
        gravity != 0,
    ) else {
        return -3;
    };
    let output = std::slice::from_raw_parts_mut(out, n * 7);
    for (i, row) in result.iter().enumerate() {
        output[i * 7..(i + 1) * 7].copy_from_slice(row);
    }
    0
}

/// # Safety
/// xy/membership span elements*2/elements*channels; uv and output span their shapes.
#[no_mangle]
pub unsafe extern "C" fn meteor_array_response(
    xy: *const f64,
    mask: *const i64,
    active: *const u8,
    elements: usize,
    channels: usize,
    uv: *const f64,
    n: usize,
    wavelength: f64,
    out: *mut f64,
) -> i32 {
    if xy.is_null()
        || mask.is_null()
        || active.is_null()
        || uv.is_null()
        || out.is_null()
        || wavelength <= 0.
    {
        return -1;
    }
    let xy = std::slice::from_raw_parts(xy, elements * 2);
    let mask = std::slice::from_raw_parts(mask, elements * channels);
    let active = std::slice::from_raw_parts(active, elements);
    let uv = std::slice::from_raw_parts(uv, n * 2);
    let result = std::slice::from_raw_parts_mut(out, n * channels * 2);
    result.fill(0.);
    for p in 0..n {
        for i in 0..elements {
            if active[i] == 0 {
                continue;
            }
            let phase =
                2. * PI / wavelength * (xy[2 * i] * uv[2 * p] + xy[2 * i + 1] * uv[2 * p + 1]);
            let (im, re) = phase.sin_cos();
            for ch in 0..channels {
                if mask[i * channels + ch] != 0 {
                    let j = (p * channels + ch) * 2;
                    result[j] += re;
                    result[j + 1] += im;
                }
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fractional_tone_and_centre_phase() {
        let fd = 23000.;
        let fs = 500000.;
        let start = 40.35;
        let rx = (0..200)
            .map(|i| Complex64::from_polar(1., 2. * PI * fd * i as f64 / fs))
            .collect::<Vec<_>>();
        let tx = vec![Complex64::new(1., 0.); 32];
        let got = matched(&rx, &tx, start, fd, fs);
        let truth = Complex64::from_polar(32f64.sqrt(), 2. * PI * fd * (start + 15.5) / fs);
        assert!((got - truth).norm() < 1e-4);
    }
    #[test]
    fn no_drag_and_shrinking_radius() {
        let initial = [0., 0., 100000., 1000., 0., -50000., 0.001];
        let result = propagate(
            initial,
            &[0., 0.1],
            &[0., 200000.],
            &[1e-30, 1e-30],
            4287.,
            3000.,
            1e-8,
            0.0005,
            false,
        )
        .unwrap();
        assert!((result[1][0] - 100.).abs() < 1e-8);
        assert!((result[1][2] - 95000.).abs() < 1e-8);
        let drag = propagate(
            initial,
            &[0., 0.1],
            &[0., 200000.],
            &[1e-7, 1e-7],
            4287.,
            3000.,
            1e-8,
            0.0005,
            false,
        )
        .unwrap();
        assert!(drag[1][6] < initial[6] && drag[1][5] > initial[5]);
    }
}
