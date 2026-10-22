//! Calibrated Capon screening of the selected extreme main peaks in XC2 only.
use anyhow::{ensure, Context, Result};
use clap::Parser;
use hdf5::{
    types::{CompoundField, CompoundType, TypeDescriptor},
    File, Group, H5Type,
};
use ndarray::{s, Array3};
use pansy_legacy_capon::{axes, covariance, validate_pairs, Capon, C};
use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
    process::Command,
    time::Duration,
};
#[derive(Parser)]
struct Args {
    #[arg(long)]
    input: PathBuf,
    #[arg(long)]
    output: PathBuf,
    #[arg(long)]
    reference: PathBuf,
    #[arg(long, required = true)]
    root: Vec<PathBuf>,
}
#[repr(C)]
#[derive(Clone, Copy, Default)]
struct Z {
    r: f32,
    i: f32,
}
unsafe impl H5Type for Z {
    fn type_descriptor() -> TypeDescriptor {
        TypeDescriptor::Compound(CompoundType {
            fields: vec![
                CompoundField::typed::<f32>("r", 0, 0),
                CompoundField::typed::<f32>("i", 4, 1),
            ],
            size: 8,
        })
    }
}
fn read<T: H5Type>(g: &Group, n: &str) -> Result<Vec<T>> {
    Ok(g.dataset(n)?.read_raw()?)
}
fn scalar(g: &Group, n: &str) -> Result<u64> {
    Ok(g.dataset(n)?.read_scalar()?)
}
fn optional(g: &Group, n: &str) -> Result<Option<Vec<f64>>> {
    if g.link_exists(n) {
        Ok(Some(read(g, n)?))
    } else {
        Ok(None)
    }
}
fn write<T: H5Type>(g: &Group, n: &str, v: &[T]) -> Result<()> {
    let b = g.new_dataset_builder().with_data(v);
    if v.len() > 256 {
        b.shuffle().deflate(4).create(n)?;
    } else {
        b.create(n)?;
    }
    Ok(())
}
fn is_xc2(p: &Path) -> bool {
    p.components().any(|c| c.as_os_str() == "xc2")
}
fn classify(theta: f64, boundary: bool) -> u8 {
    if boundary {
        4
    } else if theta <= 2. {
        1
    } else {
        2
    }
}
fn evaluate(
    source: &Group,
    v: f64,
    r: f64,
    imager: &Capon,
    phases: &[f64; 7],
) -> Result<Option<pansy_legacy_capon::Peak>> {
    let path = PathBuf::from(String::from_utf8(read(source, "source_file")?)?);
    ensure!(is_xc2(&path), "XC2 only: rejected non-XC2 source");
    let file = File::open(path)?;
    let g = file.group(&String::from_utf8(read(source, "source_group")?)?)?;
    let ds = g.dataset("xc_arr")?;
    let shape = ds.shape();
    ensure!(
        shape.len() == 4 && shape[0] == 28 && shape[1] == 5,
        "invalid XC shape"
    );
    let ax = axes(
        shape[2],
        shape[3],
        scalar(&g, "n_fft")? as usize,
        scalar(&g, "ipp")? as usize,
        1e6,
        scalar(&g, "f0")? as usize,
        scalar(&g, "f1")? as usize,
        scalar(&g, "r0")? as usize,
        scalar(&g, "r1")? as usize,
        scalar(&g, "rdec")? as usize,
        optional(&g, "fvec")?,
        optional(&g, "rvec")?,
    )?;
    let fi = (0..ax.frequency_hz.len())
        .min_by(|&a, &b| {
            (ax.frequency_hz[a] * 299792458. / 94e6 - v)
                .abs()
                .total_cmp(&(ax.frequency_hz[b] * 299792458. / 94e6 - v).abs())
        })
        .context("empty frequency axis")?;
    let ri = (0..ax.range_km.len())
        .min_by(|&a, &b| {
            (ax.range_km[a] - r)
                .abs()
                .total_cmp(&(ax.range_km[b] - r).abs())
        })
        .context("empty range axis")?;
    ensure!(
        (ax.frequency_hz[fi] * 299792458. / 94e6 - v).abs() < 1e-6
            && (ax.range_km[ri] - r).abs() < 1e-6,
        "main peak does not match native axes"
    );
    let p = g
        .dataset(if g.link_exists("ch_pairs") {
            "ch_pairs"
        } else {
            "pairs"
        })?
        .read_2d::<u64>()?;
    ensure!(p.shape() == [28, 2], "invalid pair shape");
    let pairs: Vec<_> = p
        .outer_iter()
        .map(|p| [p[0] as usize, p[1] as usize])
        .collect();
    validate_pairs(&pairs)?;
    let data: Array3<Z> = ds.read_slice(s![.., 0, fi..fi + 1, ri..ri + 1])?;
    let packed: Vec<_> = data
        .iter()
        .map(|z| C::new(z.r as f64, z.i as f64))
        .collect();
    Ok(covariance(&packed, &pairs, phases)
        .and_then(|r| imager.peak(r, 0.01))
        .ok())
}
fn process(
    date: &str,
    a: &Args,
    imager: &Capon,
    positions: &[[f64; 3]; 7],
    phases: &[f64; 7],
) -> Result<(u64, u64, u64)> {
    let dest = a.output.join(format!("{date}.capon-v2.h5"));
    if dest.exists() {
        let f = File::open(dest)?;
        return Ok((
            read::<u64>(&f, "tested_count")?[0],
            read::<u64>(&f, "detection_count")?[0],
            read::<u64>(&f, "input_errors")?[0],
        ));
    }
    let input = File::open(a.input.join(format!("{date}.full-doppler.h5")))?;
    let flags = read::<u8>(&input, "extreme_peak_flag")?;
    let v = read::<f64>(&input, "main_peak_velocity_mps")?;
    let ranges = read::<f64>(&input, "main_peak_range_km")?;
    let peak_db = read::<f32>(&input, "main_peak_total_noise_db")?;
    ensure!(
        flags.len() == v.len() && v.len() == ranges.len(),
        "main peak shape mismatch"
    );
    let mut decision = vec![0u8; flags.len()];
    let mut theta = vec![f64::NAN; flags.len()];
    let mut u = theta.clone();
    let mut w = theta.clone();
    let mut power = theta.clone();
    let mut loading = theta.clone();
    let mut errors = Vec::new();
    for i in 0..flags.len() {
        let source = input.group(&format!("sources/{i:06}"))?;
        let path = PathBuf::from(String::from_utf8(read(&source, "source_file")?)?);
        ensure!(
            is_xc2(&path) && a.root.iter().any(|p| path.starts_with(p)),
            "non-XC2 source in XC2 day"
        );
        if flags[i] == 0 {
            continue;
        }
        if peak_db[i] < 5. {
            decision[i] = 6;
            continue;
        }
        ensure!(
            v[i].abs() >= 30. && v[i].abs() <= 90. && ranges[i] >= 80. && ranges[i] <= 92.,
            "invalid extreme main peak"
        );
        match evaluate(&source, v[i], ranges[i], imager, phases) {
            Ok(Some(p)) => {
                decision[i] = classify(p.theta_deg, p.boundary);
                theta[i] = p.theta_deg;
                u[i] = p.u;
                w[i] = p.v;
                power[i] = p.power;
                loading[i] = p.loading;
            }
            Ok(None) => decision[i] = 3,
            Err(e) => {
                decision[i] = 5;
                errors.push(format!("mode {i}: {e:#}"));
            }
        }
    }
    let tested = decision.iter().filter(|&&f| (1..=5).contains(&f)).count() as u64;
    let detected = decision.iter().filter(|&&f| f == 1).count() as u64;
    let input_errors = read::<u64>(&input, "input_errors")?[0] + errors.len() as u64;
    let tmp = dest.with_extension("h5.partial");
    ensure!(!tmp.exists(), "partial output exists");
    let f = File::create_excl(&tmp)?;
    write(&f, "schema", b"pansy-xc2-main-peak-capon-v2")?;
    write(&f, "minimum_total_noise_db", &[5_f64])?;
    write(
        &f,
        "below_power_floor_count",
        &[decision.iter().filter(|&&f| f == 6).count() as u64],
    )?;
    write(&f, "date_utc", date.as_bytes())?;
    write(&f, "capon_decision", &decision)?;
    write(&f, "capon_theta_deg", &theta)?;
    write(&f, "capon_u", &u)?;
    write(&f, "capon_v", &w)?;
    write(&f, "capon_power", &power)?;
    write(&f, "diagonal_loading_power", &loading)?;
    write(&f, "tested_count", &[tested])?;
    write(&f, "detection_count", &[detected])?;
    write(&f, "input_errors", &[input_errors])?;
    write(&f, "phasecal_radians", phases)?;
    write(
        &f,
        "positions_m_xyz",
        &positions.iter().flatten().copied().collect::<Vec<_>>(),
    )?;
    write(&f,"definitions",b"XC2 only, zenith. Main peak in available +/-90 m/s; abs(v)>=30 m/s; total/noise >=5 dB (not excess SNR). Calibrated covariance Cij=<zi conj(zj)>, Rij=Cij exp(+i(phi_i-phi_j)); steering exp(-ik r.s), historical negative vertical cosine. Retain measured noise; load 0.01 trace(R)/7 I, no subtraction or eigenvalue clipping. 61x61 grid within 9 deg; detection theta<=2 deg, not boundary. Decisions: 0 not extreme,1 detection,2 off-zenith,3 invalid covariance,4 search boundary,5 source or IO error,6 below 5 dB total/noise. No complementary-code exclusion.")?;
    write(
        &f,
        "overview_source",
        a.input
            .join(format!("{date}.full-doppler.h5"))
            .to_string_lossy()
            .as_bytes(),
    )?;
    write(
        &f,
        "reference_source",
        a.reference.to_string_lossy().as_bytes(),
    )?;
    write(&f, "errors", errors.join("\n").as_bytes())?;
    f.flush()?;
    drop(f);
    fs::rename(tmp, dest)?;
    println!("CAPON_DAY {date} tested={tested} detected={detected} off_zenith={} invalid={} boundary={} errors={input_errors}",decision.iter().filter(|&&v|v==2).count(),decision.iter().filter(|&&v|v==3).count(),decision.iter().filter(|&&v|v==4).count());
    Ok((tested, detected, input_errors))
}
fn index(a: &Args, days: &[String], status: &[u8], tested: &[u64], detected: &[u64]) -> Result<()> {
    let tmp = a.output.join("capon-index.partial.h5");
    let f = File::create(&tmp)?;
    write(
        &f,
        "utc_day_yyyymmdd",
        &days
            .iter()
            .map(|d| d.replace('-', "").parse::<u32>().unwrap())
            .collect::<Vec<_>>(),
    )?;
    write(&f, "status", status)?;
    write(&f, "tested_count", tested)?;
    write(&f, "detection_count", detected)?;
    f.flush()?;
    drop(f);
    fs::rename(tmp, a.output.join("capon-index.h5"))?;
    Ok(())
}
fn main() -> Result<()> {
    let a = Args::parse();
    ensure!(
        a.root
            .iter()
            .all(|p| p.file_name().is_some_and(|n| n == "xc2")),
        "only XC2 roots allowed"
    );
    fs::create_dir_all(&a.output)?;
    let (positions, phases) = pansy_legacy_capon::calibration::load_reference(&a.reference, 0)?;
    let imager = Capon::new(&positions, 47e6, 9., 61)?;
    let mut unique = BTreeSet::new();
    for root in &a.root {
        for e in fs::read_dir(root)? {
            let e = e?;
            let n = e.file_name().to_string_lossy().into_owned();
            if n.len() < 11 || n.as_bytes()[10] != b'T' || !e.file_type()?.is_dir() {
                continue;
            }
            if chrono::NaiveDate::parse_from_str(&n[..10], "%Y-%m-%d").is_err() {
                continue;
            }
            if fs::read_dir(e.path())?.filter_map(|v| v.ok()).any(|v| {
                let n = v.file_name().to_string_lossy().into_owned();
                n.starts_with("xc@") && n.ends_with(".h5")
            }) {
                unique.insert(n[..10].to_owned());
            }
        }
    }
    let days: Vec<_> = unique.into_iter().collect();
    ensure!(!days.is_empty(), "no XC2 dates");
    let mut status = vec![0; days.len()];
    let mut tested = vec![0; days.len()];
    let mut detected = vec![0; days.len()];
    index(&a, &days, &status, &tested, &detected)?;
    loop {
        let source = File::open(a.input.join("overview-index.h5"))?;
        let source_days = read::<u32>(&source, "utc_day_yyyymmdd")?;
        let source_status = read::<u8>(&source, "status")?;
        let all_done = !source_status.contains(&0);
        drop(source);
        if all_done {
            let _ = Command::new("systemctl")
                .args([
                    "--user",
                    "set-property",
                    "--runtime",
                    "pansy-xc2-peak-capon.service",
                    "CPUQuota=100%",
                ])
                .status();
        }
        for (i, d) in days.iter().enumerate() {
            if status[i] != 0 {
                continue;
            }
            if !a.input.join(format!("{d}.full-doppler.h5")).exists() {
                let n = d.replace('-', "").parse::<u32>()?;
                if all_done
                    || source_days
                        .iter()
                        .zip(&source_status)
                        .any(|(&day, &s)| day == n && s == 2)
                {
                    // Recover the readable part of a failed overview day in
                    // this same one-core cgroup; retain explicit source errors.
                    let reader = std::env::current_exe()?.with_file_name("pansy-spectra-overviews");
                    let mut cmd = Command::new(reader);
                    for root in &a.root {
                        cmd.arg("--root").arg(root);
                    }
                    let result = cmd
                        .arg("--day")
                        .arg(d)
                        .arg("--output")
                        .arg(&a.input)
                        .status();
                    eprintln!("OVERVIEW_RECOVERY {d}: {result:?}");
                    if !a.input.join(format!("{d}.full-doppler.h5")).exists() {
                        status[i] = 2;
                        continue;
                    }
                } else {
                    continue;
                }
            }
            match process(d, &a, &imager, &positions, &phases) {
                Ok((t, n, e)) => {
                    tested[i] = t;
                    detected[i] = n;
                    status[i] = if e == 0 { 1 } else { 2 };
                }
                Err(e) => {
                    eprintln!("CAPON_DAY_ERROR {d}: {e:#}");
                    status[i] = 2;
                }
            }
            index(&a, &days, &status, &tested, &detected)?;
        }
        index(&a, &days, &status, &tested, &detected)?;
        println!(
            "CAPON_PROGRESS done={} pending={} failed={} detections={} days={}",
            status.iter().filter(|&&s| s == 1).count(),
            status.iter().filter(|&&s| s == 0).count(),
            status.iter().filter(|&&s| s == 2).count(),
            detected.iter().sum::<u64>(),
            days.len()
        );
        if !status.contains(&0) {
            break;
        }
        std::thread::sleep(Duration::from_secs(30));
    }
    ensure!(
        !status.contains(&2),
        "XC2 days with errors; inspect capon index"
    );
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn strict_xc2() {
        assert!(is_xc2(Path::new("/media/analysis/metadata/xc2/2026/x.h5")));
        assert!(!is_xc2(Path::new(
            "/media/archive/metadata/xc_rvd/2026/x.h5"
        )));
    }
    #[test]
    fn angular_screen() {
        assert_eq!(classify(1., false), 1);
        assert_eq!(classify(2., false), 1);
        assert_eq!(classify(2.1, false), 2);
        assert_eq!(classify(1., true), 4);
    }
}
