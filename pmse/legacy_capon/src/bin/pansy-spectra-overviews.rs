//! Sequential native-spectrum daily overviews and main-peak flags.
use anyhow::{ensure, Context, Result};
use clap::Parser;
use hdf5::{
    types::{CompoundField, CompoundType, TypeDescriptor},
    File, Group, H5Type,
};
use ndarray::{s, Array2};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
};
#[derive(Parser, Debug)]
struct Args {
    #[arg(long, required = true)]
    root: Vec<PathBuf>,
    #[arg(long)]
    output: PathBuf,
    #[arg(long, default_value = "2026-01-30")]
    first_day: String,
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
fn write<T: H5Type>(g: &Group, n: &str, v: &[T]) -> Result<()> {
    let b = g.new_dataset_builder().with_data(v);
    if v.len() > 256 {
        b.shuffle().deflate(4).create(n)?;
    } else {
        b.create(n)?;
    }
    Ok(())
}
fn optional(g: &Group, n: &str) -> Result<Option<Vec<f64>>> {
    if g.link_exists(n) {
        Ok(Some(read(g, n)?))
    } else {
        Ok(None)
    }
}
fn peak(v: &[f64], db: &[f32]) -> Result<usize> {
    ensure!(v.len() == db.len(), "peak shape mismatch");
    (0..v.len())
        .filter(|&i| v[i].abs() <= 90. && db[i].is_finite())
        .max_by(|&i, &j| {
            db[i]
                .total_cmp(&db[j])
                .then_with(|| v[j].abs().total_cmp(&v[i].abs()))
        })
        .context("no stored bins within 90 m/s")
}
struct Mode {
    start: u64,
    stop: u64,
    v: Vec<f64>,
    db: Vec<f32>,
    range_at_peak: Vec<f64>,
    noise: Vec<f64>,
    ranges: Vec<f64>,
    source: String,
    group: String,
}
fn process(g: &Group, path: &Path) -> Result<Mode> {
    let start = scalar(g, "i0")?;
    let stop = scalar(g, "i1")?.checked_add(1).context("stop overflow")?;
    ensure!(stop > start, "invalid mode bounds");
    let ds = g.dataset("xc_arr")?;
    let shape = ds.shape();
    ensure!(
        shape.len() == 4 && shape[0] == 28 && shape[1] == 5,
        "expected 28x5xfrequencyxrange"
    );
    let nf = shape[2];
    let nr = shape[3];
    let axes = pansy_legacy_capon::axes(
        nf,
        nr,
        scalar(g, "n_fft")? as usize,
        scalar(g, "ipp")? as usize,
        1e6,
        scalar(g, "f0")? as usize,
        scalar(g, "f1")? as usize,
        scalar(g, "r0")? as usize,
        scalar(g, "r1")? as usize,
        scalar(g, "rdec")? as usize,
        optional(g, "fvec")?,
        optional(g, "rvec")?,
    )?;
    let p = g
        .dataset(if g.link_exists("ch_pairs") {
            "ch_pairs"
        } else {
            "pairs"
        })?
        .read_2d::<u64>()?;
    ensure!(p.shape() == [28, 2], "pair shape");
    let p: Vec<_> = p
        .outer_iter()
        .map(|p| [p[0] as usize, p[1] as usize])
        .collect();
    let autos = pansy_legacy_capon::validate_pairs(&p)?;
    let mut power = vec![0.; nf * nr];
    for c in autos {
        let z: Array2<Z> = ds.read_slice(s![c, 0, .., ..])?;
        for (p, z) in power.iter_mut().zip(z.iter()) {
            ensure!(z.r.is_finite() && z.r >= 0., "invalid auto power");
            *p += z.r as f64;
        }
    }
    let maxf = axes.frequency_hz.iter().map(|v| v.abs()).fold(0., f64::max);
    let guard: Vec<_> = axes
        .frequency_hz
        .iter()
        .enumerate()
        .filter_map(|(f, v)| (v.abs() >= 0.75 * maxf).then_some(f))
        .collect();
    ensure!(guard.len() >= 8, "too few noise bins");
    let selected: Vec<_> = axes
        .range_km
        .iter()
        .enumerate()
        .filter_map(|(r, v)| (80. ..=92.).contains(v).then_some(r))
        .collect();
    ensure!(!selected.is_empty(), "80--92 km absent");
    let mut noise = Vec::new();
    let mut db = vec![f32::NEG_INFINITY; nf];
    let mut range_at_peak = vec![0.; nf];
    for &r in &selected {
        let mut scratch: Vec<_> = guard.iter().map(|&f| power[f * nr + r]).collect();
        scratch.sort_unstable_by(f64::total_cmp);
        let k = scratch.len() / 2;
        let n = if scratch.len() % 2 == 0 {
            (scratch[k - 1] + scratch[k]) / 2.
        } else {
            scratch[k]
        };
        ensure!(n > 0. && n.is_finite(), "invalid noise");
        noise.push(n);
        for f in 0..nf {
            let value = (10. * (power[f * nr + r] / n).log10()) as f32;
            if value > db[f] {
                db[f] = value;
                range_at_peak[f] = axes.range_km[r];
            }
        }
    }
    ensure!(db.iter().all(|v| v.is_finite()), "invalid spectrum");
    Ok(Mode {
        start,
        stop,
        v: axes
            .frequency_hz
            .iter()
            .map(|f| f * 299792458. / 94e6)
            .collect(),
        db,
        range_at_peak,
        noise,
        ranges: selected.iter().map(|&r| axes.range_km[r]).collect(),
        source: path.to_string_lossy().into_owned(),
        group: g.name(),
    })
}
fn index(root: &Path, days: &[String], status: &[u8], flags: &[u64]) -> Result<()> {
    let tmp = root.join("overview-index.partial.h5");
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
    write(&f, "flagged_integrations", flags)?;
    write(&f,"definitions",b"0 pending; 1 complete; 2 failed or explicit input errors. Main peak of zenith max total/noise over 80--92 km, searched in available bins within +/-90 m/s; flag abs velocity >=30 m/s. No SNR threshold, angular screening or code-ghost subtraction.")?;
    f.flush()?;
    drop(f);
    fs::rename(tmp, root.join("overview-index.h5"))?;
    Ok(())
}
fn day(root: &Path, date: &str, files: &[PathBuf]) -> Result<(u64, u64)> {
    let out = root.join(format!("{date}.full-doppler.h5"));
    if out.exists() {
        let f = File::open(out)?;
        return Ok((
            read::<u8>(&f, "extreme_peak_flag")?
                .iter()
                .filter(|&&v| v != 0)
                .count() as u64,
            read::<u64>(&f, "input_errors")?[0],
        ));
    }
    let mut modes = BTreeMap::new();
    let mut errors = Vec::new();
    for path in files {
        let file = match File::open(path) {
            Ok(f) => f,
            Err(e) => {
                errors.push(format!("{}: {e}", path.display()));
                continue;
            }
        };
        let mut groups = file.member_names()?;
        groups.sort_by_key(|n| n.parse::<u64>().unwrap_or(u64::MAX));
        for name in groups {
            if name.parse::<u64>().is_err() {
                continue;
            }
            let g = file.group(&name)?;
            if !g.link_exists("xc_arr") {
                continue;
            }
            let t = scalar(&g, "i0")?;
            if modes.contains_key(&t) {
                continue;
            }
            match process(&g, path) {
                Ok(m) => {
                    modes.insert(t, m);
                }
                Err(e) => {
                    eprintln!("BLOCK_ERROR {date} {} {name}: {e:#}", path.display());
                    errors.push(format!("{} {name}: {e:#}", path.display()));
                }
            }
        }
    }
    ensure!(!modes.is_empty(), "no readable modes");
    let modes: Vec<_> = modes.into_values().collect();
    let velocity = &modes[0].v;
    ensure!(
        modes.iter().all(|m| m.v == *velocity),
        "stored frequency axis changes within day; no data discarded or resampled"
    );
    let partial = out.with_extension("h5.partial");
    ensure!(!partial.exists(), "partial exists; inspect before retry");
    let f = File::create_excl(&partial)?;
    write(&f, "date_utc", date.as_bytes())?;
    write(&f, "schema", b"pansy-zenith-main-peak-overview-v1")?;
    write(&f, "velocity_mps", velocity)?;
    let day = modes[0].start / 86_400_000_000 * 86_400_000_000;
    write(
        &f,
        "mode_start_hours",
        &modes
            .iter()
            .map(|m| (m.start - day) as f64 / 3.6e9)
            .collect::<Vec<_>>(),
    )?;
    write(
        &f,
        "mode_stop_hours",
        &modes
            .iter()
            .map(|m| (m.stop - day) as f64 / 3.6e9)
            .collect::<Vec<_>>(),
    )?;
    write(
        &f,
        "peak_total_noise_db_mode_frequency",
        &modes
            .iter()
            .flat_map(|m| m.db.iter().copied())
            .collect::<Vec<_>>(),
    )?;
    let peaks = modes
        .iter()
        .map(|m| peak(&m.v, &m.db))
        .collect::<Result<Vec<_>>>()?;
    let pv: Vec<_> = modes.iter().zip(&peaks).map(|(m, &p)| m.v[p]).collect();
    let flags: Vec<_> = pv.iter().map(|v| u8::from(v.abs() >= 30.)).collect();
    let nflags = flags.iter().filter(|&&v| v != 0).count() as u64;
    write(&f, "main_peak_velocity_mps", &pv)?;
    write(&f, "extreme_peak_flag", &flags)?;
    write(
        &f,
        "main_peak_total_noise_db",
        &modes
            .iter()
            .zip(&peaks)
            .map(|(m, &p)| m.db[p])
            .collect::<Vec<_>>(),
    )?;
    write(
        &f,
        "main_peak_range_km",
        &modes
            .iter()
            .zip(&peaks)
            .map(|(m, &p)| m.range_at_peak[p])
            .collect::<Vec<_>>(),
    )?;
    write(&f, "input_errors", &[errors.len() as u64])?;
    write(&f,"definitions",b"ALL native stored Dopplers, zenith beam 0. Seven-channel summed auto power; max total/noise over 80--92 km. Noise per range: median outer 25 percent of stored frequency band, a proxy not a calibrated estimate. Main peak is argmax of displayed spectrum among AVAILABLE bins with abs(v)<=90 m/s; extreme flag abs(v)>=30 m/s. No minimum SNR or significance cut. Wing extent is not the main peak. Flags are mechanical inspection markers, not validated winds. xc2 may store less than +/-90 m/s. Source spectra unchanged.")?;
    let sources = f.create_group("sources")?;
    for (i, m) in modes.iter().enumerate() {
        let g = sources.create_group(&format!("{i:06}"))?;
        write(&g, "source_file", m.source.as_bytes())?;
        write(&g, "source_group", m.group.as_bytes())?;
        write(&g, "noise_range_km", &m.ranges)?;
        write(&g, "estimated_total_noise_power", &m.noise)?;
    }
    write(&f, "error_details", errors.join("\n").as_bytes())?;
    f.flush()?;
    drop(sources);
    drop(f);
    fs::rename(partial, &out)?;
    println!(
        "DAY_COMPLETE {date} modes={} flagged={nflags} errors={} bytes={}",
        modes.len(),
        errors.len(),
        fs::metadata(out)?.len()
    );
    Ok((nflags, errors.len() as u64))
}
fn main() -> Result<()> {
    let a = Args::parse();
    fs::create_dir_all(&a.output)?;
    let mut dates: BTreeMap<String, Vec<PathBuf>> = BTreeMap::new();
    for root in &a.root {
        for e in fs::read_dir(root)? {
            let e = e?;
            let name = e.file_name().to_string_lossy().into_owned();
            if name.len() < 11 || name.as_bytes()[10] != b'T' || !e.file_type()?.is_dir() {
                continue;
            }
            let date = &name[..10];
            if chrono::NaiveDate::parse_from_str(date, "%Y-%m-%d").is_err() {
                continue;
            }
            for e in fs::read_dir(e.path())? {
                let p = e?.path();
                let n = p.file_name().unwrap().to_string_lossy();
                if n.starts_with("xc@") && n.ends_with(".h5") {
                    dates.entry(date.to_string()).or_default().push(p);
                }
            }
        }
    }
    let mut days: Vec<_> = dates.keys().cloned().collect();
    if let Some(i) = days.iter().position(|d| d == &a.first_day) {
        let d = days.remove(i);
        days.insert(0, d);
    }
    ensure!(!days.is_empty(), "no populated days");
    let mut status = vec![0; days.len()];
    let mut flags = vec![0; days.len()];
    index(&a.output, &days, &status, &flags)?;
    for (i, d) in days.iter().enumerate() {
        let mut files = dates.remove(d).unwrap();
        let mut seen = BTreeSet::new();
        files.retain(|p| seen.insert(p.clone()));
        files.sort_by_key(|p| p.file_name().unwrap().to_os_string());
        println!("DAY_START {}/{} {d}", i + 1, days.len());
        match day(&a.output, d, &files) {
            Ok((n, e)) => {
                flags[i] = n;
                status[i] = if e == 0 { 1 } else { 2 };
            }
            Err(e) => {
                eprintln!("DAY_ERROR {d}: {e:#}");
                status[i] = 2;
            }
        }
        index(&a.output, &days, &status, &flags)?;
        println!(
            "PROGRESS completed={} failed={} total={}",
            status.iter().filter(|&&v| v == 1).count(),
            status.iter().filter(|&&v| v == 2).count(),
            days.len()
        );
    }
    ensure!(
        !status.contains(&2),
        "explicit failed days; inspect index and journal"
    );
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn main_peak_not_wing() {
        let v = [-100., -40., 0., 40., 100.];
        assert_eq!(peak(&v, &[50., 10., 20., 15., 60.]).unwrap(), 2);
        assert_eq!(peak(&v, &[50., 25., 20., 15., 60.]).unwrap(), 1);
    }
    #[test]
    fn inclusive_limits() {
        let v = [-90., -30., 0., 30., 90.];
        assert_eq!(peak(&v, &[1., 2., 3., 4., 5.]).unwrap(), 4);
        assert!((v[1] as f64).abs() >= 30.);
    }
    #[test]
    fn ties_nearest_zero() {
        assert_eq!(peak(&[-40., 0., 40.], &[2., 2., 2.]).unwrap(), 1);
    }
}
