//! Bin-by-bin XC2 implementation of the supplied extreme-filter procedure.
use anyhow::{ensure, Result};
use clap::Parser;
use hdf5::{
    types::{CompoundField, CompoundType, TypeDescriptor},
    File, Group, H5Type,
};
use ndarray::{s, Array3};
use pansy_legacy_capon::{axes, koki, validate_pairs, Capon, C};
use std::{
    fs,
    path::{Path, PathBuf},
    time::Instant,
};
#[derive(Parser)]
struct Args {
    #[arg(long)]
    input: PathBuf,
    #[arg(long)]
    output: PathBuf,
    #[arg(long)]
    reference: PathBuf,
    /// Discover and prepare XC2 days directly, without mixed-product overviews.
    #[arg(long)]
    root: Vec<PathBuf>,
    #[arg(long)]
    pmse_season: bool,
    #[arg(long)]
    day: Option<String>,
    #[arg(long, default_value_t = 0)]
    max_modes: usize,
    /// Aggregate bytes for this entire sweep, including partial files.
    #[arg(long, default_value_t = 400_000_000)]
    max_output_bytes: u64,
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
fn scalar(g: &Group, n: &str) -> Result<usize> {
    Ok(g.dataset(n)?.read_scalar::<u64>()? as usize)
}
fn opt(g: &Group, n: &str) -> Result<Option<Vec<f64>>> {
    if g.link_exists(n) {
        Ok(Some(read(g, n)?))
    } else {
        Ok(None)
    }
}
fn write<T: H5Type>(g: &Group, n: &str, v: &[T]) -> Result<()> {
    if v.len() > 256 {
        g.new_dataset_builder()
            .with_data(v)
            .shuffle()
            .deflate(4)
            .create(n)?;
    } else {
        g.new_dataset_builder().with_data(v).create(n)?;
    }
    Ok(())
}
fn text(g: &Group, n: &str, v: &str) -> Result<()> {
    write(g, n, v.as_bytes())
}
fn output_bytes(p:&Path)->Result<u64> {
    let mut total=0u64;
    for e in fs::read_dir(p)? {let e=e?;let m=e.metadata()?;
        total+=if m.is_dir(){output_bytes(&e.path())?}else{m.len()};
    } Ok(total)
}
fn source(g: &Group) -> Result<(String, String)> {
    Ok((
        String::from_utf8(read(g, "source_file")?)?,
        String::from_utf8(read(g, "source_group")?)?,
    ))
}
fn xc2(p: &str) -> bool {
    Path::new(p).components().any(|c| c.as_os_str() == "xc2")
}
struct Block {
    fi: Vec<u32>,
    ri: Vec<u32>,
    vel: Vec<f64>,
    theta: Vec<f64>,
    db: Vec<f64>,
    u: Vec<f64>,
    v: Vec<f64>,
    accepted: Vec<u32>,
    raw: Vec<u64>,
    clean: Vec<u64>,
    filtered: Vec<u64>,
    tested: u64,
    invalid: u64,
}
fn process_block(
    g: &Group,
    archive: &Group,
    mode: usize,
    imager: &Capon,
    phases: &[f64; 7],
    native: &[f64],
    bound: i32,
) -> Result<Block> {
    let ds = g.dataset("xc_arr")?;
    let sh = ds.shape();
    ensure!(
        sh.len() == 4 && sh[0] == 28 && sh[1] == 5,
        "invalid XC2 shape"
    );
    let (nf, nr) = (sh[2], sh[3]);
    let ax = axes(
        nf,
        nr,
        scalar(g, "n_fft")?,
        scalar(g, "ipp")?,
        1e6,
        scalar(g, "f0")?,
        scalar(g, "f1")?,
        scalar(g, "r0")?,
        scalar(g, "r1")?,
        scalar(g, "rdec")?,
        opt(g, "fvec")?,
        opt(g, "rvec")?,
    )?;
    ensure!(
        nf == native.len()
            && ax
                .frequency_hz
                .iter()
                .zip(native)
                .all(|(f, v)| (f * 299792458. / 94e6 - v).abs() < 1e-6),
        "native axis mismatch"
    );
    let pairs = g
        .dataset(if g.link_exists("ch_pairs") {
            "ch_pairs"
        } else {
            "pairs"
        })?
        .read_2d::<u64>()?;
    ensure!(pairs.shape() == [28, 2], "invalid channel pairs");
    let pairs: Vec<_> = pairs
        .outer_iter()
        .map(|p| [p[0] as usize, p[1] as usize])
        .collect();
    let autos = validate_pairs(&pairs)?;
    let z: Array3<Z> = ds.read_slice(s![.., 0, .., ..])?;
    ensure!(
        z.iter().all(|z| z.r.is_finite() && z.i.is_finite()),
        "nonfinite source spectrum"
    );
    let power: Vec<_> = (0..nf)
        .flat_map(|f| {
            let z = &z;
            let autos = &autos;
            (0..nr).map(move |r| {
                autos
                    .iter()
                    .map(|&c| {
                        let a = z[[c, f, r]];
                        (a.r as f64).hypot(a.i as f64)
                    })
                    .sum::<f64>()
            })
        })
        .collect();
    let maxf = ax.frequency_hz.iter().map(|f| f.abs()).fold(0., f64::max);
    let outer: Vec<_> = (0..nf)
        .filter(|&f| ax.frequency_hz[f].abs() >= 0.75 * maxf)
        .collect();
    ensure!(!outer.is_empty(), "no noise bins");
    let noise: Vec<_> = (0..nr)
        .map(|r| koki::median(&mut outer.iter().map(|&f| power[f * nr + r]).collect::<Vec<_>>()))
        .collect::<Result<_>>()?;
    let mut channel_noise = [0.; 7];
    for c in 0..7 {
        let mut scratch = Vec::with_capacity(outer.len() * nr);
        for &f in &outer {
            for r in 0..nr {
                let a = z[[autos[c], f, r]];
                scratch.push((a.r as f64).hypot(a.i as f64));
            }
        }
        channel_noise[c] = koki::median(&mut scratch)?;
    }
    let nh = (2 * bound + 1) as usize;
    let mut b = Block {
        fi: vec![],
        ri: vec![],
        vel: vec![],
        theta: vec![],
        db: vec![],
        u: vec![],
        v: vec![],
        accepted: vec![0; nf],
        raw: vec![0; nh],
        clean: vec![0; nh],
        filtered: vec![0; nh],
        tested: 0,
        invalid: 0,
    };
    // Match the supplied code's nominal 0.15 km sample spacing for selection.
    let r0 = scalar(g, "r0")?;
    let rd = scalar(g, "rdec")?;
    for r in 0..nr {
        let range = (r0 + r * rd) as f64 * 0.15;
        if !(80. ..=92.).contains(&range) {
            continue;
        }
        let p: Vec<_> = (0..nf).map(|f| power[f * nr + r]).collect();
        let keep = koki::leakage(&p, &ax.frequency_hz, noise[r])?;
        for f in 0..nf {
            let db = 10. * (p[f] / noise[r].max(1e-12)).max(1e-12).log10();
            if db < 5. {
                continue;
            }
            let h = ((native[f] + bound as f64 + 0.5).floor() as usize).min(nh - 1);
            b.raw[h] += 1;
            if !keep[f] {
                continue;
            }
            b.clean[h] += 1;
            b.tested += 1;
            let packed: Vec<_> = (0..28)
                .map(|c| {
                    let a = z[[c, f, r]];
                    C::new(a.r as f64, a.i as f64)
                })
                .collect();
            let peak = koki::covariance(&packed, &pairs, phases, &channel_noise)
                .and_then(|r| imager.peak(r, 0.01));
            let peak = match peak {
                Ok(p) => p,
                Err(_) => {
                    b.invalid += 1;
                    continue;
                }
            };
            if peak.theta_deg > 2. {
                continue;
            }
            b.filtered[h] += 1;
            if native[f].abs() < 15. {
                continue;
            }
            b.accepted[f] += 1;
            b.fi.push(f as u32);
            b.ri.push(r as u32);
            b.vel.push(native[f]);
            b.theta.push(peak.theta_deg);
            b.db.push(db);
            b.u.push(peak.u);
            b.v.push(peak.v);
        }
    }
    if !b.fi.is_empty() {
        let a = archive.create_group(&format!("{mode:06}"))?;
        // Native indices and source references avoid duplicating dense XC blocks.
        a.new_dataset_builder()
            .with_data(pairs_array(&pairs).view())
            .create("ch_pairs")?;
        for n in ["n_fft", "ipp", "f0", "f1", "r0", "r1", "rdec", "i0"] {
            write(&a, n, &[scalar(g, n)? as u64])?;
        }
        write(&a, "frequency_hz", &ax.frequency_hz)?;
        write(&a, "velocity_mps", native)?;
        write(&a, "range_idx", &b.ri)?;
        write(&a, "freq_idx", &b.fi)?;
        for (name,values) in [("capon_theta_deg",&b.theta),("total_noise_db",&b.db),("capon_u",&b.u),("capon_v",&b.v)] {
            write(&a,name,&values.iter().map(|&x|x as f32).collect::<Vec<_>>())?;
        }
        write(&a, "noise_power_by_range", &noise)?;
        write(&a, "channel_noise_power", &channel_noise)?;
    }
    Ok(b)
}
fn pairs_array(p: &[[usize; 2]]) -> ndarray::Array2<u64> {
    ndarray::Array2::from_shape_fn((28, 2), |(i, j)| p[i][j] as u64)
}
fn day(
    date: &str,
    a: &Args,
    imager: &Capon,
    positions: &[[f64; 3]; 7],
    phases: &[f64; 7],
) -> Result<Option<(u64, u64)>> {
    let output = a.output.join(format!("{date}.koki.h5"));
    if output.exists() {
        let f = File::open(output)?;
        return Ok(Some((
            read::<u64>(&f, "extreme_cell_count")?[0],
            read::<u64>(&f, "input_errors")?[0],
        )));
    }
    let overview = File::open(a.input.join(format!("{date}.full-doppler.h5")))?;
    let sources = overview.group("sources")?;
    let mut names = sources.member_names()?;
    names.sort();
    if names.is_empty() || !xc2(&source(&sources.group(&names[0])?)?.0) {
        return Ok(None);
    }
    if a.max_modes > 0 {
        names.truncate(a.max_modes);
    }
    let native = read::<f64>(&overview, "velocity_mps")?;
    let nf = native.len();
    let bound = native.iter().map(|x| x.abs()).fold(25., f64::max).ceil() as i32;
    let nh = (2 * bound + 1) as usize;
    let mut raw = vec![0u64; nh];
    let mut clean = raw.clone();
    let mut filtered = raw.clone();
    let mut accepted = vec![0u32; names.len() * nf];
    let mut tested = 0u64;
    let mut invalid = 0u64;
    let mut errors = Vec::new();
    let mut cells = 0u64;
    let start = Instant::now();
    let local = a.output.join("event-tables");
    fs::create_dir_all(&local)?;
    let fullpath = local.join(format!("{date}.events.h5.partial"));
    ensure!(
        !fullpath.exists(),
        "partial archive exists, inspect before rerunning: {}",
        fullpath.display()
    );
    let full = File::create_excl(&fullpath)?;
    let blocks = full.create_group("blocks")?;
    text(&full, "date_utc", date)?;
    for (m, n) in names.iter().enumerate() {
        let src = sources.group(n)?;
        let (path, key) = source(&src)?;
        let result = (|| -> Result<Block> {
            ensure!(xc2(&path), "non-XC2 source");
            let f = File::open(&path)?;
            let g=f.group(&key)?;
            let shape=g.dataset("xc_arr")?.shape();
            ensure!(shape.len()==4,"invalid shape");
            // Uncompressed upper bound plus generous HDF5 and day-summary headroom.
            let reserve=(shape[2] as u64)*(shape[3] as u64)*64+16_000_000;
            ensure!(output_bytes(&a.output)?.saturating_add(reserve)<=a.max_output_bytes,"ARCHIVE_BUDGET_STOP: aggregate limit {} bytes",a.max_output_bytes);
            process_block(&g, &blocks, m, imager, phases, &native, bound)
        })();
        match result {
            Ok(b) => {
                tested += b.tested;
                invalid += b.invalid;
                cells += b.fi.len() as u64;
                for k in 0..nh {
                    raw[k] += b.raw[k];
                    clean[k] += b.clean[k];
                    filtered[k] += b.filtered[k];
                }
                accepted[m * nf..(m + 1) * nf].copy_from_slice(&b.accepted);
                if !b.fi.is_empty() {
                    let block = blocks.group(&format!("{m:06}"))?;
                    text(&block, "source_file", &path)?;
                    text(&block, "source_group", &key)?;
                }
            }
            Err(e) => {
                if e.to_string().starts_with("ARCHIVE_BUDGET_STOP") {return Err(e);}
                errors.push(format!("{path} {key}: {e:#}"));
                eprintln!("KOKI_SOURCE_ERROR {date} mode={m}: {e:#}");
            }
        }
        full.flush()?;
        println!("KOKI_MODE {date} mode={}/{} tested={tested} extreme_cells={cells} invalid={invalid} elapsed_s={:.1}",m+1,names.len(),start.elapsed().as_secs_f64());
    }
    full.flush()?;
    drop(blocks);
    drop(full);
    fs::rename(fullpath, local.join(format!("{date}.events.h5")))?;
    let partial = output.with_extension("h5.partial");
    let f = File::create_excl(&partial)?;
    text(&f, "date_utc", date)?;
    text(&f, "schema", "pansy-koki-xc2-v1")?;
    text(&f,"definitions","Supplied bin-by-bin procedure: zenith; nominal range 80--92 km at 0.15 km/sample; total/noise >=5 dB; per-range median outer 25 percent; per-channel median outer bins over all ranges. Strong peak >30 dB: monotonic primary cluster, baseline-subtracted second moment, minimum sigma=half frequency bin, retain +/-3 sigma. Covariance normalized by max(auto-noise,noise), unit diagonal, Hermitian PSD eigenvalue clipping, 1 percent trace diagonal loading. 60x60 direction grid within 9 deg; theta<=2 deg (no additional boundary exclusion). All surviving bins enter 1 m/s histograms; extreme events abs(v)>=15 m/s. Compact event tables and source references only: no copied XC arrays. No complementary-code veto.")?;
    write(&f,"aggregate_output_limit_bytes",&[a.max_output_bytes])?;
    write(&f, "velocity_mps", &native)?;
    write(&f, "accepted_extreme_cells_mode_frequency", &accepted)?;
    write(
        &f,
        "histogram_velocity_mps",
        &(-bound..=bound).map(|v| v as f64).collect::<Vec<_>>(),
    )?;
    write(&f, "histogram_raw", &raw)?;
    write(&f, "histogram_leakage_cleaned", &clean)?;
    write(&f, "histogram_capon", &filtered)?;
    let mode_start = read::<f64>(&overview, "mode_start_hours")?;
    let mode_stop = read::<f64>(&overview, "mode_stop_hours")?;
    write(&f, "mode_start_hours", &mode_start[..names.len()])?;
    write(&f, "mode_stop_hours", &mode_stop[..names.len()])?;
    write(&f, "tested_count", &[tested])?;
    write(&f, "invalid_covariance_count", &[invalid])?;
    write(&f, "extreme_cell_count", &[cells])?;
    write(
        &f,
        "detected_mode_count",
        &[accepted
            .chunks(nf)
            .filter(|m| m.iter().any(|&x| x > 0))
            .count() as u64],
    )?;
    let input_errors = errors.len() as u64 + read::<u64>(&overview, "input_errors")?[0];
    write(&f, "input_errors", &[input_errors])?;
    text(&f, "error_details", &errors.join("\n"))?;
    let pg = f.create_group("reference")?;
    write(&pg, "positions_m", &positions.concat())?;
    write(&pg, "phases_rad", phases)?;
    f.flush()?;
    drop(pg);
    drop(f);
    fs::rename(partial, output)?;
    println!("KOKI_DAY_COMPLETE {date} tested={tested} extreme_cells={cells} errors={input_errors} elapsed_s={:.1}",start.elapsed().as_secs_f64());
    Ok(Some((cells, input_errors)))
}
fn main() -> Result<()> {
    let a = Args::parse();
    ensure!(a.max_output_bytes<=450_000_000,"maximum sweep output budget is 450 MB; leave room for mirrored summaries within 1000 MB");
    fs::create_dir_all(&a.output)?;
    fs::create_dir_all(&a.input)?;
    let (positions, phases) = pansy_legacy_capon::calibration::load_reference(&a.reference, 0)?;
    let imager = Capon::new(&positions, 47e6, 9., 60)?;
    let mut days: Vec<String> = fs::read_dir(&a.input)?
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            e.file_name()
                .to_str()
                .and_then(|s| s.strip_suffix(".full-doppler.h5"))
                .map(str::to_owned)
        })
        .collect();
    days.sort();
    if !a.root.is_empty() {
        let mut unique=std::collections::BTreeSet::new();
        for root in &a.root {
            ensure!(root.file_name().is_some_and(|n|n=="xc2"),"only XC2 roots allowed");
            for e in fs::read_dir(root)? {
                let e=e?;let n=e.file_name().to_string_lossy().into_owned();
                if n.len()>=11 && n.as_bytes()[10]==b'T' && e.file_type()?.is_dir()
                    && chrono::NaiveDate::parse_from_str(&n[..10],"%Y-%m-%d").is_ok()
                    && fs::read_dir(e.path())?.filter_map(|x|x.ok()).any(|x|x.file_name().to_string_lossy().starts_with("xc@")) {
                    unique.insert(n[..10].to_owned());
                }
            }
        }
        days=unique.into_iter().collect();
    }
    if a.pmse_season {days.retain(|d|matches!(&d[5..7],"11"|"12"|"01"|"02"));}
    days.sort_by_key(|d|(d!="2025-01-30",d.clone()));
    if let Some(d) = &a.day {
        days.retain(|v| v == d);
        ensure!(!days.is_empty(), "requested day absent");
    }
    let mut dates = Vec::new();
    let mut status = Vec::new();
    let mut counts = Vec::new();
    let mut failures = 0;
    for d in &days {
        let result=(||->Result<Option<(u64,u64)>> {
            if !a.input.join(format!("{d}.full-doppler.h5")).exists() {
                ensure!(!a.root.is_empty(),"no XC2 source roots for missing overview");
                ensure!(output_bytes(&a.input)?+5_000_000<=100_000_000,"ARCHIVE_BUDGET_STOP: overview context limit 100 MB");
                let mut cmd=std::process::Command::new(std::env::current_exe()?.with_file_name("pansy-spectra-overviews"));
                for root in &a.root {cmd.arg("--root").arg(root);}
                let status=cmd.arg("--day").arg(d).arg("--output").arg(&a.input).status()?;
                ensure!(a.input.join(format!("{d}.full-doppler.h5")).exists(),"XC2 overview preparation failed: {status}");
            }
            day(d, &a, &imager, &positions, &phases)
        })();
        match result {
            Ok(Some((n, e))) => {
                dates.push(d.replace('-', "").parse::<u32>()?);
                status.push(if e == 0 { 1u8 } else { 2 });
                counts.push(n);
                if e > 0 {
                    failures += 1;
                }
            }
            Ok(None) => continue,
            Err(e) => {
                if e.to_string().starts_with("ARCHIVE_BUDGET_STOP") {return Err(e);}
                eprintln!("KOKI_DAY_ERROR {d}: {e:#}");
                dates.push(d.replace('-', "").parse::<u32>()?);
                status.push(2);
                counts.push(0);
                failures += 1;
            }
        }
        let temp = a.output.join("koki-index.h5.partial");
        let f = File::create(&temp)?;
        write(&f, "utc_day_yyyymmdd", &dates)?;
        write(&f, "status", &status)?;
        write(&f, "extreme_cell_count", &counts)?;
        write(&f, "sweep_complete", &[0u8])?;
        f.flush()?;
        drop(f);
        fs::rename(temp, a.output.join("koki-index.h5"))?;
    }
    let f = File::open_rw(a.output.join("koki-index.h5"))?;
    f.dataset("sweep_complete")?.write_raw(&[1u8])?;
    f.flush()?;
    println!(
        "KOKI_SWEEP_COMPLETE days={} source_error_days={failures}",
        dates.len()
    );
    ensure!(failures == 0, "source errors retained in index");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn binwise_passes_a_fast_bin_without_dominant_peak_selection() {
        let path=std::env::temp_dir().join(format!("pansy-koki-fixture-{}.h5",std::process::id()));
        let f=File::create_excl(&path).unwrap();
        let g=f.create_group("input").unwrap();let archive=f.create_group("archive").unwrap();
        for (n,x) in [("n_fft",1248u64),("ipp",1600),("f0",497),("f1",753),("r0",574),("r1",575),("rdec",1),("i0",1738210500000000)] {
            g.new_dataset::<u64>().create(n).unwrap().write_scalar(&x).unwrap();
        }
        let pairs:Vec<_>=(0..7).map(|i|[i,i]).chain((0..7).flat_map(|i|(i+1..7).map(move|j|[i,j]))).collect();
        g.new_dataset_builder().with_data(pairs_array(&pairs).view()).create("ch_pairs").unwrap();
        let mut z=ndarray::Array4::<Z>::default((28,5,256,1));
        for (p,&[i,j]) in pairs.iter().enumerate() {for k in 0..256 {
            let signal=if k==127 {500.} else if k==20 {100.}else{0.};
            z[[p,0,k,0]]=Z{r:signal+if i==j{1.}else{0.},i:0.};
        }}
        g.new_dataset_builder().with_data(z.view()).create("xc_arr").unwrap();
        let native:Vec<_>=(497..753).map(|k|(k as f64-624.)/9.984*299792458./94e6).collect();
        let pos=[[0.,0.,0.],[10.,0.,0.],[0.,11.,0.],[13.,17.,0.],[-17.,12.,0.],[6.,-21.,0.],[-19.,-16.,0.]];
        let imager=Capon::new(&pos,47e6,9.,60).unwrap();
        let b=process_block(&g,&archive,0,&imager,&[0.;7],&native,90).unwrap();
        assert_eq!(b.raw.iter().sum::<u64>(),2);
        assert_eq!(b.filtered.iter().sum::<u64>(),2);
        assert_eq!(b.fi,vec![20]); // near-zero main peak did not hide this wing
        assert_eq!(b.accepted[20],1);
        assert!(!archive.group("000000").unwrap().link_exists("xc_arr"));
        assert!(archive.group("000000").unwrap().link_exists("freq_idx"));
        drop(archive);drop(g);drop(f);std::fs::remove_file(path).unwrap();
    }
}
