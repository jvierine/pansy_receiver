use anyhow::{bail, ensure, Context, Result};
use clap::Parser;
use hdf5::{
    types::{CompoundField, CompoundType, TypeDescriptor},
    File, Group, H5Type,
};
use ndarray::{s, Array3};
use pansy_legacy_capon::{axes, covariance, validate_pairs, Capon, C, LIGHT};
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs,
    path::{Path, PathBuf},
    time::Instant,
};

#[derive(Parser, Debug)]
#[command(
    about = "Stream legacy dense XC2/xc_rvd; emit reversible extreme radial-velocity candidates"
)]
struct Args {
    #[arg(long, required = true)]
    root: Vec<PathBuf>,
    #[arg(long)]
    day: String,
    #[arg(long)]
    reference: PathBuf,
    #[arg(long)]
    output: PathBuf,
    #[arg(long, default_value_t = 0)]
    beam: usize,
    #[arg(long, default_value_t = 80.)]
    range_min: f64,
    #[arg(long, default_value_t = 92.)]
    range_max: f64,
    /// Historical threshold: 10log10(total power / estimated noise), NOT excess SNR.
    #[arg(long, default_value_t = 5.)]
    snr_db: f64,
    #[arg(long, default_value_t = 15.)]
    min_abs_velocity: f64,
    #[arg(long, default_value_t = 2.)]
    zenith_max_deg: f64,
    #[arg(long, default_value_t = 9.)]
    radius_deg: f64,
    #[arg(long, default_value_t = 61)]
    grid_size: usize,
    #[arg(long, default_value_t = 0.01)]
    loading: f64,
    #[arg(long, default_value_t = 47e6)]
    radar_frequency_hz: f64,
    #[arg(long, default_value_t = 1e6)]
    sample_rate_hz: f64,
    #[arg(long, default_value_t = 0.75)]
    noise_outer_fraction: f64,
    #[arg(long, default_value_t = 600.)]
    event_gap_seconds: f64,
    #[arg(long, default_value_t = 0.7)]
    event_range_link_km: f64,
    #[arg(long, default_value_t = 2.)]
    event_velocity_link_mps: f64,
    /// Bound a pilot without changing selection thresholds. Zero scans the day.
    #[arg(long, default_value_t = 0)]
    max_blocks: usize,
}
#[repr(C)]
#[derive(Clone, Copy, Default)]
struct StoredComplex {
    r: f32,
    i: f32,
}
// Explicit layout avoids derive-time offset assumptions; matches h5py complex64.
unsafe impl H5Type for StoredComplex {
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
fn scalar(g: &Group, n: &str) -> Result<u64> {
    Ok(g.dataset(n)?.read_scalar::<u64>()?)
}
fn string(g: &Group, n: &str, s: &str) -> Result<()> {
    g.new_dataset_builder().with_data(s.as_bytes()).create(n)?;
    Ok(())
}
fn write<T: H5Type>(g: &Group, n: &str, v: &[T]) -> Result<()> {
    g.new_dataset_builder().with_data(v).create(n)?;
    Ok(())
}
fn optional_axis(g: &Group, n: &str) -> Result<Option<Vec<f64>>> {
    if g.link_exists(n) {
        Ok(Some(g.dataset(n)?.read_raw()?))
    } else {
        Ok(None)
    }
}

fn reference(root: &Path, beam: usize) -> Result<([[f64; 3]; 7], [f64; 7])> {
    let text = fs::read_to_string(root.join("cfg/antpos.csv"))?;
    let mut groups: HashMap<String, Vec<[f64; 3]>> = HashMap::new();
    for line in text.lines().skip(1).filter(|l| !l.trim().is_empty()) {
        let f: Vec<_> = line.split(',').collect();
        ensure!(f.len() >= 7, "invalid antenna position row");
        groups.entry(f[1].to_owned()).or_default().push([
            f[4].trim().parse()?,
            f[5].trim().parse()?,
            f[6].trim().parse()?,
        ]);
    }
    let connections = fs::read_to_string(root.join("cfg/connections.txt"))?;
    let names: Vec<_> = connections
        .lines()
        .filter_map(|l| l.split(',').next())
        .filter(|n| !n.is_empty() && *n != "RFTX")
        .collect();
    ensure!(names.len() == 7, "expected seven receive subarrays");
    let mut positions = [[0.; 3]; 7];
    for (i, name) in names.iter().enumerate() {
        let p = groups.get(*name).context("unknown subarray")?;
        ensure!(!p.is_empty(), "empty subarray");
        for j in 0..3 {
            positions[i][j] = p.iter().map(|r| r[j]).sum::<f64>() / p.len() as f64;
        }
    }
    let phase = File::open(root.join(format!("data/phases{beam}.h5")))?
        .dataset("phasecal")?
        .read_raw::<f64>()?;
    let phases: [f64; 7] = phase
        .try_into()
        .map_err(|_| anyhow::anyhow!("expected seven phasecal values"))?;
    Ok((positions, phases))
}

#[derive(Clone)]
struct Cell {
    start: u64,
    stop: u64,
    range: f64,
    velocity: f64,
    snr: f64,
    theta: f64,
}
struct Events {
    cells: Vec<Cell>,
    parent: Vec<usize>,
    previous: HashMap<(i64, i64), Vec<usize>>,
    last_start: Option<u64>,
}
impl Events {
    fn new() -> Self {
        Self {
            cells: Vec::new(),
            parent: Vec::new(),
            previous: HashMap::new(),
            last_start: None,
        }
    }
    fn find(&mut self, mut i: usize) -> usize {
        while self.parent[i] != i {
            self.parent[i] = self.parent[self.parent[i]];
            i = self.parent[i];
        }
        i
    }
    fn append(&mut self, cells: Vec<Cell>, start: u64, a: &Args) {
        let mut current: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        if self.last_start.is_some_and(|t| start < t) {
            self.previous.clear();
        }
        for c in cells {
            let key = (
                (c.range / a.event_range_link_km).floor() as i64,
                (c.velocity / a.event_velocity_link_mps).floor() as i64,
            );
            let i = self.cells.len();
            self.parent.push(i);
            self.cells.push(c);
            current.entry(key).or_default().push(i);
        }
        // Union every neighboring pair, including bridges between components.
        let ids: Vec<_> = current.values().flatten().copied().collect();
        for i in ids {
            let x = self.cells[i].clone();
            let key = (
                (x.range / a.event_range_link_km).floor() as i64,
                (x.velocity / a.event_velocity_link_mps).floor() as i64,
            );
            for dr in -1..=1 {
                for dv in -1..=1 {
                    let neighbors: Vec<_> = [&self.previous, &current]
                        .iter()
                        .flat_map(|m| {
                            m.get(&(key.0 + dr, key.1 + dv))
                                .into_iter()
                                .flatten()
                                .copied()
                        })
                        .collect();
                    for j in neighbors {
                        let y = &self.cells[j];
                        if x.velocity.signum() == y.velocity.signum()
                            && x.start.saturating_sub(y.stop) as f64 <= a.event_gap_seconds * 1e6
                            && (x.range - y.range).abs() <= a.event_range_link_km
                            && (x.velocity - y.velocity).abs() <= a.event_velocity_link_mps
                        {
                            let p = self.find(i);
                            let q = self.find(j);
                            self.parent[p.max(q)] = p.min(q);
                        }
                    }
                }
            }
        }
        self.previous = current;
        self.last_start = Some(start);
    }
    fn save(&mut self, out: &File) -> Result<usize> {
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for i in 0..self.cells.len() {
            let p = self.find(i);
            groups.entry(p).or_default().push(i);
        }
        let root = out.create_group("events")?;
        for (number, ids) in groups.values().enumerate() {
            let c: Vec<_> = ids.iter().map(|&i| &self.cells[i]).collect();
            let g = root.create_group(&format!("{number:06}"))?;
            let start = c.iter().map(|x| x.start).min().unwrap();
            let stop = c.iter().map(|x| x.stop).max().unwrap();
            let max = c
                .iter()
                .max_by(|x, y| x.velocity.abs().total_cmp(&y.velocity.abs()))
                .unwrap();
            let modes: HashSet<_> = c.iter().map(|x| x.start).collect();
            write(&g, "start_us", &[start])?;
            write(&g, "stop_us_exclusive", &[stop])?;
            write(&g, "cell_count", &[c.len() as u64])?;
            write(&g, "mode_count", &[modes.len() as u64])?;
            write(
                &g,
                "range_min_km",
                &[c.iter().map(|x| x.range).fold(f64::INFINITY, f64::min)],
            )?;
            write(
                &g,
                "range_max_km",
                &[c.iter().map(|x| x.range).fold(f64::NEG_INFINITY, f64::max)],
            )?;
            write(&g, "maximum_abs_velocity_mps", &[max.velocity.abs()])?;
            write(&g, "velocity_at_maximum_mps", &[max.velocity])?;
            write(
                &g,
                "peak_total_noise_db",
                &[c.iter().map(|x| x.snr).fold(f64::NEG_INFINITY, f64::max)],
            )?;
            write(
                &g,
                "maximum_theta_deg",
                &[c.iter().map(|x| x.theta).fold(0., f64::max)],
            )?;
            let utc = chrono::DateTime::from_timestamp_micros(start as i64)
                .context("bad event timestamp")?;
            println!("EVENT {number} {utc} modes={} cells={} range={:.2}..{:.2}km max_radial={:.2}m/s peak_total_noise={:.2}dB",modes.len(),c.len(),c.iter().map(|x|x.range).fold(f64::INFINITY,f64::min),c.iter().map(|x|x.range).fold(0.,f64::max),max.velocity,c.iter().map(|x|x.snr).fold(f64::NEG_INFINITY,f64::max));
        }
        Ok(groups.len())
    }
}

fn process(
    g: &Group,
    path: &Path,
    index: usize,
    out: &File,
    a: &Args,
    imager: &Capon,
    phases: &[f64; 7],
) -> Result<(Vec<Cell>, u64)> {
    let start = scalar(g, "i0")?;
    let stop = scalar(g, "i1")?.checked_add(1).context("invalid stop")?;
    ensure!(stop > start, "invalid mode bounds");
    let ds = g.dataset("xc_arr")?;
    let shape = ds.shape();
    ensure!(
        shape.len() == 4 && shape[0] == 28 && shape[1] == 5,
        "expected [28,5,frequency,range]"
    );
    let nf = shape[2];
    let nr = shape[3];
    let axes = axes(
        nf,
        nr,
        scalar(g, "n_fft")? as usize,
        scalar(g, "ipp")? as usize,
        a.sample_rate_hz,
        scalar(g, "f0")? as usize,
        scalar(g, "f1")? as usize,
        scalar(g, "r0")? as usize,
        scalar(g, "r1")? as usize,
        scalar(g, "rdec")? as usize,
        optional_axis(g, "fvec")?,
        optional_axis(g, "rvec")?,
    )?;
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
        .map(|r| [r[0] as usize, r[1] as usize])
        .collect();
    let autos = validate_pairs(&pairs)?;
    let ranges: Vec<_> = axes
        .range_km
        .iter()
        .enumerate()
        .filter_map(|(i, &r)| (r >= a.range_min && r <= a.range_max).then_some(i))
        .collect();
    ensure!(!ranges.is_empty(), "selected range interval absent");
    let lo = ranges[0];
    let hi = ranges[ranges.len() - 1] + 1;
    let n = hi - lo;
    let data: Array3<StoredComplex> = ds.read_slice(s![.., a.beam, .., lo..hi])?;
    let maxf = axes.frequency_hz.iter().map(|v| v.abs()).fold(0., f64::max);
    let guard: Vec<_> = axes
        .frequency_hz
        .iter()
        .enumerate()
        .filter_map(|(i, f)| (f.abs() >= a.noise_outer_fraction * maxf).then_some(i))
        .collect();
    ensure!(guard.len() >= 8, "insufficient noise guards");
    let power =
        |f: usize, r: usize| -> f64 { autos.iter().map(|&c| data[[c, f, r]].r as f64).sum() };
    let mut noise = vec![0.; n];
    let mut scratch = Vec::with_capacity(guard.len());
    for (r, v) in noise.iter_mut().enumerate() {
        scratch.clear();
        for &f in &guard {
            let p = power(f, r);
            ensure!(p.is_finite() && p > 0., "invalid guard power");
            scratch.push(p);
        }
        scratch.sort_unstable_by(f64::total_cmp);
        let k = scratch.len() / 2;
        *v = if scratch.len() % 2 == 0 {
            (scratch[k - 1] + scratch[k]) / 2.
        } else {
            scratch[k]
        };
    }
    let block = out.group("blocks")?.create_group(&format!("{index:06}"))?;
    string(&block, "source_file", &path.to_string_lossy())?;
    string(&block, "source_group", &g.name())?;
    string(&block, "frequency_axis_source", axes.frequency_source)?;
    write(&block, "mode_start_us", &[start])?;
    write(&block, "mode_stop_us_exclusive", &[stop])?;
    write(&block, "frequency_hz", &axes.frequency_hz)?;
    write(&block, "range_km", &axes.range_km)?;
    write(&block, "noise_range_start_index", &[lo as u64])?;
    write(&block, "estimated_total_noise_power", &noise)?;
    let mut fi = Vec::new();
    let mut ri = Vec::new();
    let mut snrs = Vec::new();
    let mut theta = Vec::new();
    let mut us = Vec::new();
    let mut vs = Vec::new();
    let mut flags = Vec::new();
    let mut loadings = Vec::new();
    let mut accepted = Vec::new();
    let mut raw_counts = vec![0u64; nf];
    let mut accepted_counts = vec![0u64; nf];
    let threshold = 10f64.powf(a.snr_db / 10.);
    let mut packed = [C::new(0., 0.); 28];
    for f in 0..nf {
        let velocity = axes.frequency_hz[f] * LIGHT / (2. * a.radar_frequency_hz);
        for r in 0..n {
            let total = power(f, r);
            ensure!(total.is_finite() && total >= 0., "invalid signal power");
            let ratio = total / noise[r];
            if ratio < threshold {
                continue;
            }
            raw_counts[f] += 1;
            if velocity.abs() < a.min_abs_velocity {
                continue;
            }
            let snr = 10. * ratio.log10();
            for c in 0..28 {
                let z = data[[c, f, r]];
                packed[c] = C::new(z.r as f64, z.i as f64);
            }
            let peak = covariance(&packed, &pairs, phases).and_then(|r| imager.peak(r, a.loading));
            let (t, u, v, flag, loading) = match peak {
                Ok(p) => {
                    let flag = if p.boundary {
                        3u8
                    } else if p.theta_deg > a.zenith_max_deg {
                        1
                    } else {
                        0
                    };
                    (p.theta_deg, p.u, p.v, flag, p.loading)
                }
                Err(_) => (f64::NAN, f64::NAN, f64::NAN, 2, f64::NAN),
            };
            fi.push(f as u32);
            ri.push((r + lo) as u32);
            snrs.push(snr);
            theta.push(t);
            us.push(u);
            vs.push(v);
            flags.push(flag);
            loadings.push(loading);
            if flag == 0 {
                accepted_counts[f] += 1;
                accepted.push(Cell {
                    start,
                    stop,
                    range: axes.range_km[r + lo],
                    velocity,
                    snr,
                    theta: t,
                });
            }
        }
    }
    write(&block, "frequency_index", &fi)?;
    write(&block, "range_index", &ri)?;
    write(&block, "total_noise_db", &snrs)?;
    write(&block, "capon_theta_deg", &theta)?;
    write(&block, "capon_u", &us)?;
    write(&block, "capon_v", &vs)?;
    write(&block, "decision", &flags)?;
    write(&block, "diagonal_loading_power", &loadings)?;
    write(&block, "raw_snr_cell_counts_by_frequency", &raw_counts)?;
    write(
        &block,
        "accepted_extreme_cell_counts_by_frequency",
        &accepted_counts,
    )?;
    println!(
        "BLOCK {start} source={} axis={} tested_extreme={} accepted={} invalid_covariance={}",
        path.display(),
        axes.frequency_source,
        fi.len(),
        accepted.len(),
        flags.iter().filter(|&&x| x == 2).count()
    );
    Ok((accepted, start))
}

fn main() -> Result<()> {
    let a = Args::parse();
    let begun = Instant::now();
    let _ = chrono::NaiveDate::parse_from_str(&a.day, "%Y-%m-%d").context("invalid UTC day")?;
    ensure!(
        a.beam < 5
            && a.range_min.is_finite()
            && a.range_max.is_finite()
            && a.range_min < a.range_max
            && a.snr_db.is_finite()
            && a.min_abs_velocity.is_finite()
            && a.min_abs_velocity >= 0.
            && a.zenith_max_deg > 0.
            && a.zenith_max_deg < a.radius_deg
            && a.loading.is_finite()
            && a.loading > 0.
            && a.sample_rate_hz.is_finite()
            && a.sample_rate_hz > 0.
            && a.noise_outer_fraction > 0.
            && a.noise_outer_fraction < 1.
            && a.event_gap_seconds.is_finite()
            && a.event_gap_seconds > 0.
            && a.event_range_link_km.is_finite()
            && a.event_range_link_km > 0.
            && a.event_velocity_link_mps.is_finite()
            && a.event_velocity_link_mps > 0.,
        "invalid options"
    );
    ensure!(
        !a.output.exists(),
        "output exists; choose a new output path"
    );
    let partial = a.output.with_extension("h5.partial");
    ensure!(
        !partial.exists(),
        "partial output exists; inspect before retrying"
    );
    let (positions, phases) = reference(&a.reference, a.beam)?;
    let imager = Capon::new(&positions, a.radar_frequency_hz, a.radius_deg, a.grid_size)?;
    let mut files = Vec::new();
    for root in &a.root {
        for entry in fs::read_dir(root).with_context(|| format!("opening {}", root.display()))? {
            let dir = entry?;
            if !dir
                .file_name()
                .to_string_lossy()
                .starts_with(&format!("{}T", a.day))
            {
                continue;
            }
            for entry in fs::read_dir(dir.path())? {
                let p = entry?.path();
                let name = p.file_name().unwrap().to_string_lossy();
                if name.starts_with("xc@") && name.ends_with(".h5") {
                    files.push(p);
                }
            }
        }
    }
    files.sort_by_key(|p| p.file_name().unwrap().to_os_string());
    ensure!(!files.is_empty(), "no matching legacy XC files");
    if let Some(parent) = a.output.parent() {
        fs::create_dir_all(parent)?;
    }
    let out = File::create_excl(&partial)?;
    out.create_group("blocks")?;
    let meta = out.create_group("metadata")?;
    string(&meta, "schema", "pansy-legacy-extreme-capon-v1")?;
    string(&meta, "options", &format!("{a:?}"))?;
    string(&meta,"definitions","Candidate cells, not confirmed vertical drafts. decision: 0 accepted angle; 1 off-zenith; 2 invalid covariance; 3 search boundary. Noise is the median of outer-band total power, an uncalibrated legacy proxy. No three-width spectral clipping. Event components join neighboring cells within a mode and the immediately previous successfully read mode; singleton events retained. Histogram counts are correlated spectral-cell occupancy, not independent events.")?;
    write(&meta, "phasecal_radians", &phases)?;
    let flat: Vec<_> = positions.iter().flatten().copied().collect();
    write(&meta, "positions_m_xyz_flat", &flat)?;
    let mut seen = HashSet::new();
    let mut events = Events::new();
    let mut count = 0usize;
    let mut errors = 0usize;
    'files: for path in files {
        let file = match File::open(&path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("FILE_ERROR {}: {e}", path.display());
                errors += 1;
                continue;
            }
        };
        let mut groups = file.member_names()?;
        groups.sort_by_key(|s| s.parse::<u64>().unwrap_or(u64::MAX));
        for name in groups {
            if name.parse::<u64>().is_err() {
                continue;
            }
            let g = file.group(&name)?;
            if !g.link_exists("xc_arr") {
                continue;
            }
            let start = scalar(&g, "i0")?;
            if !seen.insert(start) {
                eprintln!("DUPLICATE {start} {} (first copy retained)", path.display());
                continue;
            }
            match process(&g, &path, count, &out, &a, &imager, &phases) {
                Ok((cells, t)) => events.append(cells, t, &a),
                Err(e) => {
                    eprintln!("BLOCK_ERROR {} {name}: {e:#}", path.display());
                    errors += 1;
                    events.previous.clear();
                }
            }
            count += 1;
            out.flush()?;
            if a.max_blocks > 0 && count >= a.max_blocks {
                break 'files;
            }
        }
    }
    let n = events.save(&out)?;
    write(&meta, "blocks_attempted", &[count as u64])?;
    write(&meta, "input_errors", &[errors as u64])?;
    write(&meta, "event_count", &[n as u64])?;
    write(&meta, "elapsed_seconds", &[begun.elapsed().as_secs_f64()])?;
    out.flush()?;
    drop(meta);
    drop(out);
    fs::rename(&partial, &a.output)?;
    println!("SUMMARY blocks={count} input_errors={errors} grouped_candidates={n} accepted_cells={} seconds={:.3} output={}",events.cells.len(),begun.elapsed().as_secs_f64(),a.output.display());
    if errors > 0 {
        bail!("output retained with {errors} explicit input errors");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args() -> Args {
        Args::parse_from([
            "test",
            "--root",
            "unused",
            "--day",
            "2026-01-30",
            "--reference",
            "unused",
            "--output",
            "unused.h5",
            "--grid-size",
            "11",
        ])
    }
    #[test]
    fn real_hdf5_rvd_shape_and_known_source() -> Result<()> {
        let path = std::env::temp_dir().join(format!(
            "pansy-legacy-fixture-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)?
                .as_nanos()
        ));
        fs::create_dir(&path)?;
        let input = File::create(path.join("input.h5"))?;
        let g = input.create_group("1738210500000000")?;
        for (n, v) in [
            ("i0", 1738210500000000u64),
            ("i1", 1738210585999999),
            ("n_fft", 32),
            ("ipp", 1600),
            ("f0", 0),
            ("f1", 31),
            ("r0", 534),
            ("r1", 546),
            ("rdec", 4),
        ] {
            g.new_dataset::<u64>()
                .shape(())
                .create(n)?
                .write_scalar(&v)?;
        }
        let pairs: Vec<_> = (0..7)
            .map(|i| [i, i])
            .chain((0..7).flat_map(|i| ((i + 1)..7).map(move |j| [i, j])))
            .collect();
        let pp = ndarray::Array2::from_shape_vec(
            (28, 2),
            pairs.iter().flatten().map(|&x| x as u64).collect(),
        )?;
        g.new_dataset_builder()
            .with_data(pp.view())
            .create("ch_pairs")?;
        let mut z = ndarray::Array4::<StoredComplex>::default((28, 5, 32, 3));
        for (p, &[i, j]) in pairs.iter().enumerate() {
            for b in 0..5 {
                for f in 0..32 {
                    for r in 0..3 {
                        z[[p, b, f, r]].r = if i == j { 1. } else { 0. };
                        if f == 18 && r == 1 {
                            z[[p, b, f, r]].r += 100.;
                        }
                    }
                }
            }
        }
        g.new_dataset_builder()
            .with_data(z.view())
            .create("xc_arr")?;
        write(
            &g,
            "fvec",
            &(0..32)
                .map(|i| (i as f64 - 16.) / 0.256)
                .collect::<Vec<_>>(),
        )?;
        let out = File::create(path.join("output.h5"))?;
        out.create_group("blocks")?;
        let positions = std::array::from_fn(|i| [i as f64 * 9., (i * i) as f64 * 3., 0.]);
        let imager = Capon::new(&positions, 47e6, 9., 11)?;
        let (cells, _) = process(
            &g,
            &path.join("input.h5"),
            0,
            &out,
            &args(),
            &imager,
            &[0.; 7],
        )?;
        assert_eq!(cells.len(), 1);
        assert!(cells[0].theta < 1e-10);
        assert_eq!(out.dataset("blocks/000000/frequency_hz")?.size(), 32);
        assert_eq!(
            out.dataset("blocks/000000/frequency_index")?
                .read_raw::<u32>()?,
            vec![18]
        );
        assert_eq!(
            out.dataset("blocks/000000/range_index")?
                .read_raw::<u32>()?,
            vec![1]
        );
        drop(g);
        drop(input);
        drop(out);
        fs::remove_dir_all(path)?;
        Ok(())
    }
    #[test]
    fn connected_events_merge_bridges_and_keep_singletons() {
        let a = args();
        let mut e = Events::new();
        let cell = |start, range, velocity| Cell {
            start,
            stop: start + 86_000_000,
            range,
            velocity,
            snr: 10.,
            theta: 1.,
        };
        e.append(
            vec![cell(1_000_000_000, 85., 16.), cell(1_000_000_000, 85., 19.)],
            1_000_000_000,
            &a,
        );
        e.append(vec![cell(1_220_000_000, 85., 17.5)], 1_220_000_000, &a);
        assert_eq!(e.find(0), e.find(1));
        assert_eq!(e.find(0), e.find(2));
        e.append(vec![cell(2_500_000_000, 85., 17.5)], 2_500_000_000, &a);
        assert_ne!(e.find(0), e.find(3));
        e.append(vec![cell(2_720_000_000, 85., -17.5)], 2_720_000_000, &a);
        assert_ne!(e.find(3), e.find(4));
    }
}
