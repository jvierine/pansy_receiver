//! Inspect central wings without modifying spectra or relabelling ghosts.
use anyhow::{ensure, Result};
use hdf5::{
    types::{CompoundField, CompoundType, TypeDescriptor},
    File, Group, H5Type,
};
use ndarray::{s, Array2};
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
fn write<T: H5Type>(g: &Group, n: &str, v: &[T]) -> Result<()> {
    g.new_dataset_builder()
        .with_data(v)
        .shuffle()
        .deflate(4)
        .create(n)?;
    Ok(())
}
struct Block {
    name: String,
    t: f64,
    start: f64,
    stop: f64,
    range: Vec<f64>,
    v: Vec<f64>,
    db: Vec<f32>,
}
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    ensure!(a.len() == 3, "usage: pansy-wing-context CATALOGUE OUTPUT");
    ensure!(!std::path::Path::new(&a[2]).exists(), "output exists");
    let input = File::open(&a[1])?;
    let partial = format!("{}.partial", a[2]);
    let out = File::create_excl(&partial)?;
    write(&out, "catalogue_source", a[1].as_bytes())?;
    write(&out,"definitions",b"Exploratory central-wing inspection, not confirmed drafts. Selection: largest noise-normalized excess integrated over 15<|v|<70 m/s at one 80--92 km gate in each specified time window. All stored Dopplers retained in plots. No independent-trial significance claimed; noise is historical outer-band median proxy. Covariance screening copied unchanged from catalogue; absent values were not tested. Display time edges may bridge ordinary between-mode gaps, never gaps >10 min.")?;
    for (w, (lo, hi)) in [(0.75, 2.25), (9., 10.), (21.5, 22.25)]
        .into_iter()
        .enumerate()
    {
        let mut blocks = Vec::new();
        for name in input.group("blocks")?.member_names()? {
            let g = input.group(&format!("blocks/{name}"))?;
            if !g.link_exists("frequency_hz") {
                continue;
            }
            let t0 = read::<u64>(&g, "mode_start_us")?[0];
            let day = t0 / 86_400_000_000 * 86_400_000_000;
            let start = (t0 - day) as f64 / 3.6e9;
            let stop = (read::<u64>(&g, "mode_stop_us_exclusive")?[0] - day) as f64 / 3.6e9;
            let t = (start + stop) / 2.;
            if t < lo || t > hi {
                continue;
            }
            let v = read::<f64>(&g, "frequency_hz")?
                .iter()
                .map(|f| f * 299792458. / 94e6)
                .collect::<Vec<_>>();
            let range = read::<f64>(&g, "range_km")?;
            let nf = v.len();
            let nr = range.len();
            let source = File::open(String::from_utf8(read(&g, "source_file")?)?)?;
            let sg = source.group(&String::from_utf8(read(&g, "source_group")?)?)?;
            let pairs = sg
                .dataset(if sg.link_exists("ch_pairs") {
                    "ch_pairs"
                } else {
                    "pairs"
                })?
                .read_2d::<u64>()?;
            let pairs: Vec<_> = pairs
                .outer_iter()
                .map(|p| [p[0] as usize, p[1] as usize])
                .collect();
            let autos = pansy_legacy_capon::validate_pairs(&pairs)?;
            let mut power = vec![0.; nf * nr];
            for channel in autos {
                let z: Array2<Z> = sg.dataset("xc_arr")?.read_slice(s![channel, 0, .., ..])?;
                ensure!(z.shape() == [nf, nr], "shape mismatch");
                for (p, z) in power.iter_mut().zip(z.iter()) {
                    *p += z.r as f64;
                }
            }
            let noise = read::<f64>(&g, "estimated_total_noise_power")?;
            let offset = read::<u64>(&g, "noise_range_start_index")?[0] as usize;
            ensure!(offset == 0 && noise.len() == nr, "noise dimensions");
            let db = power
                .iter()
                .enumerate()
                .map(|(i, p)| (10. * (p / noise[i % nr]).log10()) as f32)
                .collect();
            blocks.push(Block {
                name,
                t,
                start,
                stop,
                range,
                v,
                db,
            });
        }
        ensure!(!blocks.is_empty(), "empty time window");
        let mut best = (f64::NEG_INFINITY, 0, 0);
        for (m, b) in blocks.iter().enumerate() {
            let nr = b.range.len();
            for (r, &range) in b.range.iter().enumerate() {
                if !(80. ..=92.).contains(&range) {
                    continue;
                }
                let score =
                    b.v.iter()
                        .enumerate()
                        .filter(|(_, v)| v.abs() > 15. && v.abs() < 70.)
                        .map(|(f, _)| (10_f64.powf(b.db[f * nr + r] as f64 / 10.) - 1.).max(0.))
                        .sum::<f64>();
                if score > best.0 {
                    best = (score, m, r);
                }
            }
        }
        let (_, m, r) = best;
        let b = &blocks[m];
        let nf = b.v.len();
        let nr = b.range.len();
        let g = out.create_group(&format!("episode{w}"))?;
        write(&g, "window_hours", &[lo, hi])?;
        write(&g, "selected_mode_hours", &[b.t])?;
        write(&g, "selected_range_km", &[b.range[r]])?;
        write(&g, "range_km", &b.range)?;
        write(&g, "velocity_mps", &b.v)?;
        write(&g, "snapshot_db_frequency_range", &b.db)?;
        write(
            &g,
            "mode_start_hours",
            &blocks.iter().map(|b| b.start).collect::<Vec<_>>(),
        )?;
        write(
            &g,
            "mode_stop_hours",
            &blocks.iter().map(|b| b.stop).collect::<Vec<_>>(),
        )?;
        let mut td = Vec::new();
        for other in &blocks {
            ensure!(other.range == b.range && other.v == b.v, "axes changed");
            for f in 0..nf {
                td.push(other.db[f * nr + r]);
            }
        }
        write(&g, "time_doppler_db_mode_frequency", &td)?;
        let original = input.group(&format!("blocks/{}", b.name))?;
        for n in ["frequency_index", "range_index"] {
            write(&g, n, &read::<u32>(&original, n)?)?;
        }
        for n in [
            "total_noise_db",
            "capon_theta_deg",
            "ghost_predicted_excess_fraction",
            "ghost_spatial_similarity",
        ] {
            write(&g, n, &read::<f64>(&original, n)?)?;
        }
        write(&g, "decision", &read::<u8>(&original, "decision")?)?;
        write(&g, "selected_range_index", &[r as u32])?;
        write(&g, "source_file", &read::<u8>(&original, "source_file")?)?;
        write(&g, "source_group", &read::<u8>(&original, "source_group")?)?;
        let fi = read::<u32>(&original, "frequency_index")?;
        let ri = read::<u32>(&original, "range_index")?;
        let theta = read::<f64>(&original, "capon_theta_deg")?;
        let flags = read::<u8>(&original, "decision")?;
        println!(
            "EPISODE {w} time={:.6} range={:.3} block={} modes={} score={:.2}",
            b.t,
            b.range[r],
            b.name,
            blocks.len(),
            best.0
        );
        for sign in [-1., 1.] {
            let mut max = 0_f64;
            let mut n = 0;
            let mut rej = [0; 5];
            for i in 0..fi.len() {
                let v = b.v[fi[i] as usize];
                if ri[i] as usize == r && v * sign > 15. && v * sign < 70. {
                    if flags[i] < 5 {
                        rej[flags[i] as usize] += 1;
                    }
                    if flags[i] == 0 && theta[i] <= 2. {
                        max = max.max(v.abs());
                        n += 1;
                    }
                }
            }
            println!(" WING sign={sign} retained_bins={n} max_abs_velocity={max:.2} flags={rej:?}");
        }
    }
    out.flush()?;
    drop(out);
    std::fs::rename(partial, &a[2])?;
    println!("OUTPUT bytes={}", std::fs::metadata(&a[2])?.len());
    Ok(())
}
