//! Small HDF5 plotting projection. Never copies dense source spectra.
use anyhow::{Context, Result};
use hdf5::{File, Group, H5Type};
use std::collections::BTreeMap;
fn read<T: H5Type>(g: &Group, n: &str) -> Result<Vec<T>> {
    Ok(g.dataset(n)?.read_raw()?)
}
fn write<T: H5Type>(g: &Group, n: &str, v: &[T]) -> Result<()> {
    if v.len() < 256 {
        g.new_dataset_builder().with_data(v).create(n)?;
    } else {
        g.new_dataset_builder()
            .with_data(v)
            .shuffle()
            .deflate(4)
            .create(n)?;
    }
    Ok(())
}
struct Block {
    time: f64,
    frequency: Vec<f64>,
    fi: Vec<u32>,
    snr: Vec<f64>,
    range: Vec<f64>,
    ri: Vec<u32>,
}
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    anyhow::ensure!(
        a.len() == 3,
        "usage: pansy-capon-plot-data INPUT.h5 SMALL_OUTPUT.h5"
    );
    let input = File::open(&a[1])?;
    let destination = std::path::Path::new(&a[2]);
    anyhow::ensure!(!destination.exists(), "plot projection exists");
    let partial = destination.with_extension("h5.partial");
    let output = File::create_excl(&partial)?;
    let mut times = Vec::new();
    let mut velocities = Vec::new();
    let mut snrs = Vec::new();
    let mut decisions = Vec::new();
    let mut blocks = BTreeMap::new();
    let mut midnight = None;
    for name in input.group("blocks")?.member_names()? {
        let g = input.group(&format!("blocks/{name}"))?;
        if !g.link_exists("decision") {
            continue;
        }
        let t0 = read::<u64>(&g, "mode_start_us")?[0];
        let t1 = read::<u64>(&g, "mode_stop_us_exclusive")?[0];
        let day = *midnight.get_or_insert(t0 / 86_400_000_000 * 86_400_000_000);
        let b = Block {
            time: ((t0 as f64 + t1 as f64) / 2. - day as f64) / 3.6e9,
            frequency: read(&g, "frequency_hz")?,
            fi: read(&g, "frequency_index")?,
            snr: read(&g, "total_noise_db")?,
            range: read(&g, "range_km")?,
            ri: read(&g, "range_index")?,
        };
        let decision = read::<u8>(&g, "decision")?;
        let mut bins: BTreeMap<(u8, i32), f64> = BTreeMap::new();
        for (i, &d) in decision.iter().enumerate() {
            if d != 0 && d != 4 {
                continue;
            }
            let v = b.frequency[b.fi[i] as usize] * 299_792_458. / (2. * 47e6);
            let bin = (v * 2.).round() as i32;
            bins.entry((d, bin))
                .and_modify(|x| *x = x.max(b.snr[i]))
                .or_insert(b.snr[i]);
        }
        for ((d, v), s) in bins {
            times.push(b.time);
            velocities.push(v as f64 / 2.);
            snrs.push(s);
            decisions.push(d);
        }
        blocks.insert(name.parse::<u64>()?, b);
    }
    let day = midnight.context("empty catalogue")?;
    let date = chrono::DateTime::from_timestamp_micros(day as i64)
        .context("bad UTC")?
        .format("%Y-%m-%d")
        .to_string();
    write(&output, "date_utc", date.as_bytes())?;
    let o = output.create_group("overview")?;
    write(&o, "time_hours", &times)?;
    write(&o, "velocity_mps", &velocities)?;
    write(&o, "peak_total_noise_db", &snrs)?;
    write(&o, "decision", &decisions)?;
    let mut events = Vec::new();
    for name in input.group("events")?.member_names()? {
        let g = input.group(&format!("events/{name}"))?;
        events.push((
            read::<u64>(&g, "mode_count")?[0],
            read::<u64>(&g, "cell_count")?[0],
            name,
        ));
    }
    events.sort_by(|a, b| b.cmp(a));
    let eroot = output.create_group("events")?;
    for (_, _, name) in events.into_iter().take(6) {
        let g = input.group(&format!("events/{name}"))?;
        let e = eroot.create_group(&name)?;
        let bi = read::<u64>(&g, "member_block_index")?;
        let ci = read::<u64>(&g, "member_candidate_index")?;
        let mut bins: BTreeMap<(u64, i32), (f64, f64, f64)> = BTreeMap::new();
        for (&block, &cell) in bi.iter().zip(&ci) {
            let b = &blocks[&block];
            let i = cell as usize;
            let v = b.frequency[b.fi[i] as usize] * 299_792_458. / (2. * 47e6);
            let entry = (b.time, b.snr[i], b.range[b.ri[i] as usize]);
            bins.entry((block, (v * 2.).round() as i32))
                .and_modify(|x| {
                    if entry.1 > x.1 {
                        *x = entry;
                    }
                })
                .or_insert(entry);
        }
        let mut t = Vec::new();
        let mut v = Vec::new();
        let mut s = Vec::new();
        let mut r = Vec::new();
        for ((_, vel), (time, snr, range)) in bins {
            t.push(time);
            v.push(vel as f64 / 2.);
            s.push(snr);
            r.push(range);
        }
        write(&e, "time_hours", &t)?;
        write(&e, "velocity_mps", &v)?;
        write(&e, "peak_total_noise_db", &s)?;
        write(&e, "range_km", &r)?;
        for n in ["mode_count", "cell_count"] {
            write(&e, n, &read::<u64>(&g, n)?)?;
        }
    }
    output.flush()?;
    drop(eroot);
    drop(o);
    drop(output);
    std::fs::rename(partial, destination)?;
    println!(
        "small_plot_projection_bytes={}",
        std::fs::metadata(&a[2])?.len()
    );
    Ok(())
}
