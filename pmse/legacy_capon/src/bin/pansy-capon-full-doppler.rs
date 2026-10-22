//! Full native Doppler context, NOT the selected candidate cells. No velocity,
//! angular or SNR cuts. A compact HDF5 projection, not a dense XC transfer.
use anyhow::{ensure, Context, Result};
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
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    ensure!(
        a.len() == 3,
        "usage: pansy-capon-full-doppler CATALOGUE.h5 SMALL_OUTPUT.h5"
    );
    let input = File::open(&a[1])?;
    let dest = std::path::Path::new(&a[2]);
    ensure!(!dest.exists(), "output exists");
    let partial = dest.with_extension("h5.partial");
    let out = File::create_excl(&partial)?;
    let mut frequency = Vec::new();
    let mut start = Vec::new();
    let mut stop = Vec::new();
    let mut values = Vec::<f32>::new();
    let mut midnight = None;
    for name in input.group("blocks")?.member_names()? {
        let g = input.group(&format!("blocks/{name}"))?;
        if !g.link_exists("frequency_hz") {
            continue;
        }
        let frequencies = read::<f64>(&g, "frequency_hz")?;
        if frequency.is_empty() {
            frequency = frequencies.clone();
        } else {
            ensure!(
                frequencies == frequency,
                "native frequency axis changed within day"
            );
        }
        let t0 = read::<u64>(&g, "mode_start_us")?[0];
        let t1 = read::<u64>(&g, "mode_stop_us_exclusive")?[0];
        let day = *midnight.get_or_insert(t0 / 86_400_000_000 * 86_400_000_000);
        start.push((t0 - day) as f64 / 3.6e9);
        stop.push((t1 - day) as f64 / 3.6e9);
        let path = String::from_utf8(read(&g, "source_file")?)?;
        let group = String::from_utf8(read(&g, "source_group")?)?;
        let source = File::open(path)?;
        let sg = source.group(&group)?;
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
        let ranges = read::<f64>(&g, "range_km")?;
        let noise = read::<f64>(&g, "estimated_total_noise_power")?;
        let offset = read::<u64>(&g, "noise_range_start_index")?[0] as usize;
        let selected: Vec<_> = ranges
            .iter()
            .enumerate()
            .filter_map(|(i, &r)| (80. ..=92.).contains(&r).then_some(i))
            .collect();
        ensure!(!selected.is_empty(), "80--92 km range interval absent");
        let nf = frequency.len();
        let nr = ranges.len();
        let mut power = vec![0f64; nf * nr];
        for channel in autos {
            let z: Array2<Z> = sg.dataset("xc_arr")?.read_slice(s![channel, 0, .., ..])?;
            ensure!(
                z.shape() == [nf, nr],
                "source shape disagrees with catalogue axes"
            );
            for (p, z) in power.iter_mut().zip(z.iter()) {
                *p += z.r as f64;
            }
        }
        for f in 0..nf {
            let mut maximum = f64::NEG_INFINITY;
            for &r in &selected {
                ensure!(
                    r >= offset && r - offset < noise.len(),
                    "noise range unavailable"
                );
                let ratio = power[f * nr + r] / noise[r - offset];
                ensure!(
                    ratio.is_finite() && ratio > 0.,
                    "invalid spectral power/noise"
                );
                maximum = maximum.max(10. * ratio.log10());
            }
            values.push(maximum as f32);
        }
        println!(
            "FULL_DOPPLER {name} bins={nf} min_hz={} max_hz={}",
            frequency[0],
            frequency[nf - 1]
        );
    }
    let day = midnight.context("empty catalogue")?;
    let date = chrono::DateTime::from_timestamp_micros(day as i64)
        .context("bad UTC")?
        .format("%Y-%m-%d")
        .to_string();
    write(&out, "date_utc", date.as_bytes())?;
    write(&out, "frequency_hz", &frequency)?;
    write(
        &out,
        "velocity_mps",
        &frequency
            .iter()
            .map(|f| f * 299_792_458. / (2. * 47e6))
            .collect::<Vec<_>>(),
    )?;
    write(&out, "mode_start_hours", &start)?;
    write(&out, "mode_stop_hours", &stop)?;
    write(&out, "peak_total_noise_db_mode_frequency", &values)?;
    write(&out,"definitions",b"ALL stored Dopplers; NO candidate, velocity, angular or SNR cuts; seven-channel summed auto power; maximum total/noise over zenith 80--92 km; source is legacy 8 ms per-beam data; original source unchanged")?;
    let eroot = out.create_group("events")?;
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
    for (_, _, name) in events.into_iter().take(6) {
        let g = input.group(&format!("events/{name}"))?;
        let e = eroot.create_group(&name)?;
        write(
            &e,
            "start_hours",
            &[(read::<u64>(&g, "start_us")?[0] - day) as f64 / 3.6e9],
        )?;
        write(
            &e,
            "stop_hours",
            &[(read::<u64>(&g, "stop_us_exclusive")?[0] - day) as f64 / 3.6e9],
        )?;
        write(&e, "mode_count", &read::<u64>(&g, "mode_count")?)?;
    }
    out.flush()?;
    drop(eroot);
    drop(out);
    std::fs::rename(partial, dest)?;
    println!(
        "COMPLETE modes={} all_native_doppler_bins={} bytes={}",
        start.len(),
        frequency.len(),
        std::fs::metadata(dest)?.len()
    );
    Ok(())
}
