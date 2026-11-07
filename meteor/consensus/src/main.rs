use anyhow::{bail, ensure, Context, Result};
use pansy_meteor_consensus::{split, Config, Point, Track};
use rdigital_rf::metadata::{Numeric, Reader};
use std::path::Path;

fn write_output(
    path: &Path,
    start: u64,
    end: u64,
    points: &[Point],
    tracks: &[Track],
    read_ms: f64,
    fit_ms: f64,
) -> Result<()> {
    ensure!(!path.exists(), "refusing to overwrite {}", path.display());
    let temporary = path.with_extension(format!("{}.tmp", std::process::id()));
    let result = (|| -> Result<()> {
        let file = hdf5::File::create(&temporary)?;
        let mut membership = vec![-1_i32; points.len()];
        for (id, track) in tracks.iter().enumerate() {
            for &member in &track.members {
                membership[member] = i32::try_from(id)?;
            }
        }
        let samples: Vec<i64> = points
            .iter()
            .map(|p| i64::try_from(p.sample_idx))
            .collect::<std::result::Result<_, _>>()?;
        let epochs: Vec<i64> = tracks
            .iter()
            .map(|t| i64::try_from(t.epoch_sample))
            .collect::<std::result::Result<_, _>>()?;
        file.new_dataset_builder()
            .with_data(&[1_i32])
            .create("schema_version")?;
        file.new_dataset_builder()
            .with_data(&[i64::try_from(start)?, i64::try_from(end)?])
            .create("window_samples")?;
        file.new_dataset_builder()
            .with_data(&samples)
            .create("sample_idx")?;
        file.new_dataset_builder()
            .with_data(&membership)
            .create("track_id")?;
        file.new_dataset_builder()
            .with_data(&epochs)
            .create("track_epoch_sample")?;
        for (name, values) in [
            (
                "range_km",
                points.iter().map(|p| p.range_km as f32).collect::<Vec<_>>(),
            ),
            (
                "radial_velocity_m_s",
                points
                    .iter()
                    .map(|p| (p.radial_velocity_km_s * 1000.) as f32)
                    .collect(),
            ),
            (
                "snr_db",
                points
                    .iter()
                    .map(|p| (10. * p.snr_linear.log10()) as f32)
                    .collect(),
            ),
            (
                "track_range_at_epoch_km",
                tracks.iter().map(|t| t.range_at_epoch_km as f32).collect(),
            ),
            (
                "track_radial_velocity_m_s",
                tracks
                    .iter()
                    .map(|t| (t.radial_velocity_km_s * 1000.) as f32)
                    .collect(),
            ),
            (
                "track_radial_acceleration_m_s2",
                tracks
                    .iter()
                    .map(|t| (t.radial_acceleration_km_s2 * 1000.) as f32)
                    .collect(),
            ),
            ("read_fit_ms", vec![read_ms as f32, fit_ms as f32]),
        ] {
            file.new_dataset_builder().with_data(&values).create(name)?;
        }
        file.flush()?;
        drop(file);
        std::fs::rename(&temporary, path)?;
        Ok(())
    })();
    if temporary.exists() {
        let _ = std::fs::remove_file(temporary);
    }
    result
}

fn floats(v: &Numeric) -> Result<Vec<f64>> {
    match v {
        Numeric::F32(x) => Ok(x.iter().map(|&v| f64::from(v)).collect()),
        Numeric::F64(x) => Ok(x.clone()),
        _ => bail!("range, Doppler and SNR must be floating-point arrays"),
    }
}
fn indices(v: &Numeric) -> Result<Vec<u64>> {
    match v {
        Numeric::U64(x) => Ok(x.clone()),
        Numeric::I64(x) => x.iter().map(|&v| Ok(u64::try_from(v)?)).collect(),
        Numeric::F64(x) => x
            .iter()
            .map(|&v| {
                ensure!(
                    v.is_finite()
                        && (0.0..9_007_199_254_740_992.0).contains(&v)
                        && v.fract() == 0.0,
                    "legacy floating pulse index is not an exact integer below 2^53"
                );
                Ok(v as u64)
            })
            .collect(),
        _ => bail!("pulse indices must be 64-bit integers"),
    }
}

fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    ensure!(
        a.len() >= 4,
        "usage: pansy-meteor-consensus METADATA_DIR START_SAMPLE END_EXCLUSIVE [--detections] [--output FILE.h5]"
    );
    let mut detections = false;
    let mut output = None;
    let mut options = a[4..].iter();
    while let Some(option) = options.next() {
        match option.as_str() {
            "--detections" => detections = true,
            "--output" => output = Some(options.next().context("missing output path")?),
            _ => bail!("unknown argument {option}"),
        }
    }
    let start: u64 = a[2].parse()?;
    let end: u64 = a[3].parse()?;
    let reader = Reader::open(&a[1])?;
    let (n, d) = reader.sample_rate();
    ensure!(
        end > start && u128::from(end - start) * u128::from(d) <= 10 * u128::from(n),
        "offline replay accepts at most 10 seconds per invocation"
    );
    let c = Config {
        sample_rate_hz: n as f64 / d as f64,
        ..Config::default()
    };
    let fields = if detections {
        vec!["tx_idx", "range", "doppler", "snr"]
    } else {
        ensure!(
            n == 1_000_000 * d,
            "PANSY MF batches assume 1 MHz sample indices"
        );
        vec!["tx_idxs", "max_range", "max_dopvel", "max_snr", "tx_pwr"]
    };
    let begin = std::time::Instant::now();
    // MF records contain 20 IPPs starting at the record key. Include the
    // preceding batch so a window boundary cannot drop its trailing pulses.
    let lower = if detections {
        start
    } else {
        start.saturating_sub(32_000)
    };
    let records = reader.read(lower..end, &fields)?;
    let mut points = vec![];
    let mut rejected_tx_records = 0;
    for record in records.values() {
        if !detections {
            let power = floats(&record["tx_pwr"].values)?;
            ensure!(!power.is_empty(), "empty TX-power array");
            if power.iter().any(|v| !v.is_finite() || *v <= 1e9) {
                rejected_tx_records += 1;
                continue;
            }
        }
        let t = indices(&record[fields[0]].values)?;
        let r = floats(&record[fields[1]].values)?;
        let v = floats(&record[fields[2]].values)?;
        let s = floats(&record[fields[3]].values)?;
        ensure!(
            t.len() == r.len() && r.len() == v.len() && v.len() == s.len(),
            "mismatched pulse-array lengths"
        );
        for i in 0..t.len() {
            if (start..end).contains(&t[i])
                && s[i].is_finite()
                && s[i] >= c.grow_snr_linear
                && r[i].is_finite()
                && v[i].is_finite()
            {
                points.push(Point {
                    sample_idx: t[i],
                    range_km: r[i],
                    radial_velocity_km_s: v[i] / 1000.,
                    snr_linear: s[i],
                });
            }
        }
    }
    let read_ms = begin.elapsed().as_secs_f64() * 1000.;
    let start_fit = std::time::Instant::now();
    let tracks = split(&points, &c).context("offline consensus")?;
    let fit_ms = start_fit.elapsed().as_secs_f64() * 1000.;
    if let Some(path) = output {
        write_output(
            Path::new(path),
            start,
            end,
            &points,
            &tracks,
            read_ms,
            fit_ms,
        )?;
    }
    println!(
        "records={} rejected_tx_records={rejected_tx_records} points={} tracks={} read_ms={read_ms:.3} fit_ms={fit_ms:.3}",
        records.len(),
        points.len(),
        tracks.len()
    );
    for (i, t) in tracks.iter().enumerate() {
        println!(
            "track={i} points={} first_sample={} last_sample={} radial_velocity_km_s={:.3} radial_acceleration_km_s2={:.3}",
            t.members.len(),
            points[t.members[0]].sample_idx,
            points[*t.members.last().unwrap()].sample_idx,
            t.radial_velocity_km_s,t.radial_acceleration_km_s2
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn legacy_float_timestamps_are_checked_not_rounded() {
        assert_eq!(
            indices(&Numeric::F64(vec![1_746_541_806_935_733.])).unwrap(),
            vec![1_746_541_806_935_733]
        );
        for x in [-1., 1.5, f64::NAN, 9_007_199_254_740_992.] {
            assert!(indices(&Numeric::F64(vec![x])).is_err());
        }
    }
}
