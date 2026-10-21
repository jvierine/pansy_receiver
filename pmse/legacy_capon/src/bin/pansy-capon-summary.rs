//! Read only compact result metadata, not source XC spectral arrays.
use anyhow::{Context, Result};
use hdf5::{File, Group};
fn value<T: hdf5::H5Type>(g: &Group, name: &str) -> Result<T> {
    g.dataset(name)?
        .read_raw::<T>()?
        .into_iter()
        .next()
        .context("empty summary field")
}
fn main() -> Result<()> {
    let path = std::env::args()
        .nth(1)
        .context("usage: pansy-capon-summary RESULT.h5")?;
    let f = File::open(path)?;
    let m = f.group("metadata")?;
    let mut counts = [0u64; 4];
    let mut raw = 0u64;
    for name in f.group("blocks")?.member_names()? {
        let g = f.group(&format!("blocks/{name}"))?;
        if !g.link_exists("decision") {
            continue;
        }
        for x in g.dataset("decision")?.read_raw::<u8>()? {
            if x < 4 {
                counts[x as usize] += 1;
            }
        }
        raw += g
            .dataset("raw_snr_cell_counts_by_frequency")?
            .read_raw::<u64>()?
            .iter()
            .sum::<u64>();
    }
    let root = f.group("events")?;
    let mut rows = Vec::new();
    let mut multi = 0usize;
    for name in root.member_names()? {
        let g = root.group(&name)?;
        let start = value::<u64>(&g, "start_us")?;
        let stop = value::<u64>(&g, "stop_us_exclusive")?;
        let modes = value::<u64>(&g, "mode_count")?;
        if modes > 1 {
            multi += 1;
        }
        rows.push((
            value::<f64>(&g, "maximum_abs_velocity_mps")?,
            name,
            start,
            stop,
            modes,
            value::<u64>(&g, "cell_count")?,
            value::<f64>(&g, "range_min_km")?,
            value::<f64>(&g, "range_max_km")?,
            value::<f64>(&g, "velocity_at_maximum_mps")?,
            value::<f64>(&g, "peak_total_noise_db")?,
            value::<f64>(&g, "maximum_theta_deg")?,
        ));
    }
    rows.sort_by(|a, b| b.0.total_cmp(&a.0));
    println!(
        "blocks={} input_errors={} seconds={:.2} raw_above_threshold_cells={raw}",
        value::<u64>(&m, "blocks_attempted")?,
        value::<u64>(&m, "input_errors")?,
        value::<f64>(&m, "elapsed_seconds")?
    );
    println!(
        "extreme_cells: accepted={} off_zenith={} invalid_covariance={} boundary={} total={}",
        counts[0],
        counts[1],
        counts[2],
        counts[3],
        counts.iter().sum::<u64>()
    );
    println!(
        "grouped_candidates={} multiple_modes={multi} single_mode={}",
        rows.len(),
        rows.len() - multi
    );
    println!("Top 20 by maximum absolute radial velocity; these are candidates, not confirmed vertical drafts:");
    for (_, id, start, stop, modes, cells, r0, r1, v, snr, theta) in rows.into_iter().take(20) {
        let utc = chrono::DateTime::from_timestamp_micros(start as i64).context("invalid UTC")?;
        println!("{id} {utc} span_s={:.2} modes={modes} cells={cells} range_km={r0:.2}..{r1:.2} radial_mps={v:+.2} peak_total_noise_db={snr:.2} max_theta_deg={theta:.2}",(stop-start)as f64/1e6);
    }
    Ok(())
}
