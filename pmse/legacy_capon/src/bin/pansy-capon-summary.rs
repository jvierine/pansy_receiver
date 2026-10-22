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
    let args: Vec<_> = std::env::args().collect();
    let path = args
        .get(1)
        .context("usage: pansy-capon-summary RESULT.h5 [MAX_EVENT_SPEED_MPS]")?;
    let max_speed: f64 = args
        .get(2)
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(f64::INFINITY);
    anyhow::ensure!(
        max_speed > 0. && !max_speed.is_nan(),
        "invalid maximum event speed"
    );
    let f = File::open(path)?;
    let m = f.group("metadata")?;
    eprintln!("WARNING: Candidate events, not validated extreme winds. Parent-compatible exclusion is conservative, not exact decoding; missing/unmodelled parents can leave contamination.");
    let mut counts = [0u64; 5];
    let mut raw = 0u64;
    let mut uncertain = 0u64;
    for name in f.group("blocks")?.member_names()? {
        let g = f.group(&format!("blocks/{name}"))?;
        if !g.link_exists("decision") {
            continue;
        }
        let decisions = g.dataset("decision")?.read_raw::<u8>()?;
        if g.link_exists("ghost_parent_search_complete") {
            uncertain += decisions
                .iter()
                .zip(
                    g.dataset("ghost_parent_search_complete")?
                        .read_raw::<u8>()?,
                )
                .filter(|(d, c)| **d == 0 && *c == 0)
                .count() as u64;
        }
        for x in decisions {
            if x < 5 {
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
        "extreme_cells: retained={} off_zenith={} invalid_covariance={} boundary={} code_ghost_excluded={} total={} retained_incomplete_parent_coverage={uncertain}",
        counts[0],
        counts[1],
        counts[2],
        counts[3],
        counts[4],
        counts.iter().sum::<u64>()
    );
    println!(
        "grouped_candidates={} multiple_modes={multi} single_mode={}",
        rows.len(),
        rows.len() - multi
    );
    let below: Vec<_> = rows.iter().filter(|r| r.0 < max_speed).collect();
    println!(
        "events_with_maximum_speed_below_{max_speed}={} multiple_modes={}",
        below.len(),
        below.iter().filter(|r| r.4 > 1).count()
    );
    println!("Top 20 within the requested event-speed bound; candidates, not confirmed vertical drafts. Bounding an event's maximum is NOT cell-level reclassification:");
    for (_, id, start, stop, modes, cells, r0, r1, v, snr, theta) in
        rows.into_iter().filter(|r| r.0 < max_speed).take(20)
    {
        let utc = chrono::DateTime::from_timestamp_micros(start as i64).context("invalid UTC")?;
        println!("{id} {utc} span_s={:.2} modes={modes} cells={cells} range_km={r0:.2}..{r1:.2} radial_mps={v:+.2} peak_total_noise_db={snr:.2} max_theta_deg={theta:.2}",(stop-start)as f64/1e6);
    }
    Ok(())
}
