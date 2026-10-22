//! Sequential, restartable sweep. The userspace systemd cgroup limits the
//! runner and every child to one CPU core in total. No source files change.
use anyhow::{ensure, Context, Result};
use clap::Parser;
use hdf5::File;
use std::{collections::BTreeSet, fs, path::PathBuf, process::Command};
#[derive(Parser)]
struct Args {
    #[arg(long, required = true)]
    root: Vec<PathBuf>,
    #[arg(long)]
    reference: PathBuf,
    #[arg(long)]
    output: PathBuf,
    #[arg(long, default_value = "2026-01-30")]
    first_day: String,
}
fn index(out: &std::path::Path, days: &[String], status: &[u8]) -> Result<()> {
    let tmp = out.join("sweep-index.partial.h5");
    let f = File::create(&tmp)?;
    let dates: Vec<u32> = days
        .iter()
        .map(|d| d.replace('-', "").parse().unwrap())
        .collect();
    f.new_dataset_builder()
        .with_data(&dates)
        .create("utc_day_yyyymmdd")?;
    f.new_dataset_builder().with_data(status).create("status")?;
    f.new_dataset_builder().with_data(b"0 pending; 1 completed; 2 failed; source data unchanged; sequential one-core userspace service".as_slice()).create("definitions")?;
    f.flush()?;
    drop(f);
    fs::rename(tmp, out.join("sweep-index.h5"))?;
    Ok(())
}
fn main() -> Result<()> {
    let a = Args::parse();
    fs::create_dir_all(&a.output)?;
    let mut unique = BTreeSet::new();
    for root in &a.root {
        for entry in fs::read_dir(root)? {
            let entry = entry?;
            let name = entry.file_name().to_string_lossy().into_owned();
            if name.len() < 11 || name.as_bytes()[10] != b'T' || !entry.file_type()?.is_dir() {
                continue;
            }
            let day = &name[..10];
            if chrono::NaiveDate::parse_from_str(day, "%Y-%m-%d").is_err() {
                continue;
            }
            let has_data = fs::read_dir(entry.path())?.filter_map(|e| e.ok()).any(|e| {
                let n = e.file_name().to_string_lossy().into_owned();
                n.starts_with("xc@") && n.ends_with(".h5")
            });
            if has_data {
                unique.insert(day.to_string());
            }
        }
    }
    let mut days: Vec<_> = unique.into_iter().collect();
    if let Some(i) = days.iter().position(|d| *d == a.first_day) {
        let first = days.remove(i);
        days.insert(0, first);
    }
    ensure!(!days.is_empty(), "no populated XC days");
    let mut status = vec![0u8; days.len()];
    index(&a.output, &days, &status)?;
    let processor = std::env::current_exe()?.with_file_name("pansy-legacy-capon");
    for (i, day) in days.iter().enumerate() {
        let output = a.output.join(format!("{day}.h5"));
        if output.exists() {
            let f = File::open(&output)?;
            let schema = f.dataset("metadata/schema")?.read_raw::<u8>()?;
            ensure!(
                schema == b"pansy-legacy-extreme-capon-v2",
                "refusing old schema output {}",
                output.display()
            );
            status[i] = if f.dataset("metadata/input_errors")?.read_raw::<u64>()?[0] == 0 {
                1
            } else {
                2
            };
        } else {
            println!("SWEEP_START {}/{} {day}", i + 1, days.len());
            let mut cmd = Command::new(&processor);
            for root in &a.root {
                cmd.arg("--root").arg(root);
            }
            let result = cmd
                .arg("--day")
                .arg(day)
                .arg("--reference")
                .arg(&a.reference)
                .arg("--output")
                .arg(&output)
                .status()
                .context("starting native processor")?;
            status[i] = if result.success() { 1 } else { 2 };
        }
        index(&a.output, &days, &status)?;
        println!(
            "SWEEP_DAY {day} status={} completed={} failed={} total={}",
            status[i],
            status.iter().filter(|&&s| s == 1).count(),
            status.iter().filter(|&&s| s == 2).count(),
            days.len()
        );
    }
    ensure!(
        !status.contains(&2),
        "sweep finished with failed days; inspect index and journal"
    );
    println!("SWEEP_COMPLETE days={}", days.len());
    Ok(())
}
