//! Export only small plot projections as completed days arrive. Run under
//! the shared <=one-core quota documented by its userspace service.
use anyhow::{ensure, Result};
use hdf5::File;
use std::{path::PathBuf, process::Command, time::Duration};
fn main() -> Result<()> {
    let root = PathBuf::from(
        std::env::args()
            .nth(1)
            .expect("usage: pansy-capon-plot-watch CATALOGUE_DIR"),
    );
    let output = root.join("plots");
    std::fs::create_dir_all(&output)?;
    let exporter = std::env::current_exe()?.with_file_name("pansy-capon-plot-data");
    loop {
        let index = File::open(root.join("sweep-index.h5"))?;
        let days = index.dataset("utc_day_yyyymmdd")?.read_raw::<u32>()?;
        let status = index.dataset("status")?.read_raw::<u8>()?;
        drop(index);
        for (&day, &state) in days.iter().zip(&status) {
            if state != 1 {
                continue;
            }
            let date = format!("{:04}-{:02}-{:02}", day / 10000, day / 100 % 100, day % 100);
            let dest = output.join(format!("{date}.plot.h5"));
            if dest.exists() {
                continue;
            }
            let result = Command::new(&exporter)
                .arg(root.join(format!("{date}.h5")))
                .arg(&dest)
                .status()?;
            ensure!(
                result.success(),
                "plot projection failed for {date}; no automatic retries"
            );
        }
        if !status.contains(&0) {
            ensure!(!status.contains(&2), "completed sweep contains failed days");
            break;
        }
        std::thread::sleep(Duration::from_secs(30));
    }
    println!("PLOT_PROJECTIONS_COMPLETE");
    Ok(())
}
