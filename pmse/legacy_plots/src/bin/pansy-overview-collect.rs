//! Finite all-date plot collection. Spectral analysis stays on Syowa.
use anyhow::{Result, ensure};
use hdf5::File;
use std::{fs, path::PathBuf, process::Command, time::Duration};
const SSH: &str = "ssh -o BatchMode=yes -o ConnectTimeout=60 -o ControlMaster=auto -o ControlPersist=600 -o ControlPath=/home/j/.ssh/rvd-zenith-plot-%C -J j@juha.no -p3131";
fn pull(name: &str, dest: &std::path::Path) -> Result<()> {
    let status = Command::new("rsync")
        .args([
            "-az",
            "--timeout=120",
            "--max-size=5m",
            "--no-perms",
            "--no-owner",
            "--no-group",
            "-e",
            SSH,
        ])
        .arg(format!(
            "radar@localhost:/media/analysis/metadata/zenith_peak_overviews_v1/{name}"
        ))
        .arg(dest)
        .status()?;
    ensure!(status.success(), "transfer failed: {name}");
    Ok(())
}
fn main() -> Result<()> {
    let out = PathBuf::from(
        std::env::args()
            .nth(1)
            .expect("usage: pansy-overview-collect OUTPUT_DIR"),
    );
    fs::create_dir_all(&out)?;
    let renderer = std::env::current_exe()?.with_file_name("pansy-full-doppler");
    loop {
        pull("overview-index.h5", &out)?;
        let f = File::open(out.join("overview-index.h5"))?;
        let days = f.dataset("utc_day_yyyymmdd")?.read_raw::<u32>()?;
        let status = f.dataset("status")?.read_raw::<u8>()?;
        let flags = f.dataset("flagged_integrations")?.read_raw::<u64>()?;
        drop(f);
        let pending = status.iter().filter(|&&v| v == 0).count();
        let failed = status.iter().filter(|&&v| v == 2).count();
        let mut html = format!(
            "<!doctype html><meta charset=utf-8><title>PANSY zenith main-peak flags</title><style>body{{font:17px system-ui;max-width:1200px;margin:2em auto;padding:0 1em;background:white;color:#172432}}img{{width:100%;height:auto}}article{{margin:2em 0}}a{{color:#075ba0}}</style><h1>PANSY zenith spectra</h1><p>Main peak within ±90 m/s; flagged when |v| ≥ 30 m/s. All stored Dopplers shown. Zenith, maximum total/noise over 80–92 km. No Capon or SNR selection. Flags are inspection markers, not validated winds.</p><p>{} dates; {} processed; {pending} pending; {failed} dates with errors.</p>",
            days.len(),
            days.len() - pending
        );
        let mut plots = 0;
        for i in 0..days.len() {
            if status[i] == 0 || flags[i] == 0 {
                continue;
            }
            let n = days[i];
            let date = format!("{:04}-{:02}-{:02}", n / 10000, (n / 100) % 100, n % 100);
            let data = format!("{date}.full-doppler.h5");
            let plot = format!("{date}-FULL-doppler-overview");
            if !out.join(&data).exists() {
                pull(&data, &out)?;
                ensure!(
                    out.join(&data).exists(),
                    "projection missing or exceeds transfer limit: {data}"
                );
            }
            if !out.join(format!("{plot}.png")).exists()
                || !out.join(format!("{plot}.pdf")).exists()
            {
                ensure!(
                    Command::new(&renderer)
                        .arg(out.join(&data))
                        .arg(&out)
                        .status()?
                        .success(),
                    "render failed: {date}"
                );
            }
            html.push_str(&format!("<article><h2>{date}: {} flagged integrations{}</h2><a href='{plot}.png'><img src='{plot}.png' loading='lazy' alt='{date} full-day Doppler spectrum'></a><p><a href='{plot}.pdf'>PDF</a></p></article>",flags[i],if status[i]==2{" (input errors recorded)"}else{""}));
            plots += 1;
        }
        fs::write(out.join("index.partial.html"), html)?;
        fs::rename(out.join("index.partial.html"), out.join("index.html"))?;
        println!(
            "PLOT_PROGRESS plots={plots} pending={pending} input_error_days={failed} total={}",
            days.len()
        );
        if pending == 0 {
            break;
        }
        std::thread::sleep(Duration::from_secs(30));
    }
    Ok(())
}
