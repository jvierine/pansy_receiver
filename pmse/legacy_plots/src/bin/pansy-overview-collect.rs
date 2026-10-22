//! Finite all-date plot collection. Spectral analysis stays on Syowa.
use anyhow::{Result, ensure};
use hdf5::File;
use std::{fs, path::PathBuf, process::Command, time::Duration};
const SSH: &str = "ssh -o BatchMode=yes -o ConnectTimeout=60 -o ControlMaster=auto -o ControlPersist=600 -o ControlPath=/home/j/.ssh/rvd-zenith-plot-%C -J j@juha.no -p3131";
// Transfers run sequentially. Reserve headroom below 10,000 byte/s for the
// existing low-rate archive transfers and SSH overhead.
const TRANSFER_LIMIT: &str = "--bwlimit=7";
fn koki_bytes(out:&std::path::Path)->Result<u64> {
    let mut bytes=0;
    for e in fs::read_dir(out)? {let e=e?;let name=e.file_name().to_string_lossy().into_owned();
        if name.ends_with(".koki.h5") || name.contains("-XC2-KOKI-") {bytes+=e.metadata()?.len();}
    } Ok(bytes)
}
fn pull(name: &str, dest: &std::path::Path) -> Result<()> {
    pull_from("zenith_peak_overviews_v1", name, dest)
}
fn pull_from(folder: &str, name: &str, dest: &std::path::Path) -> Result<()> {
    let local=match folder {
        "xc2_koki_v1"=>std::env::var_os("PANSY_KOKI_LOCAL"),
        "zenith_peak_overviews_v1"=>std::env::var_os("PANSY_OVERVIEWS_LOCAL"),
        _=>None,
    };
    if let Some(root)=local {
        let src=PathBuf::from(root).join(name);
        ensure!(fs::metadata(&src)?.len()<=5*1024*1024,"oversized local summary");
        fs::copy(src,dest.join(name))?;return Ok(());
    }
    let status = Command::new("rsync")
        .args([
            "-az",
            TRANSFER_LIMIT,
            "--timeout=120",
            "--max-size=5m",
            "--no-perms",
            "--no-owner",
            "--no-group",
            "-e",
            SSH,
        ])
        .arg(format!(
            "radar@localhost:/media/analysis/metadata/{folder}/{name}"
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
    let publish = std::env::args().any(|v| v == "--publish");
    let capon = std::env::args().any(|v| v == "--capon");
    let koki = std::env::args().any(|v| v == "--koki");
    let renderer = std::env::current_exe()?.with_file_name("pansy-full-doppler");
    let local=std::env::var_os("PANSY_KOKI_LOCAL").is_some();
    loop {
        let (days,status,flags)=if local {(vec![],vec![],vec![])} else {
            pull("overview-index.h5", &out)?;
            let f=File::open(out.join("overview-index.h5"))?;
            (f.dataset("utc_day_yyyymmdd")?.read_raw::<u32>()?,f.dataset("status")?.read_raw::<u8>()?,f.dataset("flagged_integrations")?.read_raw::<u64>()?)
        };
        let pending = status.iter().filter(|&&v| v == 0).count();
        let failed = status.iter().filter(|&&v| v == 2).count();
        let mut html = format!(
            "<!doctype html><meta charset=utf-8><title>PANSY zenith main-peak flags</title><style>body{{font:17px system-ui;max-width:1200px;margin:2em auto;padding:0 1em;background:white;color:#172432}}img{{width:100%;height:auto}}article{{margin:2em 0}}a{{color:#075ba0}}</style><h1>PANSY zenith spectra</h1><p>Main peak within ±90 m/s; flagged when |v| ≥ 30 m/s. All stored Dopplers shown. Zenith, maximum total/noise over 80–92 km.</p><p>{} dates; {} processed; {pending} pending; {failed} dates with errors.</p>",
            days.len(),
            days.len() - pending
        );
        if local {html="<!doctype html><meta charset=utf-8><title>PANSY PMSE-season XC2 events</title><style>body{font:17px system-ui;max-width:1200px;margin:2em auto;padding:0 1em;background:white;color:#172432}img{width:100%;height:auto}article{margin:2em 0}a{color:#075ba0}</style><h1>PANSY PMSE-season XC2 events</h1><p>Austral summer: November–February. Full stored Doppler spectra, event context plots and velocity histograms.</p>".to_owned();}
        let mut plots = 0;
        let mut new_raw = 0;
        let mut raw_deferred = 0;
        let mut capon_pending = 0;
        let mut koki_running = false;
        if koki {
            pull_from("xc2_koki_v1", "koki-index.h5", &out)?;
            let k = File::open(out.join("koki-index.h5"))?;
            let kd = k.dataset("utc_day_yyyymmdd")?.read_raw::<u32>()?;
            let ks = k.dataset("status")?.read_raw::<u8>()?;
            let kn = k.dataset("extreme_cell_count")?.read_raw::<u64>()?;
            koki_running = k.dataset("sweep_complete")?.read_raw::<u8>()?[0] == 0;
            html.push_str(&format!("<h2>XC2 bin-by-bin leakage and Capon screening</h2><p>Zenith, 80–92 km; total/noise ≥5 dB. Peaks above 30 dB trigger the 3σ primary-cluster mask. Capon acceptance ≤2°; extreme bins |v| ≥15 m/s are marked in cyan. Covariance normalization and PSD projection follow the supplied procedure. Full stored Doppler range remains visible.</p><p>{} dates processed; sweep {}. {} dates with source errors. Histograms count spectral cells, not independent events.</p>",kd.len(),if koki_running{"running"}else{"complete"},ks.iter().filter(|&&s|s==2).count()));
            for i in 0..kd.len() {
                let n = kd[i];
                let date = format!("{:04}-{:02}-{:02}", n / 10000, (n / 100) % 100, n % 100);
                let data = format!("{date}.full-doppler.h5");
                let sidecar = format!("{date}.koki.h5");
                let plot = format!("{date}-XC2-KOKI-overview-v1");
                let hist = format!("{date}-XC2-KOKI-histogram-v1");
                // Along with the 400 MB station limit, bound collected results
                // and their static web mirror well below 1,000 MB combined.
                ensure!(koki_bytes(&out)?+8_000_000<=200_000_000,"Koki collection budget reached; retain results at source, no further transfer or rendering");
                // A hard-error day may not have produced a transferable summary.
                if ks[i] == 2 && kn[i] == 0 {
                    continue;
                }
                if !out.join(&data).exists() {
                    pull(&data, &out)?;
                }
                if !out.join(&sidecar).exists() {
                    pull_from("xc2_koki_v1", &sidecar, &out)?;
                }
                if !out.join(format!("{plot}.png")).exists()
                    || !out.join(format!("{hist}.png")).exists()
                    || !out.join(format!("{date}-XC2-KOKI-events-v1.html")).exists()
                {
                    ensure!(
                        Command::new(&renderer)
                            .arg(out.join(&data))
                            .arg(&out)
                            .arg("--koki")
                            .arg(out.join(&sidecar))
                            .status()?
                            .success(),
                        "binwise render failed: {date}"
                    );
                }
                html.push_str(&format!("<article><h3>{date}: {} accepted extreme spectral cells</h3><a href='{plot}.png'><img src='{plot}.png' loading='lazy' alt='{date} binwise Capon detections'></a><p><a href='{date}-XC2-KOKI-events-v1.html'>Individual event plots</a> · <a href='{plot}.pdf'>Spectrum PDF</a> · <a href='{hist}.png'>Velocity histograms</a> · <a href='{hist}.pdf'>Histogram PDF</a></p></article>",kn[i]));
            }
        }
        if capon {
            pull_from("xc2_peak_capon_v2", "capon-index.h5", &out)?;
            let c = File::open(out.join("capon-index.h5"))?;
            let cd = c.dataset("utc_day_yyyymmdd")?.read_raw::<u32>()?;
            let cs = c.dataset("status")?.read_raw::<u8>()?;
            let cn = c.dataset("detection_count")?.read_raw::<u64>()?;
            capon_pending = cs.iter().filter(|&&s| s == 0).count();
            html.push_str(&format!("<h2>XC2 Capon detections</h2><p>Main peak |v| ≥ 30 m/s within ±90 m/s; total/noise ≥ 5 dB; arrival within 2° of zenith. Diagonal loading: 1%. Cyan diamonds mark detections; grey crosses show other extreme main peaks.</p><p>{} XC2 dates; {} processed; {capon_pending} pending; {} dates with errors.</p>",cd.len(),cd.len()-capon_pending,cs.iter().filter(|&&s|s==2).count()));
            for i in 0..cd.len() {
                if cs[i] == 0 || cn[i] == 0 {
                    continue;
                }
                let n = cd[i];
                let date = format!("{:04}-{:02}-{:02}", n / 10000, (n / 100) % 100, n % 100);
                let data = format!("{date}.full-doppler.h5");
                let sidecar = format!("{date}.capon-v2.h5");
                let plot = format!("{date}-XC2-CAPON-overview-v2");
                if !out.join(&data).exists() {
                    pull(&data, &out)?;
                }
                if !out.join(&sidecar).exists() {
                    pull_from("xc2_peak_capon_v2", &sidecar, &out)?;
                }
                if !out.join(format!("{plot}.png")).exists()
                    || !out.join(format!("{plot}.pdf")).exists()
                {
                    ensure!(
                        Command::new(&renderer)
                            .arg(out.join(&data))
                            .arg(&out)
                            .arg("--capon")
                            .arg(out.join(&sidecar))
                            .status()?
                            .success(),
                        "Capon render failed: {date}"
                    );
                }
                html.push_str(&format!("<article><h3>{date}: {} Capon detections</h3><a href='{plot}.png'><img src='{plot}.png' loading='lazy' alt='{date} XC2 Capon detections'></a><p><a href='{plot}.pdf'>PDF</a></p></article>",cn[i]));
            }
            html.push_str("<h2>All-data main-peak flags</h2>");
        }
        for i in 0..days.len() {
            if status[i] == 0 || flags[i] == 0 {
                continue;
            }
            let n = days[i];
            let date = format!("{:04}-{:02}-{:02}", n / 10000, (n / 100) % 100, n % 100);
            let data = format!("{date}.full-doppler.h5");
            let plot = format!("{date}-FULL-doppler-overview");
            let needs_render = !out.join(format!("{plot}.png")).exists()
                || !out.join(format!("{plot}.pdf")).exists();
            if capon && needs_render && new_raw >= 2 {
                raw_deferred += 1;
                continue;
            }
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
                new_raw += 1;
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
        if publish {
            ensure!(
                Command::new("ssh")
                    .args([
                        "-o",
                        "BatchMode=yes",
                        "-o",
                        "ConnectTimeout=30",
                        "j@juha.no",
                        "mkdir -p /mnt/shovel/share/pansy-zenith-main-peaks"
                    ])
                    .status()?
                    .success(),
                "creating plot destination failed"
            );
            ensure!(
                Command::new("rsync")
                    .args([
                        "-r",
                        TRANSFER_LIMIT,
                        "--size-only",
                        "--no-times",
                        "--omit-dir-times",
                        "--timeout=120",
                        "--no-perms",
                        "--no-owner",
                        "--no-group",
                        "--include=/*-FULL-doppler-overview.png",
                        "--include=/*-FULL-doppler-overview.pdf",
                        "--include=/*-XC2-CAPON-overview-v2.png",
                        "--include=/*-XC2-CAPON-overview-v2.pdf",
                        "--include=/*-XC2-KOKI-*-v1.png",
                        "--include=/*-XC2-KOKI-*-v1.pdf",
                        "--include=/*-XC2-KOKI-events-v1.html",
                        "--exclude=*",
                        "-e",
                        "ssh -o BatchMode=yes -o ConnectTimeout=30"
                    ])
                    .arg(format!("{}/", out.display()))
                    .arg("j@juha.no:/mnt/shovel/share/pansy-zenith-main-peaks/")
                    .status()?
                    .success(),
                "publishing static plots failed"
            );
            // The destination does not support timestamp preservation. Plot
            // files are immutable; the small changing index needs content checks.
            ensure!(
                Command::new("rsync")
                    .args([
                        "-r",
                        TRANSFER_LIMIT,
                        "--checksum",
                        "--no-times",
                        "--omit-dir-times",
                        "--no-perms",
                        "--no-owner",
                        "--no-group",
                        "--timeout=120",
                        "-e",
                        "ssh -o BatchMode=yes -o ConnectTimeout=30"
                    ])
                    .arg(out.join("index.html"))
                    .arg("j@juha.no:/mnt/shovel/share/pansy-zenith-main-peaks/index.html")
                    .status()?
                    .success(),
                "publishing gallery index failed"
            );
        }
        println!(
            "PLOT_PROGRESS plots={plots} pending={pending} input_error_days={failed} total={}",
            days.len()
        );
        if pending == 0 && capon_pending == 0 && raw_deferred == 0 && !koki_running {
            break;
        }
        std::thread::sleep(Duration::from_secs(30));
    }
    Ok(())
}
