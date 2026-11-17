mod cuts;
mod interferometry;
mod mst;
mod pipeline;
use anyhow::{ensure, Result};
use std::path::PathBuf;
use tibet_meteor::{raw_paths, Reader};

fn main() -> Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.first().map(String::as_str) == Some("repack-cuts") {
        ensure!(
            args.len() == 3,
            "Usage: tibet-meteor repack-cuts SOURCE NEW_DESTINATION"
        );
        return tibet_meteor::voltage::repack(
            PathBuf::from(&args[1]).as_path(),
            PathBuf::from(&args[2]).as_path(),
        );
    }
    if args.first().map(String::as_str) == Some("verify-cuts") {
        ensure!(
            args.len() == 2,
            "Usage: tibet-meteor verify-cuts CUT_DIRECTORY"
        );
        return cuts::verify(PathBuf::from(&args[1]).as_path());
    }
    if args.first().map(String::as_str) == Some("mst-fine") {
        ensure!(
            args.len() == 2,
            "Usage: tibet-meteor mst-fine REDUCTION_DIRECTORY"
        );
        return mst::fine(PathBuf::from(&args[1]).as_path());
    }
    if args.first().map(String::as_str) == Some("mst") {
        ensure!(
            args.len() == 4,
            "Usage: tibet-meteor mst RAW_DIRECTORY DETECTION_RUN NEW_OUTPUT"
        );
        return mst::run(
            PathBuf::from(&args[1]).as_path(),
            PathBuf::from(&args[2]).as_path(),
            PathBuf::from(&args[3]).as_path(),
        );
    }
    ensure!(args.len()>=2,"Usage: tibet-meteor inspect INPUT | process INPUT OUTPUT [--max-files N] [--max-pulses N] [--threshold-db D]");
    if args[0] == "image" {
        ensure!(args.len() == 3, "Usage: tibet-meteor image RUN GEOMETRY.h5");
        return interferometry::run(&PathBuf::from(&args[1]), &PathBuf::from(&args[2]));
    }
    if args[0] == "cut" {
        ensure!(
            args.len() == 3 || args.len() == 5,
            "Usage: tibet-meteor cut RUN OUTPUT [--pointing-schedule SCHEDULE.h5]"
        );
        let schedule = if args.len() == 5 {
            ensure!(args[3] == "--pointing-schedule", "unknown cut option");
            Some(PathBuf::from(&args[4]))
        } else {
            None
        };
        return cuts::run(
            &PathBuf::from(&args[1]),
            &PathBuf::from(&args[2]),
            schedule.as_deref(),
        );
    }
    let paths = raw_paths(&PathBuf::from(&args[1]))?;
    if args[0] == "inspect" {
        for path in paths {
            let mut r = Reader::open(&path)?;
            r.pulse(0)?;
            r.pulse(r.header.pulses - 1)?;
            let h = &r.header;
            println!("{}\n  frame={} channels={} gates={} pulses={} fs={} Hz IPP={} us waveform={} width_field={} us INACC={} padding={} bytes\n  TX alpha/beta={}/{}; RX alpha={:?}; beta={:?}\n  amplitude mask IDs={:?}; timestamp header={:?}; filename UTC hypothesis={} us",
                path.display(),h.frame,h.channels,h.gates,h.pulses,h.fs,h.ipp_us,h.waveform,h.width_us,h.inacc,h.padding_bytes,
                h.tx_alpha,h.tx_beta,h.alpha,h.beta,h.amplitude_index,&h.raw[16..32],h.epoch_us);
        }
    } else {
        ensure!(
            args[0] == "process" && args.len() >= 3,
            "unknown command or missing output"
        );
        let mut cfg = pipeline::Config::default();
        for pair in args[3..].chunks(2) {
            ensure!(pair.len() == 2, "option missing value");
            match pair[0].as_str() {
                "--max-files" => cfg.max_files = pair[1].parse()?,
                "--max-pulses" => cfg.max_pulses = pair[1].parse()?,
                "--threshold-db" => cfg.threshold_db = pair[1].parse()?,
                "--write-legacy-cuts" => cfg.write_legacy_cuts = pair[1].parse::<u8>()? != 0,
                _ => anyhow::bail!("unknown option {}", pair[0]),
            }
        }
        pipeline::run(&paths, &PathBuf::from(&args[2]), cfg)?;
    }
    Ok(())
}
