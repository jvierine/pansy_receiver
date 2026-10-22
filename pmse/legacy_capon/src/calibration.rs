//! Historical channel ordering, subarray centres and phase calibration.
use anyhow::{ensure, Context, Result};
use std::{collections::HashMap, fs, path::Path};
pub fn load_reference(root: &Path, beam: usize) -> Result<([[f64; 3]; 7], [f64; 7])> {
    let text = fs::read_to_string(root.join("cfg/antpos.csv"))?;
    let mut groups: HashMap<String, Vec<[f64; 3]>> = HashMap::new();
    for line in text.lines().skip(1).filter(|l| !l.trim().is_empty()) {
        let f: Vec<_> = line.split(',').collect();
        ensure!(f.len() >= 7, "invalid antenna row");
        groups.entry(f[1].to_owned()).or_default().push([
            f[4].trim().parse()?,
            f[5].trim().parse()?,
            f[6].trim().parse()?,
        ]);
    }
    let connections = fs::read_to_string(root.join("cfg/connections.txt"))?;
    let names: Vec<_> = connections
        .lines()
        .filter_map(|l| l.split(',').next())
        .filter(|n| !n.is_empty() && *n != "RFTX")
        .collect();
    ensure!(names.len() == 7, "expected seven receive subarrays");
    let mut positions = [[0.; 3]; 7];
    for (i, name) in names.iter().enumerate() {
        let p = groups.get(*name).context("unknown subarray")?;
        ensure!(!p.is_empty(), "empty subarray");
        for j in 0..3 {
            positions[i][j] = p.iter().map(|r| r[j]).sum::<f64>() / p.len() as f64;
        }
    }
    let phase = hdf5::File::open(root.join(format!("data/phases{beam}.h5")))?
        .dataset("phasecal")?
        .read_raw::<f64>()?;
    let phases: [f64; 7] = phase
        .try_into()
        .map_err(|_| anyhow::anyhow!("expected seven phasecal values"))?;
    Ok((positions, phases))
}
