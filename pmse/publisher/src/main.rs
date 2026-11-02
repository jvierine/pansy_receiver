use chrono::{DateTime, NaiveDate, Utc};
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{self, Read, Write},
    path::{Path, PathBuf},
    process::Command,
    time::{Duration, Instant, UNIX_EPOCH},
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
struct Stamp {
    bytes: u64,
    modified_ns: u128,
}
fn stamp(path: &Path) -> io::Result<Stamp> {
    let m = fs::metadata(path)?;
    Ok(Stamp {
        bytes: m.len(),
        modified_ns: m
            .modified()?
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos(),
    })
}
#[derive(Clone, Serialize, Deserialize)]
struct Entry {
    hint: Stamp,
    native_sha256: String,
    native_bytes: u64,
    ready: bool,
    checked: Option<Stamp>,
}
#[derive(Default, Serialize, Deserialize)]
struct State {
    directories: BTreeMap<String, Stamp>,
    files: BTreeMap<String, Entry>,
    auxiliary: BTreeMap<String, Stamp>,
    #[serde(default)]
    pending_days: BTreeSet<String>,
}
#[derive(Deserialize)]
struct Hint {
    native_sha256: String,
    native_bytes: u64,
    packets: Vec<Value>,
}
fn read_hint(path: &Path) -> Result<Hint> {
    let mut data = Vec::new();
    GzDecoder::new(File::open(path)?).read_to_end(&mut data)?;
    Ok(serde_json::from_slice(&data)?)
}
fn read_json(path: &Path) -> Result<Value> {
    Ok(serde_json::from_slice(&fs::read(path)?)?)
}
fn atomic_json(path: &Path, value: &impl Serialize) -> Result<()> {
    fs::create_dir_all(path.parent().unwrap())?;
    let tmp = path.with_extension(format!("json.{}.partial", std::process::id()));
    let result = (|| -> Result<()> {
        let mut writer = io::BufWriter::new(File::create(&tmp)?);
        serde_json::to_writer(&mut writer, value)?;
        writer.flush()?;
        drop(writer);
        fs::rename(&tmp, path)?;
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(tmp);
    }
    result
}
fn day_of(us: i64) -> String {
    DateTime::<Utc>::from_timestamp_micros(us)
        .unwrap()
        .format("%Y-%m-%d")
        .to_string()
}
fn affected_days(packets: &[Value]) -> BTreeSet<String> {
    let mut days = BTreeSet::new();
    for p in packets {
        if let (Some(a), Some(b)) = (p["start_us"].as_i64(), p["stop_us"].as_i64()) {
            let mut date = NaiveDate::parse_from_str(&day_of(a), "%Y-%m-%d").unwrap();
            let end = day_of(b.saturating_sub(1));
            while date.to_string() <= end {
                days.insert(date.to_string());
                date = date.succ_opt().unwrap();
            }
        }
    }
    days
}
fn day_bounds(day: &str) -> Result<(i64, i64)> {
    let start = NaiveDate::parse_from_str(day, "%Y-%m-%d")?
        .and_hms_opt(0, 0, 0)
        .unwrap()
        .and_utc()
        .timestamp_micros();
    Ok((start, start + 86_400_000_000))
}
fn digest_revision(day: &str, packets: &[Value]) -> String {
    let joined = packets
        .iter()
        .map(|p| {
            format!(
                "{}:{}:{}",
                p["key"].as_str().unwrap_or(""),
                p["source_revision"].as_str().unwrap_or(""),
                if p["spectra_ready"].as_bool().unwrap_or(true) {
                    "True"
                } else {
                    "False"
                }
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("{:x}", Sha256::digest(format!("16ms:{day}{joined}")))[..16].to_string()
}

fn overview(packet: &Value) -> Result<Value> {
    let first = packet["axis_start_us"]
        .as_i64()
        .ok_or("missing time axis")?;
    let stop = packet["axis_stop_us"].as_i64().ok_or("missing time axis")?;
    let mut edges = Vec::<i64>::new();
    let mut cells = Vec::<Value>::new();
    let mut column = 0;
    for frame in packet["frame_packets"].as_array().ok_or("missing frames")? {
        let left = first.max(frame["start_us"].as_i64().ok_or("missing start")?);
        let right = stop.min(frame["stop_us"].as_i64().ok_or("missing stop")?);
        if right <= left {
            continue;
        }
        if let Some(previous) = edges.last() {
            if left < *previous {
                return Err("overlapping PMSE intervals".into());
            }
            if left > *previous {
                edges.push(left);
                column += 1;
            }
        }
        if edges.is_empty() {
            edges.push(left);
        }
        let mut peak = BTreeMap::<u64, f64>::new();
        if let (Some(index), Some(snr)) = (frame["index"].as_array(), frame["snr_db"].as_array()) {
            for (index, snr) in index.iter().zip(snr) {
                if index[0].as_u64() != Some(0) {
                    continue;
                }
                if let (Some(gate), Some(snr)) = (index[2].as_u64(), snr.as_f64()) {
                    let value = peak.entry(gate).or_insert(f64::NEG_INFINITY);
                    *value = value.max(snr);
                }
            }
        }
        for (gate, snr) in peak {
            cells.push(json!([column, gate, snr]));
        }
        edges.push(right);
        column += 1;
    }
    let heights = packet["geometry"]["altitude_km"][0]
        .as_array()
        .ok_or("missing height axis")?
        .iter()
        .map(|v| v.as_f64().ok_or("invalid height"))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    if heights.len() < 2 {
        return Err("height axis requires two bins".into());
    }
    let mut y = vec![heights[0] - (heights[1] - heights[0]) / 2.];
    y.extend(heights.windows(2).map(|v| (v[0] + v[1]) / 2.));
    let n = heights.len();
    y.push(heights[n - 1] + (heights[n - 1] - heights[n - 2]) / 2.);
    Ok(
        json!({"schema":"pansy_zenith_display_v1","selected_day":packet["selected_day"],
        "available_days":packet["available_days"],"revision":packet["revision"],"frame_count":packet["frame_count"],
        "latest_mode_us":packet["latest_mode_us"],"axis_start_us":first,
        "x":edges.iter().map(|t|(*t-first) as f64/3.6e9).collect::<Vec<_>>(),"y":y,"cells":cells}),
    )
}

fn noise_display(model: &Value) -> Result<Value> {
    let origin = model["noise_model_origin_unix_s"]
        .as_f64()
        .ok_or("missing noise origin")?;
    let stop = model["noise_model_observation_stop_unix_s"]
        .as_f64()
        .ok_or("missing noise end")?;
    let mut times = Vec::new();
    let mut powers = Vec::new();
    let mut sky = Vec::new();
    for (i, t) in model["noise_model_data"]["times_h"]
        .as_array()
        .ok_or("missing noise times")?
        .iter()
        .enumerate()
    {
        if origin + t.as_f64().ok_or("invalid noise time")? * 3600. < stop - 86400. {
            continue;
        }
        times.push(t.clone());
        powers.push(json!([model["noise_model_data"]["noise_power"][i][0]]));
        sky.push(json!(
            [model["noise_model_data"]["sky_temperature_k"][i][0]]
        ));
    }
    Ok(
        json!({"schema":"pansy_monitor_noise_v1","noise_model_origin_unix_s":origin,
        "noise_model_observation_stop_unix_s":stop,"noise_zenith_model_status":model["noise_zenith_model_status"],
        "noise_model_data":{"times_h":times,"noise_power":powers,"sky_temperature_k":sky}}),
    )
}

struct Config {
    source: PathBuf,
    native: PathBuf,
    output: PathBuf,
    state: PathBuf,
    sidecar: Option<(String, String)>,
    budget: u64,
    once: bool,
}
struct Verification {
    name: String,
    stamp: Stamp,
    expected: String,
    file: File,
    hash: Sha256,
}
struct Publisher {
    cfg: Config,
    state: State,
    verification: Option<Verification>,
    seed: BTreeMap<String, BTreeMap<String, String>>,
    retry_cursor: usize,
}
impl Publisher {
    fn new(cfg: Config) -> Self {
        let state = File::open(&cfg.state)
            .ok()
            .and_then(|f| serde_json::from_reader(io::BufReader::new(f)).ok())
            .unwrap_or_default();
        Self {
            cfg,
            state,
            verification: None,
            seed: BTreeMap::new(),
            retry_cursor: 0,
        }
    }
    fn discover(&mut self) -> Result<BTreeSet<String>> {
        let mut dirty = BTreeSet::new();
        for dir in fs::read_dir(&self.cfg.source)? {
            let dir = dir?;
            let name = dir.file_name().to_string_lossy().into_owned();
            if name.len() < 13
                || !name.starts_with("20")
                || !name.contains('T')
                || !dir.file_type()?.is_dir()
            {
                continue;
            }
            let now = stamp(&dir.path())?;
            if self.state.directories.get(&name) == Some(&now) {
                continue;
            }
            let mut complete = true;
            for item in fs::read_dir(dir.path())? {
                let item = item?;
                let filename = item.file_name().to_string_lossy().into_owned();
                if !filename.starts_with("pmse@") || !filename.ends_with(".h5.json.gz") {
                    continue;
                }
                let key = format!("{name}/{}", filename.trim_end_matches(".json.gz"));
                let signature = stamp(&item.path())?;
                if self.state.files.get(&key).map(|e| &e.hint) == Some(&signature) {
                    continue;
                }
                // A partial or invalid source delivery must not hide existing good output.
                let hint = match read_hint(&item.path()) {
                    Ok(v) => v,
                    Err(e) => {
                        eprintln!("delayed {key}: {e}");
                        complete = false;
                        continue;
                    }
                };
                if stamp(&item.path())? != signature {
                    complete = false;
                    continue;
                }
                dirty.extend(affected_days(&hint.packets));
                let previous = self.state.files.get(&key).cloned();
                let unchanged = previous
                    .as_ref()
                    .map(|e| e.native_sha256 == hint.native_sha256)
                    .unwrap_or(false);
                let ready = if unchanged {
                    previous.as_ref().unwrap().ready
                } else {
                    self.seed_ready(&hint.packets)
                };
                let checked = if unchanged {
                    previous.as_ref().unwrap().checked.clone()
                } else {
                    None
                };
                self.state.files.insert(
                    key,
                    Entry {
                        hint: signature,
                        native_sha256: hint.native_sha256,
                        native_bytes: hint.native_bytes,
                        ready,
                        checked,
                    },
                );
            }
            if complete {
                self.state.directories.insert(name, now);
            }
        }
        Ok(dirty)
    }
    fn seed_ready(&mut self, packets: &[Value]) -> bool {
        // Preserve readiness that the preceding publisher already verified. A changed
        // source_revision never inherits readiness from an older published record.
        if packets.is_empty() {
            return false;
        }
        for day in affected_days(packets) {
            if self.seed.contains_key(&day) {
                continue;
            }
            let mut published = BTreeMap::new();
            if let Ok(v) = read_json(&self.cfg.output.join("days").join(format!("{day}.json"))) {
                if let Some(a) = v["frame_packets"].as_array() {
                    for p in a {
                        if p["spectra_ready"].as_bool().unwrap_or(true) {
                            published.insert(
                                p["key"].to_string(),
                                p["source_revision"].as_str().unwrap_or("").to_string(),
                            );
                        }
                    }
                }
            }
            self.seed.insert(day, published);
        }
        packets.iter().all(|p| {
            self.seed.values().any(|published| {
                published
                    .get(&p["key"].to_string())
                    .map(|old| Some(old.as_str()) == p["source_revision"].as_str())
                    .unwrap_or(false)
            })
        })
    }
    fn verify(&mut self, dirty: &mut BTreeSet<String>) -> Result<u64> {
        let mut used = 0;
        // The shared filesystem is touched only for pending deliveries. Ready files
        // are never walked, opened, or rehashed on idle cycles or process restarts.
        if self.verification.is_none() {
            let mut candidates = self
                .state
                .files
                .iter()
                .filter(|(_, e)| !e.ready)
                .map(|(n, _)| n.clone())
                .collect::<Vec<_>>();
            candidates.reverse();
            if !candidates.is_empty() {
                let offset = self.retry_cursor % candidates.len();
                candidates.rotate_left(offset);
                self.retry_cursor = (offset + 16) % candidates.len();
            }
            let mut checked = 0;
            for name in candidates {
                let e = &self.state.files[&name];
                if e.ready {
                    continue;
                }
                checked += 1;
                if checked > 16 {
                    break;
                }
                let path = self.cfg.native.join(&name);
                let metadata = match stamp(&path) {
                    Ok(v) => v,
                    Err(_) => continue,
                };
                if metadata.bytes != e.native_bytes || e.checked.as_ref() == Some(&metadata) {
                    continue;
                }
                if let Ok(file) = File::open(path) {
                    self.verification = Some(Verification {
                        name,
                        stamp: metadata,
                        expected: e.native_sha256.clone(),
                        file,
                        hash: Sha256::new(),
                    });
                    break;
                }
            }
        }
        if let Some(mut job) = self.verification.take() {
            let still_expected = self
                .state
                .files
                .get(&job.name)
                .map(|e| e.native_sha256 == job.expected && e.native_bytes == job.stamp.bytes)
                .unwrap_or(false);
            if !still_expected {
                return Ok(0);
            }
            let mut buffer = vec![0; 64 * 1024];
            let mut finished = false;
            while used < self.cfg.budget {
                let size = (self.cfg.budget - used).min(buffer.len() as u64) as usize;
                let n = match job.file.read(&mut buffer[..size]) {
                    Ok(n) => n,
                    Err(e) => {
                        eprintln!("native verification delayed: {e}");
                        return Ok(used);
                    }
                };
                if n == 0 {
                    finished = true;
                    break;
                }
                used += n as u64;
                job.hash.update(&buffer[..n]);
            }
            if finished {
                let e = self.state.files.get_mut(&job.name).unwrap();
                if stamp(&self.cfg.native.join(&job.name)).ok().as_ref() == Some(&job.stamp) {
                    e.ready = format!("{:x}", job.hash.finalize()) == e.native_sha256;
                    e.checked = Some(job.stamp);
                    if e.ready {
                        dirty.extend(affected_days(
                            &read_hint(&self.cfg.source.join(format!("{}.json.gz", job.name)))?
                                .packets,
                        ));
                    }
                }
            } else {
                self.verification = Some(job)
            }
        }
        Ok(used)
    }
    fn build_day(&self, day: &str, base: &Value) -> Result<Value> {
        let (start, end) = day_bounds(day)?;
        let mut result = base.clone();
        let old = read_json(&self.cfg.output.join("days").join(format!("{day}.json")))
            .unwrap_or(Value::Null);
        let mut packets = BTreeMap::new();
        if let Some(a) = old["frame_packets"].as_array() {
            for p in a {
                packets.insert(p["key"].to_string(), p.clone());
            }
        }
        packets.retain(|_, p| {
            let Some(key) = p["key"].as_str().and_then(|v| v.parse::<i64>().ok()) else {
                return true;
            };
            let hour = DateTime::<Utc>::from_timestamp_micros(key)
                .unwrap()
                .format("%Y-%m-%dT%H-00-00");
            !self
                .state
                .files
                .contains_key(&format!("{hour}/pmse@{}.h5", key / 600_000_000 * 600))
        });
        let previous = day_of(start - 1);
        for (name, entry) in &self.state.files {
            if !name.starts_with(day) && !name.starts_with(&format!("{previous}T23")) {
                continue;
            }
            let hint = read_hint(&self.cfg.source.join(format!("{name}.json.gz")))?;
            if hint.native_sha256 != entry.native_sha256 {
                return Err("source hint changed during publication; retrying".into());
            }
            for mut p in hint.packets {
                if p["stop_us"].as_i64().unwrap_or(0) <= start
                    || p["start_us"].as_i64().unwrap_or(end) >= end
                {
                    continue;
                }
                p["spectra_ready"] = json!(entry.ready);
                packets.insert(p["key"].to_string(), p);
            }
        }
        let mut packets = packets.into_values().collect::<Vec<_>>();
        packets.sort_by_key(|p| p["start_us"].as_i64());
        let revision = digest_revision(day, &packets);
        result["selected_day"] = json!(day);
        result["revision"] = json!(revision);
        result["generated_utc"] = json!(Utc::now().to_rfc3339());
        result["axis_start_us"] = json!(start);
        result["axis_stop_us"] = json!(end);
        result["first_mode_us"] = json!(packets
            .iter()
            .filter_map(|p| p["start_us"].as_i64())
            .min()
            .unwrap_or(start));
        result["latest_mode_us"] = json!(packets
            .iter()
            .filter_map(|p| p["stop_us"].as_i64())
            .max()
            .unwrap_or(start));
        result["frame_count"] = json!(packets.len());
        result["spectra_transfer_pending"] = json!(packets
            .iter()
            .filter(|p| p["spectra_ready"] == false)
            .count());
        result["frames"] = json!(packets
            .iter()
            .map(|p| {
                let mut frame = serde_json::Map::new();
                for key in [
                    "key",
                    "start_us",
                    "stop_us",
                    "duration_s",
                    "coverage_s",
                    "looks",
                    "quality",
                ] {
                    frame.insert(key.to_string(), p[key].clone());
                }
                frame
            })
            .collect::<Vec<_>>());
        result["frame_packets"] = json!(packets);
        result["sidecar"] = json!(format!("gui-data-{day}.h5"));
        Ok(result)
    }
    fn cycle(&mut self) -> Result<()> {
        let began = Instant::now();
        let old_state = serde_json::to_vec(&self.state)?;
        let mut dirty = self.state.pending_days.clone();
        dirty.extend(self.discover()?);
        let verified_bytes = self.verify(&mut dirty)?;
        self.state.pending_days = dirty.clone();
        let mut auxiliary_pending = false;
        for (name, path) in [
            (
                "processing-status.json",
                self.cfg.native.join("processing-status.json"),
            ),
            (
                "sky-noise/status.json",
                self.cfg.output.join("sky-noise/status.json"),
            ),
            (
                "archive-plot-status.json",
                self.cfg.output.join("archive-plot-status.json"),
            ),
            (
                "monitor-noise.json",
                self.cfg.output.join("monitor-noise.json"),
            ),
        ] {
            if let Ok(s) = stamp(&path) {
                auxiliary_pending |= self.state.auxiliary.get(name) != Some(&s);
            }
        }
        if dirty.is_empty() && !auxiliary_pending {
            if serde_json::to_vec(&self.state)? != old_state {
                atomic_json(&self.cfg.state, &self.state)?;
            }
            println!(
                "cycle_ms={} dirty_days=0 native_bytes_read={verified_bytes}",
                began.elapsed().as_millis()
            );
            return Ok(());
        }
        let mut base = read_json(&self.cfg.output.join("manifest.json"))?;
        let mut days = base["available_days"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(|v| v.as_str().map(String::from))
            .collect::<BTreeSet<_>>();
        days.extend(dirty.iter().cloned());
        base["available_days"] = json!(days);
        let mut auxiliary_changed = false;
        let noise_path = self.cfg.output.join("sky-noise/status.json");
        if let Ok(s) = stamp(&noise_path) {
            if self.state.auxiliary.get("sky-noise/status.json") != Some(&s) {
                if let Ok(value) = read_json(&noise_path).and_then(|v| noise_display(&v)) {
                    atomic_json(&self.cfg.output.join("monitor-noise.json"), &value)?;
                    self.state
                        .auxiliary
                        .insert("sky-noise/status.json".into(), s);
                }
            }
        }
        let status_path = self.cfg.native.join("processing-status.json");
        if let Ok(s) = stamp(&status_path) {
            if self.state.auxiliary.get("processing-status.json") != Some(&s) {
                if let Ok(value) = read_json(&status_path) {
                    base["source_processing"] = value;
                    self.state
                        .auxiliary
                        .insert("processing-status.json".into(), s);
                    auxiliary_changed = true;
                }
            }
        }
        for (file, field) in [
            ("archive-plot-status.json", "archive_plot"),
            ("monitor-noise.json", "monitor_noise"),
        ] {
            let path = self.cfg.output.join(file);
            if let Ok(s) = stamp(&path) {
                if self.state.auxiliary.get(file) != Some(&s) {
                    if let Ok(value) = read_json(&path) {
                        if field == "archive_plot" {
                            base[field] = value;
                        }
                        self.state.auxiliary.insert(file.to_string(), s);
                        auxiliary_changed = true;
                    }
                }
            }
        }
        let latest = days.last().cloned();
        let mut wrote_latest = false;
        for day in &dirty {
            let packet = self.build_day(day, &base)?;
            let target = self.cfg.output.join("days").join(format!("{day}.json"));
            let unchanged = read_json(&target)
                .ok()
                .map(|old| old["revision"] == packet["revision"])
                .unwrap_or(false);
            if unchanged {
                if latest.as_ref() == Some(day) && base["revision"] != packet["revision"] {
                    base = packet;
                    wrote_latest = true;
                }
                continue;
            }
            // Publish the downloadable HDF5 before exposing its new JSON revision.
            if let Some((python, script)) = &self.cfg.sidecar {
                let staging = self.cfg.output.join(format!(".sidecar-{day}.json"));
                atomic_json(&staging, &packet)?;
                let status = Command::new(python)
                    .arg(script)
                    .arg(&staging)
                    .arg(self.cfg.output.join(format!("gui-data-{day}.h5")))
                    .status()?;
                let _ = fs::remove_file(staging);
                if !status.success() {
                    return Err("HDF5 sidecar generation failed".into());
                }
            }
            atomic_json(&target, &packet)?;
            if latest.as_ref() == Some(day) {
                base = packet;
                wrote_latest = true;
            }
        }
        if wrote_latest || auxiliary_changed {
            base["available_days"] = json!(days);
            atomic_json(&self.cfg.output.join("manifest.json"), &base)?;
            if let Ok(value) = overview(&base) {
                atomic_json(&self.cfg.output.join("monitor-overview.json"), &value)?;
            }
        }
        for day in &dirty {
            self.state.pending_days.remove(day);
        }
        if serde_json::to_vec(&self.state)? != old_state {
            atomic_json(&self.cfg.state, &self.state)?;
        }
        self.seed.clear();
        println!(
            "cycle_ms={} dirty_days={} native_bytes_read={verified_bytes}",
            began.elapsed().as_millis(),
            dirty.len()
        );
        Ok(())
    }
}
fn main() -> Result<()> {
    let mut args = std::env::args().skip(1);
    let mut values = BTreeMap::new();
    let mut once = false;
    while let Some(arg) = args.next() {
        if arg == "--once" {
            once = true
        } else {
            let value = args.next().ok_or("missing argument value")?;
            values.insert(arg, value);
        }
    }
    let output = PathBuf::from(values.get("--output").ok_or("--output is required")?);
    let mut publisher = Publisher::new(Config {
        source: values
            .get("--source")
            .map(PathBuf::from)
            .unwrap_or_else(|| output.join("source")),
        native: PathBuf::from(
            values
                .get("--native-root")
                .ok_or("--native-root is required")?,
        ),
        state: values
            .get("--state")
            .map(PathBuf::from)
            .unwrap_or_else(|| output.join("publisher-state/cache.json")),
        sidecar: values.get("--sidecar-script").map(|s| {
            (
                values
                    .get("--python")
                    .cloned()
                    .unwrap_or_else(|| "python3".into()),
                s.clone(),
            )
        }),
        budget: values
            .get("--verify-bytes")
            .map(|v| v.parse())
            .transpose()?
            .unwrap_or(4 * 1024 * 1024),
        output,
        once,
    });
    loop {
        if let Err(e) = publisher.cycle() {
            eprintln!("publisher delayed: {e}");
            if publisher.cfg.once {
                return Err(e);
            }
        }
        if publisher.cfg.once {
            break;
        }
        std::thread::sleep(Duration::from_secs(60));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{write::GzEncoder, Compression};
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(0);
    struct Fixture(PathBuf);
    impl Fixture {
        fn new() -> Self {
            let p = std::env::temp_dir().join(format!(
                "pansy-publisher-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir_all(p.join("output/source")).unwrap();
            atomic_json(&p.join("output/manifest.json"), &json!({"available_days":[],"geometry":{"frequency_hz":47000000},"custom_monitor":{"preserve":true}})).unwrap();
            Self(p)
        }
        fn config(&self) -> Config {
            Config {
                source: self.0.join("output/source"),
                output: self.0.join("output"),
                native: self.0.join("native"),
                state: self.0.join("state.json"),
                sidecar: None,
                budget: 4,
                once: true,
            }
        }
        fn hint(&self, bytes: &[u8], packets: &[Value]) -> String {
            let key = packets[0]["key"].as_str().unwrap().parse::<i64>().unwrap();
            let hour = DateTime::<Utc>::from_timestamp_micros(key)
                .unwrap()
                .format("%Y-%m-%dT%H-00-00");
            let name = format!("{hour}/pmse@{}.h5", key / 600_000_000 * 600);
            let path = self.0.join("output/source").join(format!("{name}.json.gz"));
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            let temporary = path.with_extension("pending");
            let mut out = GzEncoder::new(File::create(&temporary).unwrap(), Compression::fast());
            serde_json::to_writer(&mut out,&json!({"native_bytes":bytes.len(),"native_sha256":format!("{:x}",Sha256::digest(bytes)),"packets":packets})).unwrap();
            out.finish().unwrap();
            fs::rename(temporary, path).unwrap();
            name
        }
        fn native(&self, name: &str, bytes: &[u8]) {
            let p = self.0.join("native").join(name);
            fs::create_dir_all(p.parent().unwrap()).unwrap();
            fs::write(p, bytes).unwrap();
        }
        fn day(&self, day: &str) -> Value {
            read_json(&self.0.join(format!("output/days/{day}.json"))).unwrap()
        }
    }
    impl Drop for Fixture {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    fn packet(key: i64, revision: &str) -> Value {
        json!({"key":key.to_string(),"start_us":key,"stop_us":key+90_000_000,
        "source_revision":revision,"duration_s":90,"coverage_s":90,"looks":10,"quality":"ok"})
    }
    #[test]
    fn incremental_verification_survives_restart_and_idle_does_not_rewrite() {
        let f = Fixture::new();
        let key = day_bounds("2026-11-02").unwrap().0;
        let name = f.hint(b"abcdefgh", &[packet(key, "v1")]);
        f.native(&name, b"abcdefgh");
        let mut p = Publisher::new(f.config());
        p.cycle().unwrap();
        assert_eq!(f.day("2026-11-02")["spectra_transfer_pending"], 1);
        p.cycle().unwrap();
        p.cycle().unwrap();
        assert_eq!(f.day("2026-11-02")["spectra_transfer_pending"], 0);
        assert!(p.state.files[&name].ready);
        let before = stamp(&f.0.join("output/manifest.json")).unwrap();
        fs::remove_dir_all(f.0.join("native")).unwrap();
        let mut restarted = Publisher::new(f.config());
        restarted.cycle().unwrap();
        assert_eq!(stamp(&f.0.join("output/manifest.json")).unwrap(), before);
        assert!(restarted.verification.is_none());
    }
    #[test]
    fn same_size_wrong_native_does_not_become_ready() {
        let f = Fixture::new();
        let key = day_bounds("2026-11-02").unwrap().0;
        let name = f.hint(b"good", &[packet(key, "v1")]);
        f.native(&name, b"evil");
        let mut p = Publisher::new(f.config());
        p.cycle().unwrap();
        p.cycle().unwrap();
        assert!(!p.state.files[&name].ready);
        assert_eq!(f.day("2026-11-02")["spectra_transfer_pending"], 1);
    }
    #[test]
    fn replacement_removes_old_frames_and_resets_readiness() {
        let f = Fixture::new();
        let key = day_bounds("2026-11-02").unwrap().0;
        let name = f.hint(
            b"good",
            &[packet(key, "v1"), packet(key + 100_000_000, "v1")],
        );
        f.native(&name, b"good");
        let mut p = Publisher::new(f.config());
        p.cycle().unwrap();
        p.cycle().unwrap();
        f.hint(b"next", &[packet(key, "v2")]);
        p.cycle().unwrap();
        p.cycle().unwrap();
        let day = f.day("2026-11-02");
        assert_eq!(day["frame_count"], 1);
        assert_eq!(day["spectra_transfer_pending"], 1);
        assert_eq!(day["frame_packets"][0]["source_revision"], "v2");
    }
    #[test]
    fn cross_midnight_packet_appears_on_both_days() {
        let f = Fixture::new();
        let key = day_bounds("2026-11-02").unwrap().0 - 30_000_000;
        f.hint(b"x", &[packet(key, "v1")]);
        let mut p = Publisher::new(f.config());
        p.cycle().unwrap();
        assert_eq!(f.day("2026-11-01")["frame_count"], 1);
        assert_eq!(f.day("2026-11-02")["frame_count"], 1);
        assert_eq!(f.day("2026-11-02")["custom_monitor"]["preserve"], true);
    }
    #[test]
    fn malformed_hint_preserves_existing_publication() {
        let f = Fixture::new();
        let key = day_bounds("2026-11-02").unwrap().0;
        let name = f.hint(b"good", &[packet(key, "v1")]);
        let mut p = Publisher::new(f.config());
        p.cycle().unwrap();
        let before = f.day("2026-11-02");
        let path = f.0.join("output/source").join(format!("{name}.json.gz"));
        let temp = path.with_extension("pending");
        fs::write(&temp, b"broken").unwrap();
        fs::rename(temp, path).unwrap();
        p.cycle().unwrap();
        assert_eq!(f.day("2026-11-02"), before);
    }
    #[test]
    fn overview_preserves_gaps_and_zenith_snr_maximum() {
        let mut a = packet(0, "v1");
        a["index"] = json!([[0, 1, 0], [0, 2, 0], [1, 2, 0]]);
        a["snr_db"] = json!([5., 8., 100.]);
        let b = packet(120_000_000, "v1");
        let o=overview(&json!({"axis_start_us":0,"axis_stop_us":300_000_000i64,"geometry":{"altitude_km":[[80.,81.,82.]]},"frame_packets":[a,b]})).unwrap();
        assert_eq!(o["cells"], json!([[0, 0, 8.]]));
        assert_eq!(o["x"].as_array().unwrap().len(), 4);
        assert_eq!(o["y"], json!([79.5, 80.5, 81.5, 82.5]));
    }

    #[test]
    fn json_preserves_existing_float_values_exactly() {
        let expected = 50.664925401999994_f64;
        let parsed: Value = serde_json::from_str("{\"height\":50.664925401999994}").unwrap();
        assert_eq!(
            parsed["height"].as_f64().unwrap().to_bits(),
            expected.to_bits()
        );
        let again: Value = serde_json::from_slice(&serde_json::to_vec(&parsed).unwrap()).unwrap();
        assert_eq!(
            again["height"].as_f64().unwrap().to_bits(),
            expected.to_bits()
        );
    }

    #[test]
    fn failed_sidecar_is_retried_without_another_source_change() {
        let f = Fixture::new();
        let key = day_bounds("2026-11-02").unwrap().0;
        f.hint(b"good", &[packet(key, "v1")]);
        let mut cfg = f.config();
        cfg.sidecar = Some(("/usr/bin/false".into(), "unused".into()));
        let mut p = Publisher::new(cfg);
        assert!(p.cycle().is_err());
        assert!(p.state.pending_days.contains("2026-11-02"));
        p.cfg.sidecar = None;
        p.cycle().unwrap();
        assert_eq!(f.day("2026-11-02")["frame_count"], 1);
        assert!(p.state.pending_days.is_empty());
    }
}
