use anyhow::{ensure, Result};
use hdf5::{File, Group, H5Type};
use pansy_meteor_matched_filter::Filter;
use rustfft::num_complex::Complex32;
use std::{
    path::{Path, PathBuf},
    time::Instant,
};
use tibet_meteor::{Header, Reader};

const C: f64 = 299_792_458.;
const TX_START: usize = 7;
const TX_LEN: usize = 32;
const FFT: usize = 128;
const CHANNELS: usize = 15;
const REF: usize = 8;

pub struct Config {
    pub max_files: usize,
    pub max_pulses: usize,
    pub threshold_db: f32,
    pub write_legacy_cuts: bool,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            max_files: usize::MAX,
            max_pulses: usize::MAX,
            threshold_db: 12.,
            write_legacy_cuts: true,
        }
    }
}

#[derive(Clone, Debug)]
struct Point {
    pulse: usize,
    us: u64,
    lag: usize,
    bin: usize,
    doppler: f64,
    range: f64,
    statistic_db: f32,
    beam: usize,
    amplitudes: Vec<Complex32>,
    noise: Vec<f32>,
}
#[derive(Default)]
struct Track {
    points: Vec<Point>,
}

pub(crate) fn put<T: H5Type>(g: &Group, name: &str, shape: &[usize], values: &[T]) -> Result<()> {
    let ds = g.new_dataset::<T>().shape(shape).create(name)?;
    ds.write_raw(values)?;
    Ok(())
}
pub(crate) fn scalar<T: H5Type>(g: &Group, name: &str, value: T) -> Result<()> {
    put(g, name, &[], &[value])
}
pub(crate) fn text(g: &Group, name: &str, value: &str) -> Result<()> {
    put(g, name, &[value.len()], value.as_bytes())
}
pub(crate) fn save_header(g: &Group, h: &Header) -> Result<()> {
    put(g, "raw_header", &[h.raw.len()], &h.raw)?;
    text(g, "source_path", &h.path.to_string_lossy())?;
    text(g,"time_source","filename Beijing time hypothesis; header clock invalid; filename has one-second resolution")?;
    scalar(g, "sample_rate_hz", h.fs)?;
    scalar(g, "ipp_us", h.ipp_us)?;
    scalar(g, "frequency_hz", h.frequency_hz)?;
    scalar(g, "header_pulse_width_us", h.width_us)?;
    scalar(g, "waveform", h.waveform)?;
    scalar(g, "inacc", h.inacc)?;
    scalar(g, "filename_epoch_us", h.epoch_us)?;
    scalar(g, "header_pulses", h.pulses as u64)?;
    scalar(g, "padding_bytes", h.padding_bytes)?;
    scalar(g, "tx_alpha_code", h.tx_alpha)?;
    scalar(g, "tx_beta_code", h.tx_beta)?;
    put(
        g,
        "rx_direction_uv_inferred",
        &[h.channels, 2],
        &h.inferred_rx_uv(),
    )?;
    text(g,"rx_pointing_inference","signed phase codes on 2.1 m x / 3.637 m y lattice, negative receive-steering sign; consistent with supplied S,N,Z,E,W labels, not hardware-confirmed")?;
    for (name, values) in [
        ("rx_alpha_code", &h.alpha),
        ("rx_beta_code", &h.beta),
        ("phase_weight_index", &h.phase_index),
        ("amplitude_weight_index", &h.amplitude_index),
    ] {
        put(g, name, &[values.len()], values)?;
    }
    Ok(())
}
fn median(values: &mut [f32]) -> f32 {
    let n = values.len() / 2;
    values.select_nth_unstable_by(n, f32::total_cmp);
    values[n]
}
fn noises(z: &[Complex32], gates: usize) -> Vec<f32> {
    (0..CHANNELS)
        .map(|ch| {
            let mut p = z[ch * gates + 250..(ch + 1) * gates]
                .iter()
                .map(|x| x.norm_sqr())
                .collect::<Vec<_>>();
            (median(&mut p) / std::f32::consts::LN_2).max(1.)
        })
        .collect()
}
fn phase_coherence(a: &Point, b: &Point) -> f32 {
    let mut cross = Complex32::default();
    let mut aa = 0.;
    let mut bb = 0.;
    for ch in 8..15 {
        let x = a.amplitudes[ch] / a.noise[ch].sqrt();
        let y = b.amplitudes[ch] / b.noise[ch].sqrt();
        cross += x * y.conj();
        aa += x.norm_sqr();
        bb += y.norm_sqr();
    }
    cross.norm() / (aa * bb).sqrt().max(1e-20)
}
fn predicted_code_alias(weak: &Point, strong: &Point, tx: &[Complex32], fs: f64) -> bool {
    let shift = weak.lag as isize - strong.lag as isize;
    if shift.unsigned_abs() >= tx.len() {
        return false;
    }
    let mut response = Complex32::default();
    for (k, t) in tx.iter().enumerate() {
        let j = k as isize + shift;
        if j < 0 || j >= tx.len() as isize {
            continue;
        }
        let phase =
            2. * std::f64::consts::PI * (strong.doppler * j as f64 - weak.doppler * k as f64) / fs;
        response += tx[j as usize] * t.conj() * Complex32::from_polar(1., phase as f32);
    }
    let mut residual = 0.;
    let mut observed = 0.;
    for ch in 0..CHANNELS {
        residual +=
            (weak.amplitudes[ch] - strong.amplitudes[ch] * response).norm_sqr() / weak.noise[ch];
        observed += weak.amplitudes[ch].norm_sqr() / weak.noise[ch];
    }
    // The full complex 15-channel response, not only range proximity, must agree.
    observed > 0. && residual / observed < 0.35
}
fn attach(tracks: &mut Vec<Track>, point: Point, lambda: f64) {
    let mut best = None;
    for (i, track) in tracks.iter().enumerate() {
        let previous = track.points.last().unwrap();
        if point.pulse <= previous.pulse {
            continue;
        }
        let dt = (point.us - previous.us) as f64 * 1e-6;
        if dt > 0.05 {
            continue;
        }
        let speed = -previous.doppler * lambda / 2. / 1000.;
        let error = (point.range - (previous.range + speed * dt)).abs();
        let dv = (point.doppler - previous.doppler).abs() * lambda / 2. / 1000.;
        let phase = phase_coherence(&point, previous);
        if error > 0.75 || dv > 25. || phase < 0.65 {
            continue;
        }
        let score = error + dv / 50. + (1. - phase) as f64;
        if best.is_none_or(|(_, old)| score < old) {
            best = Some((i, score));
        }
    }
    if let Some((i, _)) = best {
        tracks[i].points.push(point);
    } else {
        tracks.push(Track {
            points: vec![point],
        });
    }
}
fn point_fields(g: &Group, points: &[Point]) -> Result<()> {
    let n = points.len();
    put(
        g,
        "pulse_index",
        &[n],
        &points.iter().map(|p| p.pulse as u64).collect::<Vec<_>>(),
    )?;
    put(
        g,
        "tx_time_us",
        &[n],
        &points.iter().map(|p| p.us).collect::<Vec<_>>(),
    )?;
    put(
        g,
        "range_km_nominal",
        &[n],
        &points.iter().map(|p| p.range).collect::<Vec<_>>(),
    )?;
    put(
        g,
        "doppler_hz",
        &[n],
        &points.iter().map(|p| p.doppler).collect::<Vec<_>>(),
    )?;
    put(
        g,
        "detector_excess_over_median_db",
        &[n],
        &points.iter().map(|p| p.statistic_db).collect::<Vec<_>>(),
    )?;
    put(
        g,
        "strongest_full_array_rx",
        &[n],
        &points.iter().map(|p| p.beam as u16).collect::<Vec<_>>(),
    )?;
    put(
        g,
        "delay_samples",
        &[n],
        &points.iter().map(|p| p.lag as u32).collect::<Vec<_>>(),
    )?;
    put(
        g,
        "doppler_bin",
        &[n],
        &points.iter().map(|p| p.bin as u32).collect::<Vec<_>>(),
    )?;
    let iq = points
        .iter()
        .flat_map(|p| p.amplitudes.iter().flat_map(|z| [z.re, z.im]))
        .collect::<Vec<_>>();
    put(g, "matched_amplitude_re_im", &[n, CHANNELS, 2], &iq)?;
    put(
        g,
        "noise_power",
        &[n, CHANNELS],
        &points
            .iter()
            .flat_map(|p| p.noise.iter().copied())
            .collect::<Vec<_>>(),
    )?;
    Ok(())
}

fn range_rate(track: &Track) -> f64 {
    let t0 = track.points[0].us;
    let n = track.points.len() as f64;
    let tmean = track
        .points
        .iter()
        .map(|p| (p.us - t0) as f64 * 1e-6)
        .sum::<f64>()
        / n;
    let rmean = track.points.iter().map(|p| p.range).sum::<f64>() / n;
    let numerator = track
        .points
        .iter()
        .map(|p| ((p.us - t0) as f64 * 1e-6 - tmean) * (p.range - rmean))
        .sum::<f64>();
    let denominator = track
        .points
        .iter()
        .map(|p| ((p.us - t0) as f64 * 1e-6 - tmean).powi(2))
        .sum::<f64>();
    numerator / denominator
}

fn cut(reader: &mut Reader, track: &Track, path: &Path) -> Result<()> {
    let h = reader.header.clone();
    let padding = (0.08e6 / h.ipp_us as f64).ceil() as usize;
    let first = track.points.first().unwrap().pulse.saturating_sub(padding);
    let end = (track.points.last().unwrap().pulse + padding + 1).min(h.pulses);
    let start_gate =
        (track.points.iter().map(|p| p.lag + TX_START).min().unwrap()).saturating_sub(16);
    let end_gate = (track
        .points
        .iter()
        .map(|p| p.lag + TX_START + TX_LEN)
        .max()
        .unwrap()
        + 16)
        .min(h.gates);
    let width = end_gate - start_gate;
    let f = File::create(path)?;
    save_header(&f, &h)?;
    point_fields(&f.create_group("detections")?, &track.points)?;
    let mut raw = Vec::with_capacity((end - first) * CHANNELS * width * 2);
    let mut tx = Vec::with_capacity((end - first) * TX_LEN * 2);
    for pulse in first..end {
        let z = reader.pulse_i32(pulse)?;
        for ch in 0..CHANNELS {
            raw.extend_from_slice(
                &z[2 * (ch * h.gates + start_gate)..2 * (ch * h.gates + end_gate)],
            );
        }
        tx.extend_from_slice(
            &z[2 * (2 * h.gates + TX_START)..2 * (2 * h.gates + TX_START + TX_LEN)],
        );
    }
    put(
        &f,
        "zrx_echoes_iq_i32",
        &[end - first, CHANNELS, width, 2],
        &raw,
    )?;
    put(&f, "ztx_iq_i32", &[end - first, TX_LEN, 2], &tx)?;
    put(
        &f,
        "pulse_index",
        &[end - first],
        &(first as u64..end as u64).collect::<Vec<_>>(),
    )?;
    put(
        &f,
        "tx_time_us",
        &[end - first],
        &(first..end).map(|i| h.pulse_time_us(i)).collect::<Vec<_>>(),
    )?;
    scalar(&f, "range_start_gate", start_gate as u32)?;
    scalar(&f, "tx_start_gate", TX_START as u32)?;
    scalar(&f, "txlen", TX_LEN as u32)?;
    scalar(&f, "time_padding_clipped_at_source_boundary", 1u8)?;
    text(
        &f,
        "status",
        "meteor_candidate; no manual or angular classification yet",
    )?;
    f.flush()?;
    Ok(())
}

fn calibration(tracks: &[Track], path: &Path) -> Result<()> {
    let mut selected = Vec::new();
    for track in tracks {
        let p = track
            .points
            .iter()
            .max_by(|a, b| a.statistic_db.total_cmp(&b.statistic_db))
            .unwrap();
        let own = p.amplitudes[2].norm_sqr() / p.noise[2];
        let other = (0..5)
            .filter(|&i| i != 2)
            .map(|i| p.amplitudes[i].norm_sqr() / p.noise[i])
            .fold(0., f32::max);
        if p.beam == 2 && p.statistic_db >= 18. && own > 2. * other {
            selected.push(p);
        }
    }
    let f = File::create(path)?;
    scalar(&f, "selected_events", selected.len() as u64)?;
    scalar(&f, "reference_channel", REF as u32)?;
    text(&f,"status","PROVISIONAL boresight-ensemble estimate; not applied to data; angular bias and phase slope unresolved")?;
    text(&f,"selection","one peak per track; strongest full-array channel=2; >18 dB detector statistic; zenith/other full-array power >3 dB")?;
    let mut phase = vec![f32::NAN; CHANNELS];
    let mut coherence = phase.clone();
    for ch in 5..CHANNELS {
        let sum: Complex32 = selected
            .iter()
            .map(|p| {
                let cross = p.amplitudes[ch] * p.amplitudes[REF].conj();
                cross / cross.norm().max(1e-20)
            })
            .sum();
        if !selected.is_empty() {
            phase[ch] = sum.arg();
            coherence[ch] = sum.norm() / selected.len() as f32;
        }
    }
    put(&f, "phasecal_rad", &[CHANNELS], &phase)?;
    put(&f, "circular_concentration", &[CHANNELS], &coherence)?;
    put(
        &f,
        "between_event_phase_dispersion_rad",
        &[CHANNELS],
        &coherence
            .iter()
            .map(|r| (-2. * r.clamp(1e-12, 1.).ln()).sqrt())
            .collect::<Vec<_>>(),
    )?;
    put(
        &f,
        "voltage_correction_re_im",
        &[CHANNELS, 2],
        &phase
            .iter()
            .flat_map(|p| [p.cos(), -p.sin()])
            .collect::<Vec<_>>(),
    )?;
    put(
        &f,
        "selected_peak_time_us",
        &[selected.len()],
        &selected.iter().map(|p| p.us).collect::<Vec<_>>(),
    )?;
    println!(
        "Phase calibration: {} eligible events; provisional, not applied",
        selected.len()
    );
    Ok(())
}

pub fn run(paths: &[PathBuf], output: &Path, cfg: Config) -> Result<()> {
    ensure!(
        cfg.max_files > 0 && cfg.max_pulses > 0,
        "processing limits must be positive"
    );
    ensure!(
        cfg.threshold_db.is_finite() && cfg.threshold_db > 0.,
        "invalid threshold"
    );
    ensure!(
        !output.exists(),
        "output already exists; choose a fresh directory to prevent overwriting results"
    );
    std::fs::create_dir_all(output.join("cuts"))?;
    let catalog = File::create(output.join("catalogue.h5"))?;
    text(
        &catalog,
        "producer",
        "tibet-meteor Rust; measured common TX reference; native-sample range/Doppler",
    )?;
    text(
        &catalog,
        "time_basis",
        "filename assumed Beijing time; header year=1970; subsecond absolute timing unverified",
    )?;
    text(&catalog,"channel_map_status","order inferred from supplied 15-panel image and header amplitude indices; mask-table contents not encoded in raw header")?;
    text(&catalog,"pointing_status","raw phase increments retained; physical alpha/beta conversion and axis/sign require instrument confirmation")?;
    text(
        &catalog,
        "classification",
        "range-Doppler and phase-compatible candidates, not a validated complete meteor catalogue",
    )?;
    scalar(&catalog, "threshold_db", cfg.threshold_db)?;
    scalar(&catalog, "tx_start_gate", TX_START as u32)?;
    scalar(&catalog, "tx_reference_channel", 2u32)?;
    scalar(&catalog, "tx_samples", TX_LEN as u32)?;
    let sources = catalog.create_group("sources")?;
    let events = catalog.create_group("events")?;
    let mut all_tracks = Vec::new();
    let mut signature = None;
    let mut segment = 0u32;
    let started = Instant::now();
    for (file_id, path) in paths.iter().take(cfg.max_files).enumerate() {
        let file_start = Instant::now();
        let mut reader = Reader::open(path)?;
        let h = reader.header.clone();
        ensure!(
            h.channels == 15
                && h.gates == 640
                && h.fs == 500_000
                && h.ipp_us == 1280
                && h.iq_start == 0,
            "unexpected campaign dimensions; inspect before processing"
        );
        ensure!(
            h.waveform == 3 && (h.width_us - 8.).abs() < 0.01,
            "waveform differs from inspected 32-sample complementary-code leakage"
        );
        let expected = vec![15, 15, 15, 15, 15, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32];
        ensure!(
            h.amplitude_index == expected,
            "channel-mask indices changed; update mapping before processing"
        );
        if signature.as_ref().is_some_and(|s| s != &h.signature()) {
            segment += 1;
        }
        signature = Some(h.signature());
        let source = sources.create_group(&format!("{file_id:04}"))?;
        save_header(&source, &h)?;
        scalar(&source, "pointing_segment", segment)?;
        let gates = (TX_START + 234..=h.gates - TX_LEN).collect::<Vec<_>>();
        let mut filter = Filter::new(gates, TX_LEN, 1, FFT)?;
        let n = h.pulses.min(cfg.max_pulses);
        let mut tracks = Vec::new();
        let mut points = Vec::new();
        let mut aliases = Vec::new();
        let mut mean_power = vec![0f64; h.gates * CHANNELS];
        let mut rti = Vec::<f32>::new();
        let lambda = C / h.frequency_hz;
        for pulse in 0..n {
            let z = reader.pulse(pulse)?;
            for (i, v) in z.iter().enumerate() {
                mean_power[i] += v.norm_sqr() as f64;
            }
            let noise = noises(&z, h.gates);
            let mut tx = z[2 * h.gates + TX_START..2 * h.gates + TX_START + TX_LEN].to_vec();
            let energy: f32 = tx.iter().map(|x| x.norm_sqr()).sum();
            ensure!(
                energy > noise[2] * TX_LEN as f32 * 3.,
                "weak/absent leakage at file {file_id} pulse {pulse}"
            );
            for v in &mut tx {
                *v /= energy.sqrt();
            }
            let mut rx = Vec::with_capacity(5 * h.gates);
            for ch in 0..5 {
                let scale = noise[ch].sqrt();
                rx.extend(
                    z[ch * h.gates..(ch + 1) * h.gates]
                        .iter()
                        .map(|v| *v / scale),
                );
            }
            let result = filter.process(&rx, 5, h.gates, &tx)?;
            if pulse % 16 == 0 {
                rti.extend(
                    result
                        .power
                        .as_chunks::<FFT>()
                        .0
                        .iter()
                        .map(|row| row.iter().copied().fold(0., f32::max) / result.median_power),
                );
            }
            let peaks = result.peaks(10f32.powf(cfg.threshold_db / 10.), 3, 5, 8)?;
            let mut accepted = Vec::new();
            for peak in peaks {
                let lag = peak.range_sample - TX_START;
                let bin = peak.doppler_bin;
                let row = result
                    .gates
                    .iter()
                    .position(|&g| g == peak.range_sample)
                    .unwrap();
                let at = |b: usize| result.power[row * FFT + b].max(1e-20).ln() as f64;
                let (left, middle, right) =
                    (at((bin + FFT - 1) % FFT), at(bin), at((bin + 1) % FFT));
                let curvature = left - 2. * middle + right;
                let delta = if curvature < 0. {
                    ((left - right) / (2. * curvature)).clamp(-0.5, 0.5)
                } else {
                    0.
                };
                let frequency =
                    result.frequency_hz(bin, h.fs as f64) + delta * h.fs as f64 / FFT as f64;
                let mut amplitudes = Vec::with_capacity(CHANNELS);
                for ch in 0..CHANNELS {
                    let mut amp = Complex32::default();
                    for k in 0..TX_LEN {
                        let ph =
                            -2. * std::f32::consts::PI * frequency as f32 * k as f32 / h.fs as f32;
                        amp += z[ch * h.gates + peak.range_sample + k]
                            * tx[k].conj()
                            * Complex32::from_polar(1., ph);
                    }
                    amplitudes.push(amp);
                }
                let beam = (0..5)
                    .max_by(|&i, &j| {
                        (amplitudes[i].norm_sqr() / noise[i])
                            .total_cmp(&(amplitudes[j].norm_sqr() / noise[j]))
                    })
                    .unwrap();
                let point = Point {
                    pulse,
                    us: h.pulse_time_us(pulse),
                    lag,
                    bin,
                    doppler: frequency,
                    range: lag as f64 * C / (2. * h.fs as f64) / 1000.,
                    statistic_db: 10. * peak.excess_over_median.log10(),
                    beam,
                    amplitudes,
                    noise: noise.clone(),
                };
                if accepted
                    .iter()
                    .any(|strong| predicted_code_alias(&point, strong, &tx, h.fs as f64))
                {
                    aliases.push(point);
                    continue;
                }
                accepted.push(point.clone());
                attach(&mut tracks, point.clone(), lambda);
                points.push(point);
            }
        }
        let tracks = tracks
            .into_iter()
            .filter(|t| {
                t.points.len() >= 8 && t.points.last().unwrap().us - t.points[0].us >= 30_000
            })
            .filter(|t| {
                let mean = t.points.iter().map(|p| p.doppler).sum::<f64>() / t.points.len() as f64;
                (range_rate(t) + mean * lambda / 2. / 1000.).abs() < 10.
            })
            .collect::<Vec<_>>();
        point_fields(&source.create_group("candidate_points")?, &points)?;
        point_fields(&source.create_group("code_alias_vetoes")?, &aliases)?;
        scalar(&source, "processed_pulses", n as u64)?;
        put(
            &source,
            "mean_raw_power",
            &[CHANNELS, h.gates],
            &mean_power.iter().map(|p| p / n as f64).collect::<Vec<_>>(),
        )?;
        put(
            &source,
            "range_doppler_peak_ratio_rti",
            &[n.div_ceil(16), h.gates - TX_LEN - (TX_START + 234) + 1],
            &rti,
        )?;
        for track in &tracks {
            let id = all_tracks.len();
            let name = format!("event_{id:05}");
            let event = events.create_group(&name)?;
            point_fields(&event, &track.points)?;
            scalar(&event, "source_file_index", file_id as u32)?;
            scalar(&event, "pointing_segment", segment)?;
            scalar(&event, "range_rate_fit_km_s", range_rate(track))?;
            scalar(
                &event,
                "mean_doppler_velocity_km_s",
                -track.points.iter().map(|p| p.doppler).sum::<f64>() / track.points.len() as f64
                    * lambda
                    / 2.
                    / 1000.,
            )?;
            text(&event, "cut_file", &format!("cuts/{name}.h5"))?;
            if cfg.write_legacy_cuts {
                cut(
                    &mut reader,
                    track,
                    &output.join("cuts").join(format!("{name}.h5")),
                )?;
            }
            all_tracks.push(Track {
                points: track.points.clone(),
            });
        }
        println!("{}: {} validated records, {} peak candidates, {} code-alias vetoes, {} tracks, {:.2} s",path.file_name().unwrap().to_string_lossy(),n,points.len(),aliases.len(),tracks.len(),file_start.elapsed().as_secs_f64());
        catalog.flush()?;
    }
    if segment == 0 {
        calibration(
            &all_tracks,
            &output.join("phasecal_boresight_provisional.h5"),
        )?;
    } else {
        eprintln!("Pointing/mask changed: joint phase calibration intentionally withheld; calibrate segments separately");
    }
    println!(
        "Total: {} candidate tracks; {:.2} s",
        all_tracks.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn point(pulse: usize, range: f64, doppler: f64) -> Point {
        Point {
            pulse,
            us: 1000000 + pulse as u64 * 1280,
            lag: 300,
            bin: 64,
            doppler,
            range,
            statistic_db: 20.,
            beam: 2,
            amplitudes: vec![Complex32::new(1., 0.); 15],
            noise: vec![1.; 15],
        }
    }
    #[test]
    fn tracks_do_not_bridge_gaps_or_same_pulse() {
        let mut t = vec![];
        attach(&mut t, point(0, 100., 1000.), 6.);
        attach(&mut t, point(1, 99.99616, 1000.), 6.);
        assert_eq!(t.len(), 1);
        attach(&mut t, point(1, 120., 1000.), 6.);
        assert_eq!(t.len(), 2);
        attach(&mut t, point(100, 99.6, 1000.), 6.);
        assert_eq!(t.len(), 3);
    }
}
