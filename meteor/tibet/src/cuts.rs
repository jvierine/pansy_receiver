//! Lossless ten-subarray cuts. A common TX reference preserves receiver phases.
use crate::pipeline::{put, save_header, scalar, text};
use anyhow::{ensure, Context, Result};
use hdf5::{File, Group, H5Type};
use pansy_meteor_matched_filter::Filter;
use rustfft::num_complex::Complex32;
use std::path::Path;
use tibet_meteor::{u16at, Reader};

const START: usize = 7;
const LEN: usize = 32;
const RX: std::ops::Range<usize> = 5..15;

pub fn verify(directory: &Path) -> Result<()> {
    let mut count = 0;
    let mut largest = 0f32;
    let directory = if directory.join("cuts").is_dir() {
        directory.join("cuts")
    } else {
        directory.to_owned()
    };
    for entry in std::fs::read_dir(directory)? {
        let path = entry?.path();
        if path.extension().and_then(|s| s.to_str()) != Some("h5") {
            continue;
        }
        let h = File::open(&path)?;
        let rx = tibet_meteor::voltage::read(&h, "zrx_echoes")?;
        let tx = tibet_meteor::voltage::read(&h, "ztx_history")?;
        let shape = if h.link_exists("zrx_echoes_iq_i16") {
            h.dataset("zrx_echoes_iq_i16")?.shape()
        } else {
            h.dataset("zrx_echoes_iq_i32")?.shape()
        };
        ensure!(
            shape.len() == 4 && shape[1] == 10 && shape[3] == 2,
            "expected 10-channel interleaved cuts"
        );
        let width = shape[2];
        let pulse = h.dataset("pulse_index")?.read_raw::<u64>()?;
        let lo = h.dataset("range_start_gate")?.read_raw::<u32>()?[0] as usize;
        let fs = h.dataset("sample_rate_hz")?.read_raw::<u32>()?[0];
        let d = h.group("detections")?;
        let detected = d.dataset("pulse_index")?.read_raw::<u64>()?;
        let lag = d.dataset("delay_samples")?.read_raw::<u32>()?;
        let fd = d.dataset("doppler_hz")?.read_raw::<f64>()?;
        let reference = d.dataset("matched_amplitude_re_im")?.read_raw::<f32>()?;
        for j in 0..detected.len() {
            let row = pulse
                .binary_search(&detected[j])
                .map_err(|_|anyhow::anyhow!("detection missing from cut"))?;
            let start = lag[j] as usize + START - lo;
            ensure!(start + LEN <= width, "detection outside cut");
            let t = complex(&tx[row * 3 * LEN * 2..(row * 3 + 1) * LEN * 2]);
            for c in 0..10 {
                let k = 2 * ((row * 10 + c) * width + start);
                let a = matched(&complex(&rx[k..k + 2 * LEN]), &t, fd[j], fs);
                let b =
                    Complex32::new(reference[(j * 10 + c) * 2], reference[(j * 10 + c) * 2 + 1]);
                let error = (a - b).norm() / b.norm().max(1.);
                largest = largest.max(error);
                ensure!(
                    error < 5e-5,
                    "{}: detection {j}, channel {c}, replay error {error}",
                    path.display()
                );
                count += 1;
            }
        }
    }
    ensure!(count > 0, "no detections verified");
    println!("PASS: {count} matched amplitudes; maximum relative error {largest:.3e}");
    Ok(())
}

fn pointing_header_fields(raw: &[u8]) -> (u16, u16) {
    // IQ.docx includes the microsecond word missing from the older MATLAB reader.
    (u16at(raw, 42), u16at(raw, 44))
}

fn copy<T: H5Type>(from: &Group, to: &Group, name: &str) -> Result<()> {
    let d = from.dataset(name)?;
    put(to, name, &d.shape(), &d.read_raw::<T>()?)
}

fn complex(z: &[i32]) -> Vec<Complex32> {
    z.as_chunks::<2>()
        .0
        .iter()
        .map(|v| Complex32::new(v[0] as f32, v[1] as f32))
        .collect()
}

fn matched(rx: &[Complex32], tx: &[Complex32], fd: f64, fs: u32) -> Complex32 {
    let norm = tx.iter().map(|v| v.norm_sqr()).sum::<f32>().sqrt();
    let mut a = Complex32::default();
    for k in 0..LEN {
        let phase = -2. * std::f32::consts::PI * fd as f32 * k as f32 / fs as f32;
        a += rx[k] * (tx[k] / norm).conj() * Complex32::from_polar(1., phase);
    }
    a
}

fn schedule_field<T: H5Type + Clone>(
    schedule: &Option<Group>,
    name: &str,
    n: usize,
    default: T,
) -> Result<Vec<T>> {
    if let Some(g) = schedule {
        let d = g.dataset(name)?;
        ensure!(
            d.shape() == [n],
            "schedule {name} must cover every source pulse"
        );
        Ok(d.read_raw::<T>()?)
    } else {
        Ok(vec![default; n])
    }
}

pub fn run(run: &Path, output: &Path, schedule_path: Option<&Path>) -> Result<()> {
    ensure!(!output.exists(), "cut output exists; use a new directory");
    let partial = output.with_extension("partial");
    ensure!(
        !partial.exists(),
        "incomplete cut output exists: {}",
        partial.display()
    );
    let cat = File::open(run.join("catalogue.h5"))?;
    let schedule = schedule_path.map(File::open).transpose()?;
    std::fs::create_dir_all(partial.join("cuts"))?;
    let index = File::create(partial.join("index.h5"))?;
    text(&index, "schema", "tibet_meteor_cut_v2")?;
    text(
        &index,
        "catalogue_path",
        &run.join("catalogue.h5").to_string_lossy(),
    )?;
    text(&index, "range_alias_status", "test previous two sampled codes at observed delay/Doppler; even-IPP aliases remain unresolved; no unaliased-range claim")?;
    let events = index.create_group("events")?;
    let sources = index.create_group("sources")?;
    for name in cat.group("events")?.member_names()? {
        let event = cat.group(&format!("events/{name}"))?;
        let file_id = event.dataset("source_file_index")?.read_scalar::<u32>()?;
        let source_name = format!("{file_id:04}");
        let source = cat.group(&format!("sources/{source_name}"))?;
        let path = String::from_utf8(source.dataset("source_path")?.read_raw::<u8>()?)?;
        let mut reader = Reader::open(Path::new(&path))?;
        let h = reader.header.clone();
        ensure!(
            h.raw == source.dataset("raw_header")?.read_raw::<u8>()?,
            "source header changed"
        );
        ensure!(
            h.channels == 15
                && h.fs == 500_000
                && h.gates == 640
                && h.ipp_us == 1280
                && h.iq_start == 0
                && h.waveform == 3
                && (h.width_us - 8.).abs() < 0.01,
            "unvalidated waveform/dimensions"
        );
        ensure!(
            h.amplitude_index == [15, 15, 15, 15, 15, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32],
            "channel map changed"
        );
        let per_pulse = schedule
            .as_ref()
            .map(|s| s.group(&source_name))
            .transpose()?;
        let (beam_total, current_frame) = pointing_header_fields(&h.raw);
        ensure!(
            beam_total <= 1 || per_pulse.is_some(),
            "multi-frame pointing requires an explicit per-pulse schedule"
        );
        let alphas = schedule_field(&per_pulse, "tx_alpha_code", h.pulses, h.tx_alpha)?;
        let betas = schedule_field(&per_pulse, "tx_beta_code", h.pulses, h.tx_beta)?;
        let beam_ids = schedule_field(&per_pulse, "tx_beam_id", h.pulses, 0i32)?;
        let zenith_default = if h.tx_alpha == 0 && h.tx_beta == 0 {
            0.
        } else {
            f64::NAN
        };
        let zenith = schedule_field(&per_pulse, "tx_zenith_angle_deg", h.pulses, zenith_default)?;
        let azimuth = schedule_field(&per_pulse, "tx_azimuth_deg", h.pulses, f64::NAN)?;
        if per_pulse.is_some() {
            ensure!(
                zenith.iter().zip(&azimuth).all(|(&z, &a)| z.is_finite()
                    && (0. ..=90.).contains(&z)
                    && (z == 0. || (a.is_finite() && (0. ..360.).contains(&a)))),
                "invalid physical TX pointing schedule"
            );
        }
        ensure!(
            beam_ids.iter().all(|&v| v >= 0),
            "negative schedule beam ID"
        );
        if !sources.link_exists(&source_name) {
            let dest = sources.create_group(&source_name)?;
            save_header(&dest, &h)?;
            scalar(&dest, "header_beam_total", beam_total)?;
            scalar(&dest, "header_current_frame", current_frame)?;
            text(
                &dest,
                "pointing_source",
                if schedule.is_some() {
                    "explicit per-pulse external schedule; hardware codes retained"
                } else {
                    "per-file header; within-file switching not encoded or independently confirmed"
                },
            )?;
        }
        let pulses = event.dataset("pulse_index")?.read_raw::<u64>()?;
        let lags = event.dataset("delay_samples")?.read_raw::<u32>()?;
        let fd = event.dataset("doppler_hz")?.read_raw::<f64>()?;
        ensure!(
            !pulses.is_empty() && pulses.len() == lags.len() && pulses.len() == fd.len(),
            "invalid event arrays"
        );
        ensure!(
            pulses.windows(2).all(|p| p[1] > p[0]) && *pulses.last().unwrap() < h.pulses as u64,
            "invalid event pulse order"
        );
        let pad = (80_000usize).div_ceil(h.ipp_us as usize);
        let first = (pulses[0] as usize).saturating_sub(pad);
        let end = (*pulses.last().unwrap() as usize + pad + 1).min(h.pulses);
        let start_gate = (*lags.iter().min().unwrap() as usize + START).saturating_sub(64);
        let end_gate = (*lags.iter().max().unwrap() as usize + START + LEN + 64).min(h.gates);
        let width = end_gate - start_gate;
        let n = end - first;
        let f = File::create(partial.join("cuts").join(format!("{name}.h5")))?;
        save_header(&f, &h)?;
        text(&f, "schema", "tibet_meteor_cut_v2")?;
        text(&f, "event_id", &name)?;
        text(&f, "channel_names", "S0 S1 S2 L0 L1 L2 L3 L4 L5 L6")?;
        put(
            &f,
            "source_channel_indices",
            &[10],
            &(5u16..15).collect::<Vec<_>>(),
        )?;
        scalar(&f, "tx_reference_source_channel", 2u16)?;
        scalar(&f, "source_file_index", file_id)?;
        scalar(&f, "frame_sequence", h.frame)?;
        scalar(&f, "iq_start_header", h.iq_start)?;
        scalar(&f, "range_start_gate", start_gate as u32)?;
        scalar(&f, "tx_start_gate", START as u32)?;
        scalar(&f, "txlen", LEN as u32)?;
        scalar(&f, "samples_per_baud", 4u32)?;
        scalar(&f, "baud_count", 8u32)?;
        scalar(&f, "nominal_code_duration_us", 64u32)?;
        scalar(
            &f,
            "padding_before_clipped",
            u8::from(pulses[0] < pad as u64),
        )?;
        scalar(
            &f,
            "padding_after_clipped",
            u8::from(*pulses.last().unwrap() as usize + pad + 1 > h.pulses),
        )?;
        text(
            &f,
            "pointing_source",
            if schedule.is_some() {
                "explicit per-pulse schedule"
            } else {
                "per-file header expanded to pulse indices; no within-file pointing telemetry"
            },
        )?;
        text(&f,"tx_reference_status","sampled receive leakage, not a calibrated measurement of radiated RF power; same channel used for all echo channels")?;
        text(
            &f,
            "phase_status",
            "original integer I/Q; no phase correction applied",
        )?;
        text(&f,"alias_status","odd-IPP code comparison at detected delay and Doppler; even-IPP physical range remains ambiguous; 0/1/2 are trial orders, not chosen ranges")?;
        let mut echoes = Vec::with_capacity(n * 10 * width * 2);
        let mut tx_history = vec![0i32; n * 3 * LEN * 2];
        let mut valid = vec![0u8; n * 3];
        let mut origin_alpha = vec![0u16; n * 3];
        let mut origin_beta = vec![0u16; n * 3];
        let mut origin_beam = vec![-1i32; n * 3];
        let mut origin_zenith = vec![f64::NAN; n * 3];
        let mut origin_azimuth = vec![f64::NAN; n * 3];
        let mut decoded = vec![0i8; n * 3 * 8];
        for (i, p) in (first..end).enumerate() {
            let z = reader.pulse_i32(p)?;
            for ch in RX {
                echoes.extend_from_slice(
                    &z[2 * (ch * h.gates + start_gate)..2 * (ch * h.gates + end_gate)],
                );
            }
            for order in 0..3 {
                let Some(origin) = p.checked_sub(order) else {
                    continue;
                };
                let prior = if order == 0 {
                    z.clone()
                } else {
                    reader.pulse_i32(origin)?
                };
                let tx = &prior[2 * (2 * h.gates + START)..2 * (2 * h.gates + START + LEN)];
                let offset = (i * 3 + order) * LEN * 2;
                tx_history[offset..offset + LEN * 2].copy_from_slice(tx);
                valid[i * 3 + order] = 1;
                origin_alpha[i * 3 + order] = alphas[origin];
                origin_beta[i * 3 + order] = betas[origin];
                origin_beam[i * 3 + order] = beam_ids[origin];
                origin_zenith[i * 3 + order] = zenith[origin];
                origin_azimuth[i * 3 + order] = azimuth[origin];
                let c = complex(tx);
                let ref_baud: Complex32 = c[..4].iter().sum();
                for b in 0..8 {
                    let value: Complex32 = c[b * 4..b * 4 + 4].iter().sum();
                    decoded[(i * 3 + order) * 8 + b] = if (value * ref_baud.conj()).re >= 0. {
                        1
                    } else {
                        -1
                    };
                }
            }
        }
        tibet_meteor::voltage::write(&f, "zrx_echoes", &[n, 10, width, 2], &echoes)?;
        tibet_meteor::voltage::write(&f, "ztx_history", &[n, 3, LEN, 2], &tx_history)?;
        text(&f,"voltage_encoding","interleaved int16, HDF5 shuffle + deflate level 4; sparse zero-based C-order overflow indices and low/high signed int16 words reconstruct every original int32 exactly")?;
        let echoes = tibet_meteor::voltage::read(&f, "zrx_echoes")?;
        let tx_history = tibet_meteor::voltage::read(&f, "ztx_history")?;
        put(&f, "tx_history_valid", &[n, 3], &valid)?;
        put(&f, "tx_history_alpha_code", &[n, 3], &origin_alpha)?;
        put(&f, "tx_history_beta_code", &[n, 3], &origin_beta)?;
        put(&f, "tx_history_beam_id", &[n, 3], &origin_beam)?;
        put(&f, "tx_history_zenith_angle_deg", &[n, 3], &origin_zenith)?;
        put(&f, "tx_history_azimuth_deg", &[n, 3], &origin_azimuth)?;
        put(&f, "tx_history_decoded_baud_sign", &[n, 3, 8], &decoded)?;
        put(
            &f,
            "pulse_index",
            &[n],
            &(first as u64..end as u64).collect::<Vec<_>>(),
        )?;
        put(
            &f,
            "record_counter",
            &[n],
            &(first as u32 + 1..end as u32 + 1).collect::<Vec<_>>(),
        )?;
        put(
            &f,
            "ipp_start_time_us",
            &[n],
            &(first..end).map(|p| h.pulse_time_us(p)).collect::<Vec<_>>(),
        )?;
        put(&f, "tx_alpha_per_pulse", &[n], &alphas[first..end])?;
        put(&f, "tx_beta_per_pulse", &[n], &betas[first..end])?;
        put(&f, "tx_beam_id", &[n], &beam_ids[first..end])?;
        put(&f, "tx_zenith_angle_deg", &[n], &zenith[first..end])?;
        put(&f, "tx_azimuth_deg", &[n], &azimuth[first..end])?;
        text(&f,"tx_angle_status", "zero phase increments identified as zenith; nonzero physical angles require an external schedule; zenith azimuth is undefined")?;
        scalar(
            &f,
            "within_cut_tx_pointing_changes",
            (first + 1..end)
                .filter(|&p| {
                    alphas[p] != alphas[p - 1]
                        || betas[p] != betas[p - 1]
                        || beam_ids[p] != beam_ids[p - 1]
                })
                .count() as u64,
        )?;
        let d = f.create_group("detections")?;
        for field in ["pulse_index", "tx_time_us"] {
            copy::<u64>(&event, &d, field)?;
        }
        for field in ["delay_samples", "doppler_bin"] {
            copy::<u32>(&event, &d, field)?;
        }
        for field in ["doppler_hz", "range_km_nominal"] {
            copy::<f64>(&event, &d, field)?;
        }
        copy::<f32>(&event, &d, "detector_excess_over_median_db")?;
        let noise = event.dataset("noise_power")?.read_raw::<f32>()?;
        let amp = event
            .dataset("matched_amplitude_re_im")?
            .read_raw::<f32>()?;
        let noise10 = noise
            .as_chunks::<15>()
            .0
            .iter()
            .flat_map(|row| row[5..].iter().copied())
            .collect::<Vec<_>>();
        let amp10 = amp
            .as_chunks::<30>()
            .0
            .iter()
            .flat_map(|row| row[10..].iter().copied())
            .collect::<Vec<_>>();
        put(&d, "noise_power", &[pulses.len(), 10], &noise10)?;
        put(
            &d,
            "matched_amplitude_re_im",
            &[pulses.len(), 10, 2],
            &amp10,
        )?;
        let mut alias_power = vec![f64::NAN; pulses.len() * 3];
        let mut rd_power = vec![f32::NAN; pulses.len() * 3];
        let mut rd_delay = vec![f64::NAN; pulses.len() * 3];
        let mut rd_fd = vec![f64::NAN; pulses.len() * 3];
        let mut filter = Filter::new((0..=width - LEN).collect(), LEN, 1, 128)?;
        let mut replay = Vec::new();
        let mut max_error = 0f32;
        for j in 0..pulses.len() {
            let i = pulses[j] as usize - first;
            let offset = lags[j] as usize + START - start_gate;
            let mut rx_whitened = Vec::with_capacity(10 * width);
            for ch in 0..10 {
                let start = (i * 10 + ch) * width * 2;
                rx_whitened.extend(
                    complex(&echoes[start..start + width * 2])
                        .into_iter()
                        .map(|v| v / noise10[j * 10 + ch].sqrt()),
                );
            }
            for order in 0..3 {
                if valid[i * 3 + order] == 0 {
                    continue;
                }
                let tx =
                    complex(&tx_history[(i * 3 + order) * LEN * 2..(i * 3 + order + 1) * LEN * 2]);
                let norm = tx.iter().map(|v| v.norm_sqr()).sum::<f32>().sqrt();
                let normalized = tx.iter().map(|v| *v / norm).collect::<Vec<_>>();
                let result = filter.process(&rx_whitened, 10, width, &normalized)?;
                let (best, &power_rd) = result
                    .power
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.total_cmp(b.1))
                    .unwrap();
                let bin = best % 128;
                let row = best / 128;
                let at = |b: usize| result.power[row * 128 + b].max(1e-20).ln() as f64;
                let (left, middle, right) = (at((bin + 127) % 128), at(bin), at((bin + 1) % 128));
                let curvature = left - 2. * middle + right;
                let delta = if curvature < 0. {
                    ((left - right) / (2. * curvature)).clamp(-0.5, 0.5)
                } else {
                    0.
                };
                rd_power[j * 3 + order] = power_rd;
                rd_delay[j * 3 + order] = (start_gate + result.gates[row] - START) as f64;
                rd_fd[j * 3 + order] =
                    result.frequency_hz(bin, h.fs as f64) + delta * h.fs as f64 / 128.;
                let mut power = 0.;
                for ch in 0..10 {
                    let k = ((i * 10 + ch) * width + offset) * 2;
                    let a = matched(&complex(&echoes[k..k + LEN * 2]), &tx, fd[j], h.fs);
                    power += a.norm_sqr() as f64 / noise10[j * 10 + ch] as f64;
                    if order == 0 {
                        replay.extend([a.re, a.im]);
                        let reference =
                            Complex32::new(amp10[(j * 10 + ch) * 2], amp10[(j * 10 + ch) * 2 + 1]);
                        let error = (a - reference).norm() / reference.norm().max(1.);
                        max_error = max_error.max(error);
                        ensure!(error < 1e-5, "cut replay differs from detection amplitudes");
                    }
                }
                alias_power[j * 3 + order] = power;
            }
        }
        put(
            &d,
            "cut_replayed_amplitude_re_im",
            &[pulses.len(), 10, 2],
            &replay,
        )?;
        put(&d, "alias_order", &[3], &[0u32, 1, 2])?;
        put(
            &d,
            "alias_matched_power_at_observed_delay_doppler",
            &[pulses.len(), 3],
            &alias_power,
        )?;
        scalar(
            &d,
            "range_alias_interval_km",
            299792458. * h.ipp_us as f64 * 1e-6 / 2000.,
        )?;
        put(&d, "alias_rd_peak_power", &[pulses.len(), 3], &rd_power)?;
        put(
            &d,
            "alias_rd_peak_delay_samples_mod_ipp",
            &[pulses.len(), 3],
            &rd_delay,
        )?;
        put(&d, "alias_rd_peak_doppler_hz", &[pulses.len(), 3], &rd_fd)?;
        text(&d,"alias_rd_peak_scope","maximum across each padded echo cut and all 128 Doppler bins; can select a different target or code sidelobe; diagnostic only, not an alias decision")?;
        scalar(&d, "max_relative_replay_error", max_error)?;
        f.flush()?;
        let out = events.create_group(&name)?;
        text(&out, "cut_file", &format!("cuts/{name}.h5"))?;
        scalar(&out, "source_file_index", file_id)?;
        scalar(&out, "detection_points", pulses.len() as u64)?;
        scalar(&out, "max_relative_replay_error", max_error)?;
        println!("{name}: {n} IPPs x 10 x {width} echo samples; replay error {max_error:.2e}");
    }
    scalar(&index, "complete", 1u8)?;
    index.flush()?;
    drop(index);
    std::fs::rename(&partial, output).context("publish completed cut directory")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn pointing_offsets_include_microsecond_word() {
        let mut raw = vec![0u8; 8192];
        raw[40..42].copy_from_slice(&19u16.to_le_bytes());
        raw[42..44].copy_from_slice(&5u16.to_le_bytes());
        raw[44..46].copy_from_slice(&3u16.to_le_bytes());
        assert_eq!(pointing_header_fields(&raw), (5, 3));
    }
    #[test]
    fn earlier_pulse_alias_codes() {
        let a = [1., -1., 1., 1., 1., -1., -1., -1.];
        let b = [1., -1., 1., 1., -1., 1., 1., 1.];
        let code = |v: &[f32; 8]| {
            v.iter()
                .flat_map(|&s| vec![Complex32::new(s, 0.); 4])
                .collect::<Vec<_>>()
        };
        let ta = code(&a);
        let tb = code(&b);
        let fd = 23000.;
        let mut powers = Vec::new();
        for alias_order in 0..3 {
            let emitted = if alias_order % 2 == 0 { &ta } else { &tb };
            let rx = emitted
                .iter()
                .enumerate()
                .map(|(k, &v)| {
                    v * Complex32::from_polar(
                        1.,
                        2. * std::f32::consts::PI * fd as f32 * k as f32 / 500000.,
                    )
                })
                .collect::<Vec<_>>();
            let p0 = matched(&rx, &ta, fd, 500000).norm_sqr();
            let p1 = matched(&rx, &tb, fd, 500000).norm_sqr();
            powers.push([p0, p1]);
            if alias_order == 1 {
                assert!(p1 > 31.9 && p0 < 1e-8);
            } else {
                assert!(p0 > 31.9 && p1 < 1e-8);
            }
        }
        // The two-IPP return is exactly indistinguishable from order zero here.
        assert_eq!(powers[0], powers[2]);
    }
}
