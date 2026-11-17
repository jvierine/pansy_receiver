//! Conservative paired-code MST reduction. Masks and signed spectra are retained.
use anyhow::{ensure, Result};
use hdf5::{File, Group, H5Type};
use pansy_radar_spectra::{bartlett_spectrum as spectrum, masked_acf as acf};
use rustfft::{num_complex::Complex32 as C, FftPlanner};
use std::path::Path;
use tibet_meteor::{raw_paths, Reader};

const START: usize = 39; // Entire sampled transmit pulse has finished.
const STEP: usize = 4; // One 8-us baud, not oversampled range bins.
const LAGS: usize = 128;
const NF: usize = 512;
pub fn fine(directory: &Path) -> Result<()> {
    let mut planner = FftPlanner::new();
    let forward = planner.plan_fft_forward(8192);
    let inverse = planner.plan_fft_inverse(8192);
    let spectral = planner.plan_fft_forward(4096);
    for entry in std::fs::read_dir(directory)? {
        let path = entry?.path();
        if path.extension().and_then(|s| s.to_str()) != Some("h5")
            || path.file_name().unwrap() == "summary.h5"
        {
            continue;
        }
        let h = File::open_rw(&path)?;
        let ds = h.dataset("decoded_iq")?;
        let shape = ds.shape();
        let (nt, nr) = (shape[0], shape[2]);
        let raw = ds.read_raw::<f32>()?;
        let mask = h.dataset("mask_reason")?.read_raw::<u8>()?;
        let dt = h.dataset("pair_interval_s")?.read_scalar::<f64>()?;
        put(
            &h,
            "fine_doppler_hz",
            &[4096],
            &(0..4096)
                .map(|i| (i as f64 - 2048.) / (4096. * dt))
                .collect::<Vec<_>>(),
        )?;
        for c in 0..5 {
            let group = h.create_group(&format!("fine_channel_{c}"))?;
            let mut clean = vec![];
            let mut before = vec![];
            let mut counts = vec![];
            for r in 0..nr {
                let z = (0..nt)
                    .map(|p| {
                        let j = 2 * ((p * 5 + c) * nr + r);
                        C::new(raw[j], raw[j + 1])
                    })
                    .collect::<Vec<_>>();
                let good = (0..nt).map(|p| mask[p * nr + r] == 0).collect::<Vec<_>>();
                let count = good.iter().filter(|&&v| v).count();
                let mean = z
                    .iter()
                    .zip(&good)
                    .filter(|(_, v)| **v)
                    .map(|(z, _)| *z)
                    .sum::<C>()
                    / count.max(1) as f32;
                let centred = z.iter().map(|v| *v - mean).collect::<Vec<_>>();
                let (a, n) = pansy_radar_spectra::masked_acf_fft(
                    &centred,
                    &good,
                    1024,
                    forward.as_ref(),
                    inverse.as_ref(),
                );
                let (b, _) = pansy_radar_spectra::masked_acf_fft(
                    &z,
                    &vec![true; nt],
                    1024,
                    forward.as_ref(),
                    inverse.as_ref(),
                );
                clean.extend(spectrum(&a, spectral.as_ref()));
                before.extend(spectrum(&b, spectral.as_ref()));
                counts.extend(n);
            }
            put(&group, "spectrum_clean", &[nr, 4096], &clean)?;
            put(&group, "spectrum_before", &[nr, 4096], &before)?;
            put(&group, "valid_pair_count", &[nr, 1024], &counts)?;
        }
        println!("Fine spectra: {}", path.display());
    }
    Ok(())
}
fn median(x: &mut [f32]) -> f32 {
    let k = x.len() / 2;
    x.select_nth_unstable_by(k, |a, b| a.total_cmp(b));
    x[k]
}
fn put<T: H5Type>(g: &Group, key: &str, dims: &[usize], v: &[T]) -> Result<()> {
    g.new_dataset::<T>().shape(dims).create(key)?.write_raw(v)?;
    Ok(())
}

pub fn run(input: &Path, detections: &Path, output: &Path) -> Result<()> {
    ensure!(!output.exists(), "output already exists");
    std::fs::create_dir_all(output)?;
    let cat = File::open(detections.join("catalogue.h5"))?;
    let events = cat.group("events")?;
    let mut planner = FftPlanner::new();
    let fast_fft = planner.plan_fft_forward(1024);
    let fast_ifft = planner.plan_fft_inverse(1024);
    let slow_fft = planner.plan_fft_forward(NF);
    for path in raw_paths(input)? {
        let mut reader = Reader::open(&path)?;
        let h = reader.header.clone();
        ensure!(
            h.channels == 15 && h.gates == 640 && h.fs == 500000,
            "unvalidated mode"
        );
        ensure!(h.pulses % 2 == 0, "incomplete complementary pair");
        let nt = h.pulses / 2;
        let gates = (START..=h.gates - 32).step_by(STEP).collect::<Vec<_>>();
        let nr = gates.len();
        let mut decoded = vec![C::default(); nt * 5 * nr];
        let mut noise = vec![0.; h.pulses * 5];
        let mut raw_level = vec![0.; h.pulses];
        let mut impulse = vec![0f32; h.pulses * nr];
        let mut tx_quality = vec![0.; h.pulses];
        let mut fast_psd = vec![0f64; 1024];
        let mut ambiguity = vec![0f64; 65];
        for pulse in 0..h.pulses {
            let raw = reader.pulse(pulse)?;
            for c in 0..5 {
                noise[pulse * 5 + c] = median(
                    &mut (250..640)
                        .map(|g| raw[c * 640 + g].norm_sqr())
                        .collect::<Vec<_>>(),
                ) / 2f32.ln();
            }
            raw_level[pulse] = noise[pulse * 5 + 2];
            let tx = &raw[2 * 640 + 7..2 * 640 + 39];
            let energy = tx.iter().map(|v| v.norm_sqr()).sum::<f32>();
            tx_quality[pulse] = energy / (32. * raw_level[pulse]);
            let mut reference = vec![C::default(); 1024];
            for j in 0..32 {
                reference[j] = tx[j] / energy.sqrt();
            }
            if pulse < 2 {
                for (lag, ambiguity_value) in ambiguity.iter_mut().enumerate() {
                    let d = lag as i32 - 32;
                    let v = (0..32)
                        .filter_map(|j| {
                            let k = j as i32 + d;
                            (0..32)
                                .contains(&k)
                                .then(|| reference[j] * reference[k as usize].conj())
                        })
                        .sum::<C>();
                    *ambiguity_value += v.re as f64 / 2.;
                }
            }
            fast_fft.process(&mut reference);
            for c in 0..5 {
                let mut x = vec![C::default(); 1024];
                x[..640].copy_from_slice(&raw[c * 640..(c + 1) * 640]);
                if c == 2 {
                    let mut off = vec![C::default(); 1024];
                    off[..390].copy_from_slice(&raw[c * 640 + 250..(c + 1) * 640]);
                    fast_fft.process(&mut off);
                    for j in 0..1024 {
                        fast_psd[j] += off[j].norm_sqr() as f64;
                    }
                    for (r, &g) in gates.iter().enumerate() {
                        impulse[pulse * nr + r] = raw[c * 640 + g..c * 640 + g + 32]
                            .iter()
                            .map(|v| v.norm_sqr() / raw_level[pulse].max(1.))
                            .fold(0., f32::max);
                    }
                }
                fast_fft.process(&mut x);
                for j in 0..1024 {
                    x[j] *= reference[j].conj();
                }
                fast_ifft.process(&mut x);
                for (r, &g) in gates.iter().enumerate() {
                    decoded[(pulse / 2 * 5 + c) * nr + r] += x[g] / (1024. * 2f32.sqrt());
                }
            }
        }
        // A global pulse gate uses far-range median power, not low-altitude MST power.
        let centre = median(&mut raw_level.clone());
        let impulse_baseline = (0..nr)
            .map(|r| {
                median(
                    &mut (0..h.pulses)
                        .map(|p| impulse[p * nr + r])
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<Vec<_>>();
        let mut reason = vec![0u8; nt * nr];
        for p in 0..nt {
            let bad = raw_level[2 * p] > 3. * centre
                || raw_level[2 * p + 1] > 3. * centre
                || tx_quality[2 * p] < 3.
                || tx_quality[2 * p + 1] < 3.;
            for r in 0..nr {
                if bad {
                    reason[p * nr + r] |= 1;
                }
                let threshold = 1000f32.max(20. * impulse_baseline[r]);
                if impulse[2 * p * nr + r] > threshold || impulse[(2 * p + 1) * nr + r] > threshold
                {
                    reason[p * nr + r] |= 2;
                }
            }
        }
        // Reject detected meteor voltage windows, including the full code sidelobe footprint.
        for key in events.member_names()? {
            let e = events.group(&key)?;
            let source_index = e.dataset("source_file_index")?.read_scalar::<u32>()?;
            let source = cat
                .dataset(&format!("sources/{source_index:04}/source_path"))?
                .read_raw::<u8>()?;
            let source = String::from_utf8(source)?;
            if !source.contains(&path.file_name().unwrap().to_string_lossy().to_string()) {
                continue;
            }
            let pulse = e.dataset("pulse_index")?.read_raw::<u64>()?;
            let delay = e.dataset("delay_samples")?.read_raw::<f64>()?;
            for (&p, &d) in pulse.iter().zip(&delay) {
                for pair in (p as usize / 2).saturating_sub(3)..=(p as usize / 2 + 3).min(nt - 1) {
                    for (r, &gate) in gates.iter().enumerate() {
                        if (gate as f64 - (d + 7.)).abs() <= 40. {
                            reason[pair * nr + r] |= 4;
                        }
                    }
                }
            }
        }
        // Conservative transient-target excision; persistent narrow lines are NOT silently removed.
        let mut target = vec![false; nt * nr];
        for r in 0..nr {
            let centre = (0..5)
                .map(|c| {
                    C::new(
                        median(
                            &mut (0..nt)
                                .map(|p| decoded[(p * 5 + c) * nr + r].re)
                                .collect::<Vec<_>>(),
                        ),
                        median(
                            &mut (0..nt)
                                .map(|p| decoded[(p * 5 + c) * nr + r].im)
                                .collect::<Vec<_>>(),
                        ),
                    )
                })
                .collect::<Vec<_>>();
            let power = (0..nt)
                .map(|p| {
                    (0..5)
                        .map(|c| {
                            (decoded[(p * 5 + c) * nr + r] - centre[c]).norm_sqr()
                                / noise[2 * p * 5 + c].max(1.)
                        })
                        .sum::<f32>()
                        / 5.
                })
                .collect::<Vec<_>>();
            let baseline = median(&mut power.clone()).max(1.);
            for (p, &value) in power.iter().enumerate() {
                if value > 20. * baseline {
                    for q in p.saturating_sub(2)..=(p + 2).min(nt - 1) {
                        for s in r.saturating_sub(8)..=(r + 8).min(nr - 1) {
                            target[q * nr + s] = true;
                        }
                    }
                }
            }
        }
        for i in 0..reason.len() {
            if target[i] {
                reason[i] |= 8;
            }
        }
        let out = File::create(output.join(format!(
            "{}.h5",
            path.file_stem().unwrap().to_string_lossy()
        )))?;
        put(
            &out,
            "range_km",
            &[nr],
            &gates
                .iter()
                .map(|&g| (g - 7) as f64 * 299792458. / (2. * h.fs as f64) / 1000.)
                .collect::<Vec<_>>(),
        )?;
        put(
            &out,
            "doppler_hz",
            &[NF],
            &(0..NF)
                .map(|i| (i as f64 - NF as f64 / 2.) / (NF as f64 * 2. * h.ipp_us as f64 * 1e-6))
                .collect::<Vec<_>>(),
        )?;
        put(&out, "mask_reason", &[nt, nr], &reason)?;
        put(&out, "raw_far_range_noise", &[h.pulses, 5], &noise)?;
        put(&out, "tx_reference_snr", &[h.pulses], &tx_quality)?;
        put(&out, "fast_time_psd", &[1024], &fast_psd)?;
        put(&out, "measured_pair_range_ambiguity", &[65], &ambiguity)?;
        put(&out, "source_header", &[h.raw.len()], &h.raw)?;
        put(
            &out,
            "decoded_iq",
            &[nt, 5, nr, 2],
            &decoded
                .iter()
                .flat_map(|v| [v.re, v.im])
                .collect::<Vec<_>>(),
        )?;
        put(&out, "epoch_us", &[], &[h.epoch_us])?;
        put(&out, "pair_interval_s", &[], &[2. * h.ipp_us as f64 * 1e-6])?;
        put(&out, "frequency_hz", &[], &[h.frequency_hz])?;
        for c in 0..5 {
            let group = out.create_group(&format!("channel_{c}"))?;
            let mut clean = vec![];
            let mut before = vec![];
            let mut pairs = vec![];
            let mut real = vec![];
            let mut imag = vec![];
            let mut coherent = vec![];
            for r in 0..nr {
                let z = (0..nt)
                    .map(|p| decoded[(p * 5 + c) * nr + r])
                    .collect::<Vec<_>>();
                let good = (0..nt).map(|p| reason[p * nr + r] == 0).collect::<Vec<_>>();
                let count = good.iter().filter(|&&v| v).count();
                let mean = z
                    .iter()
                    .zip(&good)
                    .filter(|(_, v)| **v)
                    .map(|(z, _)| *z)
                    .sum::<C>()
                    / count.max(1) as f32;
                coherent.extend([mean.re, mean.im]);
                let centred = z.iter().map(|v| *v - mean).collect::<Vec<_>>();
                let (a, n) = acf(&centred, &good, LAGS);
                let (b, _) = acf(&z, &vec![true; nt], LAGS);
                clean.extend(spectrum(&a, slow_fft.as_ref()));
                before.extend(spectrum(&b, slow_fft.as_ref()));
                pairs.extend(n);
                real.extend(a.iter().map(|v| v.re));
                imag.extend(a.iter().map(|v| v.im));
            }
            put(&group, "spectrum_clean", &[nr, NF], &clean)?;
            put(&group, "spectrum_before", &[nr, NF], &before)?;
            put(&group, "acf_real", &[nr, LAGS], &real)?;
            put(&group, "acf_imag", &[nr, LAGS], &imag)?;
            put(&group, "valid_pair_count", &[nr, LAGS], &pairs)?;
            put(&group, "removed_coherent_mean_iq", &[nr, 2], &coherent)?;
        }
        println!(
            "{}: masked {:.2}% of range/pair samples",
            path.display(),
            100. * reason.iter().filter(|&&v| v != 0).count() as f64 / reason.len() as f64
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn missing_samples_preserve_tone_acf() {
        let z = (0..1000)
            .map(|i| C::from_polar(2., i as f32 * 0.13))
            .collect::<Vec<_>>();
        let good = (0..1000).map(|i| i % 4 != 0).collect::<Vec<_>>();
        let (a, n) = acf(&z, &good, 50);
        for l in 0..50 {
            assert!(n[l] > 0);
            assert!((a[l] - C::from_polar(4., l as f32 * 0.13)).norm() < 0.001);
        }
    }
}
