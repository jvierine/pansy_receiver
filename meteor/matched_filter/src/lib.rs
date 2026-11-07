//! Measured-pulse range/Doppler search. No I/O or live-pipeline side effects.
use anyhow::{ensure, Result};
use rustfft::{num_complex::Complex32, Fft, FftPlanner};
use std::sync::Arc;

pub struct Filter {
    gates: Vec<usize>,
    txlen: usize,
    decimation: usize,
    fft: Arc<dyn Fft<f32>>,
    buffer: Vec<Complex32>,
    scratch: Vec<Complex32>,
}

pub struct Output {
    pub power: Vec<f32>,   // [range, shifted frequency]
    pub median_power: f32, // Legacy noise statistic, NOT calibrated per-cell noise.
    pub channels: usize,
    pub gates: Vec<usize>,
    pub nfft: usize,
    pub decimation: usize,
    spectra: Vec<Complex32>, // [receiver, range, shifted frequency]
}

#[derive(Debug)]
pub struct Peak {
    pub range_sample: usize,
    pub doppler_bin: usize,
    pub excess_over_median: f32,
    pub receiver_amplitudes: Vec<Complex32>,
}

impl Filter {
    /// decimation=1 uses every TX sample and the full native Doppler band.
    /// decimation>1 reproduces legacy boxcar sums and discarded incomplete groups.
    pub fn new(gates: Vec<usize>, txlen: usize, decimation: usize, nfft: usize) -> Result<Self> {
        ensure!(
            !gates.is_empty() && txlen > 0 && decimation > 0,
            "empty filter"
        );
        ensure!(gates.windows(2).all(|g| g[0] < g[1]), "gates must increase");
        ensure!(
            txlen / decimation > 0 && nfft >= txlen / decimation,
            "FFT too short"
        );
        ensure!(nfft.is_multiple_of(2), "use an even FFT length");
        let fft = FftPlanner::new().plan_fft_forward(nfft);
        let scratch = vec![Complex32::default(); fft.get_inplace_scratch_len()];
        Ok(Self {
            gates,
            txlen,
            decimation,
            fft,
            buffer: vec![Complex32::default(); nfft],
            scratch,
        })
    }

    /// rx is receiver-major. Apply receiver amplitude calibration before calling.
    pub fn process(
        &mut self,
        rx: &[Complex32],
        channels: usize,
        ipp: usize,
        tx: &[Complex32],
    ) -> Result<Output> {
        ensure!(
            channels > 0 && channels.checked_mul(ipp) == Some(rx.len()),
            "RX shape"
        );
        ensure!(tx.len() == self.txlen, "TX length");
        ensure!(
            self.gates
                .last()
                .unwrap()
                .checked_add(self.txlen)
                .is_some_and(|n| n <= ipp),
            "range outside pulse"
        );
        ensure!(
            rx.iter()
                .chain(tx)
                .all(|z| z.re.is_finite() && z.im.is_finite()),
            "nonfinite voltage"
        );
        ensure!(tx.iter().any(|z| z.norm_sqr() > 0.0), "zero TX pulse");
        let nfft = self.buffer.len();
        let cells = self
            .gates
            .len()
            .checked_mul(nfft)
            .ok_or_else(|| anyhow::anyhow!("shape overflow"))?;
        let total = channels
            .checked_mul(cells)
            .ok_or_else(|| anyhow::anyhow!("shape overflow"))?;
        let mut power = vec![0.0; cells];
        let mut spectra = vec![Complex32::default(); total];
        for ch in 0..channels {
            for (r, &gate) in self.gates.iter().enumerate() {
                self.buffer.fill(Complex32::default());
                for j in 0..self.txlen / self.decimation {
                    for k in 0..self.decimation {
                        let t = j * self.decimation + k;
                        self.buffer[j] += rx[ch * ipp + gate + t] * tx[t].conj();
                    }
                }
                self.fft
                    .process_with_scratch(&mut self.buffer, &mut self.scratch);
                for k in 0..nfft {
                    let z = self.buffer[(k + nfft / 2) % nfft];
                    spectra[ch * cells + r * nfft + k] = z;
                    power[r * nfft + k] += z.norm_sqr();
                }
            }
        }
        ensure!(power.iter().all(|x| x.is_finite()), "power overflow");
        let mut ordered = power.clone();
        let mid = ordered.len() / 2;
        ordered.select_nth_unstable_by(mid, f32::total_cmp);
        let median_power = if ordered.len() % 2 == 0 {
            0.5 * ordered[mid]
                + 0.5
                    * ordered[..mid]
                        .iter()
                        .copied()
                        .max_by(f32::total_cmp)
                        .unwrap()
        } else {
            ordered[mid]
        };
        Ok(Output {
            power,
            spectra,
            median_power,
            channels,
            gates: self.gates.clone(),
            nfft,
            decimation: self.decimation,
        })
    }
}

impl Output {
    pub fn frequency_hz(&self, bin: usize, sample_rate: f64) -> f64 {
        (bin as f64 - (self.nfft / 2) as f64) * sample_rate / (self.nfft * self.decimation) as f64
    }

    /// Multiple 2D maxima, with configurable resolution-cell suppression.
    /// Guards are range *samples* and FFT bins; Doppler wraps at Nyquist.
    /// This is a candidate extractor, not a calibrated false-alarm detector.
    pub fn peaks(
        &self,
        min_excess: f32,
        max_peaks: usize,
        range_guard: usize,
        doppler_guard: usize,
    ) -> Result<Vec<Peak>> {
        ensure!(
            min_excess.is_finite() && min_excess >= 0.0,
            "invalid threshold"
        );
        ensure!(self.median_power > 0.0, "no positive noise reference");
        let mut candidates = Vec::new();
        for r in 0..self.gates.len() {
            for k in 0..self.nfft {
                let idx = r * self.nfft + k;
                if self.power[idx] / self.median_power - 1.0 < min_excess {
                    continue;
                }
                let local = (r.saturating_sub(1)..=(r + 1).min(self.gates.len() - 1)).all(|rr| {
                    [(k + self.nfft - 1) % self.nfft, k, (k + 1) % self.nfft]
                        .into_iter()
                        .all(|kk| {
                            let other = rr * self.nfft + kk;
                            self.power[other] < self.power[idx]
                                || (self.power[other] == self.power[idx] && other >= idx)
                        })
                });
                if local {
                    candidates.push(idx);
                }
            }
        }
        candidates
            .sort_unstable_by(|&a, &b| self.power[b].total_cmp(&self.power[a]).then(a.cmp(&b)));
        let mut peaks: Vec<Peak> = Vec::new();
        for idx in candidates {
            if peaks.len() == max_peaks {
                break;
            }
            let gate = self.gates[idx / self.nfft];
            let bin = idx % self.nfft;
            if peaks.iter().any(|p| {
                let dk = p.doppler_bin.abs_diff(bin);
                p.range_sample.abs_diff(gate) <= range_guard
                    && dk.min(self.nfft - dk) <= doppler_guard
            }) {
                continue;
            }
            peaks.push(Peak {
                range_sample: gate,
                doppler_bin: bin,
                excess_over_median: self.power[idx] / self.median_power - 1.0,
                receiver_amplitudes: (0..self.channels)
                    .map(|ch| self.spectra[ch * self.power.len() + idx])
                    .collect(),
            });
        }
        Ok(peaks)
    }
}
