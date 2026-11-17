//! Strict little-endian Tibet reader. Integer pulse/sample indices preserve timing.
use anyhow::{bail, ensure, Context, Result};
pub mod voltage;
use rustfft::num_complex::Complex32;
use std::{
    fs::File,
    io::{Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

pub const HEADER: usize = 8192;
pub const MAGIC: u32 = 0xaabbccdd;
pub fn u16at(b: &[u8], p: usize) -> u16 {
    u16::from_le_bytes(b[p..p + 2].try_into().unwrap())
}
pub fn u32at(b: &[u8], p: usize) -> u32 {
    u32::from_le_bytes(b[p..p + 4].try_into().unwrap())
}
pub fn u64at(b: &[u8], p: usize) -> u64 {
    u64::from_le_bytes(b[p..p + 8].try_into().unwrap())
}

#[derive(Debug, Clone)]
pub struct Header {
    pub raw: Vec<u8>,
    pub path: PathBuf,
    pub channels: usize,
    pub gates: usize,
    pub pulses: usize,
    pub ipp_us: u32,
    pub fs: u32,
    pub frequency_hz: f64,
    pub waveform: u16,
    pub width_us: f64,
    pub inacc: u16,
    pub iq_start: u32,
    pub alpha: Vec<u16>,
    pub beta: Vec<u16>,
    pub phase_index: Vec<u16>,
    pub amplitude_index: Vec<u16>,
    pub tx_alpha: u16,
    pub tx_beta: u16,
    pub frame: u64,
    pub epoch_us: u64,
    pub padding_bytes: u64,
}
impl Header {
    pub fn parse(path: &Path, b: Vec<u8>, file_size: u64) -> Result<Self> {
        ensure!(b.len() == HEADER, "short header");
        ensure!(
            u32at(&b, 0) == 0xaaaa5555 && u32at(&b, 8188) == 0xaa5555aa,
            "header guards"
        );
        let fs = match b[94] {
            0 => 8_000_000,
            1 => 4_000_000,
            2 => 2_000_000,
            3 => 1_000_000,
            4 => 500_000,
            _ => bail!("unknown sample-rate code"),
        };
        let channels = u16at(&b, 426) as usize;
        let gates = u32at(&b, 240) as usize;
        let pulses = u32at(&b, 80) as usize;
        ensure!(
            (1..=20).contains(&channels) && gates > 0 && gates < 100_000 && pulses > 0,
            "invalid dimensions"
        );
        let required = HEADER as u64 + pulses as u64 * (gates as u64 * channels as u64 * 8 + 16);
        ensure!(
            file_size >= required,
            "file shorter than declared pulse records"
        );
        let name = path
            .file_name()
            .context("filename")?
            .to_str()
            .context("filename encoding")?;
        let date = chrono::NaiveDateTime::parse_from_str(
            name.get(..14).context("missing filename timestamp")?,
            "%Y%m%d%H%M%S",
        )?;
        // The instrument clock is a placeholder in this campaign. This is a stated
        // filename-Beijing hypothesis, not a claim of GPS-quality synchronization.
        let epoch_us = u64::try_from(date.and_utc().timestamp_micros() - 8 * 3600 * 1_000_000)?;
        let array = |offset| {
            (0..channels)
                .map(|i| u16at(&b, offset + 2 * i))
                .collect::<Vec<_>>()
        };
        Ok(Self {
            channels,
            gates,
            pulses,
            ipp_us: u32at(&b, 76),
            fs,
            frequency_hz: u32at(&b, 84) as f64 / (1u64 << 32) as f64 * 120e6,
            waveform: u16at(&b, 88),
            width_us: u32at(&b, 90) as f64 / 10.,
            inacc: u16at(&b, 730),
            iq_start: u32at(&b, 236),
            alpha: array(428),
            beta: array(468),
            phase_index: array(508),
            amplitude_index: array(548),
            tx_alpha: u16at(&b, 256),
            tx_beta: u16at(&b, 258),
            frame: u64at(&b, 8),
            epoch_us,
            padding_bytes: file_size - required,
            raw: b,
            path: path.to_owned(),
        })
    }
    pub fn record_bytes(&self) -> usize {
        self.channels * self.gates * 8 + 16
    }
    /// Inference consistent with supplied S,N,Z,E,W labels; hardware sign and
    /// array-axis conventions still need confirmation for absolute astrometry.
    pub fn inferred_rx_uv(&self) -> Vec<f64> {
        self.alpha
            .iter()
            .zip(&self.beta)
            .flat_map(|(&a, &b)| {
                let signed = |v: u16| {
                    if v > 32767 {
                        v as i32 - 65536
                    } else {
                        v as i32
                    }
                };
                [
                    -signed(a) as f64 / 65536. * 299792458. / self.frequency_hz / 2.1,
                    -signed(b) as f64 / 65536. * 299792458. / self.frequency_hz / 3.637,
                ]
            })
            .collect()
    }
    pub fn signature(&self) -> Vec<u16> {
        let mut s = vec![
            self.tx_alpha,
            self.tx_beta,
            self.channels as u16,
            self.waveform,
            self.inacc,
        ];
        s.extend(&self.alpha);
        s.extend(&self.beta);
        s.extend(&self.phase_index);
        s.extend(&self.amplitude_index);
        s
    }
    pub fn pulse_time_us(&self, index: usize) -> u64 {
        self.epoch_us + index as u64 * self.ipp_us as u64
    }
}
pub struct Reader {
    pub header: Header,
    file: File,
    buffer: Vec<u8>,
}
impl Reader {
    pub fn open(path: &Path) -> Result<Self> {
        let mut file = File::open(path)?;
        let size = file.metadata()?.len();
        let mut b = vec![0; HEADER];
        file.read_exact(&mut b)?;
        let header = Header::parse(path, b, size)?;
        let buffer = vec![0; header.record_bytes()];
        Ok(Self {
            header,
            file,
            buffer,
        })
    }
    pub fn pulse_i32(&mut self, index: usize) -> Result<Vec<i32>> {
        ensure!(
            index < self.header.pulses,
            "pulse index outside declared records"
        );
        let offset = HEADER as u64 + index as u64 * self.header.record_bytes() as u64;
        self.file.seek(SeekFrom::Start(offset))?;
        self.file.read_exact(&mut self.buffer)?;
        let end = self.buffer.len() - 8;
        ensure!(
            u32at(&self.buffer, 0) == MAGIC && u32at(&self.buffer, end) == MAGIC,
            "bad pulse marker at {index}"
        );
        ensure!(
            u32at(&self.buffer, 4) == (index + 1) as u32
                && u32at(&self.buffer, end + 4) == (index + 1) as u32,
            "pulse counter mismatch at {index}"
        );
        let mut z = vec![0; 2 * self.header.channels * self.header.gates];
        for gate in 0..self.header.gates {
            for channel in 0..self.header.channels {
                let p = 8 + 8 * (gate * self.header.channels + channel);
                let re = i32::from_le_bytes(self.buffer[p..p + 4].try_into().unwrap());
                let im = i32::from_le_bytes(self.buffer[p + 4..p + 8].try_into().unwrap());
                let j = 2 * (channel * self.header.gates + gate);
                z[j] = re;
                z[j + 1] = im;
            }
        }
        Ok(z)
    }
    pub fn pulse(&mut self, index: usize) -> Result<Vec<Complex32>> {
        Ok(self
            .pulse_i32(index)?
            .as_chunks::<2>()
            .0
            .iter()
            .map(|v| Complex32::new(v[0] as f32, v[1] as f32))
            .collect())
    }
}

pub fn raw_paths(root: &Path) -> Result<Vec<PathBuf>> {
    if root.is_file() {
        ensure!(
            root.extension().is_some_and(|s| s == "raw"),
            "expected a raw voltage file"
        );
        return Ok(vec![root.to_owned()]);
    }
    let mut p = std::fs::read_dir(root)?
        .filter_map(|x| x.ok())
        .map(|x| x.path())
        .filter(|p| {
            p.extension().is_some_and(|x| x == "raw")
                && p.file_name().unwrap().to_string_lossy().starts_with("20")
        })
        .collect::<Vec<_>>();
    p.sort();
    ensure!(!p.is_empty(), "no radar raw files");
    Ok(p)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn header() -> Vec<u8> {
        let mut b = vec![0; HEADER];
        for (p, v) in [
            (0, 0xaaaa5555u32),
            (8188, 0xaa5555aa),
            (76, 1280),
            (80, 2),
            (240, 640),
        ] {
            b[p..p + 4].copy_from_slice(&v.to_le_bytes());
        }
        b[426..428].copy_from_slice(&15u16.to_le_bytes());
        b[94] = 4;
        b
    }
    #[test]
    fn timing_and_padding_are_explicit() {
        let path = Path::new("20261114055156_test.raw");
        let n = HEADER as u64 + 2 * (640 * 15 * 8 + 16) + 100;
        let h = Header::parse(path, header(), n).unwrap();
        assert_eq!(h.padding_bytes, 100);
        assert_eq!(h.pulse_time_us(1) - h.pulse_time_us(0), 1280);
        assert_eq!(h.gates as u32 * 1_000_000 / h.fs, h.ipp_us);
        assert!(Header::parse(path, header(), HEADER as u64).is_err());
    }
    #[test]
    fn iq_order_and_counter_validation() {
        let path = std::env::temp_dir().join(format!(
            "20261114055156_tibet_test_{}.raw",
            std::process::id()
        ));
        let mut b = header();
        b[80..84].copy_from_slice(&1u32.to_le_bytes());
        let mut record = vec![0u8; 640 * 15 * 8 + 16];
        for p in [0, record.len() - 8] {
            record[p..p + 4].copy_from_slice(&MAGIC.to_le_bytes());
            record[p + 4..p + 8].copy_from_slice(&1u32.to_le_bytes());
        }
        let slot = 8 + 8 * (7 * 15 + 4);
        record[slot..slot + 4].copy_from_slice(&(-31i32).to_le_bytes());
        record[slot + 4..slot + 8].copy_from_slice(&72i32.to_le_bytes());
        b.extend(record);
        std::fs::write(&path, b).unwrap();
        let mut reader = Reader::open(&path).unwrap();
        let z = reader.pulse(0).unwrap();
        assert_eq!(z[4 * 640 + 7], Complex32::new(-31., 72.));
        assert!(reader.pulse(1).is_err());
        std::fs::remove_file(path).unwrap();
    }
}
