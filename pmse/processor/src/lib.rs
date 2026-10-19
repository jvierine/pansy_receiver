//! DSP only. Digital RF/Metadata I/O remains with its reference library.
//! All complex arrays are interleaved float32, matching NumPy complex64.
use rustfft::{num_complex::Complex32, FftPlanner};
use std::{panic::{catch_unwind, AssertUnwindSafe}, slice};
pub mod imaging;

fn status(f: impl FnOnce() -> Result<(), i32>) -> i32 {
    match catch_unwind(AssertUnwindSafe(f)) { Ok(Ok(())) => 0, Ok(Err(e)) => e, Err(_) => -3 }
}

fn decode(rx: &[Complex32], tx: &[Complex32], out: &mut [Complex32],
          powers: &mut [f64], ranges: &[u32], ipp: usize, txlen: usize) -> Result<(), i32> {
    if ipp <= txlen || txlen == 0 || rx.len() % ipp != 0 ||
        tx.len() != rx.len() || powers.len() != rx.len()/ipp ||
        out.len() != powers.len()*ranges.len() || ranges.iter().any(|&r| r as usize+txlen > ipp) {
        return Err(-1);
    }
    let n=(ipp+txlen-1).next_power_of_two();
    let mut plan=FftPlanner::<f32>::new();
    let fft=plan.plan_fft_forward(n); let ifft=plan.plan_fft_inverse(n);
    let scratchlen=fft.get_inplace_scratch_len().max(ifft.get_inplace_scratch_len());
    let mut scratch=vec![Complex32::default();scratchlen];
    let mut a=vec![Complex32::default();n]; let mut b=a.clone();
    for p in 0..powers.len() {
        let reference=&tx[p*ipp..p*ipp+txlen];
        let energy: f64=reference.iter().map(|z| (z.re as f64).powi(2)+(z.im as f64).powi(2)).sum();
        if !energy.is_finite() || energy <= 0. { return Err(-2); }
        powers[p]=energy/txlen as f64;
        a.fill(Complex32::default()); b.fill(Complex32::default());
        a[..ipp].copy_from_slice(&rx[p*ipp..(p+1)*ipp]);
        for (v,z) in b.iter_mut().zip(reference) { *v=*z/(energy.sqrt() as f32); }
        fft.process_with_scratch(&mut a,&mut scratch); fft.process_with_scratch(&mut b,&mut scratch);
        for (x,y) in a.iter_mut().zip(&b) { *x=*x*y.conj(); }
        ifft.process_with_scratch(&mut a,&mut scratch);
        for (j,&r) in ranges.iter().enumerate() { out[p*ranges.len()+j]=a[r as usize]/n as f32; }
    }
    Ok(())
}

/// Caller must provide nonaliasing contiguous buffers of the documented sizes.
#[no_mangle]
pub unsafe extern "C" fn pansy_decode_selected(rx: *const Complex32, tx: *const Complex32,
    pulses: usize, ipp: usize, txlen: usize, ranges: *const u32, nr: usize,
    out: *mut Complex32, power: *mut f64) -> i32 {
    status(|| {
        if [rx as usize,tx as usize,ranges as usize,out as usize,power as usize].contains(&0) ||
            pulses==0 || ipp==0 || nr==0 || pulses.checked_mul(ipp).is_none() || pulses.checked_mul(nr).is_none() { return Err(-1); }
        decode(slice::from_raw_parts(rx,pulses*ipp),slice::from_raw_parts(tx,pulses*ipp),
            slice::from_raw_parts_mut(out,pulses*nr),slice::from_raw_parts_mut(power,pulses),
            slice::from_raw_parts(ranges,nr),ipp,txlen)
    })
}

fn spectrum(input: &[Complex32], window: &[f32], out: &mut [Complex32],
            series: usize, nr: usize, nfft: usize) -> Result<(), i32> {
    let count=window.len();
    if count<4 || count>nfft || nfft%2!=0 || input.len()!=series*count*nr || out.len()!=series*nfft*nr { return Err(-1); }
    let mut plan=FftPlanner::<f32>::new(); let fft=plan.plan_fft_forward(nfft);
    let mut scratch=vec![Complex32::default();fft.get_inplace_scratch_len()];
    let mut buf=vec![Complex32::default();nfft];
    for s in 0..series { for r in 0..nr {
        buf.fill(Complex32::default());
        for p in 0..count { buf[p]=input[(s*count+p)*nr+r]*window[p]; }
        fft.process_with_scratch(&mut buf,&mut scratch);
        for k in 0..nfft { out[(s*nfft+k)*nr+r]=buf[(k+nfft/2)%nfft]; }
    }}
    Ok(())
}

fn complementary(input: &[Complex32],window: &[f32],out: &mut [Complex32],
                 series: usize,nr: usize,nfft: usize)->Result<(),i32>{
    let pairs=window.len();let count=2*pairs;
    if pairs<4 || pairs>nfft || nfft%2!=0 || input.len()!=series*count*nr || out.len()!=series*nfft*nr{return Err(-1);}
    let mut planner=FftPlanner::<f32>::new();let fft=planner.plan_fft_forward(nfft);
    let mut scratch=vec![Complex32::default();fft.get_inplace_scratch_len()];
    let mut a=vec![Complex32::default();nfft];let mut b=a.clone();
    for s in 0..series{for r in 0..nr{
        a.fill(Complex32::default());b.fill(Complex32::default());
        for k in 0..pairs{
            a[k]=input[(s*count+2*k)*nr+r]*window[k];
            b[k]=input[(s*count+2*k+1)*nr+r]*window[k];
        }
        fft.process_with_scratch(&mut a,&mut scratch);fft.process_with_scratch(&mut b,&mut scratch);
        for k in 0..nfft{let j=(k+nfft/2)%nfft;
            out[(s*nfft+k)*nr+r]=(a[j]+b[j])/2_f32.sqrt();
        }
    }}
    Ok(())
}

/// Measured-TX decoded pulse pairs; shared window and complementary sum.
#[no_mangle]
pub unsafe extern "C" fn pansy_complementary(input:*const Complex32,window:*const f32,
    series:usize,pairs:usize,nr:usize,nfft:usize,out:*mut Complex32)->i32{
    status(||{
        if [input as usize,window as usize,out as usize].contains(&0) || series==0 || nr==0 ||
            series.checked_mul(pairs).and_then(|v|v.checked_mul(2*nr)).is_none() ||
            series.checked_mul(nfft).and_then(|v|v.checked_mul(nr)).is_none(){return Err(-1);}
        complementary(slice::from_raw_parts(input,series*2*pairs*nr),slice::from_raw_parts(window,pairs),
            slice::from_raw_parts_mut(out,series*nfft*nr),series,nr,nfft)
    })
}

/// Input [series,pulse,range], output [series,fftshift_bin,range].
#[no_mangle]
pub unsafe extern "C" fn pansy_coherent(input: *const Complex32, window: *const f32,
    series: usize, count: usize, nr: usize, nfft: usize, out: *mut Complex32) -> i32 {
    status(|| {
        if [input as usize,window as usize,out as usize].contains(&0) || series==0 || nr==0 ||
            series.checked_mul(count).and_then(|v|v.checked_mul(nr)).is_none() ||
            series.checked_mul(nfft).and_then(|v|v.checked_mul(nr)).is_none() { return Err(-1); }
        spectrum(slice::from_raw_parts(input,series*count*nr),slice::from_raw_parts(window,count),
            slice::from_raw_parts_mut(out,series*nfft*nr),series,nr,nfft)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn measured_pulse_gain_phase_and_linear_lag() {
        let ipp=32; let mut tx=vec![Complex32::default();ipp];
        tx[0]=Complex32::new(0.,2.); tx[1]=Complex32::new(2.,0.);
        let mut rx=vec![Complex32::default();ipp]; rx[8]=tx[0]; rx[9]=tx[1];
        let mut out=vec![Complex32::default();2]; let mut power=vec![0.];
        decode(&rx,&tx,&mut out,&mut power,&[8,28],ipp,2).unwrap();
        assert!((out[0].re-8_f32.sqrt()).abs()<1e-5); assert!(out[0].im.abs()<1e-5);
        assert!(out[1].norm()<1e-5); assert_eq!(power[0],4.);
    }
    #[test] fn positive_doppler_is_positive_and_padding_is_zero() {
        let n=16; let x: Vec<_>=(0..n).map(|k|Complex32::from_polar(1.,2.*std::f32::consts::PI*3.*k as f32/n as f32)).collect();
        let mut out=vec![Complex32::default();n]; spectrum(&x,&vec![1.;n],&mut out,1,1,n).unwrap();
        assert!(out[11].norm()>15.99); assert!(out[5].norm()<1e-4);
    }
}
