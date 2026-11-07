use crate::Filter;
use rustfft::num_complex::Complex32;
use std::{
    panic::{catch_unwind, AssertUnwindSafe},
    ptr, slice,
};

/// # Safety
/// gates must point to n_gates readable usize values.
#[no_mangle]
pub unsafe extern "C" fn pansy_mf_new(
    gates: *const usize,
    n_gates: usize,
    txlen: usize,
    decimation: usize,
    nfft: usize,
) -> *mut Filter {
    if gates.is_null() || n_gates == 0 {
        return ptr::null_mut();
    }
    catch_unwind(|| {
        Filter::new(
            slice::from_raw_parts(gates, n_gates).to_vec(),
            txlen,
            decimation,
            nfft,
        )
        .map(|f| Box::into_raw(Box::new(f)))
        .unwrap_or(ptr::null_mut())
    })
    .unwrap_or(ptr::null_mut())
}

/// # Safety
/// handle must be a live, exclusively accessed pointer from pansy_mf_new.
/// Input/output buffers must have their declared lengths; median holds one f32.
#[no_mangle]
pub unsafe extern "C" fn pansy_mf_process(
    handle: *mut Filter,
    rx: *const Complex32,
    channels: usize,
    ipp: usize,
    tx: *const Complex32,
    txlen: usize,
    power: *mut f32,
    power_len: usize,
    median: *mut f32,
) -> i32 {
    if handle.is_null() || rx.is_null() || tx.is_null() || power.is_null() || median.is_null() {
        return -1;
    }
    let Some(nrx) = channels.checked_mul(ipp) else {
        return -1;
    };
    catch_unwind(AssertUnwindSafe(|| {
        let f = &mut *handle;
        if txlen != f.txlen || f.gates.len().checked_mul(f.buffer.len()) != Some(power_len) {
            return -1;
        }
        match f.process(
            slice::from_raw_parts(rx, nrx),
            channels,
            ipp,
            slice::from_raw_parts(tx, txlen),
        ) {
            Ok(out) => {
                ptr::copy_nonoverlapping(out.power.as_ptr(), power, power_len);
                *median = out.median_power;
                0
            }
            Err(_) => -2,
        }
    }))
    .unwrap_or(-3)
}

/// # Safety
/// handle must be null or an unfreed pointer from pansy_mf_new, not in use.
#[no_mangle]
pub unsafe extern "C" fn pansy_mf_free(handle: *mut Filter) {
    if !handle.is_null() {
        drop(Box::from_raw(handle));
    }
}

/// Extract bounded diagnostic peaks without changing the legacy process ABI.
/// # Safety
/// handle is live and exclusive; rx/tx have the declared lengths. Each scalar
/// output has capacity entries, amplitudes has capacity*channels Complex32
/// entries, and count/median each point to one writable value.
#[no_mangle]
pub unsafe extern "C" fn pansy_mf_peaks(
    handle: *mut Filter,
    rx: *const Complex32,
    channels: usize,
    ipp: usize,
    tx: *const Complex32,
    txlen: usize,
    threshold: f32,
    capacity: usize,
    range_guard: usize,
    doppler_guard: usize,
    gates: *mut usize,
    bins: *mut usize,
    excess: *mut f32,
    amplitudes: *mut Complex32,
    count: *mut usize,
    median: *mut f32,
) -> i32 {
    if handle.is_null()
        || rx.is_null()
        || tx.is_null()
        || gates.is_null()
        || bins.is_null()
        || excess.is_null()
        || amplitudes.is_null()
        || count.is_null()
        || median.is_null()
        || capacity == 0
        || capacity > 256
        || channels.checked_mul(capacity).is_none()
    {
        return -1;
    }
    *count = 0;
    let Some(nrx) = channels.checked_mul(ipp) else {
        return -1;
    };
    catch_unwind(AssertUnwindSafe(|| {
        let f = &mut *handle;
        if txlen != f.txlen {
            return -1;
        }
        let Ok(out) = f.process(
            slice::from_raw_parts(rx, nrx),
            channels,
            ipp,
            slice::from_raw_parts(tx, txlen),
        ) else {
            return -2;
        };
        let Ok(peaks) = out.peaks(threshold, capacity, range_guard, doppler_guard) else {
            return -2;
        };
        for (i, peak) in peaks.iter().enumerate() {
            *gates.add(i) = peak.range_sample;
            *bins.add(i) = peak.doppler_bin;
            *excess.add(i) = peak.excess_over_median;
            ptr::copy_nonoverlapping(
                peak.receiver_amplitudes.as_ptr(),
                amplitudes.add(i * channels),
                channels,
            );
        }
        *count = peaks.len();
        *median = out.median_power;
        0
    }))
    .unwrap_or(-3)
}
