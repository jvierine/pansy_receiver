#!/usr/bin/env python3
"""Build and query PANSY TX cross-phase quality products."""

from __future__ import annotations

import argparse
import datetime as dt
import warnings
from pathlib import Path

import h5py
import numpy as np


DEFAULT_THRESHOLD_DEG = 30.0
DEFAULT_MAX_NEAREST_AGE_S = 7200.0
DEFAULT_MIN_VALID_CHANNELS = 6
DEFAULT_FALLBACK_SEARCH_RADIUS_S = 2.0
DEFAULT_FALLBACK_TX_SAMPLES = 120


def _sample_to_datetime64(sample: int) -> np.datetime64:
    return np.datetime64(int(sample), "us")


def _parse_time_us(value: str | None) -> int | None:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    if text.isdigit():
        return int(text)
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    parsed = dt.datetime.fromisoformat(text)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return int(parsed.timestamp() * 1e6)


def circular_median(phases_rad: np.ndarray, axis: int = 0) -> np.ndarray:
    return np.angle(np.nanmedian(np.exp(1j * phases_rad), axis=axis))


def angular_delta_rad(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.angle(np.exp(1j * (a - b)))


def read_phase_file(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    samples: list[int] = []
    phases: list[np.ndarray] = []
    amplitudes: list[np.ndarray] = []
    try:
        with h5py.File(path, "r") as h5:
            for key in h5.keys():
                if not key.isdigit():
                    continue
                group = h5[key]
                if "xphase" not in group:
                    continue
                xphase = np.asarray(group["xphase"])
                if xphase.size == 0:
                    continue
                phase = np.full(8, np.nan, dtype=np.float32)
                amp = np.full(8, np.nan, dtype=np.float32)
                n = min(8, xphase.size)
                phase[:n] = np.angle(xphase[:n]).astype(np.float32)
                amp[:n] = np.abs(xphase[:n]).astype(np.float32)
                samples.append(int(key))
                phases.append(phase)
                amplitudes.append(amp)
    except OSError:
        return (
            np.empty(0, dtype=np.int64),
            np.empty((0, 8), dtype=np.float32),
            np.empty((0, 8), dtype=np.float32),
        )
    if not samples:
        return (
            np.empty(0, dtype=np.int64),
            np.empty((0, 8), dtype=np.float32),
            np.empty((0, 8), dtype=np.float32),
        )
    return (
        np.asarray(samples, dtype=np.int64),
        np.vstack(phases).astype(np.float32),
        np.vstack(amplitudes).astype(np.float32),
    )


def _directory_day_us(path: Path) -> int | None:
    try:
        day = path.parent.name[:10]
        return int(np.datetime64(day, "D").astype("datetime64[us]").astype(np.int64))
    except Exception:
        return None


def load_phase_history(phase_dir: Path, min_sample_idx: int | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    samples = []
    phases = []
    amplitudes = []
    min_day_us = None
    if min_sample_idx is not None:
        min_day_us = int(np.datetime64(int(min_sample_idx), "us").astype("datetime64[D]").astype("datetime64[us]").astype(np.int64))
    for path in sorted(phase_dir.glob("20*/txphase@*.h5")):
        if min_day_us is not None:
            day_us = _directory_day_us(path)
            if day_us is not None and day_us < min_day_us:
                continue
        s, p, a = read_phase_file(path)
        if min_sample_idx is not None and s.size:
            keep = s >= int(min_sample_idx)
            s = s[keep]
            p = p[keep]
            a = a[keep]
        if s.size:
            samples.append(s)
            phases.append(p)
            amplitudes.append(a)
    if not samples:
        return (
            np.empty(0, dtype=np.int64),
            np.empty((0, 8), dtype=np.float32),
            np.empty((0, 8), dtype=np.float32),
        )
    sample = np.concatenate(samples)
    phase = np.concatenate(phases, axis=0)
    amplitude = np.concatenate(amplitudes, axis=0)
    order = np.argsort(sample)
    return sample[order], phase[order], amplitude[order]


def load_phase_cache(cache_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    samples = []
    phases = []
    for path in sorted(cache_dir.glob("*.npz")):
        try:
            data = np.load(path)
            t = np.asarray(data["time"], dtype=np.int64)
            p = np.asarray(data["phase"], dtype=np.float32)
        except Exception:
            continue
        if t.size == 0:
            continue
        if p.ndim != 2:
            continue
        if p.shape[1] < 8:
            padded = np.full((p.shape[0], 8), np.nan, dtype=np.float32)
            padded[:, : p.shape[1]] = p
            p = padded
        elif p.shape[1] > 8:
            p = p[:, :8]
        samples.append(t)
        phases.append(p.astype(np.float32))
    if not samples:
        return (
            np.empty(0, dtype=np.int64),
            np.empty((0, 8), dtype=np.float32),
            np.empty((0, 8), dtype=np.float32),
        )
    sample = np.concatenate(samples)
    phase = np.concatenate(phases, axis=0)
    amplitude = np.where(np.isfinite(phase), 1.0, np.nan).astype(np.float32)
    order = np.argsort(sample)
    return sample[order], phase[order], amplitude[order]


def reference_phase(
    sample_idx: np.ndarray,
    phase_rad: np.ndarray,
    reference_start_us: int | None = None,
    reference_end_us: int | None = None,
    reference_latest_days: float = 30.0,
) -> tuple[np.ndarray, np.ndarray]:
    if sample_idx.size == 0:
        raise ValueError("no TX phase samples available")
    if reference_end_us is None:
        reference_end_us = int(np.nanmax(sample_idx))
    if reference_start_us is None:
        reference_start_us = int(reference_end_us - reference_latest_days * 86400.0 * 1e6)
    mask = (sample_idx >= reference_start_us) & (sample_idx <= reference_end_us)
    if np.count_nonzero(mask) < 2:
        raise ValueError("not enough TX phase samples in the reference interval")
    return circular_median(phase_rad[mask], axis=0).astype(np.float32), mask


def build_quality_table(
    sample_idx: np.ndarray,
    phase_rad: np.ndarray,
    amplitude: np.ndarray,
    reference_rad: np.ndarray,
    threshold_deg: float = DEFAULT_THRESHOLD_DEG,
    min_valid_channels: int = DEFAULT_MIN_VALID_CHANNELS,
) -> dict[str, np.ndarray]:
    valid = np.isfinite(phase_rad) & np.isfinite(reference_rad[None, :]) & np.isfinite(amplitude)
    diff_rad = angular_delta_rad(phase_rad, reference_rad[None, :]).astype(np.float32)
    abs_diff_deg = np.abs(np.rad2deg(diff_rad))
    valid_count = np.sum(valid, axis=1).astype(np.int16)
    masked = np.where(valid, abs_diff_deg, np.nan)
    max_abs_diff_deg = np.nanmax(masked, axis=1).astype(np.float32)
    rms_diff_deg = np.sqrt(np.nanmean(masked**2, axis=1)).astype(np.float32)
    max_abs_diff_deg[valid_count == 0] = np.nan
    rms_diff_deg[valid_count == 0] = np.nan
    good = (valid_count >= int(min_valid_channels)) & (max_abs_diff_deg <= float(threshold_deg))
    return {
        "sample_idx": sample_idx.astype(np.int64),
        "phase_rad": phase_rad.astype(np.float32),
        "amplitude": amplitude.astype(np.float32),
        "phase_difference_rad": diff_rad.astype(np.float32),
        "valid_channel_count": valid_count,
        "max_abs_diff_deg": max_abs_diff_deg,
        "rms_diff_deg": rms_diff_deg,
        "good": good.astype(np.bool_),
    }


def classify_phase_vector(
    phase_rad: np.ndarray,
    amplitude: np.ndarray,
    reference_rad: np.ndarray,
    threshold_deg: float = DEFAULT_THRESHOLD_DEG,
    min_valid_channels: int = DEFAULT_MIN_VALID_CHANNELS,
) -> dict[str, float | int | bool | str]:
    phase_rad = np.asarray(phase_rad, dtype=np.float32)
    amplitude = np.asarray(amplitude, dtype=np.float32)
    reference_rad = np.asarray(reference_rad, dtype=np.float32)
    valid = np.isfinite(phase_rad) & np.isfinite(reference_rad) & np.isfinite(amplitude)
    diff_rad = angular_delta_rad(phase_rad, reference_rad).astype(np.float32)
    abs_diff_deg = np.abs(np.rad2deg(diff_rad))
    valid_count = int(np.count_nonzero(valid))
    if valid_count:
        max_abs_diff_deg = float(np.nanmax(np.where(valid, abs_diff_deg, np.nan)))
        rms_diff_deg = float(np.sqrt(np.nanmean(np.where(valid, abs_diff_deg, np.nan) ** 2)))
    else:
        max_abs_diff_deg = np.nan
        rms_diff_deg = np.nan
    good = bool(valid_count >= int(min_valid_channels) and max_abs_diff_deg <= float(threshold_deg))
    return {
        "available": True,
        "good": good,
        "raw_good": good,
        "reason": "ok" if good else "computed_tx_phase_mismatch",
        "max_abs_diff_deg": max_abs_diff_deg,
        "rms_diff_deg": rms_diff_deg,
        "valid_channel_count": valid_count,
        "threshold_deg": float(threshold_deg),
    }


def compute_event_xphase(
    sample_idx: int,
    raw_voltage_dir: Path | str,
    tx_metadata_dir: Path | str,
    search_radius_s: float = DEFAULT_FALLBACK_SEARCH_RADIUS_S,
    tx_samples: int = DEFAULT_FALLBACK_TX_SAMPLES,
) -> dict[str, object]:
    import digital_rf as drf

    raw_voltage_dir = str(raw_voltage_dir)
    tx_metadata_dir = str(tx_metadata_dir)
    reader = drf.DigitalRFReader(raw_voltage_dir)
    tx_reader = drf.DigitalMetadataReader(tx_metadata_dir)
    radius = int(round(float(search_radius_s) * 1e6))
    records = tx_reader.read(int(sample_idx) - radius, int(sample_idx) + radius, "id")
    pulse_keys = []
    for key, value in records.items():
        try:
            is_tx = int(np.asarray(value).ravel()[0]) == 1
        except Exception:
            is_tx = False
        if is_tx:
            pulse_keys.append(int(key))
    if not pulse_keys:
        return {
            "available": False,
            "reason": "no_tx_pulse_metadata_near_event",
        }
    tx_key = min(pulse_keys, key=lambda key: abs(key - int(sample_idx)))
    channels = [f"ch{j:03d}" for j in range(8)]
    xphase = np.full(8, np.nan + 1j * np.nan, dtype=np.complex64)
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="The read_vector_c81d method is deprecated.*")
            z0 = reader.read_vector_c81d(tx_key, int(tx_samples), "ch000")
            for j, channel in enumerate(channels):
                z1 = reader.read_vector_c81d(tx_key, int(tx_samples), channel)
                xphase[j] = np.mean(z0 * np.conj(z1))
    except Exception as exc:
        return {
            "available": False,
            "reason": f"tx_phase_raw_read_failed:{type(exc).__name__}",
            "nearest_sample_idx": int(tx_key),
            "age_s": float(abs(int(sample_idx) - tx_key) / 1e6),
        }
    return {
        "available": True,
        "reason": "computed_from_raw_voltage",
        "nearest_sample_idx": int(tx_key),
        "age_s": float(abs(int(sample_idx) - tx_key) / 1e6),
        "xphase": xphase,
        "phase_rad": np.angle(xphase).astype(np.float32),
        "amplitude": np.abs(xphase).astype(np.float32),
        "mean_amplitude": float(np.nanmean(np.abs(xphase))),
    }


def write_quality_h5(
    phase_dir: Path,
    output_h5: Path,
    threshold_deg: float = DEFAULT_THRESHOLD_DEG,
    max_nearest_age_s: float = DEFAULT_MAX_NEAREST_AGE_S,
    min_valid_channels: int = DEFAULT_MIN_VALID_CHANNELS,
    reference_start_us: int | None = None,
    reference_end_us: int | None = None,
    reference_latest_days: float = 30.0,
    phase_cache_dir: Path | None = None,
) -> dict[str, float | int | str]:
    if phase_cache_dir is not None:
        sample_idx, phase_rad, amplitude = load_phase_cache(phase_cache_dir)
        cache_end = int(np.nanmax(sample_idx)) if sample_idx.size else None
        if cache_end is not None:
            recent_sample, recent_phase, recent_amplitude = load_phase_history(phase_dir, min_sample_idx=cache_end + 1)
            if recent_sample.size:
                sample_idx = np.concatenate([sample_idx, recent_sample])
                phase_rad = np.concatenate([phase_rad, recent_phase], axis=0)
                amplitude = np.concatenate([amplitude, recent_amplitude], axis=0)
                order = np.argsort(sample_idx)
                sample_idx = sample_idx[order]
                phase_rad = phase_rad[order]
                amplitude = amplitude[order]
        phase_source = str(phase_cache_dir)
        phase_source_kind = "phase_history_npz_cache_plus_recent_txphase_hdf5"
    else:
        sample_idx, phase_rad, amplitude = load_phase_history(phase_dir)
        phase_source = str(phase_dir)
        phase_source_kind = "txphase_hdf5_metadata"
    reference_rad, reference_mask = reference_phase(
        sample_idx,
        phase_rad,
        reference_start_us=reference_start_us,
        reference_end_us=reference_end_us,
        reference_latest_days=reference_latest_days,
    )
    table = build_quality_table(
        sample_idx,
        phase_rad,
        amplitude,
        reference_rad,
        threshold_deg=threshold_deg,
        min_valid_channels=min_valid_channels,
    )
    output_h5.parent.mkdir(parents=True, exist_ok=True)
    tmp = output_h5.with_suffix(output_h5.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()
    with h5py.File(tmp, "w") as h5:
        h5.attrs["schema_version"] = "pansy_tx_phase_quality_v1"
        h5.attrs["source_program"] = "tx_phase_quality.py"
        h5.attrs["phase_dir"] = str(phase_dir)
        h5.attrs["phase_source"] = phase_source
        h5.attrs["phase_source_kind"] = phase_source_kind
        h5.attrs["threshold_deg"] = float(threshold_deg)
        h5.attrs["max_nearest_age_s"] = float(max_nearest_age_s)
        h5.attrs["min_valid_channels"] = int(min_valid_channels)
        h5.attrs["reference_latest_days"] = float(reference_latest_days)
        ref_samples = sample_idx[reference_mask]
        h5.attrs["reference_start_sample_idx"] = int(ref_samples[0])
        h5.attrs["reference_end_sample_idx"] = int(ref_samples[-1])
        h5.attrs["reference_sample_count"] = int(ref_samples.size)
        h5.attrs["reference_start_utc"] = str(_sample_to_datetime64(int(ref_samples[0])))
        h5.attrs["reference_end_utc"] = str(_sample_to_datetime64(int(ref_samples[-1])))
        h5.create_dataset("reference_phase_rad", data=reference_rad)
        for key, value in table.items():
            h5.create_dataset(key, data=value)
    tmp.replace(output_h5)
    return {
        "phase_samples": int(sample_idx.size),
        "reference_samples": int(np.count_nonzero(reference_mask)),
        "good_samples": int(np.count_nonzero(table["good"])),
        "bad_samples": int(np.count_nonzero(~table["good"])),
        "threshold_deg": float(threshold_deg),
    }


def quality_for_sample(
    quality_h5: Path,
    sample_idx: int,
    max_age_s: float | None = None,
    compute_if_missing: bool = False,
    raw_voltage_dir: Path | str | None = None,
    tx_metadata_dir: Path | str | None = None,
    fallback_search_radius_s: float = DEFAULT_FALLBACK_SEARCH_RADIUS_S,
    fallback_tx_samples: int = DEFAULT_FALLBACK_TX_SAMPLES,
) -> dict[str, float | int | bool | str]:
    with h5py.File(quality_h5, "r") as h5:
        samples = np.asarray(h5["sample_idx"], dtype=np.int64)
        reference = np.asarray(h5["reference_phase_rad"], dtype=np.float32)
        threshold_deg = float(h5.attrs.get("threshold_deg", DEFAULT_THRESHOLD_DEG))
        min_valid_channels = int(h5.attrs.get("min_valid_channels", DEFAULT_MIN_VALID_CHANNELS))
        if samples.size == 0:
            result = {"available": False, "good": False, "reason": "no_tx_phase_quality_samples"}
            nearest_i = None
            within_age = False
            raw_good = False
        else:
            i = int(np.searchsorted(samples, int(sample_idx)))
            if i >= samples.size:
                nearest_i = samples.size - 1
            elif i == 0:
                nearest_i = 0
            else:
                left = i - 1
                right = i
                nearest_i = left if abs(int(sample_idx) - int(samples[left])) <= abs(int(samples[right]) - int(sample_idx)) else right
            age_s = abs(int(sample_idx) - int(samples[nearest_i])) / 1e6
            limit_s = float(max_age_s if max_age_s is not None else h5.attrs.get("max_nearest_age_s", DEFAULT_MAX_NEAREST_AGE_S))
            within_age = age_s <= limit_s
            raw_good = bool(np.asarray(h5["good"])[nearest_i])
            result = {
                "available": True,
                "good": bool(within_age and raw_good),
                "raw_good": raw_good,
                "reason": "ok" if within_age and raw_good else ("nearest_tx_phase_sample_too_far" if not within_age else "tx_phase_mismatch"),
                "nearest_sample_idx": int(samples[nearest_i]),
                "age_s": float(age_s),
                "max_age_s": float(limit_s),
                "max_abs_diff_deg": float(np.asarray(h5["max_abs_diff_deg"])[nearest_i]),
                "rms_diff_deg": float(np.asarray(h5["rms_diff_deg"])[nearest_i]),
                "valid_channel_count": int(np.asarray(h5["valid_channel_count"])[nearest_i]),
                "threshold_deg": threshold_deg,
                "reference_start_sample_idx": int(h5.attrs.get("reference_start_sample_idx", -1)),
                "reference_end_sample_idx": int(h5.attrs.get("reference_end_sample_idx", -1)),
                "reference_sample_count": int(h5.attrs.get("reference_sample_count", 0)),
            }
        if result.get("good", False) or not compute_if_missing:
            return result
        if result.get("reason") == "tx_phase_mismatch":
            return result
        if raw_voltage_dir is None or tx_metadata_dir is None:
            return result
        computed = compute_event_xphase(
            sample_idx,
            raw_voltage_dir=raw_voltage_dir,
            tx_metadata_dir=tx_metadata_dir,
            search_radius_s=fallback_search_radius_s,
            tx_samples=fallback_tx_samples,
        )
        if not computed.get("available", False):
            fallback = dict(result)
            fallback["computed_available"] = False
            fallback["computed_reason"] = str(computed.get("reason", "computed_tx_phase_unavailable"))
            return fallback
        classified = classify_phase_vector(
            computed["phase_rad"],
            computed["amplitude"],
            reference,
            threshold_deg=threshold_deg,
            min_valid_channels=min_valid_channels,
        )
        classified.update(
            {
                "reason": "ok" if classified["good"] else "computed_tx_phase_mismatch",
                "computed_available": True,
                "computed_from_raw_voltage": True,
                "nearest_sample_idx": int(computed["nearest_sample_idx"]),
                "age_s": float(computed["age_s"]),
                "mean_amplitude": float(computed["mean_amplitude"]),
                "max_age_s": float(fallback_search_radius_s),
                "reference_start_sample_idx": int(h5.attrs.get("reference_start_sample_idx", -1)),
                "reference_end_sample_idx": int(h5.attrs.get("reference_end_sample_idx", -1)),
                "reference_sample_count": int(h5.attrs.get("reference_sample_count", 0)),
            }
        )
        return classified


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase-dir", type=Path, required=True)
    parser.add_argument("--phase-cache-dir", type=Path, help="Daily NPZ cache directory from backup/scripts/pansy_phase_history.py.")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threshold-deg", type=float, default=DEFAULT_THRESHOLD_DEG)
    parser.add_argument("--max-nearest-age-s", type=float, default=DEFAULT_MAX_NEAREST_AGE_S)
    parser.add_argument("--min-valid-channels", type=int, default=DEFAULT_MIN_VALID_CHANNELS)
    parser.add_argument("--reference-start", help="Reference interval start as sample index or ISO UTC time.")
    parser.add_argument("--reference-end", help="Reference interval end as sample index or ISO UTC time.")
    parser.add_argument("--reference-latest-days", type=float, default=30.0)
    args = parser.parse_args()
    summary = write_quality_h5(
        args.phase_dir,
        args.output,
        threshold_deg=args.threshold_deg,
        max_nearest_age_s=args.max_nearest_age_s,
        min_valid_channels=args.min_valid_channels,
        reference_start_us=_parse_time_us(args.reference_start),
        reference_end_us=_parse_time_us(args.reference_end),
        reference_latest_days=args.reference_latest_days,
        phase_cache_dir=args.phase_cache_dir,
    )
    print(f"wrote {args.output}")
    for key, value in summary.items():
        print(f"{key} {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
