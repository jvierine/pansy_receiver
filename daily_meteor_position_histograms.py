#!/usr/bin/env python3
"""Build daily meteor-position histograms on the interferometry search grid."""

from __future__ import annotations

import argparse
import datetime as dt
import multiprocessing as mp
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

from plot_meteor_position_histograms import (
    add_histogram_panel,
    format_direction_axis,
    interferometry_grid_offsets_deg,
)


METRIC_DTYPE = np.dtype(
    [
        ("utc_day", "S10"),
        ("diagnostic_files", "i8"),
        ("fitted_events", "i8"),
        ("skipped_events", "i8"),
        ("calibration_rejected_events", "i8"),
        ("calibration_checked_hours", "i8"),
        ("calibration_bad_hours", "i8"),
        ("total_points", "i8"),
        ("points_in_view", "i8"),
        ("median_snr", "f8"),
        ("p90_snr", "f8"),
        ("median_coherence", "f8"),
        ("p10_coherence", "f8"),
    ]
)

HOUR_QUALITY_DTYPE = np.dtype(
    [
        ("hour_start_sample_idx", "i8"),
        ("diagnostic_files", "i8"),
        ("good", "?"),
        ("available", "?"),
        ("max_abs_diff_deg", "f8"),
        ("rms_diff_deg", "f8"),
        ("valid_channel_count", "i8"),
        ("reason", "S64"),
    ]
)


def read_selected_record(path: Path, use_selection_keep: bool = True) -> dict | None:
    try:
        with h5py.File(path, "r") as h:
            if bool(h.attrs.get("tx_phase_quality_good", True)) is False:
                return {"calibration_reject": True}
            label = str(h.attrs.get("selected_hypothesis", ""))
            if not label or "hypotheses" not in h or label not in h["hypotheses"] or "candidates" not in h:
                return None
            hyp = h["hypotheses"][label]
            if bool(hyp.attrs.get("combined_reject", False)):
                return None
            if "candidate_indices" not in hyp or "beam_id" not in hyp:
                return None
            cand = h["candidates"]
            needed = ("grid_row", "grid_col", "snr", "coherence")
            if any(name not in cand for name in needed):
                return None
            idx = np.asarray(hyp["candidate_indices"], dtype=np.int64)
            beam_id = np.asarray(hyp["beam_id"], dtype=np.int64)
            n = min(len(idx), len(beam_id))
            idx = idx[:n]
            beam_id = beam_id[:n]
            keep = np.ones(n, dtype=bool)
            if use_selection_keep and "selection_keep" in hyp:
                keep &= np.asarray(hyp["selection_keep"], dtype=bool)[:n]
            grid_row_all = np.asarray(cand["grid_row"], dtype=np.int64)
            valid_idx = (idx >= 0) & (idx < len(grid_row_all))
            keep &= valid_idx
            idx = idx[keep]
            if len(idx) == 0:
                return None
            return {
                "grid_row": grid_row_all[idx],
                "grid_col": np.asarray(cand["grid_col"], dtype=np.int64)[idx],
                "beam_id": beam_id[keep],
                "snr": np.asarray(cand["snr"], dtype=np.float64)[idx],
                "coherence": np.asarray(cand["coherence"], dtype=np.float64)[idx],
            }
    except (OSError, KeyError, IndexError):
        return None


def sample_idx_from_diagnostics_path(path: Path) -> int | None:
    stem = path.stem
    prefix = "pansy_disambiguation_diagnostics_"
    if not stem.startswith(prefix):
        return None
    try:
        return int(stem[len(prefix) :])
    except ValueError:
        return None


def classify_file_calibration(path: Path, calibration_context: dict | None) -> np.ndarray:
    if calibration_context is None:
        row = np.zeros((), dtype=HOUR_QUALITY_DTYPE)
        row["good"] = True
        row["reason"] = b"not_checked"
        return row
    sample_idx = sample_idx_from_diagnostics_path(path)
    row = np.zeros((), dtype=HOUR_QUALITY_DTYPE)
    row["diagnostic_files"] = 1
    row["hour_start_sample_idx"] = -1 if sample_idx is None else int((sample_idx // 1_000_000 // 3600) * 3600 * 1_000_000)
    row["good"] = False
    row["available"] = False
    row["max_abs_diff_deg"] = np.nan
    row["rms_diff_deg"] = np.nan
    row["valid_channel_count"] = 0
    row["reason"] = b"missing_sample_idx"
    if sample_idx is None:
        return row
    try:
        from plot_interferometric_disambiguation import cut_tx_reference_phase_quality, load_cut

        cut = load_cut(calibration_context["cut_dir"], int(sample_idx))
        quality = cut_tx_reference_phase_quality(
            cut,
            calibration_context["reference_phase_rad"],
            threshold_deg=calibration_context["threshold_deg"],
            min_valid_channels=calibration_context["min_valid_channels"],
        )
        row["good"] = bool(quality.get("good", False))
        row["available"] = bool(quality.get("available", False))
        row["max_abs_diff_deg"] = float(quality.get("max_abs_diff_deg", np.nan))
        row["rms_diff_deg"] = float(quality.get("rms_diff_deg", np.nan))
        row["valid_channel_count"] = int(quality.get("valid_channel_count", 0))
        row["reason"] = str(quality.get("reason", "unknown")).encode()[:64]
    except Exception as exc:
        row["reason"] = f"calibration_check_failed:{type(exc).__name__}".encode()[:64]
    return row


def add_files_to_histogram(
    files: list[Path],
    all_counts: np.ndarray,
    grid_n: int,
    use_selection_keep: bool,
    calibration_context: dict | None = None,
) -> tuple[int, int, int, int, list[np.ndarray], list[np.ndarray]]:
    event_count = 0
    skipped = 0
    calibration_rejected = 0
    total_points = 0
    snr_chunks: list[np.ndarray] = []
    coh_chunks: list[np.ndarray] = []
    for path in files:
        if calibration_context is not None and calibration_context.get("mode") == "event":
            quality_row = classify_file_calibration(path, calibration_context)
            if not bool(quality_row["good"]):
                skipped += 1
                calibration_rejected += 1
                continue
        row = read_selected_record(path, use_selection_keep=use_selection_keep)
        if row is None:
            skipped += 1
            continue
        if row.get("calibration_reject", False):
            skipped += 1
            calibration_rejected += 1
            continue
        grid_row = np.asarray(row["grid_row"], dtype=np.int64)
        grid_col = np.asarray(row["grid_col"], dtype=np.int64)
        valid = (grid_row >= 0) & (grid_row < grid_n) & (grid_col >= 0) & (grid_col < grid_n)
        if not np.any(valid):
            skipped += 1
            continue
        event_count += 1
        grid_row = grid_row[valid]
        grid_col = grid_col[valid]
        np.add.at(all_counts, (grid_row, grid_col), 1)
        snr = np.asarray(row["snr"], dtype=np.float64)[valid]
        coh = np.asarray(row["coherence"], dtype=np.float64)[valid]
        snr_chunks.append(snr[np.isfinite(snr)])
        coh_chunks.append(coh[np.isfinite(coh)])
        total_points += len(grid_row)
    return event_count, skipped, calibration_rejected, total_points, snr_chunks, coh_chunks


def collect_day(day_dir: Path, grid_n: int, use_selection_keep: bool) -> tuple[dict, tuple[float, float, float, float]]:
    files = sorted(day_dir.glob("pansy_disambiguation_diagnostics_*.h5"))
    all_counts = np.zeros((grid_n, grid_n), dtype=np.int32)
    event_count, skipped, calibration_rejected, total_points, snr_chunks, coh_chunks = add_files_to_histogram(
        files,
        all_counts,
        grid_n,
        use_selection_keep,
        calibration_context=None,
    )
    counts = {
        "all_counts": all_counts,
        "event_count": event_count,
        "skipped_count": skipped,
        "calibration_rejected_count": calibration_rejected,
        "file_count": len(files),
        "total_points": total_points,
    }
    snr_all = np.concatenate(snr_chunks) if snr_chunks else np.asarray([], dtype=np.float64)
    coh_all = np.concatenate(coh_chunks) if coh_chunks else np.asarray([], dtype=np.float64)
    quality = (
        float(np.nanmedian(snr_all)) if len(snr_all) else np.nan,
        float(np.nanpercentile(snr_all, 90)) if len(snr_all) else np.nan,
        float(np.nanmedian(coh_all)) if len(coh_all) else np.nan,
        float(np.nanpercentile(coh_all, 10)) if len(coh_all) else np.nan,
    )
    return counts, quality


def hourly_files_for_day(day_dir: Path) -> list[list[Path]]:
    by_hour = [[] for _ in range(24)]
    for path in sorted(day_dir.glob("pansy_disambiguation_diagnostics_*.h5")):
        sample_idx = sample_idx_from_diagnostics_path(path)
        if sample_idx is None:
            continue
        hour = dt.datetime.fromtimestamp(sample_idx / 1e6, tz=dt.timezone.utc).hour
        by_hour[hour].append(path)
    return by_hour


def classify_hour_calibration(hour_files: list[Path], calibration_context: dict | None) -> np.ndarray:
    row = np.zeros((), dtype=HOUR_QUALITY_DTYPE)
    row["diagnostic_files"] = len(hour_files)
    row["hour_start_sample_idx"] = -1
    row["good"] = True
    row["available"] = False
    row["max_abs_diff_deg"] = np.nan
    row["rms_diff_deg"] = np.nan
    row["valid_channel_count"] = 0
    row["reason"] = b"not_checked"
    if not hour_files or calibration_context is None:
        return row
    sample_idx = sample_idx_from_diagnostics_path(hour_files[len(hour_files) // 2])
    if sample_idx is None:
        row["good"] = False
        row["reason"] = b"missing_sample_idx"
        return row
    row["hour_start_sample_idx"] = int((sample_idx // 1_000_000 // 3600) * 3600 * 1_000_000)
    try:
        from plot_interferometric_disambiguation import cut_tx_reference_phase_quality, load_cut

        cut = load_cut(calibration_context["cut_dir"], int(sample_idx))
        quality = cut_tx_reference_phase_quality(
            cut,
            calibration_context["reference_phase_rad"],
            threshold_deg=calibration_context["threshold_deg"],
            min_valid_channels=calibration_context["min_valid_channels"],
        )
        row["good"] = bool(quality.get("good", False))
        row["available"] = bool(quality.get("available", False))
        row["max_abs_diff_deg"] = float(quality.get("max_abs_diff_deg", np.nan))
        row["rms_diff_deg"] = float(quality.get("rms_diff_deg", np.nan))
        row["valid_channel_count"] = int(quality.get("valid_channel_count", 0))
        row["reason"] = str(quality.get("reason", "unknown")).encode()[:64]
    except Exception as exc:
        row["good"] = False
        row["available"] = False
        row["reason"] = f"calibration_check_failed:{type(exc).__name__}".encode()[:64]
    return row


def collect_day_hourly(
    day_dir: Path,
    grid_n: int,
    use_selection_keep: bool,
    calibration_context: dict | None = None,
) -> tuple[dict, tuple[float, float, float, float]]:
    all_counts = np.zeros((grid_n, grid_n), dtype=np.int32)
    snr_chunks: list[np.ndarray] = []
    coh_chunks: list[np.ndarray] = []
    event_count = 0
    skipped = 0
    total_points = 0
    calibration_rejected = 0
    file_count = 0
    hour_quality = []
    for hour_files in hourly_files_for_day(day_dir):
        file_count += len(hour_files)
        quality_row = classify_hour_calibration(hour_files, calibration_context)
        hour_quality.append(quality_row)
        if calibration_context is not None and hour_files and not bool(quality_row["good"]):
            skipped += len(hour_files)
            calibration_rejected += len(hour_files)
            continue
        h_events, h_skipped, h_calibration_rejected, h_points, h_snr, h_coh = add_files_to_histogram(
            hour_files,
            all_counts,
            grid_n,
            use_selection_keep,
            calibration_context=calibration_context,
        )
        event_count += h_events
        skipped += h_skipped
        calibration_rejected += h_calibration_rejected
        total_points += h_points
        snr_chunks.extend(h_snr)
        coh_chunks.extend(h_coh)
    counts = {
        "all_counts": all_counts,
        "event_count": event_count,
        "skipped_count": skipped,
        "calibration_rejected_count": calibration_rejected,
        "hour_quality": np.asarray(hour_quality, dtype=HOUR_QUALITY_DTYPE),
        "file_count": file_count,
        "total_points": total_points,
    }
    snr_all = np.concatenate(snr_chunks) if snr_chunks else np.asarray([], dtype=np.float64)
    coh_all = np.concatenate(coh_chunks) if coh_chunks else np.asarray([], dtype=np.float64)
    quality = (
        float(np.nanmedian(snr_all)) if len(snr_all) else np.nan,
        float(np.nanpercentile(snr_all, 90)) if len(snr_all) else np.nan,
        float(np.nanmedian(coh_all)) if len(coh_all) else np.nan,
        float(np.nanpercentile(coh_all, 10)) if len(coh_all) else np.nan,
    )
    return counts, quality


def in_view_mask(grid_n: int, extent_deg: float):
    _u, _v, east_grid, north_grid, _valid = interferometry_grid_offsets_deg(grid_n)
    in_view = (
        np.isfinite(east_grid)
        & np.isfinite(north_grid)
        & (east_grid >= -extent_deg)
        & (east_grid <= extent_deg)
        & (north_grid >= -extent_deg)
        & (north_grid <= extent_deg)
    )
    return east_grid, north_grid, in_view


def write_daily_channel(day: str, counts: dict, metrics, metadata_dir: Path, grid_n: int, extent_deg: float) -> None:
    metadata_dir.mkdir(parents=True, exist_ok=True)
    path = metadata_dir / f"meteor_position_histogram_{day}.h5"
    tmp = path.with_suffix(".tmp")
    if tmp.exists():
        tmp.unlink()
    with h5py.File(tmp, "w") as h:
        h.attrs["schema"] = "pansy_daily_meteor_position_histogram_v1"
        h.attrs["utc_day"] = day
        h.attrs["grid_n"] = int(grid_n)
        h.attrs["extent_deg"] = float(extent_deg)
        h.attrs["grid_source"] = "plot_interferometric_disambiguation.horizon_grid"
        h.attrs["histogram_source"] = "diagnostics selected-hypothesis candidate grid_row/grid_col"
        h.attrs["calibration_rejection"] = "skip diagnostics with bad tx_phase_quality_good or bad cut_tx_waveform_good"
        h.create_dataset("all_counts", data=np.asarray(counts["all_counts"], dtype=np.int32), compression="gzip", shuffle=True)
        h.create_dataset("metrics", data=np.asarray(metrics, dtype=METRIC_DTYPE))
        if "hour_quality" in counts:
            h.create_dataset("hour_quality", data=np.asarray(counts["hour_quality"], dtype=HOUR_QUALITY_DTYPE))
    tmp.replace(path)


def collect_days(events_dir: Path) -> list[str]:
    days = []
    for path in sorted(events_dir.iterdir()):
        if path.is_dir() and len(path.name) == 10 and path.name[4] == "-" and path.name[7] == "-":
            if any(path.glob("pansy_disambiguation_diagnostics_*.h5")):
                days.append(path.name)
    return days


def collect_day_task(task):
    i, day, events_dir, grid_n, use_selection_keep, calibration_context = task
    counts, quality = collect_day_hourly(Path(events_dir) / day, grid_n, use_selection_keep, calibration_context)
    return i, day, counts, quality


def build(
    events_dir: Path,
    metadata_dir: Path,
    output_h5: Path,
    grid_n: int,
    extent_deg: float,
    include_clipped: bool,
    workers: int,
    calibration_context: dict | None,
):
    days = collect_days(events_dir)
    east_grid, north_grid, view = in_view_mask(grid_n, extent_deg)
    counts_stack = np.zeros((len(days), grid_n, grid_n), dtype=np.int32)
    metrics = np.zeros(len(days), dtype=METRIC_DTYPE)
    tasks = [(i, day, str(events_dir), grid_n, not include_clipped, calibration_context) for i, day in enumerate(days)]
    if workers <= 1:
        iterator = map(collect_day_task, tasks)
    else:
        pool = mp.get_context("spawn").Pool(processes=workers)
        iterator = pool.imap_unordered(collect_day_task, tasks, chunksize=1)
    try:
        for i, day, counts, quality in iterator:
            median_snr, p90_snr, median_coh, p10_coh = quality
            all_counts = np.asarray(counts["all_counts"], dtype=np.int32)
            counts_stack[i] = all_counts
            metrics[i] = (
                day.encode(),
                counts["file_count"],
                counts["event_count"],
                counts["skipped_count"],
                counts["calibration_rejected_count"],
                int(np.count_nonzero(counts.get("hour_quality", np.zeros(0, dtype=HOUR_QUALITY_DTYPE))["available"])),
                int(np.count_nonzero(~counts.get("hour_quality", np.zeros(0, dtype=HOUR_QUALITY_DTYPE))["good"])),
                counts["total_points"],
                int(np.sum(np.where(view, all_counts, 0))),
                median_snr,
                p90_snr,
                median_coh,
                p10_coh,
            )
            write_daily_channel(day, counts, metrics[i], metadata_dir, grid_n, extent_deg)
            print(
                f"daily_position_histogram {day} files {counts['file_count']} "
                f"events {counts['event_count']} points {counts['total_points']}",
                flush=True,
            )
    finally:
        if workers > 1:
            pool.close()
            pool.join()
    order = np.argsort(metrics["utc_day"])
    days = [days[i] for i in order]
    counts_stack = counts_stack[order]
    metrics = metrics[order]
    output_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_h5, "w") as h:
        h.attrs["schema"] = "pansy_daily_meteor_position_histogram_stack_v1"
        h.attrs["grid_n"] = int(grid_n)
        h.attrs["extent_deg"] = float(extent_deg)
        h.attrs["grid_source"] = "plot_interferometric_disambiguation.horizon_grid"
        h.create_dataset("utc_day", data=np.asarray(days, dtype="S10"))
        h.create_dataset("east_grid_deg", data=east_grid.astype(np.float32), compression="gzip", shuffle=True)
        h.create_dataset("north_grid_deg", data=north_grid.astype(np.float32), compression="gzip", shuffle=True)
        h.create_dataset("in_view", data=view.astype(np.uint8), compression="gzip", shuffle=True)
        h.create_dataset("all_counts", data=counts_stack, compression="gzip", shuffle=True)
        h.create_dataset("metrics", data=metrics)
    return days, counts_stack, metrics


def build_one_day(
    day: str,
    events_dir: Path,
    metadata_dir: Path,
    grid_n: int,
    extent_deg: float,
    include_clipped: bool,
    calibration_context: dict | None,
):
    east_grid, north_grid, view = in_view_mask(grid_n, extent_deg)
    counts, quality = collect_day_hourly(events_dir / day, grid_n, use_selection_keep=not include_clipped, calibration_context=calibration_context)
    median_snr, p90_snr, median_coh, p10_coh = quality
    all_counts = np.asarray(counts["all_counts"], dtype=np.int32)
    metric = np.zeros((), dtype=METRIC_DTYPE)
    metric[...] = (
        day.encode(),
        counts["file_count"],
        counts["event_count"],
        counts["skipped_count"],
        counts["calibration_rejected_count"],
        int(np.count_nonzero(counts.get("hour_quality", np.zeros(0, dtype=HOUR_QUALITY_DTYPE))["available"])),
        int(np.count_nonzero(~counts.get("hour_quality", np.zeros(0, dtype=HOUR_QUALITY_DTYPE))["good"])),
        counts["total_points"],
        int(np.sum(np.where(view, all_counts, 0))),
        median_snr,
        p90_snr,
        median_coh,
        p10_coh,
    )
    write_daily_channel(day, counts, metric, metadata_dir, grid_n, extent_deg)
    print(f"daily_position_histogram {day} files {counts['file_count']} events {counts['event_count']} points {counts['total_points']}", flush=True)


def stack_from_daily_metadata(metadata_dir: Path, output_h5: Path, grid_n: int, extent_deg: float):
    files = sorted(metadata_dir.glob("meteor_position_histogram_*.h5"))
    east_grid, north_grid, view = in_view_mask(grid_n, extent_deg)
    days = []
    counts = []
    metrics = []
    for path in files:
        with h5py.File(path, "r") as h:
            days.append(str(h.attrs["utc_day"]))
            counts.append(np.asarray(h["all_counts"], dtype=np.int32))
            metrics.append(np.asarray(h["metrics"][()], dtype=METRIC_DTYPE))
    order = np.argsort(np.asarray(days, dtype="S10"))
    days = [days[i] for i in order]
    counts_stack = np.asarray([counts[i] for i in order], dtype=np.int32)
    metrics_arr = np.asarray([metrics[i] for i in order], dtype=METRIC_DTYPE)
    output_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_h5, "w") as h:
        h.attrs["schema"] = "pansy_daily_meteor_position_histogram_stack_v1"
        h.attrs["grid_n"] = int(grid_n)
        h.attrs["extent_deg"] = float(extent_deg)
        h.attrs["grid_source"] = "plot_interferometric_disambiguation.horizon_grid"
        h.create_dataset("utc_day", data=np.asarray(days, dtype="S10"))
        h.create_dataset("east_grid_deg", data=east_grid.astype(np.float32), compression="gzip", shuffle=True)
        h.create_dataset("north_grid_deg", data=north_grid.astype(np.float32), compression="gzip", shuffle=True)
        h.create_dataset("in_view", data=view.astype(np.uint8), compression="gzip", shuffle=True)
        h.create_dataset("all_counts", data=counts_stack, compression="gzip", shuffle=True)
        h.create_dataset("metrics", data=metrics_arr)
    return days, counts_stack, metrics_arr


def plot_frame(day: bytes, counts: np.ndarray, metric, east_grid, north_grid, view, out: Path, extent_deg: float):
    image = np.where(view, counts, 0)
    rows = np.flatnonzero(np.any(view, axis=1))
    cols = np.flatnonzero(np.any(view, axis=0))
    rs = slice(int(rows[0]), int(rows[-1]) + 1)
    cs = slice(int(cols[0]), int(cols[-1]) + 1)
    fig, ax = plt.subplots(1, 1, figsize=(5.6, 5.0), constrained_layout=True)
    im, _ = add_histogram_panel(
        ax,
        east_grid[rs, cs],
        north_grid[rs, cs],
        image[rs, cs],
        f"{day.decode()} all TX beams\nN={int(metric['points_in_view'])} points",
    )
    cb = fig.colorbar(im, ax=ax, shrink=0.82)
    cb.set_label("Meteor detections per bin")
    format_direction_axis(ax, extent_deg)
    ax.text(
        0.02,
        0.02,
        f"events={int(metric['fitted_events'])}  med SNR={float(metric['median_snr']):.1f}  med coh={float(metric['median_coherence']):.2f}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        color="white",
        bbox={"facecolor": "black", "alpha": 0.35, "edgecolor": "none", "pad": 2},
    )
    fig.savefig(out, dpi=160)
    plt.close(fig)


def make_gif(stack_h5: Path, output_gif: Path, frame_dir: Path, extent_deg: float):
    try:
        import imageio.v2 as imageio
    except ImportError:
        import imageio

    with h5py.File(stack_h5, "r") as h:
        days = h["utc_day"][()]
        counts = h["all_counts"][()]
        metrics = h["metrics"][()]
        east_grid = h["east_grid_deg"][()]
        north_grid = h["north_grid_deg"][()]
        view = h["in_view"][()].astype(bool)
    frame_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    for i, day in enumerate(days):
        frame = frame_dir / f"meteor_position_{day.decode()}.png"
        plot_frame(day, counts[i], metrics[i], east_grid, north_grid, view, frame, extent_deg)
        frames.append(imageio.imread(frame))
    output_gif.parent.mkdir(parents=True, exist_ok=True)
    imageio.mimsave(output_gif, frames, duration=0.35)


def load_calibration_context(cut_dir: Path | None, tx_phase_quality_h5: Path | None) -> dict | None:
    if cut_dir is None or tx_phase_quality_h5 is None:
        return None
    with h5py.File(tx_phase_quality_h5, "r") as h5:
        return {
            "cut_dir": cut_dir,
            "reference_phase_rad": np.asarray(h5["reference_phase_rad"], dtype=np.float32),
            "threshold_deg": float(h5.attrs.get("threshold_deg", 30.0)),
            "min_valid_channels": int(h5.attrs.get("min_valid_channels", 6)),
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events-dir", type=Path, default=Path("test_plots"))
    parser.add_argument("--metadata-dir", type=Path, default=Path("data/metadata/tracks"))
    parser.add_argument("--output-h5", type=Path, default=Path("test_plots/current_daily_meteor_position_histograms.h5"))
    parser.add_argument("--output-gif", type=Path, default=Path("test_plots/current_daily_meteor_position_histograms.gif"))
    parser.add_argument("--frame-dir", type=Path, default=Path("test_plots/meteor_position_histogram_frames"))
    parser.add_argument("--extent-deg", type=float, default=15.0)
    parser.add_argument("--grid-n", type=int, default=501)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--day", type=str, nargs="*")
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--stack-from-metadata", action="store_true")
    parser.add_argument("--cut-dir", type=Path, default=None, help="Cut metadata directory used to classify hourly TX phase calibration from cut TX pulses.")
    parser.add_argument("--tx-phase-quality-h5", type=Path, default=None, help="TX phase quality HDF5 providing the reference phase vector for hourly calibration checks.")
    parser.add_argument("--require-good-calibration-hours", action="store_true", help="Skip all events in hours whose representative cut TX phase does not match the reference phase vector.")
    parser.add_argument(
        "--calibration-check-mode",
        choices=("hour", "event"),
        default="hour",
        help="Use one representative cut per hour, or check every event cut before adding points.",
    )
    parser.add_argument("--include-clipped", action="store_true")
    args = parser.parse_args()
    calibration_context = load_calibration_context(args.cut_dir, args.tx_phase_quality_h5) if args.require_good_calibration_hours else None
    if calibration_context is not None:
        calibration_context["mode"] = args.calibration_check_mode

    if args.stack_from_metadata:
        days, counts, metrics = stack_from_daily_metadata(args.metadata_dir, args.output_h5, args.grid_n, args.extent_deg)
    elif args.day:
        for day in args.day:
            build_one_day(day, args.events_dir, args.metadata_dir, args.grid_n, args.extent_deg, args.include_clipped, calibration_context)
        if args.metadata_only:
            return
        days, counts, metrics = stack_from_daily_metadata(args.metadata_dir, args.output_h5, args.grid_n, args.extent_deg)
    else:
        days, counts, metrics = build(
            args.events_dir,
            args.metadata_dir,
            args.output_h5,
            args.grid_n,
            args.extent_deg,
            args.include_clipped,
            max(1, args.workers),
            calibration_context,
        )
    if not args.metadata_only:
        make_gif(args.output_h5, args.output_gif, args.frame_dir, args.extent_deg)
    print(args.output_h5)
    print(args.output_gif)
    print(f"daily_position_histograms_done days {len(days)} total_points {int(np.sum(metrics['total_points']))}")


if __name__ == "__main__":
    main()
