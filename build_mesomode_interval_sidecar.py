#!/usr/bin/env python3
"""Build an HDF5 sidecar of PANSY mesosphere-mode measurement intervals."""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import h5py
import numpy as np

DEFAULT_METADATA = Path("/mnt/data/juha/pansy/metadata/mesomode")
DEFAULT_OUTPUT = Path("/mnt/data/juha/pansy/metadata/mesomode_intervals.h5")
US_PER_SECOND = 1_000_000


def parse_hour_dir(name: str) -> datetime | None:
    try:
        return datetime.strptime(name, "%Y-%m-%dT%H-%M-%S").replace(tzinfo=timezone.utc)
    except ValueError:
        return None


def day_range(start: date, stop: date):
    day = start
    while day <= stop:
        yield day
        day += timedelta(days=1)


def discover_days(metadata_dir: Path) -> tuple[date, date, dict[date, list[Path]]]:
    by_day: dict[date, list[Path]] = defaultdict(list)
    for hour_dir in metadata_dir.iterdir():
        if not hour_dir.is_dir():
            continue
        hour_dt = parse_hour_dir(hour_dir.name)
        if hour_dt is None:
            continue
        files = sorted(hour_dir.glob("*.h5"))
        if files:
            by_day[hour_dt.date()].extend(files)
    if not by_day:
        raise RuntimeError(f"no hourly mesomode metadata files found under {metadata_dir}")
    return min(by_day), max(by_day), by_day


def read_file_intervals(path: Path) -> list[tuple[int, int]]:
    intervals: list[tuple[int, int]] = []
    with h5py.File(path, "r") as h:
        for key in h.keys():
            group = h[key]
            try:
                start = int(group["start"][()])
                end = int(group["end"][()])
            except Exception:
                continue
            if end > start:
                intervals.append((start, end))
    return intervals


def merge_intervals(intervals: list[tuple[int, int]], tolerance_us: int) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + tolerance_us:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return [(start, end) for start, end in merged]


def clip_interval_to_day(start: int, end: int, utc_day: date) -> tuple[int, int] | None:
    day_start = int(datetime(utc_day.year, utc_day.month, utc_day.day, tzinfo=timezone.utc).timestamp() * US_PER_SECOND)
    day_end = day_start + 86400 * US_PER_SECOND
    clipped_start = max(start, day_start)
    clipped_end = min(end, day_end)
    if clipped_end <= clipped_start:
        return None
    return clipped_start, clipped_end


def sample_to_iso_utc(samples: np.ndarray) -> np.ndarray:
    values = [
        datetime.fromtimestamp(int(sample) / US_PER_SECOND, timezone.utc).isoformat(timespec="microseconds").encode("ascii")
        for sample in samples
    ]
    return np.asarray(values, dtype="S32")


def write_sidecar(
    output: Path,
    intervals: list[tuple[int, int, str]],
    daily_rows: list[tuple[str, int, float]],
    metadata_dir: Path,
    merge_tolerance_us: int,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    starts = np.asarray([row[0] for row in intervals], dtype=np.int64)
    ends = np.asarray([row[1] for row in intervals], dtype=np.int64)
    days = np.asarray([row[2].encode("ascii") for row in intervals], dtype="S10")
    durations_s = (ends - starts).astype(np.float64) / US_PER_SECOND
    t0_unix = starts.astype(np.float64) / US_PER_SECOND
    t1_unix = ends.astype(np.float64) / US_PER_SECOND

    daily_day = np.asarray([row[0].encode("ascii") for row in daily_rows], dtype="S10")
    daily_count = np.asarray([row[1] for row in daily_rows], dtype=np.int64)
    daily_hours = np.asarray([row[2] for row in daily_rows], dtype=np.float64)

    with h5py.File(output, "w") as h:
        h.attrs["description"] = "PANSY mesosphere-mode measurement intervals from metadata/mesomode"
        h.attrs["source_metadata_dir"] = str(metadata_dir)
        h.attrs["time_scale"] = "Unix UTC"
        h.attrs["sample_rate_hz"] = US_PER_SECOND
        h.attrs["merge_tolerance_us"] = int(merge_tolerance_us)
        h.attrs["created_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")

        g = h.create_group("intervals")
        g.create_dataset("t0_sample", data=starts, compression="gzip", shuffle=True)
        g.create_dataset("t1_sample", data=ends, compression="gzip", shuffle=True)
        g.create_dataset("t0_unix", data=t0_unix, compression="gzip", shuffle=True)
        g.create_dataset("t1_unix", data=t1_unix, compression="gzip", shuffle=True)
        g.create_dataset("duration_s", data=durations_s, compression="gzip", shuffle=True)
        g.create_dataset("utc_day", data=days, compression="gzip", shuffle=True)
        g.create_dataset("t0_iso_utc", data=sample_to_iso_utc(starts), compression="gzip", shuffle=True)
        g.create_dataset("t1_iso_utc", data=sample_to_iso_utc(ends), compression="gzip", shuffle=True)

        gd = h.create_group("daily")
        gd.create_dataset("utc_day", data=daily_day, compression="gzip", shuffle=True)
        gd.create_dataset("interval_count", data=daily_count, compression="gzip", shuffle=True)
        gd.create_dataset("mesosphere_hours", data=daily_hours, compression="gzip", shuffle=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--merge-tolerance-us",
        type=int,
        default=2_000_000,
        help="Merge intervals separated by at most this many microseconds.",
    )
    args = parser.parse_args()

    first_day, last_day, files_by_day = discover_days(args.metadata_dir)
    carried: list[tuple[int, int]] = []
    all_intervals: list[tuple[int, int, str]] = []
    daily_rows: list[tuple[str, int, float]] = []

    for utc_day in day_range(first_day, last_day):
        raw = carried
        carried = []
        for file_path in files_by_day.get(utc_day, []):
            raw.extend(read_file_intervals(file_path))
        merged = merge_intervals(raw, args.merge_tolerance_us)

        day_intervals: list[tuple[int, int]] = []
        for start, end in merged:
            clipped = clip_interval_to_day(start, end, utc_day)
            if clipped is not None:
                day_intervals.append(clipped)
            end_day = datetime.fromtimestamp(end / US_PER_SECOND, timezone.utc).date()
            if end_day > utc_day:
                carried.append((start, end))

        # Merge again after clipping so overlapping open-block records do not
        # inflate daily interval counts.
        day_intervals = merge_intervals(day_intervals, args.merge_tolerance_us)
        day_text = utc_day.isoformat()
        for start, end in day_intervals:
            all_intervals.append((start, end, day_text))
        daily_hours = sum((end - start) for start, end in day_intervals) / (US_PER_SECOND * 3600.0)
        daily_rows.append((day_text, len(day_intervals), daily_hours))
        print(f"{day_text} intervals={len(day_intervals):4d} mesosphere_hours={daily_hours:6.2f}", flush=True)

    write_sidecar(args.output, all_intervals, daily_rows, args.metadata_dir, args.merge_tolerance_us)
    print(f"wrote {args.output}")
    print(f"intervals={len(all_intervals)} days={len(daily_rows)}")


if __name__ == "__main__":
    main()
