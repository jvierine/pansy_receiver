#!/usr/bin/env python3
"""Collect pathological PANSY event analyses into an HDF5 catalogue."""

from __future__ import annotations

import argparse
import datetime as dt
import re
from pathlib import Path

import h5py
import numpy as np


TIMEOUT_RE = re.compile(r"event_timeout\s+sample_idx\s+(\d+)\s+timeout_s\s+([0-9.]+)")


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def day_from_sample(sample_idx: int) -> str:
    return dt.datetime.fromtimestamp(sample_idx / 1e6, tz=dt.timezone.utc).strftime("%Y-%m-%d")


def read_existing(path: Path) -> dict[int, dict[str, str]]:
    if not path.exists():
        return {}
    out: dict[int, dict[str, str]] = {}
    with h5py.File(path, "r") as h5:
        samples = np.asarray(h5.get("sample_idx", []), dtype=np.int64)
        reasons = h5.get("reason", [])
        details = h5.get("detail", [])
        source_logs = h5.get("source_log", [])
        observed_utc = h5.get("observed_utc", [])
        for i, sample in enumerate(samples):
            out[int(sample)] = {
                "reason": decode(reasons[i]) if i < len(reasons) else "",
                "detail": decode(details[i]) if i < len(details) else "",
                "source_log": decode(source_logs[i]) if i < len(source_logs) else "",
                "observed_utc": decode(observed_utc[i]) if i < len(observed_utc) else "",
            }
    return out


def decode(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def merge_record(records: dict[int, dict[str, str]], sample_idx: int, reason: str, detail: str, source_log: str) -> None:
    observed = utc_now()
    if sample_idx in records:
        prev = records[sample_idx]
        reasons = {part for part in prev.get("reason", "").split(";") if part}
        reasons.add(reason)
        details = [part for part in prev.get("detail", "").split("\n") if part]
        if detail and detail not in details:
            details.append(detail)
        sources = {part for part in prev.get("source_log", "").split(";") if part}
        if source_log:
            sources.add(source_log)
        records[sample_idx] = {
            "reason": ";".join(sorted(reasons)),
            "detail": "\n".join(details),
            "source_log": ";".join(sorted(sources)),
            "observed_utc": prev.get("observed_utc") or observed,
        }
    else:
        records[sample_idx] = {
            "reason": reason,
            "detail": detail,
            "source_log": source_log,
            "observed_utc": observed,
        }


def collect_timeouts(records: dict[int, dict[str, str]], logs: list[Path]) -> None:
    for root in logs:
        paths = [root] if root.is_file() else sorted(root.rglob("*.log"))
        for path in paths:
            try:
                text = path.read_text(errors="replace")
            except OSError:
                continue
            for match in TIMEOUT_RE.finditer(text):
                sample_idx = int(match.group(1))
                timeout_s = float(match.group(2))
                merge_record(
                    records,
                    sample_idx,
                    "event_timeout",
                    f"event subprocess exceeded {timeout_s:g} s wall-clock timeout",
                    str(path),
                )


def write_h5(path: Path, records: dict[int, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    samples = np.asarray(sorted(records), dtype=np.int64)
    str_dtype = h5py.string_dtype("utf-8")
    tmp = path.with_suffix(path.suffix + ".tmp")
    with h5py.File(tmp, "w") as h5:
        h5.attrs["schema_version"] = "pansy_pathological_events_v1"
        h5.attrs["source_program"] = "collect_pathological_events.py"
        h5.attrs["updated_utc"] = utc_now()
        h5.attrs["event_count"] = int(len(samples))
        h5.create_dataset("sample_idx", data=samples)
        h5.create_dataset("utc_day", data=np.asarray([day_from_sample(int(s)) for s in samples], dtype=object), dtype=str_dtype)
        h5.create_dataset("reason", data=np.asarray([records[int(s)]["reason"] for s in samples], dtype=object), dtype=str_dtype)
        h5.create_dataset("detail", data=np.asarray([records[int(s)]["detail"] for s in samples], dtype=object), dtype=str_dtype)
        h5.create_dataset("source_log", data=np.asarray([records[int(s)]["source_log"] for s in samples], dtype=object), dtype=str_dtype)
        h5.create_dataset("observed_utc", data=np.asarray([records[int(s)]["observed_utc"] for s in samples], dtype=object), dtype=str_dtype)
    tmp.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-h5", type=Path, required=True)
    parser.add_argument("--log", type=Path, action="append", default=[], help="Log file or directory tree to scan for event_timeout markers.")
    parser.add_argument("--manual-sample", type=int, action="append", default=[], help="Pathological sample index to add manually.")
    parser.add_argument("--manual-reason", default="manual_pathological_case")
    parser.add_argument("--manual-detail", default="")
    parser.add_argument("--replace", action="store_true", help="Ignore existing output before writing.")
    args = parser.parse_args()

    records = {} if args.replace else read_existing(args.output_h5)
    for sample_idx in args.manual_sample:
        merge_record(records, int(sample_idx), args.manual_reason, args.manual_detail, "manual")
    collect_timeouts(records, args.log)
    write_h5(args.output_h5, records)
    print(f"pathological_events_h5 {args.output_h5}")
    print(f"pathological_events_count {len(records)}")
    for sample_idx in sorted(records):
        print(sample_idx, day_from_sample(sample_idx), records[sample_idx]["reason"])


if __name__ == "__main__":
    main()
