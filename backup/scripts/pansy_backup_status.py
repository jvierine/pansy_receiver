#!/usr/bin/env python3
import datetime as dt
import json
import os
import re
from pathlib import Path


DATE_PATTERN = re.compile(r"(20\d{2})[-_]?(\d{2})[-_]?(\d{2})")


def env(name, default):
    return os.environ.get(name, default)


def load_env_file(path):
    values = {}
    if not path.exists():
        return values
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key] = value.strip().strip('"')
    return values


def channel_summary(root):
    newest = None
    newest_date = None
    count = 0
    total_bytes = 0
    if not root.exists():
        return None, None, 0, 0
    for path in root.iterdir():
        try:
            stat = path.stat()
        except OSError:
            continue
        count += 1
        total_bytes += stat.st_size
        if newest is None or stat.st_mtime > newest:
            newest = stat.st_mtime

        match = DATE_PATTERN.search(path.name)
        if match is not None:
            entry_date = "-".join(match.groups())
            if newest_date is None or entry_date > newest_date:
                newest_date = entry_date
    return newest, newest_date, count, total_bytes


def iso_from_ts(ts):
    if ts is None:
        return None
    return dt.datetime.fromtimestamp(ts, tz=dt.timezone.utc).isoformat()


def read_rsync_status(path):
    status = {}
    if not path.exists():
        return status
    for line in path.read_text(errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) != 4:
            continue
        channel, started, finished, rc = parts
        status[channel] = {
            "last_rsync_started": started,
            "last_rsync_finished": finished,
            "last_rsync_exit_code": int(rc) if rc.isdigit() else rc,
        }
    return status


def parse_active_rsync_line(line, now):
    parts = line.split("\t")
    if len(parts) not in (3, 5):
        return None, None

    channel, started, pid = parts[:3]
    active = {
        "rsync_active": True,
        "active_rsync_started": started,
        "active_rsync_pid": int(pid) if pid.isdigit() else pid,
    }
    if len(parts) == 5:
        files_total, files_left = parts[3], parts[4]
        active["active_rsync_files_total"] = int(files_total) if files_total.isdigit() else None
        active["active_rsync_files_left"] = int(files_left) if files_left.isdigit() else None
    try:
        started_dt = dt.datetime.fromisoformat(started.replace("Z", "+00:00"))
        active["active_rsync_age_seconds"] = now - started_dt.timestamp()
    except ValueError:
        active["active_rsync_age_seconds"] = None
    return channel, active


def read_active_rsync(path, now):
    if not path.exists():
        return {}

    if path.is_dir():
        active = {}
        for active_file in path.glob("*.tsv"):
            lines = active_file.read_text(errors="replace").splitlines()
            if not lines:
                continue
            channel, record = parse_active_rsync_line(lines[-1], now)
            if channel is not None:
                active[channel] = record
        return active

    lines = path.read_text(errors="replace").splitlines()
    if not lines:
        return {}
    channel, record = parse_active_rsync_line(lines[-1], now)
    return {} if channel is None else {channel: record}


def read_local_mirror_status(path):
    if not path.exists():
        return {}
    lines = path.read_text(errors="replace").splitlines()
    if not lines:
        return {}
    parts = lines[-1].split("\t")
    if len(parts) != 5:
        return {}
    source, destination, started, finished, rc = parts
    return {
        "source": source,
        "destination": destination,
        "last_started": started,
        "last_finished": finished,
        "last_exit_code": int(rc) if rc.isdigit() else rc,
    }


def read_local_mirror_active(path, now):
    if not path.exists():
        return {"active": False}
    lines = path.read_text(errors="replace").splitlines()
    if not lines:
        return {"active": False}
    parts = lines[-1].split("\t")
    if len(parts) != 3:
        return {"active": False}
    source, destination, started = parts
    active = {
        "active": True,
        "source": source,
        "destination": destination,
        "active_started": started,
    }
    try:
        started_dt = dt.datetime.fromisoformat(started.replace("Z", "+00:00"))
        active["active_age_seconds"] = now - started_dt.timestamp()
    except ValueError:
        active["active_age_seconds"] = None
    return active


def main():
    config_path = Path(env("PANSY_BACKUP_CONFIG", str(Path.home() / ".config/pansy-backup/pansy-backup.env")))
    config = load_env_file(config_path)
    local_root = Path(config.get("LOCAL_METADATA_ROOT", "/mnt/data/juha/pansy/metadata"))
    mirror_root = Path(config.get("LOCAL_MIRROR_METADATA_ROOT", "/urdr/data/juha/pansy/metadata"))
    state_dir = Path(config.get("STATE_DIR", "/mnt/data/juha/pansy/backup_state"))
    web_dir = Path(config.get("WEB_BUILD_DIR", "/mnt/data/juha/pansy/web"))
    channels = config.get(
        "METADATA_CHANNELS",
        "clock cut detections mesomode phase simple_meteor_fit tx",
    ).split()

    web_dir.mkdir(parents=True, exist_ok=True)
    rsync_status = read_rsync_status(state_dir / "rsync_status.tsv")
    now = dt.datetime.now(dt.timezone.utc).timestamp()
    active_rsync_dir = state_dir / "rsync_active"
    active_rsync = read_active_rsync(active_rsync_dir, now)
    if not active_rsync_dir.exists():
        active_rsync = read_active_rsync(state_dir / "rsync_active.tsv", now)
    mirror_status = read_local_mirror_status(state_dir / "local_mirror_status.tsv")
    mirror_active = read_local_mirror_active(state_dir / "local_mirror_active.tsv", now)

    channel_status = []
    mirror_channels = []
    for channel in channels:
        newest, newest_date, count, total_bytes = channel_summary(local_root / channel)
        mirror_newest, mirror_newest_date, mirror_count, mirror_total_bytes = channel_summary(
            mirror_root / channel)
        record = {
            "channel": channel,
            "latest_file_date": newest_date,
            "newest_entry_utc": iso_from_ts(newest),
            "age_seconds": None if newest is None else now - newest,
            "entry_count": count,
            "total_bytes": total_bytes,
            "rsync_active": False,
        }
        record.update(rsync_status.get(channel, {}))
        record.update(active_rsync.get(channel, {}))
        channel_status.append(record)
        mirror_channels.append({
            "channel": channel,
            "source_latest_file_date": newest_date,
            "mirror_latest_file_date": mirror_newest_date,
            "source_newest_entry_utc": iso_from_ts(newest),
            "mirror_newest_entry_utc": iso_from_ts(mirror_newest),
            "source_entry_count": count,
            "mirror_entry_count": mirror_count,
            "source_total_bytes": total_bytes,
            "mirror_total_bytes": mirror_total_bytes,
            "entry_count_delta": count - mirror_count,
            "byte_delta": total_bytes - mirror_total_bytes,
            "latest_file_date_match": newest_date == mirror_newest_date,
        })

    mirror = {
        "source_root": str(local_root),
        "destination_root": str(mirror_root),
        "active": False,
        "channels": mirror_channels,
    }
    mirror.update(mirror_status)
    mirror.update(mirror_active)

    status = {
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "local_metadata_root": str(local_root),
        "local_mirror_metadata_root": str(mirror_root),
        "channels": channel_status,
        "local_mirror": mirror,
    }
    (web_dir / "status.json").write_text(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
