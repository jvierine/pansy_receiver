#!/usr/bin/env python3
import datetime as dt
import json
import os
from pathlib import Path


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
    count = 0
    total_bytes = 0
    if not root.exists():
        return None, 0, 0
    for path in root.iterdir():
        try:
            stat = path.stat()
        except OSError:
            continue
        count += 1
        total_bytes += stat.st_size
        if newest is None or stat.st_mtime > newest:
            newest = stat.st_mtime
    return newest, count, total_bytes


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


def main():
    config_path = Path(env("PANSY_BACKUP_CONFIG", str(Path.home() / ".config/pansy-backup/pansy-backup.env")))
    config = load_env_file(config_path)
    local_root = Path(config.get("LOCAL_METADATA_ROOT", "/mnt/data/juha/pansy/metadata"))
    state_dir = Path(config.get("STATE_DIR", "/mnt/data/juha/pansy/backup_state"))
    web_dir = Path(config.get("WEB_BUILD_DIR", "/mnt/data/juha/pansy/web"))
    channels = config.get(
        "METADATA_CHANNELS",
        "clock cut detections mesomode phase simple_meteor_fit tx",
    ).split()

    web_dir.mkdir(parents=True, exist_ok=True)
    rsync_status = read_rsync_status(state_dir / "rsync_status.tsv")
    now = dt.datetime.now(dt.timezone.utc).timestamp()

    channel_status = []
    for channel in channels:
        newest, count, total_bytes = channel_summary(local_root / channel)
        record = {
            "channel": channel,
            "newest_entry_utc": iso_from_ts(newest),
            "age_seconds": None if newest is None else now - newest,
            "entry_count": count,
            "total_bytes": total_bytes,
        }
        record.update(rsync_status.get(channel, {}))
        channel_status.append(record)

    status = {
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "local_metadata_root": str(local_root),
        "channels": channel_status,
    }
    (web_dir / "status.json").write_text(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
