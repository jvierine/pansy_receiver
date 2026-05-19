#!/usr/bin/env python3
import datetime as dt
import html
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


def write_status_svg(path, channels):
    width = 980
    row_h = 34
    height = 80 + row_h * max(1, len(channels))
    now = dt.datetime.now(dt.timezone.utc)
    rows = []
    for i, channel in enumerate(channels):
        y = 60 + i * row_h
        age = channel.get("age_seconds")
        rc = channel.get("last_rsync_exit_code")
        fresh = age is not None and age < 3600 and rc == 0
        color = "#1b9e77" if fresh else "#d95f02"
        newest = channel.get("newest_entry_utc") or "no local entries"
        age_txt = "n/a" if age is None else f"{age / 60:.1f} min"
        rows.append(
            f'<circle cx="24" cy="{y - 5}" r="7" fill="{color}" />'
            f'<text x="44" y="{y}" class="mono">{html.escape(channel["channel"])}</text>'
            f'<text x="230" y="{y}" class="mono">{html.escape(newest)}</text>'
            f'<text x="560" y="{y}" class="mono">{html.escape(age_txt)}</text>'
            f'<text x="690" y="{y}" class="mono">{html.escape(str(rc))}</text>'
            f'<text x="780" y="{y}" class="mono">{channel["entry_count"]}</text>'
        )
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>
  .title {{ font: 700 20px sans-serif; fill: #111; }}
  .label {{ font: 700 13px sans-serif; fill: #444; }}
  .mono {{ font: 13px monospace; fill: #111; }}
</style>
<rect width="100%" height="100%" fill="#fff" />
<text x="20" y="28" class="title">PANSY metadata backup status</text>
<text x="20" y="48" class="label">generated {html.escape(now.isoformat())}</text>
<text x="44" y="54" class="label">channel</text>
<text x="230" y="54" class="label">newest local entry UTC</text>
<text x="560" y="54" class="label">age</text>
<text x="690" y="54" class="label">rsync rc</text>
<text x="780" y="54" class="label">entries</text>
{''.join(rows)}
</svg>
"""
    path.write_text(svg)


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
    write_status_svg(web_dir / "status.svg", channel_status)


if __name__ == "__main__":
    main()
