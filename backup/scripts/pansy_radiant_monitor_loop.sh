#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_BACKUP_CONFIG:-$HOME/.config/pansy-backup/pansy-backup.env}"
source "$CONFIG"

REPO_DIR="${PANSY_RECEIVER_REPO:-$HOME/src/pansy_receiver}"
INTERVAL_SECONDS="${RADIANT_MONITOR_SLEEP_SECONDS:-1800}"
RSYNC_BWLIMIT="${RADIANT_MONITOR_RSYNC_BWLIMIT_KB:-1000}"

cd "$REPO_DIR"

exec /usr/bin/python3 "$REPO_DIR/publish_radiant_monitor.py" \
  --loop \
  --interval-s "$INTERVAL_SECONDS" \
  --rsync-bwlimit "$RSYNC_BWLIMIT"
