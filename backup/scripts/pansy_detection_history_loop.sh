#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_BACKUP_CONFIG:-$HOME/.config/pansy-backup/pansy-backup.env}"
source "$CONFIG"

REPO_DIR="${PANSY_RECEIVER_REPO:-$HOME/src/pansy_receiver}"
SLEEP_SECONDS="${DETECTION_HISTORY_SLEEP_SECONDS:-1800}"

while true; do
  "$REPO_DIR/backup/scripts/pansy_detection_history.py"
  sleep "$SLEEP_SECONDS"
done
