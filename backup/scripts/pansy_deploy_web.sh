#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_BACKUP_CONFIG:-$HOME/.config/pansy-backup/pansy-backup.env}"
source "$CONFIG"

REPO_DIR="${PANSY_RECEIVER_REPO:-$HOME/src/pansy_receiver}"

mkdir -p "$WEB_BUILD_DIR"
cp "$REPO_DIR/web/index.html" "$WEB_BUILD_DIR/index.html"
"$REPO_DIR/backup/scripts/pansy_backup_status.py"

rsync -avz \
  "$WEB_BUILD_DIR/index.html" \
  "$WEB_BUILD_DIR/status.json" \
  "$WEB_BUILD_DIR/status.svg" \
  "$WEB_HOST:$WEB_ROOT/"
