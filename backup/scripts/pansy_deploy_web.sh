#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_BACKUP_CONFIG:-$HOME/.config/pansy-backup/pansy-backup.env}"
source "$CONFIG"

REPO_DIR="${PANSY_RECEIVER_REPO:-$HOME/src/pansy_receiver}"
INCLUDE_STATIC=0
if [ "${1:-}" = "--include-static" ]; then
  INCLUDE_STATIC=1
fi

mkdir -p "$WEB_BUILD_DIR"
"$REPO_DIR/backup/scripts/pansy_backup_status.py"

DEPLOY_FILES=(
  "$WEB_BUILD_DIR/status.json"
)

if [ "$INCLUDE_STATIC" -eq 1 ]; then
  cp \
    "$REPO_DIR/web/index.html" \
    "$REPO_DIR/web/pansy-title.svg" \
    "$REPO_DIR/web/uit-logo.png" \
    "$REPO_DIR/web/iap-logo.svg" \
    "$REPO_DIR/web/u-tokyo-logo.png" \
    "$REPO_DIR/web/favicon.svg" \
    "$REPO_DIR/web/favicon.ico" \
    "$WEB_BUILD_DIR/"
  DEPLOY_FILES=(
    "$WEB_BUILD_DIR/index.html"
    "$WEB_BUILD_DIR/pansy-title.svg"
    "$WEB_BUILD_DIR/uit-logo.png"
    "$WEB_BUILD_DIR/iap-logo.svg"
    "$WEB_BUILD_DIR/u-tokyo-logo.png"
    "$WEB_BUILD_DIR/favicon.svg"
    "$WEB_BUILD_DIR/favicon.ico"
    "${DEPLOY_FILES[@]}"
  )
fi

if [ -f "$WEB_BUILD_DIR/detections_daily.png" ]; then
  DEPLOY_FILES+=("$WEB_BUILD_DIR/detections_daily.png")
fi

if [ -f "$WEB_BUILD_DIR/phase_history.png" ]; then
  DEPLOY_FILES+=("$WEB_BUILD_DIR/phase_history.png")
fi

rsync -avz "${DEPLOY_FILES[@]}" "$WEB_HOST:$WEB_ROOT/"
