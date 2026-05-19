#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_BACKUP_CONFIG:-$HOME/.config/pansy-backup/pansy-backup.env}"
source "$CONFIG"

REPO_DIR="${PANSY_RECEIVER_REPO:-$HOME/src/pansy_receiver}"

while true; do
  "$REPO_DIR/backup/scripts/pansy_deploy_web.sh"
  sleep "$WEB_DEPLOY_SLEEP_SECONDS"
done
