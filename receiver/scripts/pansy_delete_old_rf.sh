#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_RECEIVER_CONFIG:-$HOME/.config/pansy-receiver/pansy-receiver.env}"
if [ -f "$CONFIG" ]; then
  # shellcheck disable=SC1090
  source "$CONFIG"
fi

ROOT="${PANSY_DELETE_RF_ROOT:-/media/archive}"
MTIME="${PANSY_DELETE_RF_MTIME:-+10}"
INTERVAL="${PANSY_DELETE_RF_INTERVAL_SECONDS:-60}"

echo "Deleting $ROOT/ch*/rf*.h5 older than mtime $MTIME every $INTERVAL seconds"
while true; do
  find "$ROOT"/ch* -type f -name 'rf*.h5' -mtime "$MTIME" -delete 2>&1 || true
  sleep "$INTERVAL"
done
