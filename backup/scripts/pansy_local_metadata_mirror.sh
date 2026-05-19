#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_BACKUP_CONFIG:-$HOME/.config/pansy-backup/pansy-backup.env}"
source "$CONFIG"

LOCAL_MIRROR_METADATA_ROOT="${LOCAL_MIRROR_METADATA_ROOT:-/urdr/data/juha/pansy/metadata}"
LOCAL_MIRROR_SLEEP_SECONDS="${LOCAL_MIRROR_SLEEP_SECONDS:-300}"
LOCAL_MIRROR_BWLIMIT_KB="${LOCAL_MIRROR_BWLIMIT_KB:-0}"

mkdir -p "$LOCAL_METADATA_ROOT" "$LOCAL_MIRROR_METADATA_ROOT" "$STATE_DIR"

LOG_FILE="$STATE_DIR/local_mirror.log"
STATUS_FILE="$STATE_DIR/local_mirror_status.tsv"
ACTIVE_FILE="$STATE_DIR/local_mirror_active.tsv"

cleanup() {
  rm -f "$ACTIVE_FILE"
}

trap cleanup EXIT INT TERM

while true; do
  started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf "%s\t%s\t%s\n" "$LOCAL_METADATA_ROOT" "$LOCAL_MIRROR_METADATA_ROOT" "$started" >"$ACTIVE_FILE"
  echo "[$started] mirroring $LOCAL_METADATA_ROOT/ to $LOCAL_MIRROR_METADATA_ROOT/" | tee -a "$LOG_FILE"

  rsync_args=(
    -a
    --partial
    --timeout=120
  )

  if [ "$LOCAL_MIRROR_BWLIMIT_KB" != "0" ]; then
    rsync_args+=(--bwlimit="$LOCAL_MIRROR_BWLIMIT_KB")
  fi

  set +e
  rsync "${rsync_args[@]}" \
    "$LOCAL_METADATA_ROOT/" \
    "$LOCAL_MIRROR_METADATA_ROOT/" >>"$LOG_FILE" 2>&1
  rc=$?
  set -e

  finished="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  rm -f "$ACTIVE_FILE"
  printf "%s\t%s\t%s\t%s\t%s\n" "$LOCAL_METADATA_ROOT" "$LOCAL_MIRROR_METADATA_ROOT" "$started" "$finished" "$rc" >"$STATUS_FILE"

  if [ "$rc" -ne 0 ]; then
    echo "[$finished] local metadata mirror exited with status $rc" | tee -a "$LOG_FILE"
  else
    echo "[$finished] local metadata mirror completed" | tee -a "$LOG_FILE"
  fi

  sleep "$LOCAL_MIRROR_SLEEP_SECONDS"
done
