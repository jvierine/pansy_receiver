#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_BACKUP_CONFIG:-$HOME/.config/pansy-backup/pansy-backup.env}"
source "$CONFIG"

mkdir -p "$LOCAL_METADATA_ROOT" "$STATE_DIR"

LOG_FILE="$STATE_DIR/rsync.log"
STATUS_FILE="$STATE_DIR/rsync_status.tsv"

while true; do
  started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  tmp_status="$(mktemp)"

  for channel in $METADATA_CHANNELS; do
    mkdir -p "$LOCAL_METADATA_ROOT/$channel"
    echo "[$started] syncing $channel" | tee -a "$LOG_FILE"

    set +e
    rsync -avz \
      --bwlimit="$RSYNC_BWLIMIT_KB" \
      --partial \
      --timeout=120 \
      -e "ssh -p ${LOCAL_TUNNEL_PORT} -o BatchMode=yes -o ConnectTimeout=30 -o StrictHostKeyChecking=accept-new" \
      "radar@127.0.0.1:${REMOTE_METADATA_ROOT}/${channel}/20*" \
      "$LOCAL_METADATA_ROOT/$channel/" >>"$LOG_FILE" 2>&1
    rc=$?
    set -e

    finished="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf "%s\t%s\t%s\t%s\n" "$channel" "$started" "$finished" "$rc" >>"$tmp_status"

    if [ "$rc" -ne 0 ]; then
      echo "[$finished] rsync for $channel exited with status $rc" | tee -a "$LOG_FILE"
    fi
  done

  mv "$tmp_status" "$STATUS_FILE"
  sleep "$RSYNC_SLEEP_SECONDS"
done
