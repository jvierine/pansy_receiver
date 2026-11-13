#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_BACKUP_CONFIG:-$HOME/.config/pansy-backup/pansy-backup.env}"
source "$CONFIG"
source "$(dirname "${BASH_SOURCE[0]}")/pansy_storage_guard.sh"

LOCAL_MIRROR_METADATA_ROOT="${LOCAL_MIRROR_METADATA_ROOT:-/urdr/data/juha/pansy/metadata}"
LOCAL_MIRROR_SLEEP_SECONDS="${LOCAL_MIRROR_SLEEP_SECONDS:-300}"
LOCAL_MIRROR_BWLIMIT_KB="${LOCAL_MIRROR_BWLIMIT_KB:-0}"

for path in "${LOCAL_METADATA_ROOT:-}" "$LOCAL_MIRROR_METADATA_ROOT" "${STATE_DIR:-}"; do
  if [[ "$path" != /* ]] || [ "$path" = / ]; then
    echo 'Mirror and state paths must be absolute, non-root directories' >&2
    exit 1
  fi
done
source_real=$(realpath -m -- "$LOCAL_METADATA_ROOT")
mirror_real=$(realpath -m -- "$LOCAL_MIRROR_METADATA_ROOT")
case "$mirror_real/" in "$source_real/"*) echo 'Mirror must not overlap its source' >&2; exit 1 ;; esac
case "$source_real/" in "$mirror_real/"*) echo 'Mirror must not overlap its source' >&2; exit 1 ;; esac
pansy_require_mount "$LOCAL_METADATA_ROOT" "${LOCAL_METADATA_MOUNT:-}"
pansy_require_mount "$LOCAL_MIRROR_METADATA_ROOT" "${LOCAL_MIRROR_METADATA_MOUNT:-}"
mkdir -p "$LOCAL_METADATA_ROOT" "$LOCAL_MIRROR_METADATA_ROOT" "$STATE_DIR"
exec {instance_fd}>"$STATE_DIR/local_mirror.instance.lock"
if ! flock -n "$instance_fd"; then
  echo 'A local metadata mirror is already running' >&2
  exit 1
fi

LOG_FILE="$STATE_DIR/local_mirror.log"
STATUS_FILE="$STATE_DIR/local_mirror_status.tsv"
ACTIVE_FILE="$STATE_DIR/local_mirror_active.tsv"

cleanup() {
  rm -f "$ACTIVE_FILE"
}

trap cleanup EXIT
trap 'exit 0' INT TERM

while true; do
  pansy_require_mount "$LOCAL_METADATA_ROOT" "${LOCAL_METADATA_MOUNT:-}"
  pansy_require_mount "$LOCAL_MIRROR_METADATA_ROOT" "${LOCAL_MIRROR_METADATA_MOUNT:-}"
  started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf "%s\t%s\t%s\n" "$LOCAL_METADATA_ROOT" "$LOCAL_MIRROR_METADATA_ROOT" "$started" >"$ACTIVE_FILE"
  echo "[$started] mirroring $LOCAL_METADATA_ROOT/ to $LOCAL_MIRROR_METADATA_ROOT/" | tee -a "$LOG_FILE"

  rsync_args=(
    -a
    --partial
    --partial-dir=.rsync-partial
    --exclude=.rsync-partial/
    --timeout=120
  )

  if [ "$LOCAL_MIRROR_BWLIMIT_KB" != "0" ]; then
    rsync_args+=(--bwlimit="$LOCAL_MIRROR_BWLIMIT_KB")
  fi

  set +e
  rsync "${rsync_args[@]}" \
    "$LOCAL_METADATA_ROOT/" \
    "$LOCAL_MIRROR_METADATA_ROOT/" {instance_fd}>&- >>"$LOG_FILE" 2>&1
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
