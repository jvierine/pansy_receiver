#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_BACKUP_CONFIG:-$HOME/.config/pansy-backup/pansy-backup.env}"
source "$CONFIG"

mkdir -p "$LOCAL_METADATA_ROOT" "$STATE_DIR"

LOG_FILE="$STATE_DIR/rsync.log"
STATUS_FILE="$STATE_DIR/rsync_status.tsv"
ACTIVE_DIR="$STATE_DIR/rsync_active"
LOCK_DIR="$STATE_DIR/rsync_locks"
STATUS_DIR="$STATE_DIR/rsync_status"
pids=()

read -r -a channels <<< "$METADATA_CHANNELS"
total_budget="${RSYNC_TOTAL_BWLIMIT_KB:-${RSYNC_BWLIMIT_KB:-50}}"
if ! [[ "$total_budget" =~ ^[0-9]+$ ]] ||
   [ "${#channels[@]}" -eq 0 ] ||
   { [ "$total_budget" -ne 0 ] && [ "$total_budget" -lt "${#channels[@]}" ]; }; then
  echo "Transfer budget must be 0 (unlimited) or at least the number of metadata channels (KiB/s)" >&2
  exit 1
fi

channel_weight() {
  case "$1" in
    cut) echo 10 ;;
    tx) echo 4 ;;
    *) echo 1 ;;
  esac
}

total_weight=0
for channel in "${channels[@]}"; do
  total_weight=$((total_weight + $(channel_weight "$channel")))
done

channel_bwlimit_kb() {
  local channel="$1"
  local weight
  if [ "$total_budget" -eq 0 ]; then
    echo 0
    return
  fi
  weight="$(channel_weight "$channel")"
  echo "$((1 + (total_budget - ${#channels[@]}) * weight / total_weight))"
}

rsync_ssh_cmd() {
  echo "ssh -p ${LOCAL_TUNNEL_PORT} -o BatchMode=yes -o ConnectTimeout=30 -o StrictHostKeyChecking=accept-new"
}

clear_active() {
  rm -rf "$ACTIVE_DIR"
  mkdir -p "$ACTIVE_DIR" "$LOCK_DIR"
}

cleanup() {
  if [ "${#pids[@]}" -gt 0 ]; then
    pkill -TERM -P "$$" 2>/dev/null || true
    kill "${pids[@]}" 2>/dev/null || true
    sleep 2
    pkill -KILL -P "$$" 2>/dev/null || true
    wait "${pids[@]}" 2>/dev/null || true
  fi
  clear_active
}

write_channel_status() {
  local channel="$1" started="$2" finished="$3" rc="$4"
  local temporary="$STATUS_DIR/.$channel.$BASHPID.tmp"
  printf "%s\t%s\t%s\t%s\n" "$channel" "$started" "$finished" "$rc" >"$temporary"
  mv "$temporary" "$STATUS_DIR/$channel.tsv"
  (
    flock -x 9
    cat "$STATUS_DIR"/*.tsv >"$STATUS_FILE.tmp"
    mv "$STATUS_FILE.tmp" "$STATUS_FILE"
  ) 9>"$STATE_DIR/rsync_status.lock"
}

sync_channel() {
  local channel="$1"
  local channel_started
  local finished
  local active_file
  local bwlimit_kb
  local lock_file
  local rc
  local remote_source

  mkdir -p "$LOCAL_METADATA_ROOT/$channel"
  mkdir -p "$LOCK_DIR"
  channel_started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  active_file="$ACTIVE_DIR/$channel.tsv"
  lock_file="$LOCK_DIR/$channel.lock"
  bwlimit_kb="$(channel_bwlimit_kb "$channel")"

  exec {lock_fd}>"$lock_file"
  if ! flock -n "$lock_fd"; then
    echo "[$channel_started] skipping $channel because another rsync is already active" | tee -a "$LOG_FILE"
    exec {lock_fd}>&-
    return 0
  fi

  printf "%s\t%s\t%s\n" "$channel" "$channel_started" "$BASHPID" >"$active_file"
  echo "[$channel_started] syncing $channel with --bwlimit=$bwlimit_kb" | tee -a "$LOG_FILE"

  remote_source="radar@127.0.0.1:${REMOTE_METADATA_ROOT}/${channel}/"

  set +e
  rsync -avz \
    --bwlimit="$bwlimit_kb" \
    --partial \
    --timeout=120 \
    -e "$(rsync_ssh_cmd)" \
    "$remote_source" \
    "$LOCAL_METADATA_ROOT/$channel/" >>"$LOG_FILE" 2>&1
  rc=$?
  set -e

  finished="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  rm -f "$active_file"
  flock -u "$lock_fd"
  exec {lock_fd}>&-
  write_channel_status "$channel" "$channel_started" "$finished" "$rc"

  if [ "$rc" -ne 0 ]; then
    echo "[$finished] rsync for $channel exited with status $rc" | tee -a "$LOG_FILE"
  fi

  return 0
}

trap cleanup EXIT
trap 'exit 0' INT TERM

mkdir -p "$STATUS_DIR"
clear_active
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] total transfer budget $total_budget KiB/s; independent channel retries" | tee -a "$LOG_FILE"

for channel in "${channels[@]}"; do
  (
    trap 'exit 0' INT TERM
    while true; do
      sync_channel "$channel"
      sleep "$RSYNC_SLEEP_SECONDS"
    done
  ) &
  pids+=("$!")
done

# If a worker unexpectedly dies, let systemd restart the complete worker set.
wait -n "${pids[@]}"
exit 1
