#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_BACKUP_CONFIG:-$HOME/.config/pansy-backup/pansy-backup.env}"
source "$CONFIG"

mkdir -p "$LOCAL_METADATA_ROOT" "$STATE_DIR"

LOG_FILE="$STATE_DIR/rsync.log"
STATUS_FILE="$STATE_DIR/rsync_status.tsv"
ACTIVE_DIR="$STATE_DIR/rsync_active"
LOCK_DIR="$STATE_DIR/rsync_locks"
tmp_dir=""
pids=()

channel_bwlimit_kb() {
  local channel="$1"
  if [ "$channel" = "cut" ]; then
    echo "${RSYNC_BWLIMIT_CUT_KB:-50}"
  else
    echo "${RSYNC_BWLIMIT_DEFAULT_KB:-10}"
  fi
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
  if [ -n "$tmp_dir" ]; then
    rm -rf "$tmp_dir"
  fi
}

sync_channel() {
  local channel="$1"
  local channel_status_file="$2"
  local channel_started
  local finished
  local active_file
  local bwlimit_kb
  local lock_file
  local rc

  mkdir -p "$LOCAL_METADATA_ROOT/$channel"
  mkdir -p "$LOCK_DIR"
  channel_started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  active_file="$ACTIVE_DIR/$channel.tsv"
  lock_file="$LOCK_DIR/$channel.lock"
  bwlimit_kb="$(channel_bwlimit_kb "$channel")"

  exec {lock_fd}>"$lock_file"
  if ! flock -n "$lock_fd"; then
    echo "[$channel_started] skipping $channel because another rsync is already active" | tee -a "$LOG_FILE"
    printf "%s\t%s\t%s\t%s\n" "$channel" "$channel_started" "$channel_started" "locked" >"$channel_status_file"
    return 0
  fi

  printf "%s\t%s\t%s\n" "$channel" "$channel_started" "$BASHPID" >"$active_file"
  echo "[$channel_started] syncing $channel with --bwlimit=$bwlimit_kb" | tee -a "$LOG_FILE"

  set +e
  rsync -avz \
    --bwlimit="$bwlimit_kb" \
    --partial \
    --timeout=120 \
    -e "$(rsync_ssh_cmd)" \
    "radar@127.0.0.1:${REMOTE_METADATA_ROOT}/${channel}/20*" \
    "$LOCAL_METADATA_ROOT/$channel/" >>"$LOG_FILE" 2>&1
  rc=$?
  set -e

  finished="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  rm -f "$active_file"
  flock -u "$lock_fd"
  printf "%s\t%s\t%s\t%s\n" "$channel" "$channel_started" "$finished" "$rc" >"$channel_status_file"

  if [ "$rc" -ne 0 ]; then
    echo "[$finished] rsync for $channel exited with status $rc" | tee -a "$LOG_FILE"
  fi

  return 0
}

trap cleanup EXIT INT TERM

while true; do
  tmp_status="$(mktemp)"
  tmp_dir="$(mktemp -d)"
  clear_active
  pids=()

  for channel in $METADATA_CHANNELS; do
    sync_channel "$channel" "$tmp_dir/$channel.tsv" &
    pids+=("$!")
  done

  for pid in "${pids[@]}"; do
    wait "$pid" || true
  done

  cat "$tmp_dir"/*.tsv >"$tmp_status"

  mv "$tmp_status" "$STATUS_FILE"
  rm -rf "$tmp_dir"
  tmp_dir=""
  sleep "$RSYNC_SLEEP_SECONDS"
done
