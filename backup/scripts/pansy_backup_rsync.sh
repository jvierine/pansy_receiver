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
recent_minutes=${RSYNC_RECENT_MINUTES:-0}
history_seconds=${RSYNC_HISTORY_SLICE_SECONDS:-180}
if ! [[ "$recent_minutes" =~ ^[0-9]+$ ]] || ! [[ "$history_seconds" =~ ^[0-9]+$ ]] || [ "$history_seconds" -lt 1 ]; then
  echo 'Invalid recent window or history slice duration' >&2
  exit 1
fi

read -r -a channels <<< "$METADATA_CHANNELS"
total_budget="${RSYNC_TOTAL_BWLIMIT_KB:-${RSYNC_BWLIMIT_KB:-0}}"
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
  if [ "${RSYNC_DIRECT_SSH:-0}" = 1 ]; then
    echo "ssh -J ${JUMP_HOST} -p ${PANSY_SSH_PORT} -o BatchMode=yes -o ConnectTimeout=30 -o ServerAliveInterval=15 -o ServerAliveCountMax=3 -o StrictHostKeyChecking=accept-new -o ControlMaster=auto -o ControlPersist=180 -o ControlPath=${STATE_DIR}/ssh/%C"
    return
  fi
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
  local mode="${2:-all}" label="$1"
  local channel_started
  local finished
  local active_file
  local bwlimit_kb
  local lock_file
  local rc
  local remote_source
  local -a filters=() runner=()
  if [ "$mode" = history ]; then label="$channel:history";runner=(timeout --signal=TERM --kill-after=15 "${history_seconds}s");fi
  if [ "$mode" = recent ]; then
    filters+=(--include=/dmd_properties.h5)
    local minute epoch hour
    epoch=$(date -u +%s)
    for ((minute=1;minute<=recent_minutes;minute++)); do
      local stamp=$((epoch/60*60-minute*60))
      hour=$(date -u -d "@$stamp" +%Y-%m-%dT%H-00-00)
      filters+=("--include=/$hour/" "--include=/$hour/*@$stamp.h5")
    done
    filters+=(--exclude=\*)
  fi

  mkdir -p "$LOCAL_METADATA_ROOT/$channel"
  mkdir -p "$LOCK_DIR"
  channel_started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  active_file="$ACTIVE_DIR/${label//:/.}.tsv"
  lock_file="$LOCK_DIR/$channel.lock"
  bwlimit_kb="$(channel_bwlimit_kb "$channel")"

  exec {lock_fd}>"$lock_file"
  if ! flock -n "$lock_fd"; then
    echo "[$channel_started] skipping $channel because another rsync is already active" | tee -a "$LOG_FILE"
    exec {lock_fd}>&-
    return 0
  fi

  printf "%s\t%s\t%s\n" "$label" "$channel_started" "$BASHPID" >"$active_file"
  echo "[$channel_started] syncing $label ($mode) with --bwlimit=$bwlimit_kb" | tee -a "$LOG_FILE"

  remote_source="radar@127.0.0.1:${REMOTE_METADATA_ROOT}/${channel}/"
  if [ "${RSYNC_DIRECT_SSH:-0}" = 1 ]; then remote_source="${PANSY_SSH_HOST}:${REMOTE_METADATA_ROOT}/${channel}/";fi

  set +e
  "${runner[@]}" rsync -avz \
    --bwlimit="$bwlimit_kb" \
    --partial \
    --partial-dir=.rsync-partial \
    --timeout=120 \
    "${filters[@]}" \
    -e "$(rsync_ssh_cmd)" \
    "$remote_source" \
    "$LOCAL_METADATA_ROOT/$channel/" >>"$LOG_FILE" 2>&1
  rc=$?
  set -e

  finished="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  rm -f "$active_file"
  flock -u "$lock_fd"
  exec {lock_fd}>&-
  if [ "$mode" = history ] && [ "$rc" = 124 ]; then
    echo "[$finished] $channel history slice paused; refreshing live data next" | tee -a "$LOG_FILE"
    rc=paused
  fi
  write_channel_status "$label" "$channel_started" "$finished" "$rc"

  if [ "$rc" != 0 ] && [ "$rc" != paused ]; then
    echo "[$finished] rsync for $channel exited with status $rc" | tee -a "$LOG_FILE"
  fi

  return 0
}

trap cleanup EXIT
trap 'exit 0' INT TERM

mkdir -p "$STATUS_DIR"
mkdir -p "$STATE_DIR/ssh"
chmod 700 "$STATE_DIR/ssh"
clear_active
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] total transfer budget $total_budget KiB/s; independent channel retries" | tee -a "$LOG_FILE"

for channel in "${channels[@]}"; do
  (
    trap 'exit 0' INT TERM
    while true; do
      if [ "$recent_minutes" -gt 0 ]; then
        sync_channel "$channel" recent
        sync_channel "$channel" history
      else
        sync_channel "$channel"
      fi
      sleep "$RSYNC_SLEEP_SECONDS"
    done
  ) &
  pids+=("$!")
done

# If a worker unexpectedly dies, let systemd restart the complete worker set.
wait -n "${pids[@]}"
exit 1
