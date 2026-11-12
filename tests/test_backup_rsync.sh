#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/.." && pwd)
fixture=$(mktemp -d)
trap 'rm -rf "$fixture"' EXIT
mkdir -p "$fixture/bin" "$fixture/state"

cat >"$fixture/bin/rsync" <<'MOCK'
#!/usr/bin/env bash
set -eu
channel=""
budget=""
recent=0
for arg in "$@"; do
  case "$arg" in
    --bwlimit=*) budget=${arg#*=} ;;
    radar@*:*) channel=${arg%/}; channel=${channel##*/} ;;
    --include=*) touch "$FIXTURE/recent.filters"; recent=1 ;;
  esac
done
printf '%s\n' "$budget" >"$FIXTURE/$channel.budget"
if [ "$channel" = cut ]; then
  touch "$FIXTURE/cut.started"
  if [ "${TEST_RECENT_MINUTES:-0}" -gt 0 ] && [ "$recent" = 1 ]; then
    printf 'recent\n' >>"$FIXTURE/cut.recent"
    exit 0
  fi
  sleep 30
else
  if [ ! -f "$FIXTURE/$channel.failed" ]; then
    touch "$FIXTURE/$channel.failed"
    exit 255
  fi
  touch "$FIXTURE/$channel.retried"
fi
MOCK
chmod +x "$fixture/bin/rsync"
cat >"$fixture/config" <<CONFIG
LOCAL_METADATA_ROOT="$fixture/data"
STATE_DIR="$fixture/state"
METADATA_CHANNELS="cut tx phase"
RSYNC_BWLIMIT_KB=${TEST_BUDGET:-50}
RSYNC_SLEEP_SECONDS=0.2
RSYNC_RECENT_MINUTES=${TEST_RECENT_MINUTES:-0}
RSYNC_HISTORY_SLICE_SECONDS=1
RSYNC_DIRECT_SSH=${TEST_DIRECT:-0}
JUMP_HOST=user@example.org
PANSY_SSH_HOST=radar@localhost
PANSY_SSH_PORT=3131
LOCAL_TUNNEL_PORT=2222
REMOTE_METADATA_ROOT=/source
CONFIG

export FIXTURE="$fixture"
export TEST_RECENT_MINUTES=${TEST_RECENT_MINUTES:-0}
set +e
PATH="$fixture/bin:$PATH" PANSY_BACKUP_CONFIG="$fixture/config" \
  timeout --kill-after=2 5 bash "$repo/backup/scripts/pansy_backup_rsync.sh" \
  >"$fixture/output" 2>&1
rc=$?
set -e
if [ "$rc" -ne 124 ] && [ "$rc" -ne 137 ]; then
  cat "$fixture/output"
  exit 1
fi
test -f "$fixture/cut.started"
if [ "$TEST_RECENT_MINUTES" -gt 0 ]; then
  test -f "$fixture/recent.filters"
  test "$(wc -l <"$fixture/cut.recent")" -ge 2
fi
test -f "$fixture/tx.retried"
test -f "$fixture/phase.retried"
awk -F '\t' '$1 == "tx" && $4 == 0 {ok=1} END {exit !ok}' "$fixture/state/rsync_status.tsv"
awk -F '\t' '$1 == "phase" && $4 == 0 {ok=1} END {exit !ok}' "$fixture/state/rsync_status.tsv"
cut_budget=$(cat "$fixture/cut.budget")
tx_budget=$(cat "$fixture/tx.budget")
phase_budget=$(cat "$fixture/phase.budget")
if [ "${TEST_BUDGET:-50}" -eq 0 ]; then
  test "$cut_budget" -eq 0
  test "$tx_budget" -eq 0
  test "$phase_budget" -eq 0
else
  test "$cut_budget" -gt "$tx_budget"
  test "$tx_budget" -gt "$phase_budget"
  test "$((cut_budget + tx_budget + phase_budget))" -le 50
fi
echo 'PASS: failed channels retry and publish success while cut remains blocked; shared budget respected'
