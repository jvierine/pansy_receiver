#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_DIR="$HOME/.config/pansy-receiver"
SYSTEMD_DIR="$HOME/.config/systemd/user"
LOGROTATE_TEMPLATE="$REPO_DIR/receiver/logrotate/pansy-receiver"
LOGROTATE_CONF="/etc/logrotate.d/pansy-receiver"

UNITS=(
  pansy-uhd-rx.service
  pansy-uhd-watchdog.service
  pansy-delete-old-rf.service
  pansy-find-mode-starts.service
  pansy-status-plot.service
  pansy-mesomode-boundary.service
  pansy-cluster-mf.service
  pansy-quick-search-meteor.service
  pansy-meso-xc.service
  pansy-tx-xphase.service
  pansy-cut-meteors.service
  pansy-process-cut-meteor.service
)

mkdir -p "$CONFIG_DIR" "$SYSTEMD_DIR" "$REPO_DIR/logs"

if [ ! -f "$CONFIG_DIR/pansy-receiver.env" ]; then
  cp "$REPO_DIR/receiver/pansy-receiver.env" "$CONFIG_DIR/pansy-receiver.env"
else
  while IFS= read -r line; do
    if [[ "$line" =~ ^([A-Za-z_][A-Za-z0-9_]*)= ]]; then
      key="${BASH_REMATCH[1]}"
      if ! grep -q "^${key}=" "$CONFIG_DIR/pansy-receiver.env"; then
        printf '%s\n' "$line" >> "$CONFIG_DIR/pansy-receiver.env"
      fi
    fi
  done < "$REPO_DIR/receiver/pansy-receiver.env"
fi

for unit in "${UNITS[@]}"; do
  sed "s#__PANSY_RECEIVER_REPO__#$REPO_DIR#g" \
    "$REPO_DIR/receiver/systemd/$unit" > "$SYSTEMD_DIR/$unit"
done

chmod +x \
  "$REPO_DIR/receiver/scripts/pansy_uhd_rx_service.sh" \
  "$REPO_DIR/receiver/scripts/pansy_uhd_watchdog.py" \
  "$REPO_DIR/receiver/scripts/pansy_processing_service.sh" \
  "$REPO_DIR/receiver/scripts/pansy_delete_old_rf.sh"

if command -v sudo >/dev/null 2>&1; then
  tmp_logrotate="$(mktemp)"
  sed "s#__PANSY_RECEIVER_REPO__#$REPO_DIR#g" "$LOGROTATE_TEMPLATE" > "$tmp_logrotate"
  sudo install -m 0644 "$tmp_logrotate" "$LOGROTATE_CONF"
  rm -f "$tmp_logrotate"
  sudo loginctl enable-linger "$USER" >/dev/null 2>&1 || true
fi

systemctl --user daemon-reload

# Stop existing services and any manually started legacy processes before
# systemd takes ownership of the receiver and processing chain.
systemctl --user stop "${UNITS[@]}" sdrreceiver.service 2>/dev/null || true
pkill -INT -x pansy_uhd_rx 2>/dev/null || true
pkill -TERM -f 'python3 find_mode_starts.py' 2>/dev/null || true
pkill -TERM -f 'python3 status_plot.py' 2>/dev/null || true
pkill -TERM -f 'python3 mesomode_boundary.py' 2>/dev/null || true
pkill -TERM -f 'python3 cluster_mf.py' 2>/dev/null || true
pkill -TERM -f 'quick_search_meteor.py' 2>/dev/null || true
pkill -TERM -f 'python3 meso_xc.py' 2>/dev/null || true
pkill -TERM -f 'python3 tx_xphase.py' 2>/dev/null || true
pkill -TERM -f 'python3 cut_meteors.py' 2>/dev/null || true
pkill -TERM -f 'process_cut_meteor.py' 2>/dev/null || true
pkill -TERM -f 'find /media/archive/ch.*rf.*h5.*-mtime' 2>/dev/null || true
sleep 5
pkill -TERM -x pansy_uhd_rx 2>/dev/null || true

systemctl --user enable --now "${UNITS[@]}"

cat <<EOF
Installed and started the PANSY receiver user services.

Status:
  systemctl --user --no-pager --type=service --state=running | grep pansy
  systemctl --user status pansy-uhd-rx.service
  systemctl --user status pansy-quick-search-meteor.service

Logs:
  ls -lh $REPO_DIR/logs
  tail -f $REPO_DIR/logs/pansy_uhd_rx.service.log
  tail -f $REPO_DIR/logs/pansy_quick_search_meteor.log

Configuration:
  $CONFIG_DIR/pansy-receiver.env

Log rotation:
  $LOGROTATE_CONF
EOF
