#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_DIR="$HOME/.config/pansy-receiver"
SYSTEMD_DIR="$HOME/.config/systemd/user"

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

for unit in pansy-uhd-rx.service pansy-uhd-watchdog.service; do
  sed "s#__PANSY_RECEIVER_REPO__#$REPO_DIR#g" \
    "$REPO_DIR/receiver/systemd/$unit" > "$SYSTEMD_DIR/$unit"
done

chmod +x \
  "$REPO_DIR/receiver/scripts/pansy_uhd_rx_service.sh" \
  "$REPO_DIR/receiver/scripts/pansy_uhd_watchdog.py"

systemctl --user daemon-reload

# Stop the old receiver service and any manually started receiver before the
# new foreground service takes ownership of the USRPs.
systemctl --user stop pansy-uhd-watchdog.service pansy-uhd-rx.service sdrreceiver.service 2>/dev/null || true
pkill -INT -x pansy_uhd_rx 2>/dev/null || true
sleep 5
pkill -TERM -x pansy_uhd_rx 2>/dev/null || true

systemctl --user enable --now pansy-uhd-rx.service
systemctl --user enable --now pansy-uhd-watchdog.service

cat <<EOF
Installed and started the PANSY UHD receiver user services.

Status:
  systemctl --user status pansy-uhd-rx.service
  systemctl --user status pansy-uhd-watchdog.service

Logs:
  tail -f $REPO_DIR/logs/pansy_uhd_rx.service.log
  tail -f $REPO_DIR/logs/pansy_uhd_watchdog.log
  journalctl --user -u pansy-uhd-rx.service -f
  journalctl --user -u pansy-uhd-watchdog.service -f

Configuration:
  $CONFIG_DIR/pansy-receiver.env
EOF
