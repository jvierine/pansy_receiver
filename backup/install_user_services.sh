#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

mkdir -p "$HOME/.config/pansy-backup" "$HOME/.config/systemd/user"

if [ ! -f "$HOME/.config/pansy-backup/pansy-backup.env" ]; then
  cp "$REPO_DIR/backup/pansy-backup.env" "$HOME/.config/pansy-backup/pansy-backup.env"
fi

for unit in \
  pansy-backup-tunnel.service \
  pansy-backup-rsync.service \
  pansy-local-metadata-mirror.service \
  pansy-backup-web.service \
  pansy-detection-history.service \
  pansy-phase-history.service \
  pansy-radiant-monitor.service; do
  sed "s#__PANSY_RECEIVER_REPO__#$REPO_DIR#g" \
    "$REPO_DIR/backup/systemd/$unit" > "$HOME/.config/systemd/user/$unit"
done

chmod +x "$REPO_DIR"/backup/scripts/*.sh "$REPO_DIR"/backup/scripts/*.py

systemctl --user daemon-reload
systemctl --user enable --now pansy-backup-tunnel.service
systemctl --user enable --now pansy-backup-rsync.service
systemctl --user enable --now pansy-local-metadata-mirror.service
systemctl --user enable --now pansy-detection-history.service
systemctl --user enable --now pansy-phase-history.service
systemctl --user enable --now pansy-radiant-monitor.service
systemctl --user enable --now pansy-backup-web.service

cat <<EOF
Installed and started PANSY backup user services.

Check status with:
  systemctl --user status pansy-backup-tunnel.service
  systemctl --user status pansy-backup-rsync.service
  systemctl --user status pansy-local-metadata-mirror.service
  systemctl --user status pansy-detection-history.service
  systemctl --user status pansy-phase-history.service
  systemctl --user status pansy-radiant-monitor.service
  systemctl --user status pansy-backup-web.service

For boot-time user services on revontuli, enable lingering once:
  sudo loginctl enable-linger "$USER"
EOF
