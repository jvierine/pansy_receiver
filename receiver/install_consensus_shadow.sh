#!/usr/bin/env bash
# Install only the diagnostic timer; never invoke the full receiver installer.
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
test -x "$REPO_DIR/meteor/consensus/target/release/pansy-meteor-consensus"
mkdir -p /media/analysis/shadow/meteor_consensus "$HOME/.config/systemd/user"
for unit in pansy-consensus-shadow.service pansy-consensus-shadow.timer; do
  sed "s#__PANSY_RECEIVER_REPO__#$REPO_DIR#g" "$REPO_DIR/receiver/systemd/$unit" \
    > "$HOME/.config/systemd/user/$unit"
done
systemctl --user daemon-reload
systemctl --user enable --now pansy-consensus-shadow.timer
systemctl --user start pansy-consensus-shadow.service
