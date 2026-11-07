#!/usr/bin/env bash
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DROPIN="$HOME/.config/systemd/user/pansy-cluster-mf.service.d/90-rust.conf"
case "${1:-}" in
  enable)
    test -x "$REPO_DIR/meteor/consensus/target/release/libpansy_meteor_consensus.so"
    test -x "$REPO_DIR/meteor/matched_filter/target/release/libpansy_meteor_matched_filter_cut.so"
    test ! -e "$DROPIN"
    mkdir -p "$(dirname "$DROPIN")"
    sed "s#__PANSY_RECEIVER_REPO__#$REPO_DIR#g" \
      "$REPO_DIR/receiver/systemd/pansy-grouping-rust.conf" > "$DROPIN"
    ;;
  rollback) rm -f "$DROPIN" ;;
  *) echo "usage: $0 enable|rollback" >&2; exit 2 ;;
esac
systemctl --user daemon-reload
systemctl --user restart pansy-cluster-mf.service
