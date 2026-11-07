#!/usr/bin/env bash
# Call only after station preflight_cut.py passes with the selected library.
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DROPIN="$HOME/.config/systemd/user/pansy-process-cut-meteor.service.d/90-rust.conf"
case "${1:-}" in
  enable)
    test -x "$REPO_DIR/meteor/matched_filter/target/release/libpansy_meteor_matched_filter_cut.so"
    test ! -e "$DROPIN"
    mkdir -p "$(dirname "$DROPIN")"
    sed "s#__PANSY_RECEIVER_REPO__#$REPO_DIR#g" \
      "$REPO_DIR/receiver/systemd/pansy-cut-rust.conf" > "$DROPIN"
    ;;
  rollback)
    rm -f "$DROPIN"
    ;;
  *) echo "usage: $0 enable|rollback" >&2; exit 2 ;;
esac
systemctl --user daemon-reload
systemctl --user restart pansy-process-cut-meteor.service
