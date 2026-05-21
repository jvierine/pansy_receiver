#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_RECEIVER_CONFIG:-$HOME/.config/pansy-receiver/pansy-receiver.env}"
if [ -f "$CONFIG" ]; then
  # shellcheck disable=SC1090
  source "$CONFIG"
fi

REPO_DIR="${PANSY_RECEIVER_REPO:-$HOME/src/git/pansy_receiver}"
RX_BINARY="${PANSY_UHD_RX_BINARY:-$REPO_DIR/pansy_uhd_rx}"
OUTDIR="${PANSY_UHD_RX_OUTDIR:-/media/archive}"
EXTRA_ARGS="${PANSY_UHD_RX_ARGS:-}"

cd "$REPO_DIR"
mkdir -p "$REPO_DIR/logs" "$OUTDIR"

if [ ! -x "$RX_BINARY" ]; then
  echo "Receiver binary is missing or not executable: $RX_BINARY" >&2
  exit 1
fi

echo "Starting $RX_BINARY --outdir $OUTDIR $EXTRA_ARGS"
# EXTRA_ARGS is intentionally shell-split so the site config can add normal
# command-line switches without editing the service file.
# shellcheck disable=SC2086
exec "$RX_BINARY" --outdir "$OUTDIR" $EXTRA_ARGS
