#!/usr/bin/env bash
set -euo pipefail

CONFIG="${PANSY_RECEIVER_CONFIG:-$HOME/.config/pansy-receiver/pansy-receiver.env}"
if [ -f "$CONFIG" ]; then
  # shellcheck disable=SC1090
  source "$CONFIG"
fi

REPO_DIR="${PANSY_RECEIVER_REPO:-$HOME/src/git/pansy_receiver}"
VENV="${PANSY_PROCESSING_VENV:-/home/radar/testenv}"
EXTRA_LD_LIBRARY_PATH="${PANSY_UHD_RX_LD_LIBRARY_PATH:-}"

cd "$REPO_DIR"
mkdir -p "$REPO_DIR/logs"

if [ ! -d "$VENV" ]; then
  echo "Python environment is missing: $VENV" >&2
  exit 1
fi

# shellcheck disable=SC1091
source "$VENV/bin/activate"

if [ -n "$EXTRA_LD_LIBRARY_PATH" ]; then
  export LD_LIBRARY_PATH="$EXTRA_LD_LIBRARY_PATH${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

export PYTHONUNBUFFERED=1

echo "Starting from $REPO_DIR with $VENV: $*"
exec "$@"
