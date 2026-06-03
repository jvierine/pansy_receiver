#!/usr/bin/env bash
set -euo pipefail

INTERVAL="${PANSY_WAIT_INTERVAL_SECONDS:-60}"
paths=()

while [ "$#" -gt 0 ]; do
  case "$1" in
    --)
      shift
      break
      ;;
    *)
      paths+=("$1")
      shift
      ;;
  esac
done

if [ "${#paths[@]}" -eq 0 ] || [ "$#" -eq 0 ]; then
  echo "Usage: pansy_wait_for_paths.sh PATH [PATH ...] -- COMMAND [ARG ...]" >&2
  exit 2
fi

while true; do
  missing=()
  for path in "${paths[@]}"; do
    if [ ! -e "$path" ]; then
      missing+=("$path")
    fi
  done

  if [ "${#missing[@]}" -eq 0 ]; then
    break
  fi

  echo "Waiting for required metadata path(s): ${missing[*]}"
  sleep "$INTERVAL"
done

exec "$@"
