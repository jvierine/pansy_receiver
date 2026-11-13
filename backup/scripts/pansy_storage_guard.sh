#!/usr/bin/env bash
# Refuse to fall back onto the root filesystem when an archive mount disappears.
pansy_require_mount() {
  local destination="$1" expected="${2:-}" resolved_destination resolved_mount
  [ -n "$expected" ] || return 0
  if [ "$expected" = / ] || ! mountpoint -q -- "$expected"; then
    echo "Required archive mount is unavailable: $expected" >&2
    return 1
  fi
  resolved_destination=$(realpath -m -- "$destination") || return 1
  resolved_mount=$(realpath -e -- "$expected") || return 1
  case "$resolved_destination/" in
    "$resolved_mount/"*) return 0 ;;
    *) echo "Archive destination $destination is not under required mount $expected" >&2; return 1 ;;
  esac
}
