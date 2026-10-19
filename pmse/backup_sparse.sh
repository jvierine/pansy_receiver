#!/usr/bin/env bash
# Dedicated archive replication. Retain remote deletions; never remove history.
set -uo pipefail
destination="${PANSY_SPARSE_BACKUP_ROOT:-/mnt/data/juha/pansy/metadata/pmse_sparse}"
bandwidth="${PANSY_SPARSE_BACKUP_KIB:-10}"
mkdir -p "$destination"
while true; do
  rsync -rt --no-perms --no-owner --no-group --partial --delay-updates \
    --timeout=120 --bwlimit="$bandwidth" --exclude='/.writer.lock' \
    -e 'ssh -p 2222 -o BatchMode=yes -o ConnectTimeout=15' \
    radar@127.0.0.1:/media/analysis/metadata/pmse_sparse/ "$destination/"
  result=$?
  echo "Sparse PMSE replication finished with exit status $result"
  sleep 300
done
