#!/usr/bin/env bash
set -euo pipefail
# Only processing services: acquisition is deliberately excluded.
case "${1:-}" in
  mesomode-boundary|quick-search-meteor|cluster-mf|cut-meteors|process-cut-meteor) ;;
  *) echo "Specify one processing stage" >&2; exit 2 ;;
esac
service="pansy-$1.service"
directory="$HOME/.config/systemd/user/$service.d"
mkdir -p "$directory"
printf '[Service]\nTimeoutStopSec=10s\n' > "$directory/80-bounded-stop.conf"
if [ "$1" = quick-search-meteor ]; then
  printf 'Environment=OPENBLAS_NUM_THREADS=1\nEnvironment=OMP_NUM_THREADS=1\n' >> "$directory/80-bounded-stop.conf"
fi
systemctl --user daemon-reload
systemctl --user restart "$service"
systemctl --user is-active "$service"
