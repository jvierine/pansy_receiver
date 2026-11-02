#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/.." && pwd)
binary="$repo/pmse/publisher/target/release/pansy-rvd-publisher"
test -x "$binary"
state="$HOME/.cache/pansy-rvd-publisher/state.json"
unit="pansy-rvd-publisher-initialize-$(date +%s)"

systemctl --user stop pansy-rvd-web.service
if ! sudo -n systemd-run --unit="$unit" --uid="$(id -un)" --gid="$(id -gn)" \
    --property=CPUQuota=20% --property=MemoryMax=192M \
    --property=RuntimeMaxSec=600 --collect --wait \
    /usr/bin/time -f 'user_s=%U system_s=%S elapsed_s=%e max_rss_kib=%M' \
    "$binary" --once --native-root /mnt/shovel/rvd/pmse_sparse16 \
    --output /var/www/rvd-live/data --state "$state" \
    --python /home/j/venv/pansy-rvd/bin/python \
    --sidecar-script "$repo/pmse/write_gui_sidecar.py"; then
  systemctl --user start pansy-rvd-web.service
  echo "Initialization failed; previous publisher restarted. See journalctl -u $unit" >&2
  exit 1
fi

directory="$HOME/.config/systemd/user/pansy-rvd-web.service.d"
mkdir -p "$directory"
sed "s|__PANSY_RECEIVER_REPO__|$repo|g" \
  "$repo/pmse/systemd/pansy-rvd-web-rust.conf" >"$directory/30-rust-publisher.conf"
systemctl --user daemon-reload
systemctl --user start pansy-rvd-web.service
systemctl --user --no-pager status pansy-rvd-web.service
sudo -n journalctl -u "$unit" -n 12 --no-pager
