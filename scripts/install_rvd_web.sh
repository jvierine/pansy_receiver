#!/usr/bin/env bash
# Install committed static assets only; never build or reduce radar data here.
set -euo pipefail
task_repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_src="$task_repo/pmse/rvd"
task_dest="${1:-/var/www/rvd-live}"
task_release="$(cat "$task_src/runtime-release.txt")"
mkdir -p "$task_dest/assets" "$task_dest/pkg/releases/$task_release"
for task_file in "$task_src"/assets/* "$task_src/pkg/releases/$task_release/"*; do
  task_rel="${task_file#"$task_src/"}"
  if ! cmp -s "$task_file" "$task_dest/$task_rel"; then
    install -m 0644 "$task_file" "$task_dest/$task_rel.partial"
    mv "$task_dest/$task_rel.partial" "$task_dest/$task_rel"
  fi
done
# Preserve compatibility for open tabs; the new entry point uses one bundle.
for task_rel in viewer.js imaging.js diagnostics.js diagnostic_data.mjs echo_filter.mjs loading.mjs palette.mjs lib/interactive_plot.js lib/interactive_series.js; do
  install -m 0644 "$task_src/$task_rel" "$task_dest/$task_rel.partial"
  mv "$task_dest/$task_rel.partial" "$task_dest/$task_rel"
done
if ! cmp -s "$task_src/apache-static.conf" /etc/apache2/conf-enabled/pansy-rvd-static.conf; then
  sudo -n install -m 0644 "$task_src/apache-static.conf" /etc/apache2/conf-available/pansy-rvd-static.conf
  sudo -n a2enconf pansy-rvd-static
  sudo -n apache2ctl configtest
  sudo -n systemctl reload apache2
fi
cp -p "$task_dest/index.html" "$task_dest/index.html.before-browser-frames-$(date -u +%Y%m%dT%H%M%SZ)"
install -m 0644 "$task_src/index.html" "$task_dest/index.html.partial"
mv "$task_dest/index.html.partial" "$task_dest/index.html"
echo 'Installed RVD browser rendering and immutable cached assets.'
