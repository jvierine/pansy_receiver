#!/usr/bin/env bash
# Run only from a git checkout on the web host. Installs static assets; no jobs.
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
web_dir="${1:-/var/www/pansy-live}"
rvd_dir="${2:-/var/www/rvd-live}"
test -d "$web_dir" && test -d "$rvd_dir/lib"
if [ -f "$web_dir/index.html" ]; then
  cp -p "$web_dir/index.html" "$web_dir/index.html.pre-rustplotlib-$(date -u +%Y%m%dT%H%M%SZ)"
fi
for file in monitor.js monitor_data.mjs monitor.css; do
  install -m 0644 "$repo_dir/web/$file" "$web_dir/$file.partial"
  mv "$web_dir/$file.partial" "$web_dir/$file"
done
runtime_release=ca8e45af62e5fd59
mkdir -p "$web_dir/lib" "$web_dir/pkg/releases/$runtime_release"
if ! cmp -s "$repo_dir/web/lib/interactive_series.js" "$web_dir/lib/interactive_series.js"; then
  install -m 0644 "$repo_dir/web/lib/interactive_series.js" "$web_dir/lib/interactive_series.js.partial"
  mv "$web_dir/lib/interactive_series.js.partial" "$web_dir/lib/interactive_series.js"
fi
for file in rustplotlib_web.js rustplotlib_web_bg.wasm; do
  if ! cmp -s "$repo_dir/web/pkg/releases/$runtime_release/$file" "$web_dir/pkg/releases/$runtime_release/$file"; then
    install -m 0644 "$repo_dir/web/pkg/releases/$runtime_release/$file" "$web_dir/pkg/releases/$runtime_release/$file.partial"
    mv "$web_dir/pkg/releases/$runtime_release/$file.partial" "$web_dir/pkg/releases/$runtime_release/$file"
  fi
done
install -m 0644 "$repo_dir/web/index.html" "$web_dir/index.html.partial"
mv "$web_dir/index.html.partial" "$web_dir/index.html"
echo "Installed static Rustplotlib monitor; previous HTML retained as a timestamped backup."
if command -v sudo >/dev/null && sudo -n true 2>/dev/null && ! cmp -s "$repo_dir/web/pansy-monitor-compression.conf" /etc/apache2/conf-enabled/pansy-monitor-compression.conf; then
  sudo install -m 0644 "$repo_dir/web/pansy-monitor-compression.conf" /etc/apache2/conf-available/pansy-monitor-compression.conf
  sudo a2enconf pansy-monitor-compression
  sudo apache2ctl configtest
  sudo systemctl reload apache2
fi
