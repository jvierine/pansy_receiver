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
install -m 0644 "$repo_dir/pmse/rvd/lib/interactive_plot.js" "$rvd_dir/lib/interactive_plot.js.partial"
mv "$rvd_dir/lib/interactive_plot.js.partial" "$rvd_dir/lib/interactive_plot.js"
install -m 0644 "$repo_dir/web/index.html" "$web_dir/index.html.partial"
mv "$web_dir/index.html.partial" "$web_dir/index.html"
echo "Installed static Rustplotlib monitor; previous HTML retained as a timestamped backup."
