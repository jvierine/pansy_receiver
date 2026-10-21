#!/usr/bin/env bash
# Immutable URLs bind the generated JS imports and WASM binary to one build.
set -euo pipefail
cd "$(dirname "$0")/.."
task_package=pmse/rvd/pkg
task_release=$(cat "$task_package/rustplotlib_web.js" "$task_package/rustplotlib_web_bg.wasm" | sha256sum | cut -c1-16)
task_destination="$task_package/releases/$task_release"
mkdir -p "$task_destination"
install -m 644 "$task_package/rustplotlib_web.js" "$task_package/rustplotlib_web_bg.wasm" "$task_destination/"
echo "Packaged runtime: $task_release"
echo "Pin this release in viewer.js, both shared wrappers, and index.html; run tests/test_pmse_runtime_assets.py."
