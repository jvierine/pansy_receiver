#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
cargo_bin="${PANSY_CARGO:-cargo}"
bindgen_bin="${PANSY_WASM_BINDGEN:-wasm-bindgen}"
if [[ "$("$bindgen_bin" --version)" != "wasm-bindgen 0.2.129" ]]; then
  echo "This package requires wasm-bindgen-cli 0.2.129." >&2
  exit 1
fi
export CARGO_TARGET_DIR="$PWD/pmse/processor/target/browser"
"$cargo_bin" build --locked --release --target wasm32-unknown-unknown --manifest-path pmse/browser/Cargo.toml
"$bindgen_bin" --target web --out-dir pmse/rvd/imaging_pkg \
  "$CARGO_TARGET_DIR/wasm32-unknown-unknown/release/pansy_browser_imaging.wasm"
