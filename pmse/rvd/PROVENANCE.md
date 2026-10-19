Viewer renderer: https://github.com/jvierine/rustplotlib, branch codex/pansy-rvd,
commit 74eff2d. MIT license; see rustplotlib repository LICENSE.

Build the rustplotlib-web crate for wasm32-unknown-unknown, then run
wasm-bindgen --target web --out-dir pmse/rvd/pkg on rustplotlib_web.wasm.
The checked-in package is a compiled deployment artifact. JavaScript handles
data transport, controls and axis text; Rust/WASM/WebGL renders the scientific
field and performs HSV color conversion on the GPU. No Matplotlib rendering.
