# Browser plotting runtime

Generated from `jvierine/rustplotlib` commit `4a1d509` on branch
`codex/monitor-map-latency`, including the native Mollweide/HEALPix density,
categorical bar and rectilinear field features.

Build locally with the rustup-managed Rust toolchain:

```
cargo build -p rustplotlib-web --release --target wasm32-unknown-unknown
wasm-bindgen --target web --out-dir PATH_TO_PANSY/web/pkg target/wasm32-unknown-unknown/release/rustplotlib_web.wasm
esbuild web/src/interactive_series.ts --format=esm --outfile=PATH_TO_PANSY/web/lib/interactive_series.js
```

These are prebuilt static assets. Never compile the runtime on juha.no.
