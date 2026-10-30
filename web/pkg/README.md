# Browser plotting runtime

Generated from `jvierine/rustplotlib` commit `d50c095` on branch
`codex/monitor-map-latency`, including the native Mollweide/HEALPix density,
categorical bar and rectilinear field features.
The interaction wrapper includes initialization retry from commit `51c0e7b`.
Runtime release `1301344d0a86620e` is the SHA-256 prefix of the concatenated
generated JavaScript and WASM. Only content-addressed releases use immutable
HTTP caching; HTML, wrappers and telemetry are revalidated.

Build locally with the rustup-managed Rust toolchain:

```
cargo build -p rustplotlib-web --release --target wasm32-unknown-unknown
wasm-bindgen --target web --out-dir PATH_TO_PANSY/web/pkg target/wasm32-unknown-unknown/release/rustplotlib_web.wasm
esbuild web/src/interactive_series.ts --format=esm --outfile=PATH_TO_PANSY/web/lib/interactive_series.js
```

These are prebuilt static assets. Never compile the runtime on juha.no.
