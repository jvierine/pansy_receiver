# Browser plotting runtime

Generated from `jvierine/rustplotlib` commit `19b1f32` on branch
`codex/monitor-map-latency`, including the native Mollweide/HEALPix density,
categorical bar and rectilinear field features.
The interaction wrapper includes initialization retry from commit `51c0e7b`.
Runtime release `ca8e45af62e5fd59` is the SHA-256 prefix of the concatenated
generated JavaScript and WASM. Only content-addressed releases use immutable
HTTP caching; HTML, wrappers and telemetry are revalidated.

Build locally with the rustup-managed Rust toolchain:

```
bash scripts/build_monitor_runtime.sh PATH_TO_RUSTPLOTLIB_CHECKOUT
```

These are prebuilt static assets. Never compile the runtime on juha.no.
