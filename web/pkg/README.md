# Browser plotting runtime

Generated from `jvierine/rustplotlib` commit `198a030` on branch
`codex/monitor-mobile-controls`, including compact native axes, shared
scroll-safe navigation, explicit Pagella font loading, Hammer/HEALPix density,
categorical bars and rectilinear fields. Plot toolbars remain visible; grouped
controls share mode/reset actions without imposing gestures on the page.
Runtime release `7caedeeb810f57e4` is the SHA-256 prefix of the concatenated
generated JavaScript and WASM. Only content-addressed releases use immutable
HTTP caching; HTML, wrappers and telemetry are revalidated.

Build locally with the rustup-managed Rust toolchain:

```
bash scripts/build_monitor_runtime.sh PATH_TO_RUSTPLOTLIB_CHECKOUT
```

These are prebuilt static assets. Never compile the runtime on juha.no.
