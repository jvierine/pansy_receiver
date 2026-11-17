# Browser plotting runtime

Generated from `jvierine/rustplotlib` commit `1f308be` on branch
`codex/monitor-count-dates`, including compact native axes, shared
scroll-safe navigation, explicit Pagella font loading, Hammer/HEALPix density,
categorical bars and rectilinear fields. Plot toolbars remain visible; grouped
controls share mode/reset actions without imposing gestures on the page.
Runtime release `0331671d4a8ba9a1` is the SHA-256 prefix of the concatenated
generated JavaScript and WASM. Only content-addressed releases use immutable
HTTP caching; HTML, wrappers and telemetry are revalidated.

Build locally with the rustup-managed Rust toolchain:

```
bash scripts/build_monitor_runtime.sh PATH_TO_RUSTPLOTLIB_CHECKOUT
```

These are prebuilt static assets. Never compile the runtime on juha.no.
