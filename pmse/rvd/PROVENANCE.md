Viewer renderer: https://github.com/jvierine/rustplotlib. MIT license; see
rustplotlib repository LICENSE.

- Checked-in Rust/WASM/WebGL package: branch `codex/pansy-rvd`, commit `74eff2d`.
- Shared interactive controls: branch `codex/pansy-wave-theme`, commit `baa494c`.
- Native axes/colorbars: pinned separately in `pmse/plotting/Cargo.toml`.

Build the rustplotlib-web crate for wasm32-unknown-unknown, then run
wasm-bindgen --target web --out-dir pmse/rvd/pkg on rustplotlib_web.wasm.
The checked-in package is a compiled deployment artifact. JavaScript handles
data transport and the shared controls; native Rust typesets scientific axes,
labels and colorbars, while Rust/WASM/WebGL renders the scientific field and
performs HSV color conversion on the GPU. No Matplotlib rendering.

Build the shared controls from the Rustplotlib checkout with:

```
npm exec --yes --package=esbuild@0.25.0 -- esbuild web/src/interactive_plot.ts --format=esm --target=es2022 --outfile=/path/to/pansy_receiver/pmse/rvd/lib/interactive_plot.js
```

The sole local artifact adjustment is the missing-observation RGBA color
`(0.6196, 0.6745, 0.7098, 1)` to match the page's grey missing-data key.
The colorbar click targets, responsive editor, validation, live adjustment and
reset are standard shared Rustplotlib behavior, not page-specific handlers.
`npm run test:pmse-colorbars` exercises real WebGL pixels, native frame updates,
HSV preservation, keyboard access and mobile layout against read-only APIs.

Radar inversion is a separate first-party WASM package in `imaging_pkg`, built
from `pmse/browser` and `pmse/processor` with the locked dependencies and
`bash scripts/build_pmse_imaging.sh`. Source algorithms, objective functions,
noise assumptions, experimental limits and references are documented in
`pmse/processor/IMAGING.md`; `npm run test:pmse-imaging` exercises this package.
