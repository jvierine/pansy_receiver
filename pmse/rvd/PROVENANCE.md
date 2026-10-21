Viewer renderer: https://github.com/jvierine/rustplotlib. MIT license; see
rustplotlib repository LICENSE.

- Checked-in Rust/WASM package and interactive series API: branch
  `codex/rvd-colorbar-endpoints`, commit `42971b5`.
- Shared interactive controls: branch `codex/pansy-wave-theme`, commit `1db580d`.
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

Local artifact adjustments are the missing-observation RGBA color
`(0.6196, 0.6745, 0.7098, 1)` and the immutable runtime import path.
`bash scripts/package_pmse_runtime.sh` copies the generated JS/WASM pair into
`pkg/releases/<sha256-of-JS-and-WASM>/`, without changing either file. Pin the
printed release in both wrappers, `viewer.js` and the HTML preloads. Publish
the complete release directory before the entry page, and retain older public
release directories. Never replace files inside an existing release directory.
The generated JS resolves its WASM sibling within the same release, preventing
mixed cached builds such as non-callable `__wbg_*` WebAssembly imports.
The colorbar click targets, responsive editor, validation, live adjustment and
reset are standard shared Rustplotlib behavior, not page-specific handlers.
Offscreen panels release their WebGL contexts and retain their arrays in CPU
memory. On visibility or GPU-context recovery, the shared controls restore the
same scalar/HSV data, color limits and viewport. This avoids silently blank
RTIs when a browser with a small GPU-context budget opens the other panels.
`scripts/test_pmse_gpu_lifecycle.mjs` exercises scrolling and forced context
loss with an eight-context limit and compares restored scalar/HSV pixels.
`InteractiveSeriesPlot` in `web/src/interactive_series.ts` uses Rust/WASM
`SeriesView.render_plot` for axes, automatic ticks, UTC formatting, legends,
line/marker geometry and white backgrounds. Its wheel/box zoom, Shift-pan,
double-click/Home/Escape reset are shared library defaults. Auxiliary plots
receive numeric JSON arrays, not SVG plot images. Rebuild its checked-in
wrapper with esbuild, leaving `../pkg/rustplotlib_web.js` external.
`npm run test:pmse-colorbars` exercises real WebGL pixels, native frame updates,
HSV preservation, keyboard access and mobile layout against read-only APIs.

Radar inversion is a separate first-party WASM package in `imaging_pkg`, built
from `pmse/browser` and `pmse/processor` with the locked dependencies and
`bash scripts/build_pmse_imaging.sh`. Source algorithms, objective functions,
noise assumptions, experimental limits and references are documented in
`pmse/processor/IMAGING.md`; `npm run test:pmse-imaging` exercises this package.
