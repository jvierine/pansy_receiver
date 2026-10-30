#!/usr/bin/env bash
# Build on a development/analysis machine, never on juha.no.
set -euo pipefail
task_plot_source="$(cd "${1:?Pass a Rustplotlib checkout}" && pwd)"
task_repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_cargo="${PANSY_CARGO_BIN:-$HOME/.cargo/bin/cargo}"
task_bindgen="${PANSY_BINDGEN_BIN:-$HOME/.cargo/bin/wasm-bindgen}"
task_esbuild="${PANSY_ESBUILD_BIN:-$task_plot_source/node_modules/.bin/esbuild}"
test -x "$task_esbuild"
(cd "$task_plot_source" && "$task_cargo" build -p rustplotlib-web --release --target wasm32-unknown-unknown)
"$task_bindgen" --target web --out-dir "$task_repo/web/pkg" "$task_plot_source/target/wasm32-unknown-unknown/release/rustplotlib_web.wasm"
"$task_esbuild" "$task_plot_source/web/src/interactive_series.ts" --format=esm --outfile="$task_repo/web/lib/interactive_series.js"
node --input-type=module - "$task_repo" "$task_plot_source" <<'JS'
import fs from 'node:fs';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {execFileSync} from 'node:child_process';
const [root,plot]=process.argv.slice(2),pkg=path.join(root,'web/pkg');
const files=['rustplotlib_web.js','rustplotlib_web_bg.wasm'];
const hash=createHash('sha256');for(const name of files)hash.update(fs.readFileSync(path.join(pkg,name)));
const release=hash.digest('hex').slice(0,16),destination=path.join(pkg,'releases',release);
fs.mkdirSync(destination,{recursive:true});for(const name of files)fs.copyFileSync(path.join(pkg,name),path.join(destination,name));
function replace(relative,pattern,replacement){const p=path.join(root,relative);fs.writeFileSync(p,fs.readFileSync(p,'utf8').replace(pattern,replacement));}
replace('web/lib/interactive_series.js','../pkg/rustplotlib_web.js',`../pkg/releases/${release}/rustplotlib_web.js`);
replace('web/index.html',/pkg\/releases\/[0-9a-f]{16}\//g,`pkg/releases/${release}/`);
replace('scripts/install_receiver_monitor_web.sh',/runtime_release=[0-9a-f]{16}/,`runtime_release=${release}`);
replace('web/pkg/README.md',/Runtime release `[0-9a-f]{16}`/,`Runtime release \`${release}\``);
const commit=execFileSync('git',['-C',plot,'rev-parse','--short','HEAD'],{encoding:'utf8'}).trim();
replace('web/pkg/README.md',/commit `[0-9a-f]+` on branch/,`commit \`${commit}\` on branch`);
console.log('Prepared content-addressed Rustplotlib runtime '+release);
JS
"$task_esbuild" "$task_repo/web/monitor.js" --bundle --format=esm '--external:../pkg/releases/*' --outfile="$task_repo/web/monitor.bundle.js"
node --input-type=module - "$task_repo" <<'JS'
import fs from 'node:fs';
import path from 'node:path';
import {createHash} from 'node:crypto';
const root=process.argv[2],web=path.join(root,'web'),assets=path.join(web,'assets');
fs.mkdirSync(assets,{recursive:true});
function name(file,suffix){return 'monitor.'+createHash('sha256').update(fs.readFileSync(file)).digest('hex').slice(0,16)+suffix;}
const staging=path.join(web,'monitor.bundle.js'),script=name(staging,'.js'),style=name(path.join(web,'monitor.css'),'.css');
fs.renameSync(staging,path.join(assets,script));fs.copyFileSync(path.join(web,'monitor.css'),path.join(assets,style));
const index=path.join(web,'index.html');let html=fs.readFileSync(index,'utf8');
html=html.replace(/href="(?:monitor\.css(?:\?[^"]*)?|assets\/monitor\.[0-9a-f]+\.css)"/g,`href="assets/${style}"`)
 .replace(/(?:href|src)="(?:monitor\.js(?:\?[^"]*)?|assets\/monitor\.[0-9a-f]+\.js)"/g,m=>m.slice(0,m.indexOf('='))+`="assets/${script}"`)
 .replace(/^.*rel="modulepreload" href="lib\/interactive_series\.js\?[^"]*".*\n/gm,'');
fs.writeFileSync(index,html);console.log('Bundled immutable dashboard assets '+script+' and '+style);
JS
