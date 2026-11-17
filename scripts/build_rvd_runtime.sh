#!/usr/bin/env bash
# Development/analysis hosts only. Browser axes are built into Rustplotlib WASM.
set -euo pipefail
task_plot_source="$(cd "${1:?Pass a Rustplotlib checkout}" && pwd)"
task_repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_esbuild="${PANSY_ESBUILD_BIN:-$task_plot_source/node_modules/.bin/esbuild}"
task_cargo="${PANSY_CARGO_BIN:-$HOME/.cargo/bin/cargo}"
task_bindgen="${PANSY_BINDGEN_BIN:-$HOME/.cargo/bin/wasm-bindgen}"
(cd "$task_plot_source" && "$task_cargo" build -p rustplotlib-web --release --target wasm32-unknown-unknown)
"$task_bindgen" --target web --out-dir "$task_repo/pmse/rvd/pkg" "$task_plot_source/target/wasm32-unknown-unknown/release/rustplotlib_web.wasm"
node --input-type=module - "$task_repo" "$task_plot_source" <<'JS'
import fs from 'node:fs';import path from 'node:path';import{createHash}from'node:crypto';import{execFileSync}from'node:child_process';
const[root,plot]=process.argv.slice(2),web=path.join(root,'pmse/rvd'),pkg=path.join(web,'pkg');
const files=['rustplotlib_web.js','rustplotlib_web_bg.wasm'],hash=createHash('sha256');for(const f of files)hash.update(fs.readFileSync(path.join(pkg,f)));
const release=hash.digest('hex').slice(0,16),dest=path.join(pkg,'releases',release);fs.mkdirSync(dest,{recursive:true});for(const f of files)fs.copyFileSync(path.join(pkg,f),path.join(dest,f));
for(const file of ['viewer.js','lib/interactive_plot.js','lib/interactive_series.js','index.html']){const p=path.join(web,file);let text=fs.readFileSync(p,'utf8');text=text.replace(/(?:\.\.\/|\.\/|\/rvd\/)?pkg\/releases\/[0-9a-f]{16}\/rustplotlib_web\.js/g,`/rvd/pkg/releases/${release}/rustplotlib_web.js`).replace(/(?:\/rvd\/)?pkg\/releases\/[0-9a-f]{16}\/rustplotlib_web_bg\.wasm/g,`/rvd/pkg/releases/${release}/rustplotlib_web_bg.wasm`);fs.writeFileSync(p,text);}
fs.writeFileSync(path.join(web,'runtime-release.txt'),release+'\n');
fs.writeFileSync(path.join(pkg,'README.md'),'Browser runtime built from rustplotlib '+execFileSync('git',['-C',plot,'rev-parse','HEAD'],{encoding:'utf8'}).trim()+'.\nRelease '+release+' is the concatenated JS/WASM SHA256 prefix. Native field decorations run in the browser.\n');
JS
"$task_esbuild" "$task_repo/pmse/rvd/viewer.js" --bundle --format=esm '--external:/rvd/pkg/releases/*' --outfile="$task_repo/pmse/rvd/viewer.bundle.js"
node --input-type=module - "$task_repo" <<'JS'
import fs from'node:fs';import path from'node:path';import{createHash}from'node:crypto';
const web=path.join(process.argv[2],'pmse/rvd'),assets=path.join(web,'assets');fs.mkdirSync(assets,{recursive:true});
const name=(file,suffix)=>'viewer.'+createHash('sha256').update(fs.readFileSync(file)).digest('hex').slice(0,16)+suffix;
const bundle=path.join(web,'viewer.bundle.js'),js=name(bundle,'.js'),css=name(path.join(web,'theme.css'),'.css');fs.renameSync(bundle,path.join(assets,js));fs.copyFileSync(path.join(web,'theme.css'),path.join(assets,css));
let html=fs.readFileSync(path.join(web,'index.html'),'utf8');
html=html.replace(/^.*rel="modulepreload".*\n/gm,'').replace(/href="(?:theme\.css(?:\?[^"]*)?|assets\/viewer\.[0-9a-f]+\.css)"/g,`href="assets/${css}"`).replace(/src="(?:viewer\.js(?:\?[^"]*)?|assets\/viewer\.[0-9a-f]+\.js)"/g,`src="assets/${js}"`);
html=html.replace(/import\(['"]\.\/(?:viewer\.js(?:\?[^'"]*)?|assets\/viewer\.[0-9a-f]+\.js)['"]\)/g,`import('./assets/${js}')`);
const release=fs.readFileSync(path.join(web,'runtime-release.txt'),'utf8').trim();
html=html.replace('</head>',`<link rel="modulepreload" href="assets/${js}">\n<link rel="modulepreload" href="/rvd/pkg/releases/${release}/rustplotlib_web.js">\n</head>`);fs.writeFileSync(path.join(web,'index.html'),html);console.log('Built immutable RVD browser assets: '+js);
JS
