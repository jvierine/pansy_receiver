import assert from 'node:assert/strict';
import {test} from 'node:test';
import {readFile} from 'node:fs/promises';
import {persistentEchoes,filterImageCells} from '../pmse/rvd/echo_filter.mjs';
const frame=(i,index,snr_db=index.map(()=>5))=>({key:String(i),start_us:i*220e6,stop_us:i*220e6+86e6,index,snr_db});

test('single-mode suppression is checked on initial page load',async()=>{
 const html=await readFile(new URL('../pmse/rvd/index.html',import.meta.url),'utf8');
 assert.match(html,/<input\b[^>]*id="hideSingleMode"[^>]*\bchecked\b[^>]*>/);
});

test('page and native decorations remain white independently of the colormap',async()=>{
 const css=await readFile(new URL('../pmse/rvd/theme.css',import.meta.url),'utf8');
 const renderer=await readFile(new URL('../pmse/plotting/src/main.rs',import.meta.url),'utf8');
 assert.match(css,/--paper:\s*#ffffff;/);assert.match(css,/--foam:\s*#ffffff;/);
 assert.match(renderer,/const WAVE_PAPER:\[f64;3\]=\[1\.,1\.,1\.\];/);
});

test('same beam and exact range survive Doppler changes; isolated gates disappear',()=>{
 const source=[frame(0,[[0,1,10],[0,1,20]]),frame(1,[[0,8,10],[1,1,20],[0,8,30]]),frame(2,[[0,2,30]])];
 const before=JSON.stringify(source),out=persistentEchoes(source);
 assert.deepEqual(out.frames.map(f=>f.index),[[[0,1,10]],[[0,8,10],[0,8,30]],[[0,2,30]]]);
 assert.equal(out.hiddenCells,2);assert.equal(out.hiddenGates,2);
 assert.equal(JSON.stringify(source),before,'Original archive metadata stays unchanged');
});
test('multiple bins in one mode do not establish persistence',()=>{
 const out=persistentEchoes([frame(0,[[0,0,5],[0,1,5],[0,2,5]])]);
 assert.equal(out.frames[0].index.length,0);assert.equal(out.hiddenGates,1);
});
test('gaps, different beams and invalid values cannot support an echo',()=>{
 assert.equal(persistentEchoes([frame(0,[[0,1,5]]),frame(10,[[0,1,5]])]).hiddenCells,2);
 assert.equal(persistentEchoes([frame(0,[[0,1,5]]),frame(1,[[1,1,5]])]).hiddenCells,2);
 assert.equal(persistentEchoes([frame(0,[[0,1,5]]),frame(1,[[0,1,5]],[NaN])]).frames[0].index.length,0);
});
test('newest isolated echo becomes visible on confirmation; all metadata arrays stay aligned',()=>{
 const a=frame(0,[[0,1,5],[0,2,6]],[4,9]);a.look_dominance=[.1,.2];
 assert.equal(persistentEchoes([a]).frames[0].index.length,0);
 const out=persistentEchoes([a,frame(1,[[0,9,6]])]);
 assert.deepEqual(out.frames[0].snr_db,[9]);assert.deepEqual(out.frames[0].look_dominance,[.2]);
});
test('image filtering preserves exact complex samples in every look/channel row',()=>{
 const header={looks:2,index:[[0,1,5],[1,1,5],[4,9,5]]};
 const values=Float32Array.from({length:2*7*3*2},(_,i)=>i+.25),saved=values.slice();
 const out=filterImageCells(header,values,['1:5','4:5']);
 assert.deepEqual(out.header.index,header.index.slice(1));assert.equal(out.header.hidden_single_mode_cells,1);
 for(let row=0;row<14;row++)assert.deepEqual(out.values.slice(row*4,row*4+4),values.slice(row*6+2,row*6+6));
 assert.deepEqual(values,saved);assert.strictEqual(filterImageCells(header,values,undefined).values,values);
 assert.equal(filterImageCells(header,values,[]).values.length,0);
});
test('wave palette retains endpoints, adds hues and increases luminance monotonically',async()=>{
 const palette=JSON.parse(await readFile(new URL('../pmse/rvd/palette.json',import.meta.url)));
 assert(palette.length>=9);assert.deepEqual(palette[0],[.03137255,.15294118,.27843137]);
 const linear=c=>c<=.04045?c/12.92:((c+.055)/1.055)**2.4;
 const luminance=palette.map(c=>c.map(linear).reduce((s,v,i)=>s+v*[.2126,.7152,.0722][i],0));
 for(let i=1;i<luminance.length;i++)assert(luminance[i]>luminance[i-1]);
});
