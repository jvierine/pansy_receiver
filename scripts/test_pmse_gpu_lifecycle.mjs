// Exercise real WASM/WebGL under the smaller context budget used by some
// mobile browsers. Only local assets and read-only deployed APIs are used.
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE||'playwright');
const base=new URL(process.env.PMSE_URL||'https://juha.no/rvd/');
const deployed=process.env.PMSE_TEST_DEPLOYED==='1';
const browser=await chromium.launch({...(process.env.CHROME?{executablePath:process.env.CHROME}:{}),
  args:['--no-sandbox','--max-active-webgl-contexts=8']});
try{
 const page=await browser.newPage({viewport:{width:1200,height:900}}),errors=[],warnings=[];
 let legacyRequests=0;
 page.on('pageerror',e=>errors.push(String(e)));
 page.on('console',m=>{if(m.type()==='warning'&&m.text().includes('Too many active WebGL'))warnings.push(m.text());});
 if(!deployed)await page.route(`${base.origin}${base.pathname}**`,async route=>{
  const name=new URL(route.request().url()).pathname.slice(base.pathname.length)||'index.html';
  // Poison the mutable old URLs to simulate a phone retaining incompatible
  // glue/binary cache entries. No entry point may request these any more.
  if(['pkg/rustplotlib_web.js','pkg/rustplotlib_web_bg.wasm'].includes(name)){
   legacyRequests++;return route.fulfill({status:500,body:'Stale mutable runtime must not be requested'});
  }
  if(['index.html','viewer.js','imaging.js','diagnostics.js','lib/interactive_plot.js','lib/interactive_series.js'].includes(name)||/^pkg\/releases\/[a-f0-9]+\/rustplotlib_web(?:_bg\.wasm|\.js)$/.test(name))
   return route.fulfill({body:await readFile(new URL(`../pmse/rvd/${name}`,import.meta.url)),contentType:name.endsWith('.wasm')?'application/wasm':name.endsWith('.html')?'text/html':'text/javascript'});
  return route.continue();
 });
 await page.goto(base.href);
 await page.waitForFunction(()=>document.querySelectorAll('#rti canvas[data-gpu-ready]').length>=3,null,{timeout:120000});
 const visibleReady=()=>page.waitForFunction(()=>[...document.querySelectorAll('#rti canvas,#rd canvas,#timeDoppler canvas')]
  .filter(c=>{const r=c.getBoundingClientRect();return r.top<innerHeight-10&&r.bottom>10;})
  .every(c=>c.dataset.gpuReady==='true'),null,{timeout:30000});
 for(let pass=0;pass<2;pass++){
  await page.locator('#timeDoppler').scrollIntoViewIfNeeded();
  await page.waitForFunction(()=>document.querySelectorAll('#timeDoppler canvas').length===5,null,{timeout:120000});
  await visibleReady();
  await page.locator('#rti').scrollIntoViewIfNeeded();await visibleReady();
 }
 // Force a real context loss on a visible RTI and prove its retained arrays
 // are re-uploaded. ReadPixels must agree before and after restoration.
 await page.evaluate(async base=>{
  const {InteractiveFieldPlot}=await import(`${base}lib/interactive_plot.js?v=gpu-lifecycle-test`);
  const root=document.createElement('div');root.id='gpu-lifecycle-test';document.querySelector('main').prepend(root);
  root.scrollIntoView();
  window.gpuTest=await InteractiveFieldPlot.mount(root,{id:'rd',frameParameters:{day:Number(document.querySelector('#day').value.replaceAll('-',''))},
   spec:{bounds:[-1,0,1,1],xlabel:'Velocity (m/s)',ylabel:'Height (km)',colorbar_label:'SNR (dB)',panels:['Test']},colorLimits:[3,6]});
  window.gpuTest.panels[0].set_rectilinear_data(new Float32Array([-1,1]),new Float32Array([0,1]),new Float32Array([4.5]));
 },base.href);
 const fixture=page.locator('#gpu-lifecycle-test canvas');
 await fixture.scrollIntoViewIfNeeded();
 await fixture.waitFor({state:'visible'});
 await page.waitForFunction(()=>document.querySelector('#gpu-lifecycle-test canvas').dataset.gpuReady==='true');
 const pixel=()=>page.evaluate(()=>{const p=window.gpuTest.panels[0];return [...p.field.pixel(Math.floor(p.canvas.width/2),Math.floor(p.canvas.height/2))];});
 const before=await pixel();
 assert.equal(before[3],255);assert(before.slice(0,3).some(v=>v>0),'Fixture must contain visible rendered data');
 await page.evaluate(()=>document.querySelector('#gpu-lifecycle-test canvas').getContext('webgl2').getExtension('WEBGL_lose_context').loseContext());
 await page.waitForFunction(()=>document.querySelector('#gpu-lifecycle-test canvas').dataset.gpuReady!=='true');
 await page.waitForFunction(()=>document.querySelector('#gpu-lifecycle-test canvas').dataset.gpuReady==='true');
 assert.deepEqual(await pixel(),before,'Context recovery must restore identical data pixels');
 await page.locator('#timeDoppler').scrollIntoViewIfNeeded();
 await page.waitForFunction(()=>document.querySelector('#gpu-lifecycle-test canvas').dataset.gpuReady!=='true');
 await fixture.scrollIntoViewIfNeeded();
 await page.waitForFunction(()=>document.querySelector('#gpu-lifecycle-test canvas').dataset.gpuReady==='true');
 assert.deepEqual(await pixel(),before,'Offscreen suspension must preserve data, scale and viewport');
 await page.evaluate(()=>window.gpuTest.panels[0].set_moment_data(new Float32Array([4.5,10,2,0]),30,10));
 const hsv=await pixel();assert.notDeepEqual(hsv,before,'HSV must render a different encoding');
 await page.locator('#timeDoppler').scrollIntoViewIfNeeded();
 await page.waitForFunction(()=>document.querySelector('#gpu-lifecycle-test canvas').dataset.gpuReady!=='true');
 await fixture.scrollIntoViewIfNeeded();
 await page.waitForFunction(()=>document.querySelector('#gpu-lifecycle-test canvas').dataset.gpuReady==='true');
 assert.deepEqual(await pixel(),hsv,'Suspension must preserve HSV data and limits too');
 await page.evaluate(()=>window.gpuTest.dispose());
 assert.deepEqual(errors,[]);assert.deepEqual(warnings,[]);assert.equal(legacyRequests,0);
 assert.equal(await page.locator('#error').innerText(),'');
 console.log(JSON.stringify({passed:true,deployed,contextBudget:8,pixel:before,checks:'RTI/TD scrolling, offscreen resource release, forced GPU loss, identical restored pixels'}));
}finally{await browser.close();}
