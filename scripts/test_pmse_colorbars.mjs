// Exercises the shared Rustplotlib controls with real WASM/WebGL and native
// frames. Only HTTP GET/POST to read-only viewer APIs; no receiver changes.
// PMSE_URL may point at a staging deployment. By default, serve local assets
// over the live viewer's read-only data/API, without deploying those assets.
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE||'playwright');
const base=new URL(process.env.PMSE_URL||'https://juha.no/rvd/');
const browser=await chromium.launch({
  ...(process.env.CHROME?{executablePath:process.env.CHROME}:{}),
  args:['--no-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader'],
});
const page=await browser.newPage({viewport:{width:1200,height:900}}),errors=[];
page.on('pageerror',e=>errors.push(String(e)));
const close=(a,b)=>a.every((v,i)=>Math.abs(v-b[i])<=2);
try{
  if(process.env.PMSE_TEST_DEPLOYED!=='1')await page.route(`${base.origin}${base.pathname}**`,async route=>{
    const url=new URL(route.request().url()),name=url.pathname.slice(base.pathname.length)||'index.html';
    if(!['index.html','viewer.js','imaging.js','theme.css','lib/interactive_plot.js'].includes(name))return route.continue();
    const body=await readFile(new URL(`../pmse/rvd/${name}`,import.meta.url));
    await route.fulfill({body,contentType:name.endsWith('.js')?'text/javascript':name.endsWith('.css')?'text/css':'text/html'});
  });
  await page.goto(base.href);
  await page.waitForFunction(()=>document.querySelector('#frameStatus').textContent.includes('retained cells'),null,{timeout:120000});
  // Each default plot exposes a native-derived colorbar target; no page popup.
  for(const id of ['rti','rd','timeDoppler']){
    await page.locator(`#${id} button[aria-label="Adjust color limits"]`).first().click();
    assert.match(await page.locator('dialog').innerText(),/SNR \(dB\)/);
    await page.getByRole('button',{name:'Done',exact:true}).click();
  }
  // Controlled constant field makes pixel checks independent of live echoes.
  await page.evaluate(async base=>{
    const {InteractiveFieldPlot}=await import(`${base}lib/interactive_plot.js?v=colorbar-test`);
    const root=document.createElement('div');root.id='colorbar-test';document.querySelector('main').prepend(root);
    const day=Number(document.querySelector('#day').value.replaceAll('-',''));
    window.colorbarTest=await InteractiveFieldPlot.mount(root,{id:'rd',frameParameters:{day},
      spec:{bounds:[-1,0,1,1],xlabel:'Velocity (m/s)',ylabel:'Height (km)',colorbar_label:'SNR (dB)',panels:['Test']},colorLimits:[3,6]});
    window.colorbarTest.panels[0].set_rectilinear_data(new Float32Array([-1,1]),new Float32Array([0,1]),new Float32Array([4.5]));
    window.colorbarPixel=()=>{const p=window.colorbarTest.panels[0];return Array.from(p.field.pixel(Math.floor(p.canvas.width/2),Math.floor(p.canvas.height/2)));};
  },base.href);
  const before=await page.evaluate(()=>window.colorbarPixel());
  const open=()=>page.locator('#colorbar-test button[aria-label="Adjust color limits"]').first().click();
  await open();
  await page.getByRole('spinbutton',{name:'Maximum value',exact:true}).fill('12');
  assert.deepEqual(await page.evaluate(()=>window.colorbarTest.getColorLimits()),[3,12]);
  const after=await page.evaluate(()=>window.colorbarPixel());assert(!close(before,after),'Limits must change actual WebGL pixels');
  await page.waitForFunction(()=>window.colorbarTest.frame?.svg&&window.colorbarTest.frameURL,null,{timeout:120000});
  // Invalid and incomplete input must not silently turn into zero or invert.
  const low=page.getByRole('spinbutton',{name:'Minimum value',exact:true});
  await low.fill('20');assert.deepEqual(await page.evaluate(()=>window.colorbarTest.getColorLimits()),[3,12]);
  assert.equal(await low.getAttribute('aria-invalid'),'true');
  await low.fill('');assert.deepEqual(await page.evaluate(()=>window.colorbarTest.getColorLimits()),[3,12]);
  await low.fill('-100');assert.deepEqual(await page.evaluate(()=>window.colorbarTest.getColorLimits()),[-100,12]);
  assert.equal(await page.getByRole('slider',{name:'Minimum',exact:true}).getAttribute('min'),'-100');
  await page.getByRole('button',{name:'Reset scale',exact:true}).click();
  assert.deepEqual(await page.evaluate(()=>window.colorbarTest.getColorLimits()),[3,6]);
  assert(close(before,await page.evaluate(()=>window.colorbarPixel())),'Reset must restore original WebGL colors');
  await page.keyboard.press('Escape');await page.locator('dialog').waitFor({state:'detached'});
  // HSV state must survive scale edits, native frame refresh and resize.
  await page.evaluate(()=>window.colorbarTest.panels[0].set_moment_data(new Float32Array([4.5,10,2,0]),30,10));
  await page.waitForTimeout(400);
  await open();assert.equal(await page.getByRole('spinbutton').count(),4);
  const hsvBefore=await page.evaluate(()=>window.colorbarPixel());
  await page.getByRole('spinbutton',{name:'Velocity +/- m/s value',exact:true}).fill('60');
  assert.deepEqual(await page.evaluate(()=>window.colorbarTest.getHsvLimits()),[60,10]);
  assert(!close(hsvBefore,await page.evaluate(()=>window.colorbarPixel())),'HSV control must update hue');
  const nativeUpdate=page.waitForResponse(r=>r.url().endsWith('/api/frame')&&r.request().postDataJSON()?.velocity_limit===60&&r.request().postDataJSON()?.color_limits?.[1]===9,{timeout:120000});
  await page.getByRole('spinbutton',{name:'Maximum value',exact:true}).fill('9');
  const adjusted=await page.evaluate(()=>window.colorbarPixel());
  assert((await nativeUpdate).ok(),'Native key must accept the same limits');
  await page.waitForTimeout(200);
  assert(close(adjusted,await page.evaluate(()=>window.colorbarPixel())),'Native redraw must preserve HSV shading');
  await page.setViewportSize({width:390,height:844});
  const size=await page.locator('dialog').evaluate(d=>({client:d.clientWidth,scroll:d.scrollWidth,left:d.getBoundingClientRect().left,right:d.getBoundingClientRect().right}));
  assert(size.scroll<=size.client&&size.left>=0&&size.right<=390,'Editor must fit mobile without overflow: '+JSON.stringify(size));
  await page.screenshot({path:'/tmp/pansy-colorbar-mobile.png'});
  await page.getByRole('button',{name:'Reset scale',exact:true}).click();
  assert.deepEqual(await page.evaluate(()=>window.colorbarTest.getHsvLimits()),[30,10]);
  await page.getByRole('button',{name:'Done',exact:true}).click();
  // Colorbar is also an ordinary keyboard-operable button.
  await page.locator('#colorbar-test button[aria-label="Adjust color limits"]').first().focus();
  await page.keyboard.press('Enter');await page.locator('dialog').waitFor({state:'visible'});
  await page.keyboard.press('Escape');await page.locator('dialog').waitFor({state:'detached'});
  await page.evaluate(()=>{window.colorbarTest.dispose();document.querySelector('#colorbar-test').remove();});
  assert.deepEqual(errors,[]);
  console.log(JSON.stringify({passed:true,scalarPixels:{before,after},hsvPixels:{before:hsvBefore,adjusted},mobile:size,errors}));
}finally{await browser.close();}
