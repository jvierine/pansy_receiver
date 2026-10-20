import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE||'playwright');
const base=new URL(process.env.PMSE_URL||'https://juha.no/rvd/');
const deployed=process.env.PMSE_TEST_DEPLOYED==='1';
const browser=await chromium.launch({...(process.env.CHROME?{executablePath:process.env.CHROME}:{}),args:['--no-sandbox']});
try{
  for(const width of [1200,390]){
    const page=await browser.newPage({viewport:{width,height:850}}),errors=[];
    page.on('pageerror',e=>errors.push(String(e)));
    let failed=false,requests=0;const images=[];
    page.on('request',r=>{if(/noise-module-.*\.svg/.test(r.url()))images.push(r.url());});
    if(!deployed)await page.route(`${base.origin}${base.pathname}**`,async route=>{
      const name=new URL(route.request().url()).pathname.slice(base.pathname.length)||'index.html';
      if(['index.html','viewer.js','theme.css'].includes(name))return route.fulfill({body:await readFile(new URL(`../pmse/rvd/${name}`,import.meta.url)),contentType:name.endsWith('.js')?'text/javascript':name.endsWith('.css')?'text/css':'text/html'});
      if(name==='data/manifest.json')return route.fulfill({json:{available_days:['2026-10-20'],selected_day:'2026-10-20',frame_packets:[],frame_count:0,latest_mode_us:Date.now()*1000,revision:'test',plot_generation_status:'unavailable'}});
      if(name==='data/sky-noise/status.json'){
        requests++;
        if(failed)return route.fulfill({status:503,body:'Temporarily unavailable'});
        return route.fulfill({json:{noise_model_status:Array.from({length:7},()=>({receiver_temp_k:null,status:'insufficient coverage or sky modulation'})),noise_model_observation_stop_unix_s:Date.now()/1000,
          noise_model_origin_unix_s:0,noise_model_data:{times_h:[0,1,2],sky_temperature_k:[[9000,10000,11000,12000,13000],[10000,11000,12000,13000,14000],[11000,12000,13000,14000,15000]],noise_power:Array.from({length:3},(_,i)=>Array.from({length:5},()=>Array.from({length:7},(_,ch)=>100+i+ch*20)))},
          noise_plot:{status:'ready',revision:'test',format:'rustplotlib-series-data-v1',generated_utc:new Date().toISOString()}}});
      }
      return route.continue();
    });
    await page.goto(base.href);
    // Wait for the asynchronous WebAssembly initialization and initial refresh.
    await page.waitForFunction(()=>document.querySelector('#day').options.length>0);
    if(!deployed)assert.equal(requests,0,'Sky diagnostics should load only when opened');
    await page.getByText('Sky noise & receiver temperature',{exact:true}).click();
    const image=page.locator('#noisePowerPlot canvas');await image.waitFor({state:'visible'});
    await page.waitForFunction(()=>document.querySelector('#noisePowerPlot').dataset.renderer==='rustplotlib-wasm-series');
    assert.match(await page.locator('#noiseStatus').innerText(),/PyGDSM GSM2008/);
    for(const ch of ['0','6']){
      await page.locator('#noiseModule').selectOption(ch);
      assert.equal(await page.locator('#noisePowerPlot').getAttribute('data-series-count'),'5');
      assert.equal(await page.locator('#skyTemperaturePlot canvas').count(),1);
    }
    if(!deployed){
      failed=true;
      await page.getByText('Sky noise & receiver temperature',{exact:true}).click();
      await page.getByText('Sky noise & receiver temperature',{exact:true}).click();
      await page.waitForFunction(()=>document.querySelector('#noiseStatus').textContent.includes('Refresh delayed'));
      assert(await image.isVisible(),'A failed refresh must preserve the cached plot');
    }else await image.screenshot({path:`/tmp/pansy-sky-noise-${width}.png`});
    assert.deepEqual(errors,[]);assert.deepEqual(images,[],'Diagnostic SVGs must never be requested');await page.close();
  }
  console.log(JSON.stringify({passed:true,deployed,checks:'sky model visible independently of RTI decorations; modules 0/6, desktop/mobile, cached fallback'}));
}finally{await browser.close();}
