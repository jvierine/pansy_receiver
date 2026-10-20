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
    let failed=false,requests=0;
    if(!deployed)await page.route(`${base.origin}${base.pathname}**`,async route=>{
      const name=new URL(route.request().url()).pathname.slice(base.pathname.length)||'index.html';
      if(['index.html','viewer.js','theme.css'].includes(name))return route.fulfill({body:await readFile(new URL(`../pmse/rvd/${name}`,import.meta.url)),contentType:name.endsWith('.js')?'text/javascript':name.endsWith('.css')?'text/css':'text/html'});
      if(name==='data/manifest.json')return route.fulfill({json:{available_days:['2026-10-20'],selected_day:'2026-10-20',frame_packets:[],frame_count:0,latest_mode_us:Date.now()*1000,revision:'test',plot_generation_status:'unavailable'}});
      if(name==='data/sky-noise/status.json'){
        requests++;
        if(failed)return route.fulfill({status:503,body:'Temporarily unavailable'});
        return route.fulfill({json:{noise_model_status:Array.from({length:7},()=>({receiver_temp_k:null,status:'insufficient coverage or sky modulation'})),noise_model_observation_stop_unix_s:Date.now()/1000,
          noise_plot:{status:'ready',revision:'test',generated_utc:new Date().toISOString(),files:Object.fromEntries(Array.from({length:7},(_,ch)=>[ch,{640:`noise-module-${ch}-640.svg`,960:`noise-module-${ch}-960.svg`}]))}}});
      }
      if(name.startsWith('data/sky-noise/noise-module-'))return route.fulfill({contentType:'image/svg+xml',body:'<svg xmlns="http://www.w3.org/2000/svg" width="640" height="550"><rect width="640" height="550" fill="#faf7ed"/></svg>'});
      return route.continue();
    });
    await page.goto(base.href);
    // Wait for the asynchronous WebAssembly initialization and initial refresh.
    await page.waitForFunction(()=>document.querySelector('#day').options.length>0);
    if(!deployed)assert.equal(requests,0,'Sky diagnostics should load only when opened');
    await page.getByText('Sky noise & receiver temperature',{exact:true}).click();
    const image=page.locator('#noisePlot');await image.waitFor({state:'visible'});
    await page.waitForFunction(()=>{const p=document.querySelector('#noisePlot');return p.complete&&p.naturalWidth>0;});
    assert.match(await page.locator('#noiseStatus').innerText(),/PyGDSM GSM2008.*Plots updated/);
    for(const ch of ['0','6']){
      await page.locator('#noiseModule').selectOption(ch);
      assert.match(await image.getAttribute('src'),new RegExp(`sky-noise/noise-module-${ch}-${width<900?640:960}\\.svg`));
      await page.waitForFunction(()=>{const p=document.querySelector('#noisePlot');return p.complete&&p.naturalWidth>0;});
    }
    if(!deployed){
      failed=true;
      await page.getByText('Sky noise & receiver temperature',{exact:true}).click();
      await page.getByText('Sky noise & receiver temperature',{exact:true}).click();
      await page.waitForFunction(()=>document.querySelector('#noiseStatus').textContent.includes('Refresh delayed'));
      assert(await image.isVisible(),'A failed refresh must preserve the cached plot');
    }else await image.screenshot({path:`/tmp/pansy-sky-noise-${width}.png`});
    assert.deepEqual(errors,[]);await page.close();
  }
  console.log(JSON.stringify({passed:true,deployed,checks:'sky model visible independently of RTI decorations; modules 0/6, desktop/mobile, cached fallback'}));
}finally{await browser.close();}
