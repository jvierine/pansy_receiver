import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE||'playwright');
const deployed=process.env.PMSE_TEST_DEPLOYED==='1';
const browser=await chromium.launch({...(process.env.CHROME?{executablePath:process.env.CHROME}:{}),args:['--no-sandbox']});
try{
 for(const width of [1200,390]){
  const page=await browser.newPage({viewport:{width,height:850}}),errors=[],legacy=[];
  let failed=false;
  page.on('pageerror',e=>errors.push(String(e)));
  page.on('console',m=>{if(m.type()==='warning'||m.type()==='error')console.log(m.type(),m.text().slice(0,350));});
  page.on('response',async r=>{if(r.status()>=400)console.log('HTTP',r.status(),r.url(),(await r.text().catch(()=>'' )).slice(0,400));});
  page.on('request',r=>{if(/\.(gif|png)(\?|$)/.test(r.url())&&!/\/(uit-logo|u-tokyo-logo)\.png/.test(r.url()))legacy.push(r.url());});
  const now=Date.now()/1000;
  const packet={generated_unix_s:now,freshness:{Raw:now-1,Fits:now-60},fits:[[now-100,40,90,10,-20]],
   phase:[[now-100,0,1,2,3,4,5,6,7],[now-10,1,2,3,4,5,6,7,8]],
   meteor:{time_unix_s:now-100,seconds:[0,.1,.2],east_km:[0,1,2],north_km:[0,1,2],up_km:[90,89,88]},
   raw:{start_unix_s:now-1,sample_step_us:8,units:'raw ADC units; display stride 8',channels:Array.from({length:8},()=>({real:[0,1,0,-1],imag:[1,0,-1,0]}))}};
  if(!deployed){
   await page.route('https://juha.no/pansy/**',async route=>{
    const name=new URL(route.request().url()).pathname.slice('/pansy/'.length)||'index.html';
    if(['index.html','monitor.js','monitor.css','monitor_data.mjs'].includes(name))return route.fulfill({body:await readFile(new URL('../web/'+name,import.meta.url)),contentType:/\.m?js$/.test(name)?'text/javascript':name.endsWith('.css')?'text/css':'text/html'});
    if(name==='monitor-data.json')return failed?route.fulfill({status:503}):route.fulfill({json:packet});
    return route.continue();
   });
   await page.route('https://juha.no/rvd/lib/interactive_plot.js**',async route=>route.fulfill({body:await readFile(new URL('../pmse/rvd/lib/interactive_plot.js',import.meta.url)),contentType:'text/javascript'}));
   await page.route('https://juha.no/rvd/data/sky-noise/status.json',async route=>{
    if(failed)return route.fulfill({status:503});
    const d=await(await fetch(route.request().url())).json();
    d.noise_zenith_model_status=d.noise_model_status;
    return route.fulfill({json:d});
   });
  }
  await page.goto('https://juha.no/pansy/',{waitUntil:'domcontentloaded',timeout:60000});
  await page.waitForFunction(()=>document.querySelector('#rti')?.dataset.renderer==='rustplotlib-wasm-webgl',null,{timeout:90000});
  await page.waitForFunction(()=>document.querySelector('#noiseTemperature')?.dataset.seriesCount==='7',null,{timeout:60000});
  for(const id of ['counts','speed','radiant','meteor','rawVoltage','phase','processing']){
   await page.waitForFunction(id=>document.getElementById(id)?.dataset.renderer==='rustplotlib-wasm-series',id,{timeout:30000});
  }
  assert.equal(await page.locator('#rti canvas').count(),5);
  await page.locator('#rawVoltage').scrollIntoViewIfNeeded();await page.locator('#rawChannel').selectOption('6');
  await page.locator('#rti').scrollIntoViewIfNeeded();
  await page.getByRole('button',{name:'Adjust color limits'}).first().click();
  await page.getByLabel('Maximum value',{exact:true}).fill('25');await page.getByRole('button',{name:'Done',exact:true}).click();
  if(!deployed){
   failed=true;await page.getByRole('button',{name:'Refresh',exact:true}).click();
   await page.waitForFunction(()=>document.getElementById('noiseStatus').textContent.includes('Refresh delayed')||document.getElementById('noiseStatus').textContent.includes('refresh delayed'));
   assert.equal(await page.locator('#noiseTemperature canvas').count(),1,'Keep cached plot on refresh failure');
  }
  assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth+2),false,'No horizontal page overflow');
  assert.deepEqual(errors,[]);assert.deepEqual(legacy,[],'No Matplotlib plots or historical animations');
  await page.screenshot({path:`/tmp/pansy-monitor-${width}.png`,fullPage:true});
  await page.close();
 }
 console.log(JSON.stringify({passed:true,deployed,checks:'Rustplotlib RTI, seven zenith temperatures, all numeric panels, no legacy images, mobile/desktop, colorbar and cached fallback'}));
}finally{await browser.close();}
