import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE||'playwright');
const base=new URL(process.env.PMSE_URL||'https://juha.no/rvd/');
const deployed=process.env.PMSE_TEST_DEPLOYED==='1';
const browser=await chromium.launch({...(process.env.CHROME?{executablePath:process.env.CHROME}:{}),args:['--no-sandbox']});
try{
  const page=await browser.newPage({viewport:{width:1200,height:800}});
  if(!deployed)await page.route(`${base.origin}${base.pathname}**`,async route=>{
    const name=new URL(route.request().url()).pathname.slice(base.pathname.length)||'index.html';
    if(['index.html','viewer.js'].includes(name))return route.fulfill({body:await readFile(new URL(`../pmse/rvd/${name}`,import.meta.url)),contentType:name.endsWith('.js')?'text/javascript':'text/html'});
    if(name==='data/manifest.json'){
      const response=await route.fetch(),m=await response.json();
      m.source_processing={generated_unix_s:Date.now()/1000,live:{state:'active',last_written_stop_us:1792479684000000},backfill:{state:'active',current_index:254,total_modes:2374,current_start_us:1791992495000000}};
      m.transfer.archive.days={'2026-10-12':{latest_mode_us:1791780549000000}};
      m.plot_generation_status='delayed';
      return route.fulfill({json:m});
    }
    return route.continue();
  });
  await page.goto(base.href);
  const status=page.locator('#processingStatus');await status.waitFor({state:'visible'});
  const message=await status.innerText();
  assert.match(message,/Source live processing/);
  assert.match(message,/Rectangular backfill.*\d+ \/ \d+/);
  assert.match(message,/Historical voltages delivered/);
  assert.match(message,/processing and transfer are separate/);
  if(!deployed){
    assert.match(message,/Optional plot refresh delayed/);
    assert(await page.locator('#archivePlot').isHidden());
    assert(await page.locator('#noisePlot').isHidden());
  }
  await page.locator('#day option').first().waitFor({state:'attached'});
  const days=await page.locator('#day option').allTextContents();
  assert(days.length>1);
  await page.waitForFunction(()=>document.querySelector('#frameStatus').textContent.includes('retained cells'),null,{timeout:120000});
  assert.equal(await page.locator('#rti canvas').count(),5);
  assert.equal(await page.locator('#error').innerText(),'');
  console.log(JSON.stringify({passed:true,deployed,message,days,observationStatus:await page.locator('#status').innerText()}));
}finally{await browser.close();}
