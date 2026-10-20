import assert from 'node:assert/strict';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE||'playwright');
const base=new URL(process.env.PMSE_URL||'https://juha.no/rvd/');
const browser=await chromium.launch({...(process.env.CHROME?{executablePath:process.env.CHROME}:{}),args:['--no-sandbox']});
try{
 for(const width of [1200,390]){
  const page=await browser.newPage({viewport:{width,height:850}}),errors=[],svg=[];
  page.on('pageerror',e=>errors.push(String(e)));
  page.on('request',r=>{if(/archive-growth-.*\.svg/.test(r.url()))svg.push(r.url());});
  await page.goto(base.href);await page.getByText('Archive & transfer',{exact:true}).click();
  for(const id of ['hourlyDataPlot','cumulativeDataPlot','retentionPlot']){
   await page.waitForFunction(id=>document.getElementById(id).dataset.renderer==='rustplotlib-wasm-series',id);
   assert(Number(await page.locator('#'+id).getAttribute('data-sample-count'))>0);
  }
  const plot=page.locator('#hourlyDataPlot');
  await plot.getByRole('button',{name:'Reset view'}).click();
  assert.equal(await plot.getAttribute('data-view-bounds'),'[0,0,1,1]');
  assert.deepEqual(svg,[]);assert.deepEqual(errors,[]);
  await page.screenshot({path:`/tmp/pansy-storage-plots-${width}.png`,fullPage:false});await page.close();
 }
 console.log(JSON.stringify({passed:true,checks:'numeric inventory, browser Rustplotlib axes/series, all three panels, no SVG image requests'}));
}finally{await browser.close();}
