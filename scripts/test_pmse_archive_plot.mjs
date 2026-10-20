import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE||'playwright');
const base=new URL(process.env.PMSE_URL||'https://juha.no/rvd/');
const deployed=process.env.PMSE_TEST_DEPLOYED==='1';
const browser=await chromium.launch({...(process.env.CHROME?{executablePath:process.env.CHROME}:{}),args:['--no-sandbox']});
try{
  for(const width of [1200,390]){
    const page=await browser.newPage({viewport:{width,height:850}}),errors=[];
    page.on('pageerror',error=>errors.push(String(error)));
    let mode='ready';
    if(!deployed)await page.route(`${base.origin}${base.pathname}**`,async route=>{
      const name=new URL(route.request().url()).pathname.slice(base.pathname.length)||'index.html';
      if(['index.html','viewer.js','theme.css'].includes(name))return route.fulfill({body:await readFile(new URL(`../pmse/rvd/${name}`,import.meta.url)),contentType:name.endsWith('.js')?'text/javascript':name.endsWith('.css')?'text/css':'text/html'});
      if(name==='data/manifest.json')return route.fulfill({json:{
        available_days:['2026-10-20'],selected_day:'2026-10-20',frame_packets:[],frame_count:0,
        latest_mode_us:Date.now()*1000,revision:'test',plot_generation_status:'unavailable',
        transfer:{archive:{bytes:1000,records:1,files:1,projected_bytes_per_day:null,hourly_unix_s:[0]}},
        archive_plot:mode==='unavailable'?{status:'unavailable'}:{status:mode,revision:mode,
          generated_utc:'2026-10-20T18:00:00+00:00',files:{640:'axes/archive-growth-640.svg',960:'axes/archive-growth-960.svg'}}
      }});
      if(name.startsWith('data/axes/archive-growth-'))return mode==='broken'?route.fulfill({status:404,body:'Not found'}):route.fulfill({contentType:'image/svg+xml',body:'<svg xmlns="http://www.w3.org/2000/svg" width="640" height="810"><rect width="640" height="810" fill="#faf7ed"/></svg>'});
      return route.continue();
    });
    const modes=deployed?['ready']:['ready','delayed','unavailable','broken'];
    for(mode of modes){
      await page.goto(base.href);
      await page.getByText('Archive & transfer',{exact:true}).click();
      const image=page.locator('#archivePlot'),status=page.locator('#archivePlotStatus');
      if(mode==='unavailable'||mode==='broken'){
        await page.waitForFunction(text=>document.querySelector('#archivePlotStatus').textContent.includes(text),mode==='broken'?'could not be loaded':'being prepared');
        assert(await image.isHidden());
      }else{
        await image.waitFor({state:'visible'});
        await page.waitForFunction(()=>{const img=document.querySelector('#archivePlot');return img.complete&&img.naturalWidth>0;});
        assert.match(await image.getAttribute('src'),new RegExp(`archive-growth-${width<900?640:960}\\.svg`));
        assert.match(await status.innerText(),mode==='delayed'?/last successful storage plots/:/Storage plots updated/);
        if(deployed)await page.screenshot({path:`/tmp/pansy-storage-plots-${width}.png`,fullPage:false});
      }
      assert.deepEqual(errors,[]);
    }
    await page.close();
  }
  console.log(JSON.stringify({passed:true,deployed,checks:'desktop/mobile native storage plots; independent readiness and stale fallback'}));
}finally{await browser.close();}
