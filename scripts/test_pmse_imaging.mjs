// Local browser assets over read-only deployed APIs; no acquisition changes.
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {createHash} from 'node:crypto';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE||'playwright');
const base=new URL(process.env.PMSE_URL||'https://juha.no/rvd/');
const browser=await chromium.launch({...(process.env.CHROME?{executablePath:process.env.CHROME}:{}),
  args:['--no-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
const context=await browser.newContext({viewport:{width:1200,height:1000}});
const page=await context.newPage(),errors=[],results=[];
page.on('pageerror',e=>errors.push(String(e)));
let spectraRequests=0;
page.on('request',r=>{if(r.url().includes('/api/spectra?'))spectraRequests++;});
try{
  if(process.env.PMSE_TEST_DEPLOYED!=='1')await context.route(`${base.origin}${base.pathname}**`,async route=>{
    const url=new URL(route.request().url()),name=url.pathname.slice(base.pathname.length)||'index.html';
    if(!['index.html','viewer.js','imaging.js','imaging_worker.js','theme.css'].includes(name)
      && !/^imaging_pkg\/[\w.]+$/.test(name))return route.continue();
    const body=await readFile(new URL(`../pmse/rvd/${name}`,import.meta.url));
    await route.fulfill({body,contentType:name.endsWith('.wasm')?'application/wasm':name.endsWith('.js')?'text/javascript':name.endsWith('.css')?'text/css':'text/html'});
  });
  await page.goto(base.href);
  await page.waitForFunction(()=>document.querySelector('#frameStatus').textContent.includes('retained cells'),null,{timeout:120000});
  assert(await page.locator('#sparseImaging').isChecked());
  assert.equal(await page.locator('#imageMethod').inputValue(),'wavelet');
  // The scientific fixture has exact orthogonal look/noise sequences, so its
  // sample covariance is source covariance + identity, without random luck.
  const synthetic=await page.evaluate(async base=>{
    const g=await(await fetch(`${base}api/geometry`)).json();
    const {default:init,BrowserImager}=await import(`${base}imaging_pkg/pansy_browser_imaging.js?v=sparse-1`);await init();
    const r=85,grid=61,looks=64,k=2*Math.PI*g.frequency_hz/299792458,limit=Math.sin(20*Math.PI/180);
    const target=[34,36],u=-limit+target[0]*2*limit/60,v=-limit+target[1]*2*limit/60;
    const direction=[u,v,-Math.sqrt(1-u*u-v*v)];
    const header={index:[[0,1,0]],ranges_km:[85,86],velocity_mps:[-10,0,10],
      noise_power:Array.from({length:7},()=>Array.from({length:5},()=>[1,1])),looks,nyquist_mps:20};
    const values=new Float32Array(looks*7*2);
    function inject(noSignal=false,offGrid=false){
      const d=direction.slice();if(offGrid){d[0]+=.0037;d[1]-=.0041;d[2]=-Math.sqrt(1-d[0]**2-d[1]**2);}
      for(let l=0;l<looks;l++)for(let ch=0;ch<7;ch++){
        const phase=-k*g.positions_m[ch].reduce((s,x,j)=>s+x*d[j],0)-g.phasecal_rad[0][ch];
        const noise=2*Math.PI*(ch+2)*l/looks;
        values[2*(l*7+ch)]=(noSignal?0:Math.sqrt(10)*Math.cos(phase))+Math.cos(noise);
        values[2*(l*7+ch)+1]=(noSignal?0:Math.sqrt(10)*Math.sin(phase))+Math.sin(noise);
      }
      return d;
    }
    const im=new BrowserImager(JSON.stringify(g),r,grid,true),out=[];
    for(const offGrid of [false,true])for(const method of ['sparse','wavelet','tsvd']){
      const d=inject(false,offGrid),started=performance.now();
      const m=im.image_with_options(JSON.stringify(header),JSON.stringify(g),values,JSON.stringify({method,strength:.5,reweight:method==='sparse'}));
      let mass=0,outside=0,peak=-Infinity,best=0;
      for(let p=0;p<grid*grid;p++)if(Number.isFinite(m[4*p])&&m[4*p]>-119){
        const power=10**(m[4*p]/10);mass+=power;
        const x=-limit+(p%grid)*2*limit/60,y=-limit+Math.floor(p/grid)*2*limit/60;
        const dot=x*d[0]+y*d[1]-Math.sqrt(1-x*x-y*y)*d[2];
        if(Math.acos(Math.min(1,Math.max(-1,dot)))>2*Math.PI/180)outside+=power;
        if(m[4*p]>peak){peak=m[4*p];best=p;}
      }
      const x=-limit+(best%grid)*2*limit/60,y=-limit+Math.floor(best/grid)*2*limit/60;
      const error=Math.acos(Math.min(1,Math.max(-1,x*d[0]+y*d[1]-Math.sqrt(1-x*x-y*y)*d[2])))*180/Math.PI;
      out.push({method,offGrid,outside:outside/mass,error,ms:performance.now()-started,diagnostics:JSON.parse(im.last_diagnostics())});
    }
    inject(true);
    const noise=im.image_with_options(JSON.stringify(header),JSON.stringify(g),values,JSON.stringify({method:'wavelet',strength:1,reweight:false}));
    const noisePower=Array.from({length:grid*grid},(_,p)=>noise[4*p+3]).some(v=>v>0);
    im.free();return {cases:out,noisePower};
  },base.href);
  console.log('synthetic',JSON.stringify(synthetic));
  assert(!synthetic.noisePower,'Noise-only covariance must not manufacture emission');
  for(const row of synthetic.cases){assert(Number.isFinite(row.outside));assert(row.error<1.1,'Wrong source direction: '+JSON.stringify(row));}
  for(const offGrid of [false,true]){
    const s=synthetic.cases.find(x=>x.method==='sparse'&&x.offGrid===offGrid),t=synthetic.cases.find(x=>x.method==='tsvd'&&x.offGrid===offGrid);
    assert(s.outside<t.outside*.5,'Sparse must materially suppress this point-source PSF');
  }
  // Choose a real archived echo, preserving all available channels/looks.
  const selection=await page.evaluate(async()=>{
    const day=document.querySelector('#day').value,m=await(await fetch(`./api/day?day=${day}`)).json();
    let best;
    for(const f of m.frame_packets)for(let c=0;c<f.index.length;c++){
      const r=m.geometry.slant_range_km[f.index[c][2]],snr=f.snr_db[c];
      if(r>=70&&r<=110&&(!best||snr>best.snr))best={key:f.key,range:r,snr};
    }
    if(!best)throw Error('No real retained echo available for this day');return best;
  });
  await page.evaluate(async({base,selection})=>{const {imageSelection}=await import(`${base}imaging.js?v=sparse-1`);await imageSelection(selection.key,selection.range);},{base:base.href,selection});
  const ready=method=>page.waitForFunction(method=>document.querySelector('#imagePlot').dataset.method===method&&document.querySelector('#imageStatus').textContent.includes('Rust/WASM'),method,{timeout:120000});
  const hash=async()=>createHash('sha256').update(await page.locator('#imagePlot canvas').screenshot()).digest('hex');
  await ready('wavelet');results.push({method:'wavelet',status:await page.locator('#imageStatus').innerText(),diagnostics:await page.locator('#imageDiagnostics').innerText(),hash:await hash()});
  await page.locator('#imageMethod').selectOption('sparse');await ready('sparse');
  results.push({method:'sparse',status:await page.locator('#imageStatus').innerText(),diagnostics:await page.locator('#imageDiagnostics').innerText(),hash:await hash()});
  await page.locator('#sparseImaging').uncheck();await ready('tsvd');
  results.push({method:'tsvd',hash:await hash()});
  assert.equal(spectraRequests,1,'Method changes must reuse the existing packet');
  assert.notEqual(results[0].hash,results[1].hash,'Methods must change actual WebGL pixels');
  assert.notEqual(results[1].hash,results[2].hash);
  await page.locator('#sparseImaging').check();
  // Rapidly supersede a running solve; stale workers must not win the race.
  await page.locator('#imageMethod').selectOption('wavelet');
  await page.locator('#imageMethod').selectOption('sparse');await ready('sparse');
  assert.equal(spectraRequests,1);
  await page.locator('#display').selectOption('hsv');
  await page.locator('#imagePlot button[aria-label="Adjust color limits"]').first().click();
  assert.equal(await page.locator('dialog').getByRole('spinbutton').count(),4);
  await page.getByRole('button',{name:'Done',exact:true}).click();
  await page.locator('.imaging-references').evaluate(d=>d.open=true);
  assert.equal(await page.locator('.imaging-references a').count(),3);
  await page.locator('#imagePlot').scrollIntoViewIfNeeded();await page.screenshot({path:'/tmp/pansy-sparse-imaging.png'});
  await page.setViewportSize({width:390,height:844});
  assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),'Mobile page overflow');
  assert.equal(await page.locator('[data-plot-error]').count(),0);
  assert.deepEqual(errors,[]);
  console.log(JSON.stringify({passed:true,selection,results,spectraRequests,errors}));
}finally{await browser.close();}
