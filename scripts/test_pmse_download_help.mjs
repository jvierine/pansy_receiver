import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {createHash} from 'node:crypto';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE||'playwright');
const base=new URL(process.env.PMSE_URL||'https://juha.no/rvd/');
const browser=await chromium.launch({...(process.env.CHROME?{executablePath:process.env.CHROME}:{}),args:['--no-sandbox']});
const context=await browser.newContext({viewport:{width:1200,height:700}});
const page=await context.newPage();
try{
  if(process.env.PMSE_TEST_DEPLOYED!=='1')await context.route(`${base.origin}${base.pathname}**`,async route=>{
    const name=new URL(route.request().url()).pathname.slice(base.pathname.length)||'index.html';
    if(!['index.html','theme.css','pansy-sparse-spectra-guide.pdf'].includes(name))return route.continue();
    await route.fulfill({body:await readFile(new URL(`../pmse/rvd/${name}`,import.meta.url)),
      contentType:name.endsWith('.pdf')?'application/pdf':name.endsWith('.css')?'text/css':'text/html'});
  });
  await page.goto(base.href);
  const link=page.locator('#downloadHelp');await link.waitFor({state:'visible'});
  assert.match(await link.getAttribute('aria-label'),/file contents and data layout.*PDF/);
  assert.equal(await link.getAttribute('target'),'_blank');
  assert.equal(await link.getAttribute('rel'),'noopener');
  const url=await link.getAttribute('href');
  const pdf=await page.evaluate(async url=>{
    const r=await fetch(url),data=await r.arrayBuffer();
    return {ok:r.ok,type:r.headers.get('content-type'),bytes:data.byteLength,
      header:new TextDecoder().decode(data.slice(0,5)),
      hash:Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',data)),b=>b.toString(16).padStart(2,'0')).join('')};
  },url);
  const expected=createHash('sha256').update(await readFile(new URL('../pmse/rvd/pansy-sparse-spectra-guide.pdf',import.meta.url))).digest('hex');
  assert(pdf.ok&&pdf.header==='%PDF-'&&pdf.type.startsWith('application/pdf'));
  assert.equal(pdf.hash,expected,'Help must serve the visually verified UiT memo');
  const layouts=[];
  for(const width of [1200,390,320]){
    await page.setViewportSize({width,height:700});
    const a=await page.locator('#downloadSpectra').boundingBox(),b=await link.boundingBox();
    assert(b.x>=a.x+a.width,'Help must stay next to Download spectra');
    assert(Math.abs((a.y+a.height/2)-(b.y+b.height/2))<2,'Help must share the download row');
    assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),'Horizontal page overflow');
    layouts.push({width,download:a,help:b});
    if(width===390)await page.screenshot({path:'/tmp/pansy-download-help-mobile.png'});
  }
  console.log(JSON.stringify({passed:true,pdf,layouts}));
}finally{await browser.close();}
