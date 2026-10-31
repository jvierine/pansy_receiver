// Real Chromium gestures against the deployed runtime, not synthetic event-only tests.
import assert from 'node:assert/strict';
import {mkdir} from 'node:fs/promises';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE||'playwright');
const url=process.env.MONITOR_URL||'http://127.0.0.1:8279/pansy/';
const artifacts=process.env.ARTIFACT_DIR||'/tmp/pansy-monitor-ui-validation';await mkdir(artifacts,{recursive:true});
const browser=await chromium.launch({executablePath:process.env.CHROME,args:['--no-sandbox']});
try{
 for(const mobile of [false,true]){
  const context=await browser.newContext({viewport:mobile?{width:390,height:844}:{width:1280,height:1000},isMobile:mobile,hasTouch:mobile,deviceScaleFactor:1});
  const page=await context.newPage(),errors=[];page.on('pageerror',e=>errors.push(String(e)));
  await page.goto(url);await page.waitForFunction(()=>document.querySelectorAll('.plot canvas').length===21,null,{timeout:90000});
  await page.waitForTimeout(150);
  assert.equal(await page.locator('details').count(),0,'No hidden/collapsed sections');
  assert.equal(await page.locator('section.card').count(),await page.locator('.maximize').count());
  assert.ok(await page.evaluate(()=>document.fonts.check('11px RustplotlibPagella')),'Explicit scientific font loaded');
  assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),'No horizontal overflow');
  const logos=page.locator('.logos img');assert.equal(await logos.count(),4);
  for(const logo of await logos.all()){
   assert.ok(await logo.isVisible(),`${await logo.getAttribute('alt')} logo visible`);
   assert.ok(await logo.evaluate(img=>img.complete&&img.naturalWidth>0),'Logo asset loaded');
   const rect=await logo.boundingBox();assert.ok(rect.width>0&&rect.x>=0&&rect.x+rect.width<=(mobile?390:1280),'Logo fits viewport');
  }
  await page.screenshot({path:`${artifacts}/${mobile?'mobile':'desktop'}-top.png`});
  const raw=page.locator('section.card').filter({has:page.locator('#rawVoltage')});
  await raw.scrollIntoViewIfNeeded();
  const sizes=await page.locator('#rawVoltage').evaluate(el=>({height:el.clientHeight,columns:getComputedStyle(el).gridTemplateColumns.split(' ').length,rects:[...el.querySelectorAll('canvas')].map(c=>({w:c.clientWidth,h:c.clientHeight}))}));
  assert.equal(sizes.columns,2);assert.ok(sizes.height<480);assert.equal(sizes.rects.length,8);assert.ok(sizes.rects.every(r=>r.w>=150&&r.h<=132));
  await raw.screenshot({path:`${artifacts}/${mobile?'mobile':'desktop'}-raw.png`});
  const before=await page.locator('#raw-0').getAttribute('data-view-bounds');
  await raw.getByRole('button',{name:'Zoom in',exact:true}).click();
  const zoomed=await page.locator('#raw-0').getAttribute('data-view-bounds');assert.notEqual(zoomed,before);
  const originalScroll=await page.evaluate(()=>scrollY);
  await raw.locator('.maximize').click();await page.waitForTimeout(100);
  assert.equal(await raw.locator('.maximize').textContent(),'Restore');assert.equal(await page.locator('#raw-0').getAttribute('data-view-bounds'),zoomed);
  assert.ok(await raw.evaluate(el=>el.getBoundingClientRect().top<10));
  await raw.locator('.maximize').click();await page.waitForTimeout(100);
  assert.ok(Math.abs(await page.evaluate(()=>scrollY)-originalScroll)<3,'Restore scroll location');
  assert.equal(await page.locator('#raw-0').getAttribute('data-view-bounds'),zoomed);
  await raw.getByRole('button',{name:'Reset view',exact:true}).click();assert.equal(await page.locator('#raw-0').getAttribute('data-view-bounds'),before);
  const plot=page.locator('#rti'),canvas=plot.locator('canvas'),card=page.locator('section.card').filter({has:plot});
  await canvas.scrollIntoViewIfNeeded();const b=await canvas.boundingBox(),initial=await plot.getAttribute('data-view-bounds'),scroll=await page.evaluate(()=>scrollY);
  await page.mouse.move(b.x+b.width*.5,b.y+b.height*.5);await page.mouse.wheel(0,150);await page.waitForTimeout(250);
  assert.equal(await plot.getAttribute('data-view-bounds'),initial);assert.ok(await page.evaluate(()=>scrollY)>scroll+30,'Wheel scrolls page, not plot');
  if(mobile){
   await canvas.scrollIntoViewIfNeeded();const box=await canvas.boundingBox(),y=Math.min(780,box.y+box.height-25),x=box.x+box.width*.5,s=await page.evaluate(()=>scrollY),cdp=await context.newCDPSession(page);
   await cdp.send('Input.dispatchTouchEvent',{type:'touchStart',touchPoints:[{x,y}]});
   for(let i=1;i<=8;i++){await cdp.send('Input.dispatchTouchEvent',{type:'touchMove',touchPoints:[{x,y:y-i*18}]});await page.waitForTimeout(20);}
   await cdp.send('Input.dispatchTouchEvent',{type:'touchEnd',touchPoints:[]});await page.waitForTimeout(300);
   assert.ok(await page.evaluate(()=>scrollY)>s+40,'Real finger swipe scrolls page');assert.equal(await plot.getAttribute('data-view-bounds'),initial);
  }
  await card.getByRole('button',{name:'Zoom',exact:true}).click();await canvas.scrollIntoViewIfNeeded();const box=await canvas.boundingBox();
  await page.mouse.move(box.x+box.width*.35,box.y+box.height*.3);await page.mouse.down();await page.mouse.move(box.x+box.width*.7,box.y+box.height*.65,{steps:5});await page.mouse.up();
  assert.notEqual(await plot.getAttribute('data-view-bounds'),initial,'Explicit box zoom works');
  const zoomBounds=await plot.getAttribute('data-view-bounds');await card.getByRole('button',{name:'Pan',exact:true}).click();
  await page.mouse.move(box.x+box.width*.5,box.y+box.height*.4);await page.mouse.down();await page.mouse.move(box.x+box.width*.58,box.y+box.height*.45,{steps:5});await page.mouse.up();assert.notEqual(await plot.getAttribute('data-view-bounds'),zoomBounds,'Explicit pan works');
  await card.getByRole('button',{name:'Reset view',exact:true}).click();assert.equal(await plot.getAttribute('data-view-bounds'),initial);
  assert.equal(await card.getByRole('button',{name:'Scroll',exact:true}).getAttribute('aria-pressed'),'true');
  // Colorbar editing must survive the shared-control refactor.
  await plot.getByRole('button',{name:'Adjust color limits'}).click();assert.equal(await plot.getByRole('slider').count(),2);await plot.getByRole('button',{name:'Done',exact:true}).click();
  assert.deepEqual(errors,[]);console.log(`PASS ${mobile?'phone touch':'desktop'}: 21 native plots, compact raw panels, visible sections, wheel/touch scroll, zoom, pan, reset, maximize/restore and colorbar controls`);
  await context.close();
 }
}finally{await browser.close();}
