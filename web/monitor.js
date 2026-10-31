import {InteractiveSeriesPlot,PlotControlGroup} from './lib/interactive_series.js?v=mobile-1';
import {monitorHealth} from './health.mjs';
import {zenithTemperature,fitCharts,phaseChart,rawChart,meteorCharts,rtiChart,latencyChart} from './monitor_data.mjs?v=restore-2';
const $=id=>document.getElementById(id),plots=new Map(),pending=new Map();
let telemetry,operations,rti,chosenDay='',manifest,followLive=true;
const cardControls=new Map();
const plotHeight=id=>id.startsWith('raw-')?(Number(id.slice(4))<6?105:132):id==='meteorGround'?280:id.startsWith('meteor')?220:['rti','processing'].includes(id)?255:id==='speed'?275:300;
let rtiRetry;
let checkingTelemetry=false,telemetryError='';
const utc=t=>new Date(t*1000).toISOString().replace('T',' ').slice(0,19)+' UTC';
function message(id,text){$(id).textContent=text;}
function freshness(){
 if(!telemetry||!Number.isFinite(telemetry.generated_unix_s))return;
 const age=Date.now()/1000-telemetry.generated_unix_s,raw=telemetry.freshness?.Raw,rawAge=Date.now()/1000-raw;
 message('freshness',`Source telemetry ${utc(telemetry.generated_unix_s)} · raw voltage ${Number.isFinite(rawAge)?Math.max(0,rawAge).toFixed(0)+' s old':'unavailable'}${age>1800?' · telemetry delayed':''}${checkingTelemetry?' · checking for updates…':''}${telemetryError?' · '+telemetryError:''}`);
 $('freshness').classList.toggle('stale',!!telemetryError||age>1800||!Number.isFinite(raw)||rawAge>600);
}
function snapshot(name){try{return JSON.parse(sessionStorage.getItem('pansy-monitor-v2:'+name)||'null');}catch{return null;}}
function saveSnapshot(name,data){try{sessionStorage.setItem('pansy-monitor-v2:'+name,JSON.stringify(data));}catch{/* Storage availability must never block the monitor. */}}
async function json(url){const controller=new AbortController(),timer=setTimeout(()=>controller.abort(),30000);try{const response=await fetch(url,{cache:'no-cache',signal:controller.signal});if(!response.ok)throw Object.assign(Error(`HTTP ${response.status}`),{status:response.status});return await response.json();}finally{clearTimeout(timer);}}
async function draw(id,spec){
 const height=plotHeight(id),controls=cardControls.get($(id).closest('.card'));
 if(!plots.has(id)){
  if(!pending.has(id))pending.set(id,InteractiveSeriesPlot.mount($(id),spec,{height,controls}).then(p=>{plots.set(id,p);return p;}).finally(()=>pending.delete(id)));
  await pending.get(id);
 }
 plots.get(id).setSpec(spec);
}
async function attempt(id,fn){try{await fn();$(id).classList.remove('error');}catch(error){if(!plots.has(id)&&!$(id).querySelector('canvas'))message(id,`Awaiting data: ${error.message}`);console.warn(id,error);}}
function table(root,headers,rows){
 const t=document.createElement('table'),head=document.createElement('tr');
 headers.forEach(s=>{const th=document.createElement('th');th.textContent=s;head.append(th);});t.append(head);
 for(const row of rows){const tr=document.createElement('tr');row.forEach(s=>{const td=document.createElement('td');td.textContent=s;tr.append(td);});t.append(tr);}
 root.replaceChildren(t);
}
async function refreshNoise(){
 try{
  message('noiseStatus','Loading receiver temperatures…');
  let model;
  try{model=await json('/rvd/data/monitor-noise.json');}
  catch(error){if(error.status!==404)throw error;model=await json('/rvd/data/sky-noise/status.json');}
  await draw('noiseTemperature',zenithTemperature(model));
  const fits=model.noise_zenith_model_status;
  table($('noiseFits'),['Receiver','Zenith fitted Trec (K)','Fit RMS','Calibration status'],fits.map((f,ch)=>[
   `RX ${ch}`,Number.isFinite(f.receiver_temp_k)?f.receiver_temp_k.toFixed(0):'Unavailable',
   Number.isFinite(f.fractional_rms)?(100*f.fractional_rms).toFixed(1)+'%':'—',f.status]));
  const age=Date.now()/1000-model.noise_model_observation_stop_unix_s;
  message('noiseStatus',`Zenith only · measured system temperature, not receiver offset · observations through ${utc(model.noise_model_observation_stop_unix_s)}${age>3600?' · source data delayed':''}`);
 }catch(error){message('noiseStatus',`Temperature refresh delayed: ${error.message}. Last successful plot is retained.`);}
}
async function refreshRTI(){
 try{
  if(!rti){const cached=snapshot('overview');if(cached){try{await draw('rti',rtiChart(cached));rti=plots.get('rti');message('pmseStatus',`${cached.selected_day} · source spectra through ${utc(cached.latest_mode_us/1e6)} · checking for updates…`);}catch{/* Retry from the live feed if a cached packet is unusable. */}}}
  let latest;
  try{latest=await json('/rvd/data/monitor-overview.json');}catch{latest=await json('/rvd/data/manifest.json');}
  const days=latest.available_days||[];
  if(followLive||!chosenDay||!days.includes(chosenDay))chosenDay=days.at(-1)||latest.selected_day;
  $('day').replaceChildren(...days.map(day=>{const option=document.createElement('option');option.value=day;option.textContent=day;return option;}));$('day').value=chosenDay;
  const data=latest.selected_day===chosenDay?latest:await json(`/rvd/data/days/${chosenDay}.json`);
  if(!(data.frame_count||data.frame_packets?.length))throw Error('No completed mesosphere modes available');
  manifest=data;
  if(data.schema==='pansy_zenith_display_v1')saveSnapshot('overview',data);
  await draw('rti',rtiChart(data));rti=plots.get('rti');
  clearTimeout(rtiRetry);
  message('pmseStatus',`${chosenDay} · ${data.frame_count??data.frame_packets.length} completed modes · source spectra through ${utc(data.latest_mode_us/1e6)}`);
 }catch(error){message('pmseStatus',`PMSE refresh delayed: ${error.message}. Last successful plot is retained.`);
  if(!rti){clearTimeout(rtiRetry);rtiRetry=setTimeout(()=>void refreshRTI(),10000);}}
}
function renderRaw(){if(!telemetry?.raw)return;return Promise.all(telemetry.raw.channels.map((_,ch)=>attempt(`raw-${ch}`,()=>draw(`raw-${ch}`,rawChart(telemetry.raw,ch)))));}
async function renderTelemetry(data,checking=false){
  telemetry=data;checkingTelemetry=checking;freshness();
  const fits=fitCharts(data);
  for(const [id,select] of [['counts','counts'],['speed','speed'],['radiant','radiant']])await attempt(id,()=>draw(id,fits[select]));
  await attempt('phase',()=>draw('phase',phaseChart(data)));
  if(data.meteor){const m=data.meteor;message('meteorStatus',utc(m.time_unix_s));
   const views=meteorCharts(m),ids={east:'meteorEast',north:'meteorNorth',up:'meteorUp',ground:'meteorGround',signal:'meteorSignal',doppler:'meteorDoppler'};
   for(const [key,spec] of Object.entries(views))await attempt(ids[key],()=>draw(ids[key],spec));
   const speed=m.v0_km_s?Math.hypot(...m.v0_km_s).toFixed(1)+' km/s · ':'';
   message('meteorDetail',speed+(m.residual_std_km?'Blue: measured points, error bars ±2 residual RMS · red: saved straight-line fit. ':'')+
    (m.range_time?'':m.detection_trace?'Range and Doppler show saved upstream detected peaks; the full range–time field appears after the next fitted event.':'This legacy fit lacks range–time and measured-Doppler arrays; showing the available SNR and match-function diagnostics.'));
  }
  if(data.raw){message('rawStatus',`${utc(data.raw.start_unix_s)} · ${data.raw.units}`);await renderRaw();}
  else message('rawStatus',`Raw trace unavailable: ${data.errors?.raw||'no sample yet'}`);
  await attempt('processing',()=>draw('processing',latencyChart(data)));
  if(!document.body.dataset.primaryReadyMs)document.body.dataset.primaryReadyMs=String(Math.round(performance.now()));
}
async function refreshTelemetry(){
 checkingTelemetry=true;freshness();
 const incoming=json('./monitor-data.json').then(data=>({data})).catch(error=>({error}));
 try{
  if(!telemetry){const cached=snapshot('telemetry');if(cached){try{await renderTelemetry(cached,true);}catch{/* Do not let an invalid cache block fresh observations. */}}}
  const {data,error}=await incoming;if(error)throw error;
  telemetryError='';await renderTelemetry(data);saveSnapshot('telemetry',data);
 }catch(error){telemetryError=`Refresh delayed: ${error.message}`;checkingTelemetry=false;if(telemetry)freshness();else message('freshness',`Telemetry refresh delayed: ${error.message}.`);$('freshness').classList.add('stale');}
}
async function refreshOperations(){
 try{const ops=await json('./ops_status.json');operations=ops;renderHealth();message('operationsStatus',`Source service report ${ops.generated_utc}`);
  const services=ops.services.map(s=>{const node=document.createElement('span');node.className='service'+(s.active_state==='active'?'':' bad');node.textContent=s.name.replace('pansy-','').replace('.service','')+': '+s.active_state;return node;});
  $('serviceStatus').replaceChildren(...services);
  const inactive=ops.services.filter(s=>s.active_state!=='active').length;
  message('healthSummary',inactive?`${inactive} services not active`:'all services active');$('healthSummary').className=inactive?'bad':'ok';
  table($('diskStatus'),['Disk','Used','Free (GiB)','Health'],ops.disks.map(d=>[d.path,d.used_percent?.toFixed(1)+'%',d.free_bytes==null?'Unavailable':(d.free_bytes/2**30).toFixed(1),d.error?'Unknown':d.critical?'CRITICAL ≥95%':d.warning?'Warning ≥80%':'OK']));
 }catch(error){message('operationsStatus',`Operations refresh delayed: ${error.message}`);}
 try{const backup=await json('./status.json');
  const description=document.createElement('p');description.textContent=`Backup report ${backup.generated_utc||'undated'}`;
  const remote=document.createElement('div'),mirror=document.createElement('div');
  $('backupStatus').replaceChildren(description,remote,mirror);
  const rows=channels=>(channels||[]).map(c=>[c.channel,c.latest_file_date||'Unavailable',c.rsync_active?'Transferring':String(c.last_rsync_exit_code??'Unknown')]);
  table(remote,['Source metadata channel','Latest UTC day','Transfer exit/state'],rows(backup.channels));
  table(mirror,['Mirror metadata channel','Latest UTC day','Transfer exit/state'],rows(backup.local_mirror?.channels));
 }catch(error){message('backupStatus',`Backup report unavailable: ${error.message}`);}
}
for(let ch=0;ch<8;ch++){const root=document.createElement('div');root.id=`raw-${ch}`;root.className='plot';$('rawVoltage').append(root);}
let lastHealth='';
function renderHealth(){
 const health=monitorHealth(operations,telemetry);
 const key=JSON.stringify(health);if(key===lastHealth)return;lastHealth=key;
 $('systemHealth').dataset.state=health.state;
 message('overallHealth',({ok:'SYSTEM OK',bad:'ACTION REQUIRED',warn:'CHECK DISK HEADROOM',unknown:'HEALTH NOT FULLY CONFIRMED'})[health.state]);
 $('healthItems').replaceChildren(...health.items.map(item=>{const node=document.createElement('div');node.className='health-item';node.dataset.state=item.state;const label=document.createElement('strong'),detail=document.createElement('span');label.textContent=`${item.label}: ${item.state==='bad'?'FAULT':item.state==='warn'?'WARNING':item.state==='unknown'?'UNKNOWN':'OK'}`;detail.textContent=item.detail;node.append(label,detail);return node;}));
}
let maximized;
function restoreCard(){if(!maximized)return;const {card,button,scroll,title}=maximized;card.classList.remove('maximized');button.textContent='Maximize';button.setAttribute('aria-label',`Maximize ${title}`);button.setAttribute('aria-expanded','false');document.body.classList.remove('has-maximized');for(const [id,p]of plots)if(card.contains($(id)))p.setHeight(plotHeight(id));maximized=undefined;window.scrollTo(0,scroll);button.focus({preventScroll:true});}
for(const card of document.querySelectorAll('section.card')){
 const top=document.createElement('div');top.className='card-top';const heading=card.querySelector(':scope > .heading,:scope > h2'),title=card.querySelector('h2')?.textContent||'panel';if(heading)top.append(heading);
 const actions=document.createElement('div');actions.className='card-actions';
 if(card.querySelector('.plot,.raw-grid')){const controls=new PlotControlGroup();cardControls.set(card,controls);actions.append(controls.element);}
 const button=document.createElement('button');button.type='button';button.className='maximize';button.textContent='Maximize';button.setAttribute('aria-expanded','false');button.setAttribute('aria-label',`Maximize ${title}`);
 button.onclick=()=>{if(maximized?.card===card){restoreCard();return;}restoreCard();maximized={card,button,scroll:window.scrollY,title};card.classList.add('maximized');document.body.classList.add('has-maximized');button.textContent='Restore';button.setAttribute('aria-label',`Restore ${title}`);button.setAttribute('aria-expanded','true');for(const [id,p]of plots)if(card.contains($(id)))p.setHeight(id.startsWith('raw-')?Math.max(plotHeight(id),Math.floor((innerHeight-200)/4)):Math.max(plotHeight(id),Math.min(650,innerHeight-220)));card.scrollTop=0;button.focus({preventScroll:true});};
 actions.append(button);top.append(actions);card.prepend(top);
}
document.addEventListener('keydown',event=>{if(event.key==='Escape')restoreCard();});
$('day').onchange=()=>{chosenDay=$('day').value;followLive=chosenDay===$('day').options[$('day').options.length-1]?.value;void refreshRTI();};
function clock(){message('utc',utc(Date.now()/1000));freshness();renderHealth();}clock();setInterval(clock,1000);
// Failures in one feed must not prevent the remaining independent panels.
let refreshing=false;
async function refresh(){if(refreshing)return;refreshing=true;$('refresh').disabled=true;
 try{await Promise.allSettled([refreshRTI(),refreshNoise(),refreshTelemetry(),refreshOperations()]);}
 finally{refreshing=false;$('refresh').disabled=false;if(!document.body.dataset.initialReadyMs)document.body.dataset.initialReadyMs=String(Math.round(performance.now()));}}
$('refresh').onclick=()=>void refresh();
void refresh();setInterval(()=>void refresh(),300_000);
