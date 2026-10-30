import {InteractiveSeriesPlot} from './lib/interactive_series.js?v=restore-2';
import {zenithTemperature,fitCharts,phaseChart,rawChart,meteorCharts,rtiChart,latencyChart} from './monitor_data.mjs?v=restore-2';
const $=id=>document.getElementById(id),plots=new Map(),pending=new Map();
let telemetry,rti,chosenDay='',manifest,followLive=true;
let rtiRetry;
const utc=t=>new Date(t*1000).toISOString().replace('T',' ').slice(0,19)+' UTC';
function message(id,text){$(id).textContent=text;}
async function json(url){const response=await fetch(url,{cache:'no-cache'});if(!response.ok)throw Error(`HTTP ${response.status}`);return response.json();}
async function draw(id,spec){
 const height=id.startsWith('raw-')?170:id.startsWith('meteor')?220:['rti','processing'].includes(id)?255:id==='speed'?275:300;
 if(!plots.has(id)){
  if(!pending.has(id))pending.set(id,InteractiveSeriesPlot.mount($(id),spec,{height,compact:true}).then(p=>{plots.set(id,p);return p;}).finally(()=>pending.delete(id)));
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
  const model=await json('/rvd/data/sky-noise/status.json');
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
  let latest;
  try{latest=await json('/rvd/data/monitor-overview.json');}catch{latest=await json('/rvd/data/manifest.json');}
  const days=latest.available_days||[];
  if(followLive||!chosenDay||!days.includes(chosenDay))chosenDay=days.at(-1)||latest.selected_day;
  $('day').replaceChildren(...days.map(day=>{const option=document.createElement('option');option.value=day;option.textContent=day;return option;}));$('day').value=chosenDay;
  const data=latest.selected_day===chosenDay?latest:await json(`/rvd/data/days/${chosenDay}.json`);
  if(!(data.frame_count||data.frame_packets?.length))throw Error('No completed mesosphere modes available');
  manifest=data;
  await draw('rti',rtiChart(data));rti=plots.get('rti');
  clearTimeout(rtiRetry);
  message('pmseStatus',`${chosenDay} · ${data.frame_count??data.frame_packets.length} completed modes · source spectra through ${utc(data.latest_mode_us/1e6)}`);
 }catch(error){message('pmseStatus',`PMSE refresh delayed: ${error.message}. Last successful plot is retained.`);
  if(!rti){clearTimeout(rtiRetry);rtiRetry=setTimeout(()=>void refreshRTI(),10000);}}
}
function renderRaw(){if(!telemetry?.raw)return;return Promise.all(telemetry.raw.channels.map((_,ch)=>attempt(`raw-${ch}`,()=>draw(`raw-${ch}`,rawChart(telemetry.raw,ch)))));}
async function refreshTelemetry(){
 try{
  const data=await json('./monitor-data.json');telemetry=data;
  const age=Date.now()/1000-data.generated_unix_s,rawAge=Date.now()/1000-data.freshness.Raw;
  message('freshness',`Source telemetry ${utc(data.generated_unix_s)} · raw voltage ${Number.isFinite(rawAge)?Math.max(0,rawAge).toFixed(0)+' s old':'unavailable'}${age>1800?' · telemetry delayed':''}`);
  $('freshness').classList.toggle('stale',age>1800||!data.freshness.Raw||rawAge>600);
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
 }catch(error){message('freshness',`Telemetry refresh delayed: ${error.message}. Last successful data are retained.`);$('freshness').classList.add('stale');}
}
async function refreshOperations(){
 try{const ops=await json('./ops_status.json');message('operationsStatus',`Source service report ${ops.generated_utc}`);
  const services=ops.services.map(s=>{const node=document.createElement('span');node.className='service'+(s.active_state==='active'?'':' bad');node.textContent=s.name.replace('pansy-','').replace('.service','')+': '+s.active_state;return node;});
  $('serviceStatus').replaceChildren(...services);
  const inactive=ops.services.filter(s=>s.active_state!=='active').length;
  message('healthSummary',inactive?`${inactive} services not active`:'all services active');$('healthSummary').className=inactive?'bad':'ok';
  table($('diskStatus'),['Disk','Used','Free (GiB)'],ops.disks.map(d=>[d.path,d.used_percent?.toFixed(1)+'%',d.free_bytes==null?'Unavailable':(d.free_bytes/2**30).toFixed(1)]));
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
$('noiseSection').ontoggle=()=>{if($('noiseSection').open)void refreshNoise();};
$('day').onchange=()=>{chosenDay=$('day').value;followLive=chosenDay===$('day').options[$('day').options.length-1]?.value;void refreshRTI();};
function clock(){message('utc',utc(Date.now()/1000));}clock();setInterval(clock,1000);
// Failures in one feed must not prevent the remaining independent panels.
let refreshing=false;
async function refresh(){if(refreshing)return;refreshing=true;$('refresh').disabled=true;
 try{await Promise.allSettled([refreshRTI(),...($('noiseSection').open?[refreshNoise()]:[]),refreshTelemetry(),refreshOperations()]);}
 finally{refreshing=false;$('refresh').disabled=false;}}
$('refresh').onclick=()=>void refresh();
void refresh();setInterval(()=>void refresh(),300_000);
