import {InteractiveSeriesPlot} from '/rvd/lib/interactive_series.js?v=runtime-c70022f0';
import {InteractiveFieldPlot} from '/rvd/lib/interactive_plot.js?v=monitor-missing-1';
import {chart,colors,zenithTemperature,fitCharts,phaseChart,rawChart,meteorChart,rtiData} from './monitor_data.mjs?v=1';
const $=id=>document.getElementById(id),plots=new Map(),pending=new Map();
let telemetry,rti,chosenDay='',manifest,followLive=true;
let rtiRetry;
const utc=t=>new Date(t*1000).toISOString().replace('T',' ').slice(0,19)+' UTC';
function message(id,text){$(id).textContent=text;}
async function json(url){const response=await fetch(url,{cache:'no-cache'});if(!response.ok)throw Error(`HTTP ${response.status}`);return response.json();}
async function draw(id,spec){
 if(!plots.has(id)){
  if(!pending.has(id))pending.set(id,InteractiveSeriesPlot.mount($(id),spec,{height:320}).then(p=>{plots.set(id,p);return p;}).finally(()=>pending.delete(id)));
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
  const latest=await json('/rvd/data/manifest.json'),days=latest.available_days||[];
  if(followLive||!chosenDay||!days.includes(chosenDay))chosenDay=days.at(-1)||latest.selected_day;
  $('day').replaceChildren(...days.map(day=>{const option=document.createElement('option');option.value=day;option.textContent=day;return option;}));$('day').value=chosenDay;
  const data=latest.selected_day===chosenDay?latest:await json(`/rvd/data/days/${chosenDay}.json`);
  if(!data.frame_packets?.length)throw Error('No completed mesosphere modes available');
  manifest=data;
  const spec={bounds:[0,50,24,120],xlabel:'Time (UTC)',ylabel:'Nominal altitude (km)',colorbar_label:'SNR (dB)',panels:['Zenith','North','East','South','West'],utc_origin_seconds:data.axis_start_us/1e6};
  const parameters={day:Number(chosenDay.replaceAll('-',''))};
  if(!rti)rti=await InteractiveFieldPlot.mount($('rti'),{id:'rti',spec,frameEndpoint:'/rvd/api/frame',frameParameters:parameters,
   colorLimits:[0,30],missingColor:[8/255,39/255,71/255,1]});
  else{rti.setFrameParameters(parameters);rti.setSpec(spec);}
  const fields=rtiData(data);
  rti.panels.forEach((panel,b)=>{const f=fields[b];panel.set_rectilinear_data(f.x,f.y,f.z);});
  $('rti').dataset.renderer='rustplotlib-wasm-webgl';
  clearTimeout(rtiRetry);
  message('pmseStatus',`${chosenDay} · ${data.frame_packets.length} completed modes · source spectra through ${utc(data.latest_mode_us/1e6)}`);
 }catch(error){message('pmseStatus',`PMSE refresh delayed: ${error.message}. Last successful plot is retained.`);
  if(!rti){clearTimeout(rtiRetry);rtiRetry=setTimeout(()=>void refreshRTI(),10000);}}
}
function renderRaw(){if(!telemetry?.raw)return;return attempt('rawVoltage',()=>draw('rawVoltage',rawChart(telemetry.raw,+$('rawChannel').value)));}
async function refreshTelemetry(){
 try{
  const data=await json('./monitor-data.json');telemetry=data;
  const age=Date.now()/1000-data.generated_unix_s,rawAge=Date.now()/1000-data.freshness.Raw;
  message('freshness',`Source telemetry ${utc(data.generated_unix_s)} · raw voltage ${Number.isFinite(rawAge)?Math.max(0,rawAge).toFixed(0)+' s old':'unavailable'}${age>1800?' · telemetry delayed':''}`);
  $('freshness').classList.toggle('stale',age>1800||!data.freshness.Raw||rawAge>600);
  for(const [id,select] of [['counts','counts'],['speed','speed'],['radiant','radiant']])await attempt(id,()=>draw(id,fitCharts(data)[select]));
  await attempt('phase',()=>draw('phase',phaseChart(data)));
  if(data.meteor){message('meteorStatus',utc(data.meteor.time_unix_s));await attempt('meteor',()=>draw('meteor',meteorChart(data.meteor)));}
  if(data.raw){message('rawStatus',`${utc(data.raw.start_unix_s)} · ${data.raw.units}`);await renderRaw();}
  else message('rawStatus',`Raw trace unavailable: ${data.errors?.raw||'no sample yet'}`);
  const origin=data.generated_unix_s-24*3600;
  const series=Object.entries(data.freshness).filter(([name,t])=>Number.isFinite(t)).map(([name,t],i)=>({name,color:colors[i%colors.length],line:false,radius:5,points:[[(t-origin)/3600,i]]}));
  await attempt('processing',()=>draw('processing',chart(series,'Latest processed sample (UTC)','Pipeline stage (legend)',{origin})));
 }catch(error){message('freshness',`Telemetry refresh delayed: ${error.message}. Last successful data are retained.`);$('freshness').classList.add('stale');}
}
async function refreshOperations(){
 try{const ops=await json('./ops_status.json');message('operationsStatus',`Source service report ${ops.generated_utc}`);
  table($('serviceStatus'),['Service','State'],ops.services.map(s=>[s.name.replace('pansy-','').replace('.service',''),s.active_state]));
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
for(let ch=0;ch<8;ch++){const o=document.createElement('option');o.value=String(ch);o.textContent=ch===7?'TX reference (7)':`RX ${ch}`;$('rawChannel').append(o);}
$('rawChannel').onchange=renderRaw;
$('day').onchange=()=>{chosenDay=$('day').value;followLive=chosenDay===$('day').options[$('day').options.length-1]?.value;void refreshRTI();};
function clock(){message('utc',utc(Date.now()/1000));}clock();setInterval(clock,1000);
// Failures in one feed must not prevent the remaining independent panels.
let refreshing=false;
async function refresh(){if(refreshing)return;refreshing=true;$('refresh').disabled=true;
 try{await Promise.allSettled([refreshRTI(),refreshNoise(),refreshTelemetry(),refreshOperations()]);}
 finally{refreshing=false;$('refresh').disabled=false;}}
$('refresh').onclick=()=>void refresh();
void refresh();setInterval(()=>void refresh(),300_000);
