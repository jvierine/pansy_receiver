// Scientific provenance: rustplotlib interactive API d3837b8; numerical imaging is Rust/WASM.
import init from '/rvd/pkg/releases/e7774d827ba5b8e0/rustplotlib_web.js';
import {imageSelection,setEchoFilter} from './imaging.js?v=loading-3';
import {persistentEchoes} from './echo_filter.mjs?v=1';
import {renderNoiseDiagnostics,renderArchiveDiagnostics} from './diagnostics.js?v=runtime-c70022f0';
import {InteractiveFieldPlot} from './lib/interactive_plot.js?v=loading-3';
import {readJSON} from './loading.mjs?v=1';
const $=id=>document.getElementById(id), names=['Zenith','North','East','South','West'];
// Start data I/O while WebAssembly initializes.
const fetchData=(url)=>fetch(url,{cache:'no-cache',signal:AbortSignal.timeout(30000)});
let initialResponse=fetchData('./data/manifest.json').catch(error=>({error}));
let manifest,frames=[],plots=[],rd,playing=null,selected=0,revision='',rtiPlot,rdPlot,tdPlot;
let chosenDay='',refreshSequence=0;
let originalFrames=[];
let summaryCache=new WeakMap();
let catalogue;
let showLoading=true;
const snapshotKey='pansy-rvd-day-view-v1';
function saveView(packet){try{const saved=JSON.stringify(packet);if(saved.length<2500000)sessionStorage.setItem(snapshotKey,saved);}catch{/* Storage must never block current observations. */}}
async function restoreView(){
 try{
  const packet=JSON.parse(sessionStorage.getItem(snapshotKey)||'null');
  if(!packet?.frame_packets?.length||!packet.geometry||!packet.selected_day)return;
  chosenDay=packet.selected_day;manifest=packet;originalFrames=packet.frame_packets;applyEchoFilter();
  const option=document.createElement('option');option.value=chosenDay;option.textContent=chosenDay;$('day').replaceChildren(option);
  selected=frames.length-1;$('frame').max=String(selected);$('frame').value=String(selected);
  await mountPlots();renderRTI();renderRD();revision=packet.revision;
  $('status').textContent=`Saved view · ${chosenDay} UTC · observations through ${utc(packet.latest_mode_us)} · checking for updates…`;
  document.body.dataset.plotsReadyMs=String(Math.round(performance.now()));loaded();
 }catch{revision='';try{sessionStorage.removeItem(snapshotKey);}catch{/* A fresh fetch remains available. */}}
}
function loading(stage,progress){
 if(!showLoading)return;
 $('loadProgress').hidden=false;$('loadStage').textContent=stage;
 const bar=$('loadBar');
 if(progress?.total){bar.max=progress.total;bar.value=progress.received;}
 else bar.removeAttribute('value');
 $('loadDetail').textContent=progress?.received?`${(progress.received/2**20).toFixed(2)} MiB received (decoded)` : 'Please wait…';
}
function loaded(){showLoading=false;$('loadProgress').hidden=true;}
function applyEchoFilter(){
 summaryCache=new WeakMap();
 const enabled=$('hideSingleMode').checked;
 const filtered=enabled?persistentEchoes(originalFrames):null;
 frames=filtered?filtered.frames:originalFrames;
 setEchoFilter(filtered?key=>[...(filtered.allowed.get(String(key))||[])]:undefined);
 $('filterStatus').textContent=filtered?`Single-mode filter: ${filtered.hiddenGates.toLocaleString()} beam/range detections hidden. Requires the same range gate in this beam in a neighboring mode (gap ≤10 min). New echoes wait for confirmation; short real echoes can also be hidden.`:'';
}
let skyNoise,noiseRequest;
let archiveInventory;
const secondaryMounts=new Map();
const visiblePlots=new IntersectionObserver(entries=>{
 for(const entry of entries)if(entry.isIntersecting){
  visiblePlots.unobserve(entry.target);
  secondaryMounts.get(entry.target)?.().then(()=>{renderRTI();renderRD();}).catch(failure);
 }
},{rootMargin:'200px'});
const axisStart=()=>manifest.axis_start_us??manifest.first_mode_us;
const axisStop=()=>manifest.axis_stop_us??manifest.latest_mode_us;
const dayNumber=()=>Number(manifest.selected_day.replaceAll('-',''));
const utc=us=>new Date(us/1000).toISOString().replace('T',' ').slice(0,19)+' UTC';
const edges=centers=>Float32Array.from([...centers.map((x,i)=>i?(centers[i-1]+x)/2:x-(centers[1]-x)/2),centers.at(-1)+(centers.at(-1)-centers.at(-2))/2]);
async function mountPlots(){
 const hours=(axisStop()-axisStart())/3.6e9,origin=axisStart()/1e6;
 const spec={bounds:[0,50,Math.max(hours,.001),120],xlabel:'Time (UTC)',ylabel:'Nominal altitude (km)',colorbar_label:'SNR (dB)',panels:names,utc_origin_seconds:origin};
 if(!rtiPlot){
  rtiPlot=await InteractiveFieldPlot.mount($('rti'),{id:'rti',spec,frameParameters:{day:dayNumber()},colorLimits:[+$('low').value,+$('high').value],onSelect:({panel,x,y})=>{
   const time=axisStart()+x*3.6e9;
   selected=frames.reduce((best,f,i)=>Math.abs((f.start_us+f.stop_us)/2-time)<Math.abs((frames[best].start_us+frames[best].stop_us)/2-time)?i:best,0);
   $('rdbeam').value=String(panel);$('frame').value=String(selected);renderRD();
   const range=y/Math.sin(manifest.geometry.beam_elevation_deg[panel]*Math.PI/180);
   const offset=Math.max(0,Math.floor((time-frames[selected].start_us)/120e6)*120);
   imageSelection(frames[selected].key,range,Math.min(offset,Math.max(0,Math.floor((frames[selected].duration_s-1)/120)*120)),frames[selected].spectra_ready!==false);
  }});
  plots=rtiPlot.panels.map(p=>({canvas:p.canvas,view:p,owner:rtiPlot}));
 }else rtiPlot.setSpec(spec);
 const mountRD=async()=>{
 const [v0,v1]=manifest.velocity_limits_mps;
 if(!rdPlot){rdPlot=await InteractiveFieldPlot.mount($('rd'),{id:'rd',frameParameters:{day:dayNumber()},spec:{bounds:[v0,50,v1,120],xlabel:'Radial velocity (m/s)',ylabel:'Nominal altitude (km)',colorbar_label:'SNR (dB)',panels:['Spectrum']},colorLimits:[+$('low').value,+$('high').value],onSelect:({y})=>imageSelection(frames[selected].key,y/Math.sin(manifest.geometry.beam_elevation_deg[+$('rdbeam').value]*Math.PI/180),0,frames[selected].spectra_ready!==false)});rd={view:rdPlot.panels[0],owner:rdPlot};}
 rdPlot.setFrameParameters({day:dayNumber()});
 };
 const mountTD=async()=>{
 const [v0,v1]=manifest.velocity_limits_mps;
 const tdSpec={bounds:[0,v0,Math.max(hours,.001),v1],xlabel:'Time (UTC)',ylabel:'Radial velocity (m/s)',colorbar_label:'SNR (dB)',panels:names,utc_origin_seconds:origin};
 if(!tdPlot){tdPlot=await InteractiveFieldPlot.mount($('timeDoppler'),{id:'td',frameParameters:{day:dayNumber()},spec:tdSpec,colorLimits:[+$('low').value,+$('high').value],onSelect:({panel,x})=>{
  const time=axisStart()+x*3.6e9;
  selected=frames.reduce((best,f,i)=>Math.abs((f.start_us+f.stop_us)/2-time)<Math.abs((frames[best].start_us+frames[best].stop_us)/2-time)?i:best,0);
  $('rdbeam').value=String(panel);$('frame').value=String(selected);renderRD();
 }});}else tdPlot.setSpec(tdSpec);
 tdPlot.setFrameParameters({day:dayNumber()});
 };
 // Paint the visible day overview first. Below-the-fold spectra mount only
 // when approached, rather than delaying the whole page behind their frames.
 for(const [root,mount,ready] of [[$('rd'),mountRD,rdPlot],[$('timeDoppler'),mountTD,tdPlot]]){
  secondaryMounts.set(root,mount);
  if(ready)await mount();else visiblePlots.observe(root);
 }
}
function scale(p){const [low,high]=p.owner.getColorLimits();p.view.set_scale(low,high,0,.1,2);}
function summaries(frame,beam,nr){
 const cached=summaryCache.get(frame);if(cached?.[beam])return cached[beam];
  const peak=new Float32Array(nr),sum=new Float64Array(nr),vsum=new Float64Array(nr),v2sum=new Float64Array(nr),vel=manifest.geometry.velocity_mps;
  const nyquist=manifest.velocity_limits_mps[1];
  frame.index.forEach(([b,f,r],i)=>{if(b!==beam)return;const w=10**(frame.snr_db[i]/10),angle=Math.PI*vel[f]/nyquist;peak[r]=Math.max(peak[r],frame.snr_db[i]);sum[r]+=w;vsum[r]+=w*Math.cos(angle);v2sum[r]+=w*Math.sin(angle);});
  const moments=new Float32Array(nr*4),variance=new Float64Array(nr);for(let r=0;r<nr;r++){moments[4*r]=peak[r];if(sum[r])moments[4*r+1]=Math.atan2(v2sum[r],vsum[r])*nyquist/Math.PI;}
  frame.index.forEach(([b,f,r],i)=>{if(b!==beam)return;const w=10**(frame.snr_db[i]/10),period=2*nyquist,dv=((vel[f]-moments[4*r+1]+nyquist)%period+period)%period-nyquist;variance[r]+=w*dv*dv;});
  for(let r=0;r<nr;r++)if(sum[r])moments[4*r+2]=Math.sqrt(variance[r]/sum[r]);
  const result={peak,moments},cache=cached||[];cache[beam]=result;summaryCache.set(frame,cache);return result;
}
function renderRTI(){
  if(!rtiPlot||!frames.length)return;
  const style={day:dayNumber(),hsv:$('display').value==='hsv'?1:0,velocity_limit:+$('vlim').value,width_limit:+$('wlim').value};
  rtiPlot.setFrameParameters(style);tdPlot?.setFrameParameters(style);
  const geom=manifest.geometry,nr=geom.slant_range_km.length,first=axisStart(),columns=[],xb=[];
  for(let i=0;i<frames.length;i++){
    const f=frames[i],previous=frames[i-1],next=frames[i+1];
    const left=previous&&f.start_us-previous.stop_us<=600e6?(previous.stop_us+f.start_us)/2:f.start_us;
    const right=next&&next.start_us-f.stop_us<=600e6?(f.stop_us+next.start_us)/2:f.stop_us;
    if(xb.length&&left>xb.at(-1)+1){columns.push(null);xb.push(left);}
    if(!xb.length)xb.push(left);columns.push(f);xb.push(right);
  }
  const x=Float32Array.from(xb,t=>(t-first)/3600e6),nx=columns.length;
  const nativeNf=geom.velocity_mps.length,nf=nativeNf+1;
  tdPlot?.panels.forEach((panel,beam)=>{
   const value=new Float32Array(nx*nf),moment=new Float32Array(nx*nf*4);
   columns.forEach((frame,c)=>{
    if(!frame){for(let f=0;f<nf;f++){value[f*nx+c]=NaN;moment[4*(f*nx+c)]=NaN;}return;}
    const summary=summaries(frame,beam,nr);
    frame.index.forEach(([b,f,r],i)=>{if(b===beam&&geom.altitude_km[beam][r]>=50&&geom.altitude_km[beam][r]<=120){
      const cell=f*nx+c,snr=frame.snr_db[i];if(snr>value[cell]){
       value[cell]=snr;moment[4*cell]=snr;moment[4*cell+1]=geom.velocity_mps[f];moment[4*cell+2]=summary.moments[4*r+2];
       if(f===0){const last=nativeNf*nx+c;value[last]=snr;moment.set(moment.subarray(4*cell,4*cell+4),4*last);moment[4*last+1]=manifest.velocity_limits_mps[1];}
      }
    }});
   });
   panel.set_rectilinear_data(x,edges([...geom.velocity_mps,manifest.velocity_limits_mps[1]]),value);panel.set_scale(...tdPlot.getColorLimits());
   if(style.hsv)panel.set_moment_data(moment,...(tdPlot.getHsvLimits()||[+$('vlim').value,+$('wlim').value]));
  });
  plots.forEach((p,beam)=>{
    const val=new Float32Array(nx*nr),moment=new Float32Array(nx*nr*4);
    columns.forEach((frame,c)=>{const s=frame?summaries(frame,beam,nr):null;
      for(let r=0;r<nr;r++){val[r*nx+c]=s?s.peak[r]:NaN;for(let k=0;k<4;k++)moment[4*(r*nx+c)+k]=s?s.moments[4*r+k]:(k===0?NaN:0);}});
    p.view.set_rectilinear_data(x,edges(geom.altitude_km[beam]),val);scale(p);
    if($('display').value==='hsv')p.view.set_moment_data(moment,...(rtiPlot.getHsvLimits()||[+$('vlim').value,+$('wlim').value]));
  });
  const hsv=$('display').value==='hsv';$('colorLabel').textContent=hsv?'Hue: radial velocity; saturation: spectral width; brightness: dB SNR. Click the HSV key to adjust.':'Peak excess signal/noise (dB). Click a colorbar to adjust its limits.';
}
function renderRD(){
  if(!rd||!frames.length)return;
  const frame=frames[selected],beam=+$('rdbeam').value,g=manifest.geometry,nr=g.slant_range_km.length;
  const bins=g.velocity_mps.map((v,i)=>i),nx=bins.length+1;
  const lookup=new Map(bins.map((f,i)=>[f,i])),value=new Float32Array(nx*nr),moment=new Float32Array(nx*nr*4),summary=summaries(frame,beam,nr);let count=0,transients=0;
  frame.index.forEach(([b,f,r],i)=>{if(b===beam&&lookup.has(f)){
    const cell=r*nx+lookup.get(f);value[cell]=frame.snr_db[i];moment[4*cell]=frame.snr_db[i];moment[4*cell+1]=g.velocity_mps[f];moment[4*cell+2]=summary.moments[4*r+2];
    if(f===0){const last=r*nx+nx-1;value[last]=value[cell];moment.set(moment.subarray(4*cell,4*cell+4),4*last);moment[4*last+1]=manifest.velocity_limits_mps[1];}
    count++;if(frame.look_dominance?.[i]>.5)transients++;
  }});
  rd.view.set_rectilinear_data(edges([...g.velocity_mps,manifest.velocity_limits_mps[1]]),edges(g.altitude_km[beam]),value);scale(rd);
  const hsv=$('display').value==='hsv';rdPlot.setFrameParameters({day:dayNumber(),hsv:hsv?1:0,velocity_limit:+$('vlim').value,width_limit:+$('wlim').value});
  if(hsv)rd.view.set_moment_data(moment,...(rdPlot.getHsvLimits()||[+$('vlim').value,+$('wlim').value]));
  $('frameStatus').textContent=`${utc(frame.start_us)} · ${names[beam]} · ${frame.duration_s.toFixed(2)} s mode, ${frame.coverage_s.toFixed(2)} s valid support · ${frame.looks} looks · ${count} retained cells · 16 ms complementary pairs · ${frame.quality}`;
}
async function refresh({useIndex=false}={}){
  const request=++refreshSequence;
  const update=(stage,progress)=>{if(request===refreshSequence)loading(stage,progress);};
  update('Downloading observations…');
  let m=catalogue;
  if(!useIndex||!m){
   const response=initialResponse?await initialResponse:await fetchData('./data/manifest.json');initialResponse=null;
   if(response.error)throw response.error;
   m=await readJSON(response,p=>update('Downloading observations…',p));
   catalogue=m;
  }
  if(request!==refreshSequence)return;
  renderProcessingStatus(m);
  const days=m.available_days||[];
  if(!chosenDay||!days.includes(chosenDay))chosenDay=days.at(-1)||m.selected_day;
  const control=$('day');control.replaceChildren(...days.map(day=>{const option=document.createElement('option');option.value=day;option.textContent=day;return option;}));control.value=chosenDay;
  if(skyNoise)renderNoise();else if(m.noise_model_status)renderNoise(m);
  if($('skyNoiseSection').open)void refreshSkyNoise();
  if(m.transfer){const a=m.transfer.archive,t=m.transfer.syowa_to_revontuli;
    $('transfer').textContent=`Delivered native archive: ${(a.bytes/2**20).toFixed(2)} MiB · ${a.records} native records · ${a.files} HDF5 files. `+
      (a.dense_complex_sample_retained_fraction!=null?`Complex spectra retain ${(100*a.dense_complex_sample_retained_fraction).toPrecision(3)}% of the full dense spectra (${a.retained_complex_samples.toLocaleString()} / ${a.dense_complex_samples.toLocaleString()} complex samples; all beams, Doppler bins, range gates, receiver channels and coherent looks, before altitude, velocity and detection masks). `:'Complex-spectra retained fraction: awaiting complete dense-spectrum counts. ')+
      'Storage: complex64 spectra (8 bytes/sample); float32 noise power and uncertainty (4 bytes/value), with lossless compression. '+
      (a.first_mode_us?`${utc(a.first_mode_us)} through ${utc(a.latest_mode_us)}. `:'')+
      (a.projected_bytes_per_day===null?'Daily growth: interval too short for extrapolation. ':`Provisional growth: ${(a.projected_bytes_per_day/2**20).toFixed(1)} MiB/day (${(a.projected_bytes_per_day/86400/1024).toFixed(2)} KiB/s average). `)+
      (t?`Latest Syowa replication: ${(t.literal_bytes/2**20).toFixed(2)} MiB in ${t.elapsed_seconds.toFixed(1)} s (${(t.payload_bytes_per_second/1024).toFixed(2)} KiB/s payload; cap ${t.bandwidth_cap_kib_s} KiB/s). `:'Transfer measurement pending. ')+
      'Web copy: seven-day retention; Revontuli archive retained.';}
  renderArchivePlot(m);
  if(m.selected_day!==chosenDay||!m.frame_packets){
   update(`Downloading ${chosenDay}…`);
   // Prepared static packets support gzip/conditional requests without involving
   // the API worker. Fall back only when publication has not prepared this day.
   let daily=await fetchData(`./data/days/${encodeURIComponent(chosenDay)}.json`);
   if(daily.status===404)daily=await fetchData(`./api/day?day=${encodeURIComponent(chosenDay)}`);
   m=await readJSON(daily,p=>update(`Downloading ${chosenDay}…`,p));
  }
  if(request!==refreshSequence)return;
  saveView(m);
  const age=(Date.now()-m.latest_mode_us/1000)/60000;
  $('status').textContent=m.frame_count?`${m.frame_count} observing modes · latest ${utc(m.latest_mode_us)} · ${Math.max(0,age).toFixed(0)} min behind now`:`${chosenDay} UTC · no observations delivered yet`;
  if(m.spectra_transfer_pending)$('status').textContent+=` · complex spectra transferring for ${m.spectra_transfer_pending} modes`;
  const download=$('downloadSpectra');
  download.href=`./api/download?day=${encodeURIComponent(chosenDay)}`;
  download.setAttribute('download',`pansy-sparse-spectra-${chosenDay}.tar`);
  download.setAttribute('aria-disabled',m.frame_count?'false':'true');
  download.title=`Native complex spectra and metadata for ${chosenDay} UTC (.tar)`;
  $('status').className=age>10?'status warn':'status';$('error').textContent='';
  if(m.revision===revision){loaded();return;}manifest=m;
  update('Preparing echo data…');
  originalFrames=m.frame_packets;applyEchoFilter();revision=m.revision;selected=frames.length-1;$('frame').max=String(Math.max(0,selected));$('frame').value=String(Math.max(0,selected));
  if(!frames.length){$('status').textContent=`${chosenDay} UTC · no complete 16 ms records yet`;loaded();return;}
  update('Drawing Rustplotlib axes…');
  await new Promise(resolve=>requestAnimationFrame(resolve));
  await mountPlots();
  if(request!==refreshSequence)return;
  update('Drawing echoes…');
  $('sidecar').href='./data/'+m.sidecar;renderRTI();renderRD();
  await rtiPlot.requestFrame();
  if(request!==refreshSequence)return;
  requestAnimationFrame(()=>{if(request===refreshSequence){const now=String(Math.round(performance.now()));if(!document.body.dataset.plotsReadyMs)document.body.dataset.plotsReadyMs=now;document.body.dataset.latestRefreshMs=now;loaded();}});
}
function renderArchivePlot(m){
  archiveInventory=m.transfer?.archive;
  if($('archiveSection').open)void drawArchive();
}
async function drawArchive(){
  try{await renderArchiveDiagnostics(archiveInventory);$('archivePlotStatus').textContent='Numeric inventory rendered in this browser by Rustplotlib. Zoom and pan each plot; no plot images are downloaded.';}
  catch(error){$('archivePlotStatus').textContent=String(error);}
}
function renderProcessingStatus(m){
  const s=m.source_processing,node=$('processingStatus');node.hidden=!s;
  if(!s)return;
  const b=s.backfill||{},live=s.live||{},parts=[];
  if(live.last_written_stop_us)parts.push(`Source live processing (${live.state}): ${utc(live.last_written_stop_us)}.`);
  if(b.current_index&&b.total_modes){
    parts.push(`Rectangular backfill (${b.complete?'complete':b.state}): ${b.current_index} / ${b.total_modes} modes; working on ${utc(b.current_start_us)}.`);
    const frontier=new Date(b.current_start_us/1000).toISOString().slice(0,10);
    const delivered=Object.entries(m.transfer?.archive?.days||{}).filter(([d,v])=>d<=frontier&&v.latest_mode_us).sort(([a],[b])=>a.localeCompare(b)).at(-1);
    parts.push(delivered?`Historical voltages delivered to viewer through ${utc(delivered[1].latest_mode_us)}.`:'Historical voltage delivery pending.');
  }
  const age=Math.max(0,(Date.now()/1000-s.generated_unix_s)/60);
  parts.push(`Processing status updated ${age.toFixed(0)} min ago; processing and transfer are separate.`);
  if(m.plot_generation_status==='delayed')parts.push('Optional plot refresh delayed; native observations continue updating.');
  node.textContent=parts.join(' ');
}
async function refreshSkyNoise(){
  if(noiseRequest)return noiseRequest;
  noiseRequest=(async()=>{
    try{
      const response=await fetchData('./data/sky-noise/status.json');
      if(!response.ok)throw Error(`Sky-noise diagnostics unavailable (${response.status})`);
      const value=await response.json();
      if(value.noise_model_status?.length!==7||!value.noise_model_data?.times_h?.length)throw Error('Numeric sky-noise data have not arrived yet');
      skyNoise=value;renderNoise();
    }catch(error){
      if(skyNoise){skyNoise.noise_plot={...skyNoise.noise_plot,status:'delayed'};renderNoise();}
      else{$('noiseStatus').textContent='Numeric sky-noise data are not available yet; the independent analysis job will retry. Observations remain available.';}
    }
  })();
  try{await noiseRequest;}finally{noiseRequest=null;}
}
function renderNoise(model=skyNoise||manifest){
  const module=+$('noiseModule').value,fit=model?.noise_model_status?.[module];
  if(!fit){void refreshSkyNoise();return;}
  $('noiseStatus').textContent=`PyGDSM GSM2008 · one-way receive pattern · 90% crossed-Yagi / 10% isotropic element. Module ${module}: `+
    (fit.receiver_temp_k===null?fit.status:`T_receiver = ${fit.receiver_temp_k.toFixed(0)} K; fractional RMS ${(fit.fractional_rms*100).toFixed(1)}%; ${fit.status}`)+
    '. Top: measured noise power; bottom: modeled sky temperature. Fits update hourly and require ≥12 hours of data.';
  if(model.noise_model_observation_stop_unix_s)$('noiseStatus').textContent+=` Measurements through ${utc(model.noise_model_observation_stop_unix_s*1e6)}.`;
  if(model.noise_plot?.status==='delayed')$('noiseStatus').textContent+=' Refresh delayed; using cached numeric data.';
  if($('skyNoiseSection').open)void renderNoiseDiagnostics(model,module).catch(error=>{$('noiseStatus').textContent+=' '+String(error);});
}
function failure(e){revision='';loading('Loading failed — retrying…');$('error').textContent=e.message;console.error(e);}
try{
 loading('Loading plotting engine…');
 await init();if(new URLSearchParams(location.search).has('embed'))document.body.classList.add('embedded');
 document.body.dataset.display=$('display').value;
 for(const id of ['display','low','high','vlim','wlim'])$(id).onchange=()=>{try{
  document.body.dataset.display=$('display').value;
  if(id==='low'||id==='high')for(const plot of [rtiPlot,rdPlot,tdPlot])plot?.setColorLimits(+$('low').value,+$('high').value);
  if(id==='vlim'||id==='wlim')for(const plot of [rtiPlot,rdPlot,tdPlot])plot?.setHsvLimits(+$('vlim').value,+$('wlim').value);
  renderRTI();renderRD();}catch(e){failure(e);}};
 $('rdbeam').onchange=renderRD;$('frame').oninput=()=>{selected=+$('frame').value;renderRD();};
 $('hideSingleMode').onchange=()=>{applyEchoFilter();renderRTI();renderRD();};
 $('noiseModule').onchange=()=>renderNoise();
 $('skyNoiseSection').ontoggle=()=>{if($('skyNoiseSection').open)void refreshSkyNoise();};
 $('archiveSection').ontoggle=()=>{if($('archiveSection').open)void drawArchive();};
 $('day').onchange=()=>{showLoading=true;loading('Loading selected day…');chosenDay=$('day').value;revision='';$('downloadSpectra').setAttribute('aria-disabled','true');$('downloadSpectra').removeAttribute('href');if(playing){clearInterval(playing);playing=null;$('play').textContent='Play';}refresh({useIndex:true}).catch(failure);};
 $('play').onclick=()=>{if(playing){clearInterval(playing);playing=null;$('play').textContent='Play';}else{$('play').textContent='Pause';playing=setInterval(()=>{selected=(selected+1)%frames.length;$('frame').value=selected;renderRD();},500);}};
 let resizeTimer;window.addEventListener('resize',()=>{clearTimeout(resizeTimer);resizeTimer=setTimeout(()=>{if(frames.length){renderRTI();renderRD();}},150);});
 // A failed first request must not prevent all future refreshes. Schedule after
 // completion so slow requests cannot accumulate overlapping polling work.
 async function pollObservations(){
  let delay=60000;
  try{await refresh();}catch(e){failure(e);delay=5000;}
  finally{setTimeout(pollObservations,delay);}
 }
 await restoreView();
 await pollObservations();
}catch(e){failure(e);}
