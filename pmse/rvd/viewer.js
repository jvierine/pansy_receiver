// Scientific provenance: rustplotlib interactive API d3837b8; numerical imaging is Rust/WASM.
import init from './pkg/rustplotlib_web.js';
import {imageSelection} from './imaging.js?v=sparse-2';
import {InteractiveFieldPlot} from './lib/interactive_plot.js?v=fast-start-2';
const $=id=>document.getElementById(id), names=['Zenith','North','East','South','West'];
// Start data I/O while WebAssembly initializes.
let initialResponse=fetch('./data/manifest.json',{cache:'no-cache'}).catch(error=>({error}));
let manifest,frames=[],plots=[],rd,playing=null,selected=0,revision='',rtiPlot,rdPlot,tdPlot;
let chosenDay='',refreshSequence=0;
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
   imageSelection(frames[selected].key,range,Math.min(offset,Math.max(0,Math.floor((frames[selected].duration_s-1)/120)*120)));
  }});
  plots=rtiPlot.panels.map(p=>({canvas:p.canvas,view:p,owner:rtiPlot}));
 }else rtiPlot.setSpec(spec);
 const mountRD=async()=>{
 const [v0,v1]=manifest.velocity_limits_mps;
 if(!rdPlot){rdPlot=await InteractiveFieldPlot.mount($('rd'),{id:'rd',frameParameters:{day:dayNumber()},spec:{bounds:[v0,50,v1,120],xlabel:'Radial velocity (m/s)',ylabel:'Nominal altitude (km)',colorbar_label:'SNR (dB)',panels:['Spectrum']},colorLimits:[+$('low').value,+$('high').value],onSelect:({y})=>imageSelection(frames[selected].key,y/Math.sin(manifest.geometry.beam_elevation_deg[+$('rdbeam').value]*Math.PI/180))});rd={view:rdPlot.panels[0],owner:rdPlot};}
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
  const peak=new Float32Array(nr),sum=new Float64Array(nr),vsum=new Float64Array(nr),v2sum=new Float64Array(nr),vel=manifest.geometry.velocity_mps;
  const nyquist=manifest.velocity_limits_mps[1];
  frame.index.forEach(([b,f,r],i)=>{if(b!==beam)return;const w=10**(frame.snr_db[i]/10),angle=Math.PI*vel[f]/nyquist;peak[r]=Math.max(peak[r],frame.snr_db[i]);sum[r]+=w;vsum[r]+=w*Math.cos(angle);v2sum[r]+=w*Math.sin(angle);});
  const moments=new Float32Array(nr*4),variance=new Float64Array(nr);for(let r=0;r<nr;r++){moments[4*r]=peak[r];if(sum[r])moments[4*r+1]=Math.atan2(v2sum[r],vsum[r])*nyquist/Math.PI;}
  frame.index.forEach(([b,f,r],i)=>{if(b!==beam)return;const w=10**(frame.snr_db[i]/10),period=2*nyquist,dv=((vel[f]-moments[4*r+1]+nyquist)%period+period)%period-nyquist;variance[r]+=w*dv*dv;});
  for(let r=0;r<nr;r++)if(sum[r])moments[4*r+2]=Math.sqrt(variance[r]/sum[r]);
  return {peak,moments};
}
function renderRTI(){
  if(!rtiPlot||!frames.length)return;
  rtiPlot.setFrameParameters({day:dayNumber(),hsv:$('display').value==='hsv'?1:0,velocity_limit:+$('vlim').value,width_limit:+$('wlim').value});
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
   const value=new Float32Array(nx*nf);
   columns.forEach((frame,c)=>{
    if(!frame){for(let f=0;f<nf;f++)value[f*nx+c]=NaN;return;}
    frame.index.forEach(([b,f,r],i)=>{if(b===beam&&geom.altitude_km[beam][r]>=50&&geom.altitude_km[beam][r]<=120){value[f*nx+c]=Math.max(value[f*nx+c],frame.snr_db[i]);if(f===0)value[nativeNf*nx+c]=value[c];}});
   });
   panel.set_rectilinear_data(x,edges([...geom.velocity_mps,manifest.velocity_limits_mps[1]]),value);panel.set_scale(...tdPlot.getColorLimits());
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
  const lookup=new Map(bins.map((f,i)=>[f,i])),value=new Float32Array(nx*nr);let count=0,transients=0;
  frame.index.forEach(([b,f,r],i)=>{if(b===beam&&lookup.has(f)){value[r*nx+lookup.get(f)]=frame.snr_db[i];if(f===0)value[r*nx+nx-1]=frame.snr_db[i];count++;if(frame.look_dominance?.[i]>.5)transients++;}});
  rd.view.set_rectilinear_data(edges([...g.velocity_mps,manifest.velocity_limits_mps[1]]),edges(g.altitude_km[beam]),value);scale(rd);
  $('frameStatus').textContent=`${utc(frame.start_us)} · ${names[beam]} · ${frame.duration_s.toFixed(2)} s mode, ${frame.coverage_s.toFixed(2)} s valid support · ${frame.looks} looks · ${count} retained cells · 16 ms complementary pairs · ${frame.quality}`;
}
async function refresh(){
  const request=++refreshSequence;
  const response=initialResponse?await initialResponse:await fetch('./data/manifest.json',{cache:'no-cache'});initialResponse=null;
  if(response.error)throw response.error;
  if(!response.ok)throw Error(`Data unavailable (${response.status})`);
  let m=await response.json();
  renderProcessingStatus(m);
  const days=m.available_days||[];
  if(!chosenDay||!days.includes(chosenDay))chosenDay=days.at(-1)||m.selected_day;
  const control=$('day');control.replaceChildren(...days.map(day=>{const option=document.createElement('option');option.value=day;option.textContent=day;return option;}));control.value=chosenDay;
  if(m.noise_model_status)renderNoise(m);
  if(m.transfer){const a=m.transfer.archive,t=m.transfer.syowa_to_revontuli;
    $('transfer').textContent=`Revontuli native archive: ${(a.bytes/2**20).toFixed(2)} MiB · ${a.records} native records · ${a.files} HDF5 files. `+
      (a.dense_complex_sample_retained_fraction!=null?`Complex spectra retain ${(100*a.dense_complex_sample_retained_fraction).toPrecision(3)}% of the full dense spectra (${a.retained_complex_samples.toLocaleString()} / ${a.dense_complex_samples.toLocaleString()} complex samples; all beams, Doppler bins, range gates, receiver channels and coherent looks, before altitude, velocity and detection masks). `:'Complex-spectra retained fraction: awaiting complete dense-spectrum counts. ')+
      'Storage: complex64 spectra (8 bytes/sample); float32 noise power and uncertainty (4 bytes/value), with lossless compression. '+
      (a.first_mode_us?`${utc(a.first_mode_us)} through ${utc(a.latest_mode_us)}. `:'')+
      (a.projected_bytes_per_day===null?'Daily growth: interval too short for extrapolation. ':`Provisional growth: ${(a.projected_bytes_per_day/2**20).toFixed(1)} MiB/day (${(a.projected_bytes_per_day/86400/1024).toFixed(2)} KiB/s average). `)+
      (t?`Latest Syowa replication: ${(t.literal_bytes/2**20).toFixed(2)} MiB in ${t.elapsed_seconds.toFixed(1)} s (${(t.payload_bytes_per_second/1024).toFixed(2)} KiB/s payload; cap ${t.bandwidth_cap_kib_s} KiB/s). `:'Transfer measurement pending. ')+
      'Web copy: seven-day retention; Revontuli archive retained.';}
  $('archivePlot').hidden=m.plot_generation_status&&m.plot_generation_status!=='ready';
  if(!($('archivePlot').hidden)&&m.transfer?.archive?.hourly_unix_s?.length){const width=document.body.clientWidth<900?640:960;
    $('archivePlot').src=`./data/axes/archive-growth-${width}.svg?v=${m.generated_utc}`;$('archivePlot').style.width=width+'px';$('archivePlot').hidden=false;}
  if(m.selected_day!==chosenDay||!m.frame_packets){
   const daily=await fetch(`./api/day?day=${encodeURIComponent(chosenDay)}`,{cache:'no-cache'});if(!daily.ok)throw Error('Selected UTC day is unavailable');
   m=await daily.json();
  }
  if(request!==refreshSequence)return;
  const age=(Date.now()-m.latest_mode_us/1000)/60000;
  $('status').textContent=m.frame_count?`${m.frame_count} observing modes · latest ${utc(m.latest_mode_us)} · ${Math.max(0,age).toFixed(0)} min behind now`:`${chosenDay} UTC · no observations delivered yet`;
  const download=$('downloadSpectra');
  download.href=`./api/download?day=${encodeURIComponent(chosenDay)}`;
  download.setAttribute('download',`pansy-sparse-spectra-${chosenDay}.tar`);
  download.setAttribute('aria-disabled',m.frame_count?'false':'true');
  download.title=`Native complex spectra and metadata for ${chosenDay} UTC (.tar)`;
  $('status').className=age>10?'status warn':'status';$('error').textContent='';
  if(m.revision===revision)return;manifest=m;
  frames=m.frame_packets;revision=m.revision;selected=frames.length-1;$('frame').max=String(Math.max(0,selected));$('frame').value=String(Math.max(0,selected));
  if(!frames.length){$('status').textContent=`${chosenDay} UTC · no complete 16 ms records yet`;return;}
  await mountPlots();
  $('sidecar').href='./data/'+m.sidecar;renderRTI();renderRD();
  requestAnimationFrame(()=>{document.body.dataset.plotsReadyMs=String(Math.round(performance.now()));});
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
function renderNoise(model=manifest){
  const module=+$('noiseModule').value,fit=model.noise_model_status[module];
  $('noiseStatus').textContent=`PyGDSM GSM2008 · one-way receive pattern · 90% crossed-Yagi / 10% isotropic element. Module ${module}: `+
    (fit.receiver_temp_k===null?fit.status:`T_receiver = ${fit.receiver_temp_k.toFixed(0)} K; fractional RMS ${(fit.fractional_rms*100).toFixed(1)}%; ${fit.status}`)+
    '. Top: measured noise power; bottom: modeled sky temperature. Fits update hourly and require ≥12 hours of data.';
  if(model.plot_generation_status&&model.plot_generation_status!=='ready'){
    $('noisePlot').hidden=true;$('noiseStatus').textContent+=' Plot refresh delayed; native observations remain current.';return;
  }
  const width=document.body.clientWidth<900?640:960;
  $('noisePlot').src=`./data/axes/noise-module-${module}-${width}.svg?v=${model.noise_model_updated_unix_s}`;
  $('noisePlot').style.width=width+'px';$('noisePlot').hidden=false;
}
function failure(e){revision='';$('error').textContent=e.message;console.error(e);}
try{
 await init();if(new URLSearchParams(location.search).has('embed'))document.body.classList.add('embedded');
 document.body.dataset.display=$('display').value;
 for(const id of ['display','low','high','vlim','wlim'])$(id).onchange=()=>{try{
  document.body.dataset.display=$('display').value;
  if(id==='low'||id==='high')for(const plot of [rtiPlot,rdPlot,tdPlot])plot?.setColorLimits(+$('low').value,+$('high').value);
  if(id==='vlim'||id==='wlim')rtiPlot?.setHsvLimits(+$('vlim').value,+$('wlim').value);
  renderRTI();renderRD();}catch(e){failure(e);}};
 $('rdbeam').onchange=renderRD;$('frame').oninput=()=>{selected=+$('frame').value;renderRD();};
 $('noiseModule').onchange=()=>renderNoise();
 $('day').onchange=()=>{chosenDay=$('day').value;revision='';$('downloadSpectra').setAttribute('aria-disabled','true');$('downloadSpectra').removeAttribute('href');if(playing){clearInterval(playing);playing=null;$('play').textContent='Play';}refresh().catch(failure);};
 $('play').onclick=()=>{if(playing){clearInterval(playing);playing=null;$('play').textContent='Play';}else{$('play').textContent='Pause';playing=setInterval(()=>{selected=(selected+1)%frames.length;$('frame').value=selected;renderRD();},500);}};
 window.addEventListener('resize',()=>{if(frames.length){renderRTI();renderRD();}});
 await refresh();setInterval(()=>refresh().catch(failure),60000);
}catch(e){failure(e);}
