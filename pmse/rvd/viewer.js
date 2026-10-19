// Scientific provenance: rustplotlib interactive API d3837b8; numerical imaging is Rust/WASM.
import init from './pkg/rustplotlib_web.js';
import {imageSelection} from './imaging.js';
import {InteractiveFieldPlot} from './lib/interactive_plot.js';
const $=id=>document.getElementById(id), names=['Zenith','North','East','South','West'];
let manifest,frames=[],plots=[],rd,playing=null,selected=0,revision='',rtiPlot,rdPlot;
const utc=us=>new Date(us/1000).toISOString().replace('T',' ').slice(0,19)+' UTC';
const edges=centers=>Float32Array.from([...centers.map((x,i)=>i?(centers[i-1]+x)/2:x-(centers[1]-x)/2),centers.at(-1)+(centers.at(-1)-centers.at(-2))/2]);
async function mountPlots(){
 const hours=(manifest.latest_mode_us-manifest.first_mode_us)/3.6e9,origin=manifest.first_mode_us/1e6;
 const spec={bounds:[0,50,Math.max(hours,.001),120],xlabel:'Time (UTC)',ylabel:'Nominal altitude (km)',colorbar_label:'SNR (dB)',panels:names,utc_origin_seconds:origin};
 if(!rtiPlot){
  rtiPlot=await InteractiveFieldPlot.mount($('rti'),{id:'rti',spec,colorLimits:[+$('low').value,+$('high').value],onSelect:({panel,x,y})=>{
   const time=manifest.first_mode_us+x*3.6e9;
   selected=frames.reduce((best,f,i)=>Math.abs((f.start_us+f.stop_us)/2-time)<Math.abs((frames[best].start_us+frames[best].stop_us)/2-time)?i:best,0);
   $('rdbeam').value=String(panel);$('frame').value=String(selected);renderRD();
   const range=y/Math.sin(manifest.geometry.beam_elevation_deg[panel]*Math.PI/180);
   const offset=Math.max(0,Math.floor((time-frames[selected].start_us)/120e6)*120);
   imageSelection(frames[selected].key,range,Math.min(offset,Math.max(0,Math.floor((frames[selected].duration_s-1)/120)*120)));
  }});
  plots=rtiPlot.panels.map(p=>({canvas:p.canvas,view:p}));
 }else rtiPlot.setSpec(spec);
 const ve=edges(manifest.geometry.velocity_mps);
 if(!rdPlot){rdPlot=await InteractiveFieldPlot.mount($('rd'),{id:'rd',spec:{bounds:[ve[0],50,ve.at(-1),120],xlabel:'Radial velocity (m/s)',ylabel:'Nominal altitude (km)',colorbar_label:'SNR (dB)',panels:['Spectrum']},colorLimits:[+$('low').value,+$('high').value],onSelect:({y})=>imageSelection(frames[selected].key,y/Math.sin(manifest.geometry.beam_elevation_deg[+$('rdbeam').value]*Math.PI/180))});rd={view:rdPlot.panels[0]};}
}
function scale(p){const low=+$('low').value,high=+$('high').value;if(!(high>low))throw Error('SNR high must exceed low');p.view.set_scale(low,high,0,.1,2);}
function summaries(frame,beam,nr){
  const peak=new Float32Array(nr),sum=new Float64Array(nr),vsum=new Float64Array(nr),v2sum=new Float64Array(nr),vel=manifest.geometry.velocity_mps;
  frame.index.forEach(([b,f,r],i)=>{if(b!==beam)return;const w=10**(frame.snr_db[i]/10),v=vel[f];peak[r]=Math.max(peak[r],frame.snr_db[i]);sum[r]+=w;vsum[r]+=w*v;v2sum[r]+=w*v*v;});
  const moments=new Float32Array(nr*4);for(let r=0;r<nr;r++){moments[4*r]=peak[r];if(sum[r]){moments[4*r+1]=vsum[r]/sum[r];moments[4*r+2]=Math.sqrt(Math.max(0,v2sum[r]/sum[r]-(vsum[r]/sum[r])**2));}}
  return {peak,moments};
}
function renderRTI(){
  const geom=manifest.geometry,nr=geom.slant_range_km.length,first=frames[0].start_us,columns=[],xb=[];
  for(let i=0;i<frames.length;i++){
    const f=frames[i],previous=frames[i-1],next=frames[i+1];
    const left=previous&&f.start_us-previous.stop_us<=600e6?(previous.stop_us+f.start_us)/2:f.start_us;
    const right=next&&next.start_us-f.stop_us<=600e6?(f.stop_us+next.start_us)/2:f.stop_us;
    if(xb.length&&left>xb.at(-1)+1){columns.push(null);xb.push(left);}
    if(!xb.length)xb.push(left);columns.push(f);xb.push(right);
  }
  const x=Float32Array.from(xb,t=>(t-first)/3600e6),nx=columns.length;
  plots.forEach((p,beam)=>{
    const val=new Float32Array(nx*nr),moment=new Float32Array(nx*nr*4);
    columns.forEach((frame,c)=>{const s=frame?summaries(frame,beam,nr):null;
      for(let r=0;r<nr;r++){val[r*nx+c]=s?s.peak[r]:NaN;for(let k=0;k<4;k++)moment[4*(r*nx+c)+k]=s?s.moments[4*r+k]:(k===0?NaN:0);}});
    p.view.set_rectilinear_data(x,edges(geom.altitude_km[beam]),val);scale(p);
    if($('display').value==='hsv')p.view.set_moment_data(moment,+$('vlim').value,+$('wlim').value);
  });
  const hsv=$('display').value==='hsv';$('colorLabel').textContent=hsv?`Hue −${$('vlim').value} to +${$('vlim').value} m/s; saturation decreases over 0–${$('wlim').value} m/s width; intensity ${$('low').value}–${$('high').value} dB`:`Peak excess signal/noise, ${$('low').value} to ${$('high').value} dB`;
}
function renderRD(){
  const frame=frames[selected],beam=+$('rdbeam').value,g=manifest.geometry,nr=g.slant_range_km.length;
  const bins=g.velocity_mps.map((v,i)=>i),nx=bins.length;
  const lookup=new Map(bins.map((f,i)=>[f,i])),value=new Float32Array(nx*nr);let count=0,transients=0;
  frame.index.forEach(([b,f,r],i)=>{if(b===beam&&lookup.has(f)){value[r*nx+lookup.get(f)]=frame.snr_db[i];count++;if(frame.look_dominance[i]>.5)transients++;}});
  rd.view.set_rectilinear_data(edges(bins.map(f=>g.velocity_mps[f])),edges(g.altitude_km[beam]),value);scale(rd);
  $('frameStatus').textContent=`${utc(frame.start_us)} · ${names[beam]} · ${frame.duration_s.toFixed(2)} s mode, ${frame.coverage_s.toFixed(2)} s valid support · ${frame.looks} looks · ${count} retained cells (${transients} dominated by one look) · ${frame.quality}`;
}
async function refresh(){
  const response=await fetch('./data/manifest.json',{cache:'no-cache'});if(!response.ok)throw Error(`Data unavailable (${response.status})`);
  const m=await response.json();
  if(m.noise_model_status){manifest=m;renderNoise();}
  if(m.transfer){const a=m.transfer.archive,t=m.transfer.syowa_to_revontuli;
    $('transfer').textContent=`Revontuli native archive: ${(a.bytes/2**20).toFixed(2)} MiB · ${a.records} native records · ${a.files} HDF5 files. `+
      (a.first_mode_us?`${utc(a.first_mode_us)} through ${utc(a.latest_mode_us)}. `:'')+
      (a.projected_bytes_per_day===null?'Daily growth: interval too short for extrapolation. ':`Provisional growth: ${(a.projected_bytes_per_day/2**20).toFixed(1)} MiB/day (${(a.projected_bytes_per_day/86400/1024).toFixed(2)} KiB/s average). `)+
      (t?`Latest Syowa replication: ${(t.literal_bytes/2**20).toFixed(2)} MiB in ${t.elapsed_seconds.toFixed(1)} s (${(t.payload_bytes_per_second/1024).toFixed(2)} KiB/s payload; cap ${t.bandwidth_cap_kib_s} KiB/s). `:'Transfer measurement pending. ')+
      'Web copy: seven-day retention; Revontuli archive retained.';}
  if(m.transfer?.archive?.hourly_unix_s?.length){const width=document.body.clientWidth<900?640:960;
    $('archivePlot').src=`./data/axes/archive-growth-${width}.svg?v=${m.generated_utc}`;$('archivePlot').style.width=width+'px';$('archivePlot').hidden=false;}
  if(m.revision===revision)return;
  const cache=new Map(frames.map(f=>[f.key,f]));manifest=m;
  // A bounded pool avoids issuing hundreds of simultaneous requests.
  const result=new Array(m.frames.length);let next=0;
  await Promise.all(Array.from({length:6},async()=>{for(;;){const i=next++;if(i>=result.length)break;const key=m.frames[i].key;if(cache.has(key)){result[i]=cache.get(key);continue;}const r=await fetch(`./data/frames/${key}.json`);if(!r.ok)throw Error('Incomplete frame transfer');result[i]=await r.json();}}));
  frames=result;revision=m.revision;selected=frames.length-1;$('frame').max=String(selected);$('frame').value=String(selected);
  await mountPlots();
  const age=(Date.now()-m.latest_mode_us/1000)/60000;$('status').textContent=`${frames.length} observing modes · latest completed ${utc(m.latest_mode_us)} · ${Math.max(0,age).toFixed(1)} min old`;$('status').className=age>30?'status warn':'status';
  $('sidecar').href='./data/'+m.sidecar;renderRTI();renderRD();
}
function renderNoise(){
  const module=+$('noiseModule').value,fit=manifest.noise_model_status[module];
  $('noiseStatus').textContent=`PyGDSM GSM2008 · one-way receive pattern · 90% crossed-Yagi / 10% isotropic element. Module ${module}: `+
    (fit.receiver_temp_k===null?fit.status:`T_receiver = ${fit.receiver_temp_k.toFixed(0)} K; fractional RMS ${(fit.fractional_rms*100).toFixed(1)}%; ${fit.status}`)+
    '. Top: measured noise power; bottom: modeled sky temperature. Fits update hourly and require ≥12 hours of data.';
  const width=document.body.clientWidth<900?640:960;
  $('noisePlot').src=`./data/axes/noise-module-${module}-${width}.svg?v=${manifest.noise_model_updated_unix_s}`;
  $('noisePlot').style.width=width+'px';$('noisePlot').hidden=false;
}
function failure(e){revision='';$('error').textContent=e.message;console.error(e);}
try{
 await init();if(new URLSearchParams(location.search).has('embed'))document.body.classList.add('embedded');
 for(const id of ['display','beam','low','high','vlim','wlim'])$(id).onchange=()=>{try{renderRTI();renderRD();}catch(e){failure(e);}};
 $('rdbeam').onchange=renderRD;$('frame').oninput=()=>{selected=+$('frame').value;renderRD();};
 $('noiseModule').onchange=renderNoise;
 $('play').onclick=()=>{if(playing){clearInterval(playing);playing=null;$('play').textContent='Play';}else{$('play').textContent='Pause';playing=setInterval(()=>{selected=(selected+1)%frames.length;$('frame').value=selected;renderRD();},500);}};
 window.addEventListener('resize',()=>{if(frames.length){renderRTI();renderRD();}});
 await refresh();setInterval(()=>refresh().catch(failure),60000);
}catch(e){failure(e);}
