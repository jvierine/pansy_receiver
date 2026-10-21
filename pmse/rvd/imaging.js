import {InteractiveFieldPlot} from './lib/interactive_plot.js?v=gpu-lifecycle-1';
const $=id=>document.getElementById(id);
let worker,geometry,plot,plotMount,sequence=0,previousRange,lastMoments,lastPacket,lastKey,echoFilter,busy=false,fetchAbort,timeout;
export function setEchoFilter(resolver){echoFilter=resolver;if(lastPacket)reconstruct();}
const grid=61;
function edges(extent){const d=2*extent/(grid-1);return Float32Array.from({length:grid+1},(_,i)=>-extent-d/2+i*d);}
function styleParameters(range){return {range_km:range,hsv:$('display').value==='hsv'?1:0,velocity_limit:+$('imageVelocity').value,width_limit:+$('imageWidth').value};}
function repaint(){if(plot&&lastMoments){
 const extent=previousRange*Math.sin(20*Math.PI/180);
 plot.panels[0].set_rectilinear_data(edges(extent),edges(extent),Float32Array.from({length:grid*grid},(_,i)=>lastMoments[4*i]));
 if($('display').value==='hsv')plot.panels[0].set_moment_data(lastMoments,...(plot.getHsvLimits()||[+$('imageVelocity').value,+$('imageWidth').value]));
 plot.setFrameParameters(styleParameters(previousRange));
}}
function options(){return {method:$('sparseImaging').checked?$('imageMethod').value:'tsvd',strength:$('sparseImaging').checked?+$('imageStrength').value:1,reweight:$('imageReweight').checked};}
function methodNote(){
 const o=options();$('sparseOptions').hidden=o.method==='tsvd';$('reweightControl').hidden=o.method!=='sparse';
 $('imageMethodNote').textContent=o.method==='wavelet'?'Positive source power with 2-D Haar-wavelet sparsity. Uncheck for the legacy TSVD comparison.':
  o.method==='sparse'?'Positive pixel ℓ₁ for compact echoes; optional reweighting and selected-support power refit.':'Legacy 18-mode TSVD comparison: squared inverse magnitude, not sparse source-power reconstruction.';
}
function makeWorker(){
 const active=worker=new Worker(new URL('./imaging_worker.js?v=persistence-2',import.meta.url),{type:'module'});
   worker.onmessage=async({data})=>{
    if(data.id!==sequence)return;
    busy=false;clearTimeout(timeout);
    $('imagePlot').setAttribute('aria-busy','false');
    if(data.error){$('imageStatus').textContent=data.error;return;}
    const {header,moments}=data,range=header.selected_range_km,extent=range*Math.sin(20*Math.PI/180);
    try {
     if(!plot){
      if(!plotMount){previousRange=range;plotMount=InteractiveFieldPlot.mount($('imagePlot'),{id:'image',frameParameters:styleParameters(range),
       spec:{bounds:[-extent,-extent,extent,extent],xlabel:'East (km)',ylabel:'North (km)',colorbar_label:'Relative power (dB)',panels:['']},colorLimits:[-30,0]});}
      try{plot=await plotMount;}finally{plotMount=undefined;}
     }
     if(data.id!==sequence)return;
     if(previousRange!==range){previousRange=range;plot.setFrameParameters(styleParameters(range));plot.setSpec({bounds:[-extent,-extent,extent,extent],xlabel:'East (km)',ylabel:'North (km)',colorbar_label:'Relative power (dB)',panels:['']});}
     lastMoments=moments;
     plot.panels[0].set_rectilinear_data(edges(extent),edges(extent),Float32Array.from({length:grid*grid},(_,i)=>moments[4*i]));repaint();
     const signal=Array.from({length:grid*grid},(_,i)=>moments[4*i+3]).some(x=>x>0);
     const d=data.diagnostics,method=d.method==='wavelet'?'PMSE wavelet ℓ₁':d.method==='sparse'?'Compact pixel ℓ₁':'Legacy TSVD';
     $('imageStatus').textContent=`${new Date(header.start_us/1000).toISOString()} — ${new Date(header.stop_us/1000).toISOString().slice(11,19)} UTC · slant range ${range.toFixed(2)} km · all five beams · ${header.looks} looks · ${header.index.length} retained spectral cells · Rust/WASM ${data.compute_ms.toFixed(0)} ms. `+
       `${method}. `+(header.hidden_single_mode_cells!=null?`${header.hidden_single_mode_cells} single-mode spectral cells hidden. `:'')+(signal?'':'No reconstructed power at this selection. ')+`Historical calibration: not independently validated for this observation.`;
     $('imageDiagnostics').textContent=d.cells?`Covariance fit: mean relative residual ${(100*d.relative_residual).toFixed(1)}% · weighted RMS ${d.residual_rms.toFixed(2)} · ${d.cells} cells · ${(d.iterations/d.cells).toFixed(0)} mean iterations. `+
      (d.unconverged?`${d.unconverged} cells reached the iteration limit; reconstruction is approximate. `:'Iterate tolerance reached. ')+
      'Noise weights assume independent Gaussian looks; correlated noise and detection censoring are not modeled.':'Legacy TSVD: no sparse-fit residual available.';
     if(d.cells&&d.relative_residual>.5)$('imageDiagnostics').textContent+=' Large residual: check calibration, angular support and source model.';
     $('imageDiagnostics').classList.toggle('warn',d.unconverged>0||d.relative_residual>.5);
     $('imagePlot').dataset.method=d.method;
    }catch(e){if(data.id===sequence)$('imageStatus').textContent=String(e);}
   };
 worker.onerror=event=>{if(worker!==active)return;busy=false;clearTimeout(timeout);$('imagePlot').setAttribute('aria-busy','false');$('imageStatus').textContent=`Imaging worker failed: ${event.message}`;active.terminate();worker=undefined;};
}
function reconstruct(){
 methodNote();if(!lastPacket)return;
 if($('sparseImaging').checked&&!$('imageStrength').checkValidity()){$('imageStrength').reportValidity();return;}
 const id=++sequence;
 if(busy){worker.terminate();worker=undefined;busy=false;}
 clearTimeout(timeout);if(!worker)makeWorker();
 const packet=lastPacket.slice(0);busy=true;
 $('imagePlot').setAttribute('aria-busy','true');
 $('imageDiagnostics').textContent='';
  $('imageStatus').textContent='Reconstructing all five beams in Rust/WebAssembly…';
 worker.postMessage({id,packet,geometry,grid,options:options(),allowedGates:echoFilter?.(lastKey)},[packet]);
 timeout=setTimeout(()=>{if(id===sequence&&busy){worker.terminate();worker=undefined;busy=false;$('imagePlot').setAttribute('aria-busy','false');$('imageStatus').textContent='Reconstruction exceeded 60 seconds. Select a different method or echo to retry.';}},60000);
}
export async function imageSelection(key,slantRange,offsetSeconds=0,spectraReady=true){
 lastKey=key;
 const id=++sequence;fetchAbort?.abort();fetchAbort=new AbortController();const signal=fetchAbort.signal;lastPacket=undefined;
 if(busy){worker.terminate();worker=undefined;busy=false;}clearTimeout(timeout);
 if(!spectraReady){$('imagePlot').setAttribute('aria-busy','false');$('imageStatus').textContent='These observations are available in the RTI; their complex spectra are still transferring from Syowa. Retry after the next update.';$('imageDiagnostics').textContent='';return;}
 $('imageStatus').textContent='Fetching existing complex spectra…';$('imageDiagnostics').textContent='';
 try {
  geometry??=await (await fetch('./api/geometry')).json();
  if(id!==sequence)return;
  const response=await fetch(`./api/spectra?key=${encodeURIComponent(key)}&range_km=${slantRange}&offset_s=${offsetSeconds}`,{signal});
  if(!response.ok)throw Error((await response.json()).error||'Spectra unavailable');
  const packet=await response.arrayBuffer();if(id!==sequence)return;
  if(packet.byteLength>64*2**20)throw Error('Selected spectra exceed the 64 MiB browser buffer limit');
  lastPacket=packet;reconstruct();
 }catch(e){if(id===sequence&&e.name!=='AbortError')$('imageStatus').textContent=String(e);}
}
for(const id of ['imageVelocity','imageWidth','display'])$(id).addEventListener('change',()=>{if(id!=='display')plot?.setHsvLimits(+$('imageVelocity').value,+$('imageWidth').value);repaint();});
for(const id of ['sparseImaging','imageMethod','imageStrength','imageReweight'])$(id).addEventListener('change',reconstruct);
methodNote();
