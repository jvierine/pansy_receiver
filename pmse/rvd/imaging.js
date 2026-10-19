import {InteractiveFieldPlot} from './lib/interactive_plot.js';
const $=id=>document.getElementById(id);
let worker,geometry,plot,sequence=0,previousRange,lastMoments;
const grid=61;
function edges(extent){const d=2*extent/(grid-1);return Float32Array.from({length:grid+1},(_,i)=>-extent-d/2+i*d);}
function repaint(){if(plot&&lastMoments){plot.panels[0].set_moment_data(lastMoments,+$('imageVelocity').value,+$('imageWidth').value);plot.setFrameParameters({range_km:previousRange,velocity_limit:+$('imageVelocity').value,width_limit:+$('imageWidth').value});}}
export async function imageSelection(key,slantRange,offsetSeconds=0){
 const id=++sequence;$('imageStatus').textContent='Fetching existing complex spectra…';
 try {
  geometry??=await (await fetch('./api/geometry')).json();
  const response=await fetch(`./api/spectra?key=${encodeURIComponent(key)}&range_km=${slantRange}&offset_s=${offsetSeconds}`);
  if(!response.ok)throw Error((await response.json()).error||'Spectra unavailable');
  const packet=await response.arrayBuffer();if(id!==sequence)return;
  if(!worker){worker=new Worker(new URL('./imaging_worker.js',import.meta.url),{type:'module'});
   worker.onmessage=async({data})=>{
    if(data.id!==sequence)return;
    if(data.error){$('imageStatus').textContent=data.error;return;}
    const {header,moments}=data,range=header.selected_range_km,extent=range*Math.sin(20*Math.PI/180);
    try {
     if(!plot||previousRange!==range){plot?.dispose();$('imagePlot').replaceChildren();
      plot=await InteractiveFieldPlot.mount($('imagePlot'),{id:'image',frameParameters:{range_km:range,velocity_limit:+$('imageVelocity').value,width_limit:+$('imageWidth').value},
       spec:{bounds:[-extent,-extent,extent,extent],xlabel:'East (km)',ylabel:'North (km)',colorbar_label:'Relative power (dB)',panels:['']},colorLimits:[-30,0]});previousRange=range;
     }
     if(data.id!==sequence)return;
     lastMoments=moments;
     plot.panels[0].set_rectilinear_data(edges(extent),edges(extent),Float32Array.from({length:grid*grid},(_,i)=>moments[4*i]));repaint();
     const signal=Array.from({length:grid*grid},(_,i)=>moments[4*i+3]).some(x=>x>0);
     $('imageStatus').textContent=`${new Date(header.start_us/1000).toISOString()} — ${new Date(header.stop_us/1000).toISOString().slice(11,19)} UTC · slant range ${range.toFixed(2)} km · all five beams · ${header.looks} looks · ${header.index.length} retained spectral cells · Rust/WASM ${data.compute_ms.toFixed(0)} ms. `+
       (signal?'':'No retained echo at this selection. ')+`Historical calibration: not independently validated for this observation.`;
    }catch(e){$('imageStatus').textContent=String(e);}
   };
  }
  $('imageStatus').textContent='Reconstructing all five beams in Rust/WebAssembly…';
  worker.postMessage({id,packet,geometry,grid},[packet]);
 }catch(e){if(id===sequence)$('imageStatus').textContent=String(e);}
}
for(const id of ['imageVelocity','imageWidth'])$(id).addEventListener('change',repaint);
