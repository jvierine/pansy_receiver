// Scientific provenance: rustplotlib viewer 74eff2d; no figure generation in JS.
import init,{FieldView} from './pkg/rustplotlib_web.js';
const $=id=>document.getElementById(id), names=['Zenith','North','East','South','West'];
let manifest,frames=[],plots=[],rd,playing=null,selected=0,revision='',layouts={};
const utc=us=>new Date(us/1000).toISOString().replace('T',' ').slice(0,19)+' UTC';
const edges=centers=>Float32Array.from([...centers.map((x,i)=>i?(centers[i-1]+x)/2:x-(centers[1]-x)/2),centers.at(-1)+(centers.at(-1)-centers.at(-2))/2]);
const resize=canvas=>{const rect=canvas.getBoundingClientRect(),dpr=Math.min(devicePixelRatio||1,2);canvas.width=Math.max(2,Math.round(rect.width*dpr));canvas.height=Math.max(2,Math.round(rect.height*dpr));};
function makePlot(parent,title,rdMode=false){
  const box=document.createElement('div');box.className='plot'+(rdMode?' rd':'');
  const h=document.createElement('h2');h.textContent=title;box.append(h);
  const field=document.createElement('div');field.className='native-frame';
  const axes=document.createElement('img');axes.alt='';axes.className='native-axes';field.append(axes);
  const canvas=document.createElement('canvas');canvas.setAttribute('aria-label',title);canvas.className='native-data';field.append(canvas);
  box.append(field);parent.append(box);canvas.width=10;canvas.height=10;
  const view=new FieldView(canvas);view.set_missing_color(.18,.23,.30,1);
  return {box,field,canvas,axes,view,kind:rdMode?'rd':'rti'};
}
function nativeFrame(p){
  const width=p.box.clientWidth<500?360:p.box.clientWidth<1000?640:960;
  const key=`${p.kind}-${width}-${$('low').value}-${$('high').value}`,layout=layouts[key];
  if(!layout)throw Error('Native rustplotlib axes are not available yet for this scale');
  p.field.style.width=layout.width+'px';p.field.style.height=layout.height+'px';
  p.axes.src='./data/axes/'+layout.file+'?v='+revision;
  const [x,y,w,h]=layout.panel;p.canvas.style.left=(x*100)+'%';p.canvas.style.top=(y*100)+'%';
  p.canvas.style.width=(w*100)+'%';p.canvas.style.height=(h*100)+'%';resize(p.canvas);
}
function scale(p){const low=+$('low').value,high=+$('high').value;if(!(high>low))throw Error('SNR high must exceed low');p.view.set_scale(low,high,0,.1,2);}
function crop(p,y){const ye=edges(y);p.view.set_view(0,(50-ye[0])/(ye.at(-1)-ye[0]),1,(120-ye[0])/(ye.at(-1)-ye[0]));}
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
    p.box.hidden=$('beam').value!=='all'&&+$('beam').value!==beam;if(p.box.hidden)return;nativeFrame(p);
    const val=new Float32Array(nx*nr),moment=new Float32Array(nx*nr*4);
    columns.forEach((frame,c)=>{const s=frame?summaries(frame,beam,nr):null;
      for(let r=0;r<nr;r++){val[r*nx+c]=s?s.peak[r]:NaN;for(let k=0;k<4;k++)moment[4*(r*nx+c)+k]=s?s.moments[4*r+k]:(k===0?NaN:0);}});
    p.view.set_rectilinear_data(x,edges(geom.altitude_km[beam]),val);scale(p);
    if($('display').value==='hsv')p.view.set_moment_data(moment,+$('vlim').value,+$('wlim').value);crop(p,geom.altitude_km[beam]);
    p.canvas.onclick=e=>{const t=xb[0]+(e.clientX-p.canvas.getBoundingClientRect().left)/p.canvas.clientWidth*(xb.at(-1)-xb[0]);selected=frames.reduce((best,f,i)=>Math.abs(f.start_us-t)<Math.abs(frames[best].start_us-t)?i:best,0);$('rdbeam').value=String(beam);$('frame').value=String(selected);renderRD();};
  });
  const hsv=$('display').value==='hsv';$('colorLabel').textContent=hsv?`Hue −${$('vlim').value} to +${$('vlim').value} m/s; saturation decreases over 0–${$('wlim').value} m/s width; intensity ${$('low').value}–${$('high').value} dB`:`Peak excess signal/noise, ${$('low').value} to ${$('high').value} dB`;
}
function renderRD(){
  const frame=frames[selected],beam=+$('rdbeam').value,g=manifest.geometry,nr=g.slant_range_km.length;
  const bins=g.velocity_mps.map((v,i)=>Math.abs(v)<=90.322?i:-1).filter(i=>i>=0),nx=bins.length;
  const lookup=new Map(bins.map((f,i)=>[f,i])),value=new Float32Array(nx*nr);let count=0,transients=0;
  frame.index.forEach(([b,f,r],i)=>{if(b===beam&&lookup.has(f)){value[r*nx+lookup.get(f)]=frame.snr_db[i];count++;if(frame.look_dominance[i]>.5)transients++;}});
  nativeFrame(rd);rd.view.set_rectilinear_data(edges(bins.map(f=>g.velocity_mps[f])),edges(g.altitude_km[beam]),value);scale(rd);crop(rd,g.altitude_km[beam]);
  $('frameStatus').textContent=`${utc(frame.start_us)} · ${names[beam]} · ${frame.duration_s.toFixed(2)} s mode, ${frame.coverage_s.toFixed(2)} s valid support · ${frame.looks} looks · ${count} retained cells (${transients} dominated by one look) · ${frame.quality}`;
}
async function refresh(){
  const response=await fetch('./data/manifest.json',{cache:'no-cache'});if(!response.ok)throw Error(`Data unavailable (${response.status})`);
  const m=await response.json();if(m.revision===revision)return;
  const cache=new Map(frames.map(f=>[f.key,f]));manifest=m;
  // A bounded pool avoids issuing hundreds of simultaneous requests.
  const result=new Array(m.frames.length);let next=0;
  await Promise.all(Array.from({length:6},async()=>{for(;;){const i=next++;if(i>=result.length)break;const key=m.frames[i].key;if(cache.has(key)){result[i]=cache.get(key);continue;}const r=await fetch(`./data/frames/${key}.json`);if(!r.ok)throw Error('Incomplete frame transfer');result[i]=await r.json();}}));
  frames=result;revision=m.revision;selected=frames.length-1;$('frame').max=String(selected);$('frame').value=String(selected);
  if(!m.plot_layouts)throw Error('Native rustplotlib plot frames are being built; observations are retained.');
  layouts=await (await fetch('./data/'+m.plot_layouts+'?v='+revision)).json();
  const age=(Date.now()-m.latest_mode_us/1000)/60000;$('status').textContent=`${frames.length} observing modes · latest completed ${utc(m.latest_mode_us)} · ${Math.max(0,age).toFixed(1)} min old`;$('status').className=age>30?'status warn':'status';
  $('sidecar').href='./data/'+m.sidecar;renderRTI();renderRD();
}
function failure(e){$('error').textContent=e.message;console.error(e);}
try{
 await init();if(new URLSearchParams(location.search).has('embed'))document.body.classList.add('embedded');
 for(let beam=0;beam<5;beam++)plots.push(makePlot($('rti'),names[beam]));rd=makePlot($('rd'),'Retained spectral bins',true);
 for(const id of ['display','beam','low','high','vlim','wlim'])$(id).onchange=()=>{try{renderRTI();renderRD();}catch(e){failure(e);}};
 $('rdbeam').onchange=renderRD;$('frame').oninput=()=>{selected=+$('frame').value;renderRD();};
 $('play').onclick=()=>{if(playing){clearInterval(playing);playing=null;$('play').textContent='Play';}else{$('play').textContent='Pause';playing=setInterval(()=>{selected=(selected+1)%frames.length;$('frame').value=selected;renderRD();},500);}};
 window.addEventListener('resize',()=>{if(frames.length){renderRTI();renderRD();}});
 await refresh();setInterval(()=>refresh().catch(failure),60000);
}catch(e){failure(e);}
