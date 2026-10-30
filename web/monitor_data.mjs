export const colors=['#0057d9','#e52521','#00872f','#9b19be','#ef8a00','#009eaa','#a44d00','#202020'];
export function chart(series,xlabel,ylabel,{origin,bounds}={}){
 const points=series.flatMap(s=>s.points).filter(p=>p.every(v=>v!==null&&Number.isFinite(v)));
 if(!points.length&&!bounds)throw Error(`No numeric observations for ${ylabel}`);
 if(!bounds){let x0=Infinity,x1=-Infinity,y0=Infinity,y1=-Infinity;
  for(const [x,y]of points){x0=Math.min(x0,x);x1=Math.max(x1,x);y0=Math.min(y0,y);y1=Math.max(y1,y);}
  const dx=Math.max((x1-x0)*.025,1e-5),dy=Math.max((y1-y0)*.1,Math.abs(y1)*.01,1e-5);
  bounds=[x0-dx,y0-dy,x1+dx,y1+dy];
 }
 return {bounds,xlabel,ylabel,series,...origin!==undefined?{utc_origin_seconds:origin}:{}};
}
export function zenithTemperature(model){
 const d=model.noise_model_data,fits=model.noise_zenith_model_status;
 if(!d?.times_h?.length||!fits||fits.length!==7)throw Error('Waiting for the zenith-only temperature calibration');
 const series=[];
 fits.forEach((fit,ch)=>{
  const gain=fit.power_per_kelvin;
  if(!(gain>0)||!Number.isFinite(gain)||!Number.isFinite(fit.receiver_temp_k))return;
  series.push({name:`RX ${ch}`,color:colors[ch],line:true,radius:1,
    points:d.times_h.flatMap((t,i)=>{if(model.noise_model_origin_unix_s+t*3600<model.noise_model_observation_stop_unix_s-86400)return [];
     const p=d.noise_power[i]?.[0]?.[ch];return [[t,p>0&&Number.isFinite(p)?p/gain:null]];})});
 });
 return chart(series,'Time (UTC)','Zenith system noise temperature (K)',{origin:model.noise_model_origin_unix_s});
}
export function fitCharts(packet){
 const fits=packet.fits||[],origin=packet.generated_unix_s-24*3600;
 const bins=new Array(24).fill(0);
 for(const row of fits){const i=Math.floor((row[0]-origin)/3600);if(i>=0&&i<24)bins[i]++;}
 const counts=chart([{color:colors[0],name:'Fitted meteors / hour',line:true,radius:2,
  points:bins.map((n,i)=>[i+.5,n])}],'Time (UTC)','Fitted meteors per hour',{origin});
 const speed=chart([{color:colors[0],line:false,radius:1.4,
  points:fits.map(r=>[(r[0]-origin)/3600,r[1]])}],'Time (UTC)','Meteor speed (km/s)',{origin,...!fits.length?{bounds:[0,0,24,80]}:{}});
 const radiant={...chart([], 'Solar-centered ecliptic longitude (deg)','Ecliptic latitude (deg)',{bounds:[-180,-90,180,90]}),
  sky:{lonlat:fits.filter(r=>Number.isFinite(r[3])&&Number.isFinite(r[4])).map(r=>[r[3],r[4]]),nside:32,central_longitude:270,west_positive:true,paper_style:true}};
 return {counts,speed,radiant};
}
export function phaseChart(packet){
 const rows=packet.phase||[],origin=packet.generated_unix_s-86400;
 const spec=chart(Array.from({length:8},(_,ch)=>({name:ch===7?'TX ref':`RX ${ch}`,color:colors[ch],line:true,radius:2.2,line_width:2.2,
  points:rows.flatMap((row,i)=>{const p=[(row[0]-origin)/3600,row[ch+1]??null],previous=rows[i-1];
   return previous&&Math.abs(row[ch+1]-previous[ch+1])>180?[[p[0],null],p]:[p];})})),
  'Time (UTC)','Relative phase (deg)',{origin});
 spec.bounds=[0,spec.bounds[1]-2,24,spec.bounds[3]+2];return spec;
}
export function rawChart(raw,ch){
 const channel=raw.channels[ch];
 return chart(['real','imag'].map((component,i)=>({name:component==='real'?'I':'Q',color:colors[i],line:true,radius:0,
  points:channel[component].map((v,n)=>[n*raw.sample_step_us,v])})),
  'Time (µs)',ch===7?'TX reference · ADC':`RX ${ch} · ADC`);
}
export function latencyChart(packet,now=packet.generated_unix_s){
 const names=['Raw','TX','Matched filter','Detections','Cuts','Modes','Fits'];
 const labels=['Raw','TX','MF','Events','Cuts','Modes','Fits'];
 const series=names.flatMap((name,i)=>{const t=packet.freshness?.[name];if(!Number.isFinite(t))return [];
  const age=Math.max(0,now-t)/60,healthy=age<(name==='Raw'?0.5:10);
  return [{color:healthy?'#138a38':'#d32626',line:false,points:[[i,age]]}];});
 const max=Math.max(1,...series.map(s=>s.points[0][1]*1.15));
 return {...chart(series,'Processing stage','Latest sample age (min)',{bounds:[-0.6,0,6.6,max]}),bars:true,categories:labels};
}
export function meteorCharts(m){
 const series=m.seconds,mean=series.reduce((a,b)=>a+b,0)/series.length,out={};
 ['east_km','north_km','up_km'].forEach((name,i)=>{
  const points=series.map((t,j)=>[t,m[name][j]]),err=m.residual_std_km?.[i];
  const traces=[{color:colors[0],line:false,radius:2.5,points,errors:Number.isFinite(err)?points.map(()=>2*err):[]}];
  if(m.r0_km&&m.v0_km_s)traces.push({color:'#e32626',line:true,radius:0,points:series.map(t=>[t,m.r0_km[i]+m.v0_km_s[i]*(t-mean)])});
  const spec=chart(traces,'Time since first echo (s)',['East (km)','North (km)','Up (km)'][i]);
  if(Number.isFinite(err)){spec.bounds[1]-=2*err;spec.bounds[3]+=2*err;}out[['east','north','up'][i]]=spec;
 });
 // Ten-degree beam cones intersect the event's actual mean-altitude plane.
 const height=m.up_km.reduce((a,b)=>a+b,0)/m.up_km.length,extent=height*Math.tan(20*Math.PI/180)*1.04;
 const footprints=[[0,0],[0,10],[90,10],[180,10],[270,10]].map(([az,za])=>{
  const a=az*Math.PI/180,z=za*Math.PI/180,b=10*Math.PI/180;
  const d=[Math.sin(z)*Math.sin(a),Math.sin(z)*Math.cos(a),Math.cos(z)],e=[Math.cos(a),-Math.sin(a),0],q=[Math.cos(z)*Math.sin(a),Math.cos(z)*Math.cos(a),-Math.sin(z)];
  return {color:'#8b969f',line:true,radius:0,line_width:.6,fill:'#00000013',points:Array.from({length:97},(_,i)=>{const p=i*2*Math.PI/96,u=d.map((v,j)=>Math.cos(b)*v+Math.sin(b)*(Math.cos(p)*e[j]+Math.sin(p)*q[j]));return [height*u[0]/u[2],height*u[1]/u[2]];})};
 });
 out.ground={...chart([...footprints,{color:colors[0],line:false,radius:2.5,points:series.map((t,i)=>[m.east_km[i],m.north_km[i]])}], 'East (km)','North (km)',{bounds:[-extent,-extent,extent,extent]}),equal_aspect:true};
 out.signal=chart([{color:colors[0],line:false,radius:2,points:series.map((t,i)=>[t,m.snr_db[i]])}], 'Time since first echo (s)','SNR (dB)');
 if(m.doppler_km_s)out.doppler=chart([{color:'#148d28',line:false,radius:3,points:series.map((t,i)=>[t,m.doppler_km_s[i]])}], 'Time since first echo (s)','Measured Doppler (km/s)');
 else if(m.match_function)out.doppler=chart([{color:'#148d28',line:false,radius:3,points:series.map((t,i)=>[t,m.match_function[i]])}], 'Time since first echo (s)','Match function');
 if(m.detection_trace){const d=m.detection_trace;
  out.signal=chart([{color:colors[0],line:false,radius:2.5,points:d.seconds.map((t,i)=>[t,d.range_km[i]])}],'Event time (s)','Detected range (km)');
  if(!m.doppler_km_s)out.doppler=chart([{color:'#148d28',line:false,radius:2.5,points:d.seconds.map((t,i)=>[t,d.doppler_km_s[i]])}],'Event time (s)','Measured Doppler (km/s)');
 }
 if(m.range_time){const r=m.range_time;out.signal={...chart([], 'Time since first echo (s)','Range (km)',{bounds:[Math.min(...r.seconds),Math.min(...r.range_km),Math.max(...r.seconds),Math.max(...r.range_km)]}),field:{x:axisEdges(r.seconds),y:axisEdges(r.range_km),values:r.snr_db.flat(),limits:[0,20],label:'SNR (dB)',palette:'turbo'}};}
 return out;
}
export const axisEdges=values=>[values[0]-(values[1]-values[0])/2,...values.slice(1).map((v,i)=>(v+values[i])/2),values.at(-1)+(values.at(-1)-values.at(-2))/2];
export function rtiChart(manifest){
 if(manifest.schema==='pansy_zenith_display_v1'){
  const nx=manifest.x.length-1,values=new Array(nx*(manifest.y.length-1)).fill(null);
  for(const [c,r,v]of manifest.cells)values[r*nx+c]=v;
  return {...chart([], 'Time (UTC)','Zenith altitude (km)',{origin:manifest.axis_start_us/1e6,bounds:[0,50,24,120]}),field:{x:manifest.x,y:manifest.y,values,limits:[0,13],palette:'plasma',label:'SNR (dB)',missing_value:0}};
 }
 const data=rtiData(manifest)[0];
 return {...chart([], 'Time (UTC)','Zenith altitude (km)',{origin:manifest.axis_start_us/1e6,bounds:[0,50,24,120]}),field:{x:Array.from(data.x),y:Array.from(data.y),values:Array.from(data.z,v=>Number.isFinite(v)?v:null),limits:[0,13],palette:'plasma',label:'SNR (dB)',missing_value:0}};
}
export function meteorChart(meteor){
 return chart(['east_km','north_km','up_km'].map((name,i)=>({name:['East','North','Up'][i],color:colors[i],line:true,radius:2,
  points:meteor.seconds.map((t,j)=>[t,meteor[name][j]])})),
  'Time since first fitted echo (s)','Meteor position (km)');
}
export function rtiData(manifest){
 const frames=[...(manifest.frame_packets||[])].sort((a,b)=>a.start_us-b.start_us);
 const geom=manifest.geometry,first=manifest.axis_start_us,nr=geom.slant_range_km.length;
 if(!frames.length)throw Error('No PMSE mode summaries have arrived for this UTC day');
 // Exact mode intervals, with explicit NaN columns. No interpolation over gaps.
 const columns=[],edges=[];
 for(const f of frames){
  const left=Math.max(first,f.start_us),right=Math.min(manifest.axis_stop_us,f.stop_us);
  if(right<=left)continue;
  if(edges.length&&left<edges.at(-1))throw Error('Overlapping PMSE intervals');
  if(edges.length&&left>edges.at(-1)){columns.push(null);edges.push(left);}
  if(!edges.length)edges.push(left);
  columns.push(f);edges.push(right);
 }
 const x=Float32Array.from(edges,t=>(t-first)/3.6e9),nx=columns.length;
 const axisEdges=values=>Float32Array.from([values[0]-(values[1]-values[0])/2,
  ...values.slice(1).map((v,i)=>(v+values[i])/2),values.at(-1)+(values.at(-1)-values.at(-2))/2]);
 return Array.from({length:5},(_,beam)=>{
  const z=new Float32Array(nx*nr).fill(NaN);
  columns.forEach((f,c)=>{if(!f)return;f.index.forEach(([b,bin,r],i)=>{
   if(b!==beam||r<0||r>=nr||!Number.isFinite(f.snr_db[i]))return;
   const cell=r*nx+c;z[cell]=Number.isNaN(z[cell])?f.snr_db[i]:Math.max(z[cell],f.snr_db[i]);
  });});
  return {x,y:axisEdges(geom.altitude_km[beam]),z};
 });
}
