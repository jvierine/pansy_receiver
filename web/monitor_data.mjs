export const colors=['#176594','#bd5b26','#328455','#9b4d98','#887718','#b24655','#5349a5','#5c737c'];
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
 const radiant=chart([{color:colors[0],line:false,radius:1.2,points:fits.map(r=>[r[3],r[4]])}],
  'Solar-centered ecliptic longitude (deg)','Ecliptic latitude (deg)',{bounds:[-180,-90,180,90]});
 return {counts,speed,radiant};
}
export function phaseChart(packet){
 const rows=packet.phase||[],origin=rows[0]?.[0];
 return chart(Array.from({length:8},(_,ch)=>({name:`RX ${ch}`,color:colors[ch],line:false,radius:1,
  points:rows.map(row=>[(row[0]-origin)/3600,row[ch+1]??null])})),
  'Time (UTC)','Transmit-pulse relative phase (deg)',{origin});
}
export function rawChart(raw,ch){
 const channel=raw.channels[ch];
 return chart(['real','imag'].map((component,i)=>({name:component==='real'?'I':'Q',color:colors[i],line:true,radius:0,
  points:channel[component].map((v,n)=>[n*raw.sample_step_us,v])})),
  'Time since trace start (microseconds)','Raw complex voltage (ADC units)');
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
