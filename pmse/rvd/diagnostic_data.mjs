const colors=['#202020','#2474a8','#b65c18','#22854a','#b52e48'];
const names=['Zenith','North','East','South','West'];
function chart(points,xlabel,ylabel,origin,{zero=false,series}={}){
 const finite=points.filter(p=>Number.isFinite(p[0])&&Number.isFinite(p[1]));
 if(!finite.length)throw Error(`No numeric data for ${ylabel}`);
 const xs=finite.map(p=>p[0]),ys=finite.map(p=>p[1]);
 let x0=Math.min(...xs),x1=Math.max(...xs),y0=zero?0:Math.min(...ys),y1=Math.max(...ys);
 if(x1<=x0)x1=x0+1;
 const pad=Math.max((y1-y0)*.12,Math.abs(y1)*.001,.00001);
 y0=zero?-pad*.15:y0-pad;y1+=pad;
 return {bounds:[x0,y0,x1,y1],xlabel,ylabel,utc_origin_seconds:origin,series:series||[{color:'#215d86',points,line:true,radius:2}]};
}
export function noiseSeries(model,module){
 const data=model.noise_model_data;
 if(!data?.times_h?.length||!Number.isFinite(model.noise_model_origin_unix_s))throw Error('Numeric sky-noise data have not arrived yet');
 return ['Receiver noise power (dB)','Modeled sky temperature (K)'].map((ylabel,panel)=>{
  const series=names.map((name,beam)=>({name,color:colors[beam],line:true,radius:1.7,points:data.times_h.map((t,i)=>[t,panel===0?10*Math.log10(data.noise_power[i][beam][module]):data.sky_temperature_k[i][beam]])}));
  return chart(series.flatMap(s=>s.points),'Time (UTC)',ylabel,model.noise_model_origin_unix_s,{series});
 });
}
export function archiveSeries(archive){
 if(!archive?.hourly_unix_s?.length)throw Error('Archive inventory data have not arrived yet');
 const origin=archive.hourly_unix_s[0],x=archive.hourly_unix_s.map(t=>(t-origin)/3600);
 const hourly=archive.hourly_bytes.map(b=>b/2**20);let total=0;
 const cumulative=[[x[0],0],...x.map((t,i)=>[t+1,total+=hourly[i]])];
 const retention=x.map((t,i)=>[t+.5,archive.hourly_dense_retained_fraction?.[i]==null?null:100*archive.hourly_dense_retained_fraction[i]]);
 const specs=[chart(x.map((t,i)=>[t+.5,hourly[i]]),'Time (UTC)','Data per UTC hour (MiB)',origin,{zero:true}),
  chart(cumulative,'Time (UTC)','Accumulated data (MiB)',origin,{zero:true})];
 // Missing retention is unknown, never a fabricated zero. Keep gaps in lines.
 if(retention.some(p=>p[1]!=null))specs.push(chart(retention,'Time (UTC)','Full spectra retained (%)',origin,{zero:true}));
 for(const spec of specs){spec.bounds[0]=0;spec.bounds[2]=Math.max(1,x.at(-1)+1);}
 return specs;
}
