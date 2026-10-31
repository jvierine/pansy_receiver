// Source health is distinct from delivery age and event-driven processing lag.
export function monitorHealth(ops,packet,now=Date.now()/1000){
 const unknown=(label,detail)=>({label,state:'unknown',detail});
 const opsTime=Date.parse(ops?.generated_utc)/1000,opsFresh=Number.isFinite(opsTime)&&now-opsTime>=-60&&now-opsTime<1800;
 let disk=unknown('Disks','Awaiting a current disk report'),receiver=unknown('Receiver','Awaiting a current service report'),phase=unknown('Phases','Awaiting a current watchdog check');
 if(opsFresh){
  const disks=ops.disks||[],critical=disks.filter(d=>d.critical||d.used_percent>=95),warning=disks.filter(d=>d.warning||d.used_percent>=80);
  const worst=[...disks].filter(d=>Number.isFinite(d.used_percent)).sort((a,b)=>b.used_percent-a.used_percent)[0];
  disk={label:'Disks',state:critical.length?'bad':disks.some(d=>d.error)||!disks.length?'unknown':warning.length?'warn':'ok',detail:worst?`${worst.path.replace('/media/','')}: ${worst.used_percent.toFixed(1)}% · ${(worst.free_bytes/2**30).toFixed(0)} GiB free`:'Disk measurements unavailable'};
  const services=ops.services||[],failed=services.filter(s=>s.active_state!=='active');
  receiver={label:'Receiver',state:failed.length?'bad':services.length?'ok':'unknown',detail:failed.length?failed.map(s=>s.name.replace('pansy-','').replace('.service','')+': '+s.active_state).join(', '):'Acquisition and processing services active'};
  const p=ops.phase_health;
  if(p&&Number.isFinite(p.checked_unix_s)&&now-p.checked_unix_s>=-60&&now-p.checked_unix_s<900&&typeof p.ok==='boolean')phase={label:'Phases',state:p.ok?'ok':'bad',detail:p.reason};
  else if(p)phase=unknown('Phases',p.ok==null?p.reason:'Last watchdog phase result is over 15 minutes old');
 }
 const source=packet?.generated_unix_s,raw=packet?.freshness?.Raw;
 if(Number.isFinite(source)&&now-source<1800&&now-source>=-60&&Number.isFinite(raw)&&source-raw>600)receiver={label:'Receiver',state:'bad',detail:`Raw recording stalled: ${Math.round((source-raw)/60)} min behind at source`};
 else if(receiver.state==='ok'&&(!Number.isFinite(source)||now-source>=1800||!Number.isFinite(raw)))receiver=unknown('Receiver','Services active; current raw-recording health unconfirmed');
 const ages=Object.values(packet?.freshness||{}).filter(Number.isFinite).map(t=>Math.max(0,source-t));
 const processing=Number.isFinite(source)&&now-source<1800&&ages.length?{label:'Processing',state:Math.max(...ages)>=43200?'bad':'ok',detail:`Largest source lag: ${(Math.max(...ages)/3600).toFixed(1)} h · red at 12 h`}:unknown('Processing','Awaiting current source-stage ages');
 const items=[disk,receiver,phase,processing];
 return {items,state:items.some(s=>s.state==='bad')?'bad':items.some(s=>s.state==='unknown')?'unknown':items.some(s=>s.state==='warn')?'warn':'ok'};
}
