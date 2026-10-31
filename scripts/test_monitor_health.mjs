import assert from 'node:assert/strict';
import {monitorHealth} from '../web/health.mjs';
import {latencyChart,rawChart} from '../web/monitor_data.mjs';
const now=1793424000;
const packet={generated_unix_s:now,freshness:{Raw:now-1,TX:now-60,Fits:now-3600}};
const ops={generated_utc:new Date(now*1000).toISOString(),services:[{name:'pansy-uhd-rx.service',active_state:'active'}],disks:[{path:'/media/analysis',used_percent:50,free_bytes:2**40}],phase_health:{ok:true,checked_unix_s:now-60,reason:'0.6 deg from reference'}};
assert.equal(monitorHealth(ops,packet,now).state,'ok');
const item=(o,p,label,t=now)=>monitorHealth(o,p,t).items.find(s=>s.label===label);
assert.equal(item({...ops,disks:[{...ops.disks[0],used_percent:95}]},packet,'Disks').state,'bad');
assert.equal(item({...ops,disks:[{...ops.disks[0],used_percent:85}]},packet,'Disks').state,'warn');
assert.equal(item({...ops,services:[{name:'pansy-uhd-rx.service',active_state:'failed'}]},packet,'Receiver').state,'bad');
assert.equal(item(ops,{...packet,freshness:{Raw:now-601}},'Receiver').state,'bad');
assert.equal(item({...ops,phase_health:{...ops.phase_health,ok:false}},packet,'Phases').state,'bad');
assert.equal(item({...ops,phase_health:{...ops.phase_health,checked_unix_s:now-900}},packet,'Phases').state,'unknown');
assert.equal(item({...ops,phase_health:undefined},packet,'Phases').state,'unknown');
assert.equal(monitorHealth(ops,packet,now+1801).state,'unknown');
for(const seconds of [0,30,600,3600,43199,43200,43201]){
 const p={...packet,freshness:{Raw:now,TX:now-seconds}};
 assert.equal(latencyChart(p).series[1].color,seconds<43200?'#138a38':'#d32626');
 assert.equal(item(ops,p,'Processing').state,seconds<43200?'ok':'bad');
}
const raw={sample_step_us:2,channels:Array.from({length:8},()=>({real:[1,2],imag:[2,3]}))};
for(let ch=0;ch<8;ch++){const spec=rawChart(raw,ch);assert.ok(spec.series.every(s=>!s.name));assert.equal(spec.compact_axes,true);assert.equal(spec.hide_x_labels,ch<6);}
console.log('PASS: disk, receiver, phase, stale/unknown health, 12-hour boundaries and compact raw traces');
