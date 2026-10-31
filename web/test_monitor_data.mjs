import assert from 'node:assert/strict';
import {zenithTemperature,rtiData,fitCharts,phaseChart,rawChart,latencyChart,meteorCharts} from './monitor_data.mjs';
const model={noise_model_origin_unix_s:0,noise_model_observation_stop_unix_s:172800,
 noise_zenith_model_status:Array.from({length:7},(_,ch)=>({power_per_kelvin:ch+1,receiver_temp_k:300+ch})),
 noise_model_data:{times_h:[0,24,36,48],sky_temperature_k:[[500,9999],[600,9999],[700,9999],[null,9999]],noise_power:Array.from({length:4},()=>Array.from({length:5},(_,b)=>Array.from({length:7},(_,ch)=>b?1e9:1000*(ch+1))))}};
const spec=zenithTemperature(model);
assert.equal(spec.series.length,14);
for(const s of spec.series.slice(0,7)){assert.deepEqual(s.points,[[24,1000],[36,1000],[48,1000]]);assert.equal(s.line,false);}
for(let ch=0;ch<7;ch++){
 assert.deepEqual(spec.series[7+ch].points,[[24,900+ch],[36,1000+ch],[48,null]],'Add fitted Trec to zenith sky only, preserving missing predictions');
 assert.equal(spec.series[7+ch].color,spec.series[ch].color);
 assert.equal(spec.series[7+ch].line,true);
}
model.noise_zenith_model_status[6].power_per_kelvin=null;
assert.equal(zenithTemperature(model).series.length,12,'No fabricated measurement or prediction for an uncalibrated receiver');
const geometry={slant_range_km:[80,81],altitude_km:Array.from({length:5},()=>[80,81])};
const data={geometry,axis_start_us:0,axis_stop_us:86400e6,frame_packets:[
 {start_us:0,stop_us:80e6,index:[[0,1,0],[0,2,0],[1,1,1]],snr_db:[5,9,12]},
 {start_us:100e6,stop_us:180e6,index:[],snr_db:[]}]};
const fields=rtiData(data);
assert.equal(fields.length,5);assert.equal(fields[0].z[0],9);
assert(Number.isNaN(fields[0].z[1]),'A time gap must remain NaN');
assert(Number.isNaN(fields[0].z[2]),'A noise-only mode must not become measured zero');
assert(Number.isNaN(fields[0].z[3]),'An unretained range cell must remain NaN');
assert.equal(fitCharts({generated_unix_s:172800,fits:[]}).counts.series[0].points.length,24);
console.log('Monitor numerical checks passed: zenith only, all receivers, P/G not Trec, 24 h, RTI gaps and empty modes.');
const map=fitCharts({generated_unix_s:172800,fits:[[172700,30,90,-90,-20]]}).radiant;
assert.equal(map.sky.central_longitude,270);assert.equal(map.sky.nside,32);
assert.deepEqual(map.sky.lonlat,[[-90,-20]],'Scientific coordinates must not be pre-rotated');
const phase=phaseChart({generated_unix_s:86400,phase:[[0,179],[600,-179]]});
assert.equal(phase.series[0].points[1][1],null,'Break the line across a phase wrap');
assert.equal(phase.series[0].radius,2.2);
const latency=latencyChart({generated_unix_s:1000,freshness:{Raw:999,TX:0}});
assert.equal(latency.bars,true);assert.equal(latency.series[0].color,'#138a38');assert.equal(latency.series[1].color,'#138a38');
assert.equal(latencyChart({generated_unix_s:43201,freshness:{TX:0}}).series[0].color,'#d32626');
assert.equal(rawChart({channels:[{real:[1,2],imag:[3,4]}],sample_step_us:8},0).series.length,2);
const meteor=meteorCharts({seconds:[0,2],east_km:[1,3],north_km:[2,4],up_km:[90,89],snr_db:[10,20],r0_km:[2,3,89.5],v0_km_s:[1,1,-.5],residual_std_km:[.1,.2,.3]});
assert.deepEqual(meteor.east.series[1].points,[[0,1],[2,3]],'Saved r0 is at mean time, not first echo');
assert.deepEqual(meteor.east.series[0].errors,[.2,.2]);
console.log('Sky density, bright phase wraps, health bars, complex raw I/Q and saved-fit timing checks passed.');
