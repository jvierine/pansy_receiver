import assert from 'node:assert/strict';
import {zenithTemperature,rtiData,fitCharts} from './monitor_data.mjs';
const model={noise_model_origin_unix_s:0,noise_model_observation_stop_unix_s:172800,
 noise_zenith_model_status:Array.from({length:7},(_,ch)=>({power_per_kelvin:ch+1,receiver_temp_k:300+ch})),
 noise_model_data:{times_h:[0,24,36,48],noise_power:Array.from({length:4},()=>Array.from({length:5},(_,b)=>Array.from({length:7},(_,ch)=>b?1e9:1000*(ch+1))))}};
const spec=zenithTemperature(model);
assert.equal(spec.series.length,7);
for(const s of spec.series){assert.deepEqual(s.points,[[24,1000],[36,1000],[48,1000]]);}
model.noise_zenith_model_status[6].power_per_kelvin=null;
assert.equal(zenithTemperature(model).series.length,6,'No fabricated temperature for an uncalibrated receiver');
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
