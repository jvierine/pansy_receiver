import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {archiveSeries,noiseSeries} from '../pmse/rvd/diagnostic_data.mjs';
test('archive numeric plots preserve hourly totals, cumulative totals and retention gaps',()=>{
 const plots=archiveSeries({hourly_unix_s:[1000,4600,8200],hourly_bytes:[1048576,2097152,0],hourly_dense_retained_fraction:[0,.001,null]});
 assert.equal(plots.length,3);assert.deepEqual(plots[0].series[0].points,[[.5,1],[1.5,2],[2.5,0]]);
 assert.deepEqual(plots[1].series[0].points,[[0,0],[1,1],[2,3],[3,3]]);
 assert.deepEqual(plots[2].series[0].points,[[.5,0],[1.5,.1],[2.5,null]]);
 for(const spec of plots){assert.equal(spec.utc_origin_seconds,1000);assert.equal(spec.bounds[0],0);assert.equal(spec.bounds[2],3);}
});
test('noise arrays use the selected module and preserve physical UTC and dB units',()=>{
 const data={noise_model_origin_unix_s:1000,noise_model_data:{times_h:[0,1],sky_temperature_k:[[1000,2000,3000,4000,5000],[1100,2100,3100,4100,5100]],noise_power:Array.from({length:2},()=>Array.from({length:5},()=>[1,10,100,1000,10000,100000,1000000]))}};
 const [power,sky]=noiseSeries(data,2);assert.equal(power.series.length,5);assert.deepEqual(power.series[0].points,[[0,20],[1,20]]);assert.deepEqual(sky.series[4].points,[[0,5000],[1,5100]]);
});
test('auxiliary plots no longer embed or fetch server SVGs',async()=>{
 const html=await readFile(new URL('../pmse/rvd/index.html',import.meta.url),'utf8');
 const viewer=await readFile(new URL('../pmse/rvd/viewer.js',import.meta.url),'utf8');
 assert(!html.includes('id="noisePlot"'));assert(!html.includes('id="archivePlot"'));
 assert(!viewer.includes('noise-module-'));assert(!viewer.includes('archive-growth-'));assert(!viewer.includes('.svg'));
});
