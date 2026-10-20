import {InteractiveSeriesPlot} from './lib/interactive_series.js?v=1';
import {noiseSeries,archiveSeries} from './diagnostic_data.mjs?v=1';
const plots=new Map(),mounting=new Map();
async function draw(id,spec){
 const root=document.getElementById(id);if(!root)return;
 if(!plots.has(id)){
  if(!mounting.has(id))mounting.set(id,InteractiveSeriesPlot.mount(root,spec,{height:310}).then(p=>{plots.set(id,p);return p;}).finally(()=>mounting.delete(id)));
  await mounting.get(id);
 }
 plots.get(id).setSpec(spec);
}
export async function renderNoiseDiagnostics(model,module){
 const specs=noiseSeries(model,module);
 await Promise.all(specs.map((s,i)=>draw(['noisePowerPlot','skyTemperaturePlot'][i],s)));
}
export async function renderArchiveDiagnostics(archive){
 const specs=archiveSeries(archive);
 await Promise.all(specs.map((s,i)=>draw(['hourlyDataPlot','cumulativeDataPlot','retentionPlot'][i],s)));
}
