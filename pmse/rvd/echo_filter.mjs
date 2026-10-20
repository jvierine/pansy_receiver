// Display-only temporal support. Native sparse spectra are never modified.
const gateKey=([beam,,range])=>`${beam}:${range}`;
export function persistentEchoes(frames,maxGapUs=600e6){
  const present=frames.map(frame=>new Set(frame.index.filter((_,i)=>Number.isFinite(frame.snr_db[i])&&frame.snr_db[i]>0).map(gateKey)));
  let hiddenCells=0,hiddenGates=0;
  const allowed=new Map();
  const filtered=frames.map((frame,i)=>{
    const previous=i>0&&frame.start_us-frames[i-1].stop_us<=maxGapUs?present[i-1]:new Set();
    const next=i+1<frames.length&&frames[i+1].start_us-frame.stop_us<=maxGapUs?present[i+1]:new Set();
    const gates=new Set([...present[i]].filter(key=>previous.has(key)||next.has(key)));
    allowed.set(String(frame.key),gates);hiddenGates+=present[i].size-gates.size;
    const keep=frame.index.map((cell,j)=>gates.has(gateKey(cell))?j:-1).filter(j=>j>=0);
    hiddenCells+=frame.index.length-keep.length;
    const result={...frame,index:keep.map(j=>frame.index[j]),snr_db:keep.map(j=>frame.snr_db[j])};
    for(const name of ['look_dominance','look_persistence'])if(frame[name])result[name]=keep.map(j=>frame[name][j]);
    return result;
  });
  return {frames:filtered,allowed,hiddenCells,hiddenGates};
}

// Keep complex sample rows and sparse indices aligned in the browser worker.
export function filterImageCells(header,values,allowedGates){
  if(allowedGates==null)return {header,values};
  const allowed=new Set(allowedGates),count=header.index.length;
  const keep=header.index.map((cell,i)=>allowed.has(gateKey(cell))?i:-1).filter(i=>i>=0);
  if(values.length!==header.looks*7*count*2)throw Error('Invalid complex sample layout');
  const filtered=new Float32Array(header.looks*7*keep.length*2);
  for(let row=0;row<header.looks*7;row++)for(let j=0;j<keep.length;j++){
    const from=2*(row*count+keep[j]),to=2*(row*keep.length+j);
    filtered[to]=values[from];filtered[to+1]=values[from+1];
  }
  return {header:{...header,index:keep.map(i=>header.index[i]),hidden_single_mode_cells:count-keep.length},values:filtered};
}
