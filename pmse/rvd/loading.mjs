// Fetch streams expose decoded bytes: a gzip Content-Length is not their total.
export async function readJSON(response, onProgress=()=>{}) {
  if(!response.ok)throw Error(`Data unavailable (${response.status})`);
  if(!response.body?.getReader)return response.json();
  const size=Number(response.headers.get('Content-Length'));
  const total=!response.headers.get('Content-Encoding')&&size>0?size:null;
  const reader=response.body.getReader(),decoder=new TextDecoder(),parts=[];
  let received=0;
  onProgress({received,total});
  try {
    while(true){
      const {done,value}=await reader.read();if(done)break;
      received+=value.byteLength;parts.push(decoder.decode(value,{stream:true}));
      onProgress({received,total});
    }
    parts.push(decoder.decode());
    return JSON.parse(parts.join(''));
  }finally{reader.releaseLock();}
}
