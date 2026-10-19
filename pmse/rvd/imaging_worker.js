// All numerical imaging runs in Rust/WASM, off the UI thread.
import init,{BrowserImager} from './imaging_pkg/pansy_browser_imaging.js';
let ready=init(),imager,identity='';
self.onmessage=async({data})=>{
 try {
  await ready;
  const {id,packet,geometry,grid}=data;
  const view=new DataView(packet),length=view.getUint32(0,true);
  const headerJSON=new TextDecoder().decode(new Uint8Array(packet,4,length));
  const header=JSON.parse(headerJSON),offset=4+Math.ceil(length/4)*4;
  const values=new Float32Array(packet,offset);
  const key=JSON.stringify([geometry,header.selected_range_km,grid]);
  if(identity!==key){imager?.free();imager=new BrowserImager(JSON.stringify(geometry),header.selected_range_km,grid,true);identity=key;}
  const started=performance.now(),moments=imager.image(headerJSON,JSON.stringify(geometry),values);
  self.postMessage({id,moments,header,compute_ms:performance.now()-started},[moments.buffer]);
 }catch(e){self.postMessage({id:data.id,error:String(e)});}
};
