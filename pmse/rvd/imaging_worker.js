// All numerical imaging runs in Rust/WASM, off the UI thread.
import init,{BrowserImager} from './imaging_pkg/pansy_browser_imaging.js?v=sparse-1';
let ready=init({module_or_path:new URL('./imaging_pkg/pansy_browser_imaging_bg.wasm?v=sparse-1',import.meta.url)}),imager,identity='';
self.onmessage=async({data})=>{
 try {
  await ready;
  const {id,packet,geometry,grid,options}=data;
  const view=new DataView(packet),length=view.getUint32(0,true);
  const headerJSON=new TextDecoder().decode(new Uint8Array(packet,4,length));
  const header=JSON.parse(headerJSON),offset=4+Math.ceil(length/4)*4;
  const values=new Float32Array(packet,offset);
  const key=JSON.stringify([geometry,header.selected_range_km,grid]);
  if(identity!==key){imager?.free();imager=new BrowserImager(JSON.stringify(geometry),header.selected_range_km,grid,true);identity=key;}
  const started=performance.now(),moments=imager.image_with_options(headerJSON,JSON.stringify(geometry),values,JSON.stringify(options));
  const diagnostics=JSON.parse(imager.last_diagnostics());
  self.postMessage({id,moments,header,diagnostics,compute_ms:performance.now()-started},[moments.buffer]);
 }catch(e){self.postMessage({id:data.id,error:String(e)});}
};
