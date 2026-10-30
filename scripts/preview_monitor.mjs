// Read-only local deployment preview; no radar data processing.
import {createServer} from 'node:http';
import {readFile} from 'node:fs/promises';
import {resolve,extname} from 'node:path';
const root=resolve(new URL('../web/',import.meta.url).pathname);
const cache=new Map();
createServer(async(req,res)=>{
 try{
  const path=new URL(req.url,'http://localhost').pathname;
  const isData=path.startsWith('/rvd/data/')||['monitor-data.json','ops_status.json','status.json'].some(n=>path==='/pansy/'+n);
  if(isData){
   let bytes=cache.get(path);
   if(!bytes){const response=await fetch('https://juha.no'+path);if(!response.ok)throw Error('Upstream HTTP '+response.status);bytes=Buffer.from(await response.arrayBuffer());cache.set(path,bytes);}
   res.setHeader('Content-Type','application/json');res.end(bytes);return;
  }
  const relative=path.replace(/^\/pansy\//,'')||'index.html';
  const file=resolve(root,relative);if(!file.startsWith(root+'/'))throw Error('Invalid preview path');
  res.setHeader('Content-Type',({'.js':'text/javascript','.mjs':'text/javascript','.html':'text/html','.css':'text/css','.wasm':'application/wasm','.svg':'image/svg+xml','.png':'image/png'})[extname(file)]||'application/octet-stream');
  res.end(await readFile(file));
 }catch(error){res.writeHead(404);res.end(String(error));}
}).listen(8279,'127.0.0.1',()=>console.log('Preview: http://127.0.0.1:8279/pansy/'));
