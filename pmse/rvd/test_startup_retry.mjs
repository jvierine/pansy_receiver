import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import vm from 'node:vm';

const source=readFileSync(new URL('./viewer.js',import.meta.url),'utf8');
const polling=source.match(/async function pollObservations\(\)\{[\s\S]*?\n \}/)?.[0];
assert.ok(polling,'startup polling function must exist');
let attempts=0,failures=0,timers=[];
const context=vm.createContext({
 refresh:async()=>{if(++attempts===1)throw Error('Transient HTTP 503');},
 failure:()=>{failures++;},
 setTimeout:(callback,delay)=>timers.push({callback,delay}),
});
vm.runInContext(polling,context);
await context.pollObservations();
assert.equal(failures,1);
assert.equal(timers.length,1);
assert.equal(timers[0].delay,5000,'failed initial load must schedule recovery');
await timers.shift().callback();
assert.equal(attempts,2);
assert.equal(timers.length,1);
assert.equal(timers[0].delay,60000,'successful load resumes normal cadence');
console.log('Startup failure retries; successful recovery resumes normal polling.');
