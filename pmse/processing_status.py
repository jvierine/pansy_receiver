"""Small, read-only source-worker status; never exposes access paths or hosts."""
import json
import re
import subprocess
import time


def parse_worker(lines, backfill=False):
    result={}
    for line in lines:
        try:
            value=json.loads(line)
        except (ValueError,TypeError):
            value={}
        if not isinstance(value,dict):value={}
        if 'start' in value and ('result' in value or 'processing_seconds' in value):
            result['last_written_start_us']=int(value['start'])
            if 'seconds' in value and not backfill:
                result['last_written_stop_us']=int(value['start'])+round(value['seconds']*1e6)
            result['last_result']=value.get('result','created')
        match=re.match(r'Rectangular backfill (\d+)/(\d+): (\d+)\.\.(\d+)',line)
        if backfill and match:
            result.update(current_index=int(match[1]),total_modes=int(match[2]),
                          current_start_us=int(match[3]))
        if line.startswith('Rectangular backfill complete'):
            result['complete']=True
    return result


def status():
    output=dict(generated_unix_s=time.time())
    for name,unit,backfill in [('live','pansy-pmse-sparse.service',False),
                               ('backfill','pansy-pmse-rectangular-backfill.service',True)]:
        try:
            journal=subprocess.run(['journalctl','--user','-u',unit,'-n','60','-o','cat','--no-pager'],
                                   text=True,capture_output=True,check=True,timeout=10)
            worker=parse_worker(journal.stdout.splitlines(),backfill)
            active=subprocess.run(['systemctl','--user','is-active',unit],
                                  text=True,capture_output=True,timeout=5)
            worker['state']=active.stdout.strip() or 'unknown'
            output[name]=worker
        except (OSError,subprocess.SubprocessError,ValueError,TypeError,KeyError):
            output[name]=dict(state='unknown')
    return output


if __name__=='__main__':
    print(json.dumps(status()))
