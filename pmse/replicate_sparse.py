"""Replicate native HDF5 without ownership changes; persist transfer evidence."""
import argparse
import os
from pathlib import Path
import re
import signal
import subprocess
import time
import h5py
from .web_retention import sha256


def rsync(source,destination,bandwidth=None,remote_shell=None,files_from=None,preserve_times=True,force=False,timeout=900,partial_ok=False):
    cmd=['rsync','-rt','--no-perms','--no-owner','--no-group','--partial','--delay-updates',
         '--timeout=120','--stats','--exclude=/.writer.lock']
    if bandwidth: cmd.append(f'--bwlimit={bandwidth}')
    if remote_shell: cmd+=['-e',remote_shell]
    if files_from is not None:cmd+=['--files-from=-']
    if not preserve_times:cmd+=['--no-times','--omit-dir-times','--size-only']
    if force:cmd+=['--ignore-times']
    began=time.monotonic()
    process=subprocess.Popen(cmd+[source,destination],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
        stdin=subprocess.PIPE,env=dict(os.environ,LC_ALL='C'),start_new_session=True)
    interrupted=False
    try:
        stdout,stderr=process.communicate(input=None if files_from is None else '\n'.join(files_from)+'\n',timeout=timeout)
    except subprocess.TimeoutExpired:
        interrupted=True;os.killpg(process.pid,signal.SIGINT)
        try:stdout,stderr=process.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGKILL);stdout,stderr=process.communicate()
    print(stdout,flush=True)
    if process.returncode and not(interrupted and partial_ok):raise RuntimeError(f'rsync {process.returncode}: {stderr}')
    elapsed=time.monotonic()-began
    match=re.search(r'Literal data:\s*([\d,]+)',stdout)
    payload=int(match[1].replace(',','')) if match else 0
    return dict(finished_unix_s=time.time(),elapsed_seconds=elapsed,literal_bytes=payload,
                payload_bytes_per_second=payload/max(elapsed,1e-9),bandwidth_cap_kib_s=bandwidth or 0,partial=interrupted)


def save_metrics(path,leg,values):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with h5py.File(path,'a') as h:
        g=h.require_group('runs').create_group(f'{time.time_ns()}-{leg}')
        g.attrs.update(leg=leg,generator='pmse/replicate_sparse.py')
        for key,value in values.items(): g[key]=value
        h.flush()


def catchup_files(files,root,now):
    """Do not strand earlier live chunks behind the multi-day backfill."""
    pending=[]
    for p,size,mtime in files[1:]:
        local=Path(root)/p
        if not local.exists() or local.stat().st_size!=size or abs(local.stat().st_mtime-mtime)>1e-6:
            pending.append(p)
    recent=[p for p in pending if int(Path(p).stem.split('@')[1])>=now-3*3600]
    historical=[p for p in pending if p not in recent]
    return recent[:1],historical[-1:]


def changed_files(root,receipt,relative_paths):
    """Force changed chunks even when a Hann replacement has the same size."""
    if not Path(receipt).exists():return list(relative_paths)
    changed=[]
    with h5py.File(receipt,'r') as h:
        for relative in relative_paths:
            if relative=='dmd_properties.h5':continue
            st=(Path(root)/relative).stat()
            g=h.get('files/'+relative)
            if g is None or int(g.attrs.get('mtime_ns',-1))!=st.st_mtime_ns or int(g.attrs.get('bytes',-1))!=st.st_size:
                changed.append(relative)
    return changed


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/mnt/data/juha/pansy/metadata/pmse_sparse16'))
    p.add_argument('--metrics',type=Path,default=Path('/mnt/data/juha/pansy/rvd16-transfer.h5'))
    p.add_argument('--syowa-kib-s',type=int,default=10)
    p.add_argument('--archive-host',default='j@juha.no')
    p.add_argument('--archive-root',default='/mnt/shovel/rvd/pmse_sparse16')
    p.add_argument('--loop',action='store_true');a=p.parse_args();a.root.mkdir(parents=True,exist_ok=True)
    while True:
        try:
            ssh=['ssh','-p','2222','-o','BatchMode=yes','-o','ConnectTimeout=30',
                 '-o','ControlMaster=auto','-o','ControlPersist=300','-o','ControlPath=/tmp/pansy16-syowa-control']
            subprocess.run(ssh+['radar@127.0.0.1',
                'cd /home/radar/src/git/pansy_pmse_sparse && OPENBLAS_NUM_THREADS=1 /home/radar/testenv/bin/python -m pmse.snapshot_archive'],check=True,timeout=120)
            listing=subprocess.run(ssh+['radar@127.0.0.1',
                "find /media/analysis/metadata/pmse_sparse16_snapshot -maxdepth 2 -type f -name '*.h5' -printf '%P %s %T@\\n'"],
                check=True,capture_output=True,text=True,timeout=60)
            files=[]
            for line in listing.stdout.splitlines():
                path,size,mtime=line.rsplit(' ',2)
                if re.fullmatch(r'20\d\d-\d\d-\d\dT\d\d-00-00/pmse@\d+\.h5',path):files.append((path,int(size),float(mtime)))
            files.sort(reverse=True)
            source='radar@127.0.0.1:/media/analysis/metadata/pmse_sparse16_snapshot/'
            shell=' '.join(ssh)
            # Latest chunk first. Bound historical transfer work so a large
            # backfill can never monopolize the link for fifteen minutes.
            hot=['dmd_properties.h5']+[p for p,_,_ in files[:1]]
            values=rsync(source,str(a.root)+'/',a.syowa_kib_s,shell,files_from=hot)
            save_metrics(a.metrics,'syowa_to_revontuli',values)
            recent,cold=catchup_files(files,a.root,time.time())
            if recent:
                values=rsync(source,str(a.root)+'/',a.syowa_kib_s,shell,files_from=recent)
                save_metrics(a.metrics,'syowa_to_revontuli_catchup',values)
            if cold:
                values=rsync(source,str(a.root)+'/',a.syowa_kib_s,shell,files_from=[cold[-1]],timeout=45,partial_ok=True)
                save_metrics(a.metrics,'syowa_to_revontuli_backfill',values)
            # Fail closed if shovel is unmounted; never fill the server root disk.
            subprocess.run(['ssh',a.archive_host,'mountpoint -q /mnt/shovel && test -d /mnt/shovel/rvd/pmse_sparse16'],check=True,timeout=30)
            cutoff=time.time()-7*86400-600
            keep=['dmd_properties.h5']+[str(path.relative_to(a.root)) for path in a.root.glob('*T*/pmse@*.h5')
                if int(path.stem.split('@')[1])>=cutoff]
            receipt=a.metrics.with_name('rvd-revontuli-receipt.h5')
            changed=changed_files(a.root,receipt,keep)
            values=rsync(str(a.root)+'/',f'{a.archive_host}:{a.archive_root}/',files_from=keep,preserve_times=False)
            if changed:
                rsync(str(a.root)+'/',f'{a.archive_host}:{a.archive_root}/',files_from=changed,preserve_times=False,force=True)
            save_metrics(a.metrics,'revontuli_to_shovel',values)
            # Only changed/new source files are rehashed. The receipt proves
            # preservation on Revontuli before any secondary expiration.
            with h5py.File(receipt,'a') as h:
                files=h.require_group('files')
                for path in a.root.glob('*T*/pmse@*.h5'):
                    relative=str(path.relative_to(a.root));g=files.require_group(relative);st=path.stat()
                    if int(g.attrs.get('mtime_ns',-1))!=st.st_mtime_ns or int(g.attrs.get('bytes',-1))!=st.st_size:
                        g.attrs.update(sha256=sha256(path),bytes=st.st_size,mtime_ns=st.st_mtime_ns)
                h.attrs.update(source=str(a.root),checked_unix_s=time.time())
            rsync(str(receipt),f'{a.archive_host}:/mnt/shovel/rvd/revontuli-receipt.h5',preserve_times=False,force=True)
            # Source deployment of this command is through git on juha.no.
            result=subprocess.run(['ssh',a.archive_host,
                'cd /home/j/src/pansy_pmse_sparse && /home/j/venv/pansy-rvd/bin/python -m pmse.web_retention '
                '--root /mnt/shovel/rvd/pmse_sparse16 --receipt /mnt/shovel/rvd/revontuli-receipt.h5 --apply'],
                capture_output=True,text=True,cwd=Path(__file__).resolve().parents[1],timeout=120)
            if result.returncode: raise RuntimeError(result.stderr)
            print(result.stdout,flush=True)
        except Exception:
            import traceback;traceback.print_exc()
            if not a.loop: raise
        if not a.loop: break
        time.sleep(30)


if __name__=='__main__': main()
