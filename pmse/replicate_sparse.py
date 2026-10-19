"""Replicate native HDF5 without ownership changes; persist transfer evidence."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import time
import h5py
from .web_retention import sha256


def rsync(source,destination,bandwidth=None,remote_shell=None,files_from=None,preserve_times=True,force=False):
    cmd=['rsync','-rt','--no-perms','--no-owner','--no-group','--partial','--delay-updates',
         '--timeout=120','--stats','--exclude=/.writer.lock']
    if bandwidth: cmd.append(f'--bwlimit={bandwidth}')
    if remote_shell: cmd+=['-e',remote_shell]
    if files_from is not None:cmd+=['--files-from=-']
    if not preserve_times:cmd+=['--no-times','--omit-dir-times','--size-only']
    if force:cmd+=['--ignore-times']
    began=time.monotonic()
    result=subprocess.run(cmd+[source,destination],text=True,capture_output=True,
                          input=None if files_from is None else '\n'.join(files_from)+'\n',
                          env=dict(os.environ,LC_ALL='C'),timeout=900)
    print(result.stdout,flush=True)
    if result.returncode: raise RuntimeError(f'rsync {result.returncode}: {result.stderr}')
    elapsed=time.monotonic()-began
    match=re.search(r'Literal data:\s*([\d,]+)',result.stdout)
    payload=int(match[1].replace(',','')) if match else 0
    return dict(finished_unix_s=time.time(),elapsed_seconds=elapsed,literal_bytes=payload,
                payload_bytes_per_second=payload/max(elapsed,1e-9),bandwidth_cap_kib_s=bandwidth or 0)


def save_metrics(path,leg,values):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with h5py.File(path,'a') as h:
        g=h.require_group('runs').create_group(f'{time.time_ns()}-{leg}')
        g.attrs.update(leg=leg,generator='pmse/replicate_sparse.py')
        for key,value in values.items(): g[key]=value
        h.flush()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/mnt/data/juha/pansy/metadata/pmse_sparse'))
    p.add_argument('--metrics',type=Path,default=Path('/mnt/data/juha/pansy/rvd-transfer.h5'))
    p.add_argument('--syowa-kib-s',type=int,default=10)
    p.add_argument('--archive-host',default='j@juha.no')
    p.add_argument('--archive-root',default='/mnt/shovel/rvd/pmse_sparse')
    p.add_argument('--loop',action='store_true');a=p.parse_args();a.root.mkdir(parents=True,exist_ok=True)
    while True:
        try:
            subprocess.run(['ssh','-p','2222','-o','BatchMode=yes','-o','ConnectTimeout=15','radar@127.0.0.1',
                'cd /home/radar/src/git/pansy_pmse_sparse && /home/radar/testenv/bin/python -m pmse.snapshot_archive'],check=True,timeout=120)
            values=rsync('radar@127.0.0.1:/media/analysis/metadata/pmse_sparse_snapshot/',str(a.root)+'/',
                         a.syowa_kib_s,'ssh -p 2222 -o BatchMode=yes -o ConnectTimeout=15')
            save_metrics(a.metrics,'syowa_to_revontuli',values)
            # Fail closed if shovel is unmounted; never fill the server root disk.
            subprocess.run(['ssh',a.archive_host,'mountpoint -q /mnt/shovel && test -d /mnt/shovel/rvd/pmse_sparse'],check=True,timeout=30)
            cutoff=time.time()-7*86400-600
            keep=['dmd_properties.h5']+[str(path.relative_to(a.root)) for path in a.root.glob('*T*/pmse@*.h5')
                if int(path.stem.split('@')[1])>=cutoff]
            values=rsync(str(a.root)+'/',f'{a.archive_host}:{a.archive_root}/',files_from=keep,preserve_times=False)
            save_metrics(a.metrics,'revontuli_to_shovel',values)
            # Only changed/new source files are rehashed. The receipt proves
            # preservation on Revontuli before any secondary expiration.
            receipt=a.metrics.with_name('rvd-revontuli-receipt.h5')
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
                '--receipt /mnt/shovel/rvd/revontuli-receipt.h5 --apply'],
                capture_output=True,text=True,cwd=Path(__file__).resolve().parents[1],timeout=120)
            if result.returncode: raise RuntimeError(result.stderr)
            print(result.stdout,flush=True)
        except Exception:
            import traceback;traceback.print_exc()
            if not a.loop: raise
        if not a.loop: break
        time.sleep(300)


if __name__=='__main__': main()
