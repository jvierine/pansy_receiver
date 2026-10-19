"""Verify Revontuli receipts before expiring the seven-day secondary copy.

Runs on juha.no. The source receipt contains SHA256 values computed from the
retained Revontuli files, not merely from the secondary copy. Files spanning the
cutoff get an additional 600-second guard. No raw or historical XC is in scope.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import h5py


def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def expire(root,receipt,now=None,days=7,dry_run=True):
    root=Path(root).resolve();now=time.time() if now is None else now
    if root!=Path('/mnt/shovel/rvd/pmse_sparse16') and not dry_run:
        raise ValueError('Refusing expiration outside the explicitly scoped secondary archive')
    cutoff=now-days*86400-600;removed=[];blocked=[]
    with h5py.File(receipt,'r') as h:
        for path in root.glob('*T*/pmse@*.h5'):
            sec=int(path.stem.split('@')[1])
            if sec>=cutoff:continue
            relative=str(path.relative_to(root))
            if relative not in h['files']:
                blocked.append(relative);continue
            expected=h['files'][relative].attrs['sha256']
            if sha256(path)!=expected:
                blocked.append(relative);continue
            removed.append(relative)
            if not dry_run:path.unlink()
    return dict(expired_files=removed,unverified_files_preserved=blocked,dry_run=dry_run,
                bytes=sum(p.stat().st_size for p in root.glob('*T*/pmse@*.h5')),
                retention_days=days,checked_unix_s=now)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/mnt/shovel/rvd/pmse_sparse16'))
    p.add_argument('--receipt',type=Path,required=True);p.add_argument('--apply',action='store_true')
    a=p.parse_args();print(json.dumps(expire(a.root,a.receipt,dry_run=not a.apply)))
