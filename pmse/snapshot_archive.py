"""Make stable transfer files under the native writer's advisory lock."""
import argparse
import fcntl
from pathlib import Path
import shutil


def snapshot(source,destination,reserve_gib=200):
    source=Path(source);destination=Path(destination)
    destination.mkdir(parents=True,exist_ok=True);copied=0
    with (source/'.writer.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_SH)
        for path in source.rglob('*.h5'):
            target=destination/path.relative_to(source);st=path.stat()
            if target.exists() and target.stat().st_size==st.st_size and target.stat().st_mtime_ns==st.st_mtime_ns:continue
            if shutil.disk_usage(destination).free<reserve_gib*2**30+st.st_size:
                raise RuntimeError('Snapshot free-space reserve reached; native source remains untouched')
            target.parent.mkdir(parents=True,exist_ok=True);tmp=target.with_suffix('.h5.partial')
            shutil.copy2(path,tmp);tmp.replace(target);copied+=1
    return copied


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path('/media/analysis/metadata/pmse_sparse16'))
    p.add_argument('--destination',type=Path,default=Path('/media/analysis/metadata/pmse_sparse16_snapshot'))
    a=p.parse_args();print(snapshot(a.source,a.destination))
