"""Make stable transfer files under the native writer's advisory lock."""
import argparse
import fcntl
from pathlib import Path
import shutil
import os
import h5py


def compact_copy(source,target):
    """Losslessly copy live HDF5 objects, omitting abandoned dataset storage."""
    with h5py.File(source,'r') as src,h5py.File(target,'w') as dst:
        for name,value in src.attrs.items():dst.attrs[name]=value
        for name in src:src.copy(name,dst,name=name)


def snapshot(source,destination,reserve_gib=200,force=False):
    source=Path(source);destination=Path(destination)
    destination.mkdir(parents=True,exist_ok=True);copied=0
    with (source/'.writer.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_SH)
        for path in source.rglob('*.h5'):
            target=destination/path.relative_to(source);st=path.stat()
            if not force and target.exists() and target.stat().st_mtime_ns==st.st_mtime_ns:continue
            if shutil.disk_usage(destination).free<reserve_gib*2**30+st.st_size:
                raise RuntimeError('Snapshot free-space reserve reached; native source remains untouched')
            target.parent.mkdir(parents=True,exist_ok=True);tmp=target.with_suffix('.h5.partial')
            compact_copy(path,tmp)
            os.utime(tmp,ns=(st.st_atime_ns,st.st_mtime_ns))
            tmp.replace(target);copied+=1
    return copied


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path('/media/analysis/metadata/pmse_sparse16'))
    p.add_argument('--destination',type=Path,default=Path('/media/analysis/metadata/pmse_sparse16_snapshot'))
    p.add_argument('--compact-all',action='store_true',help='Repack existing snapshots once; native archive is untouched')
    a=p.parse_args();print(snapshot(a.source,a.destination,force=a.compact_all))
