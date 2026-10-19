"""Explicitly scoped retirement of this RVD experiment's obsolete 8ms data.

Never targets raw, XC, historical PMSE imaging, meteor or calibration archives.
All files are audited before any deletion. A v2 record aborts the entire action.
"""
import argparse
import json
from pathlib import Path
import h5py

ALLOWED={
    '/media/analysis/metadata/pmse_sparse',
    '/media/analysis/metadata/pmse_sparse_snapshot',
    '/media/analysis/metadata/pmse_sparse_trial',
    '/media/analysis/metadata/pmse_sparse_v1',
    '/media/analysis/metadata/pmse_sparse_rust_validation',
    '/mnt/data/juha/pansy/metadata/pmse_sparse',
    '/mnt/shovel/rvd/pmse_sparse',
}


def retire(roots,apply=False):
    files=[];records=0
    for value in roots:
        root=Path(value)
        if str(root) not in ALLOWED or root.is_symlink():raise ValueError('Unapproved retirement root')
        for path in sorted(root.rglob('*.h5')):
            if path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):raise ValueError('Unexpected symlink')
            if path.name!='dmd_properties.h5' and not path.name.startswith('pmse@'):raise ValueError('Unexpected HDF5 file')
            if path.name!='dmd_properties.h5':
                with h5py.File(path,'r') as h:
                    for key,g in h.items():
                        if not key.isdigit():continue
                        schema=g['schema_version'][()]
                        if isinstance(schema,bytes):schema=schema.decode()
                        if str(schema)!='pansy_pmse_sparse_v1':raise ValueError(f'Refusing non-v1 record in {path}')
                        records+=1
            files.append((path,path.stat().st_size))
    result=dict(dry_run=not apply,files=[str(p) for p,_ in files],bytes=sum(n for _,n in files),records=records)
    if apply:
        for path,_ in files:path.unlink()
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,action='append',required=True);p.add_argument('--apply',action='store_true')
    a=p.parse_args();print(json.dumps(retire(a.root,a.apply)))
