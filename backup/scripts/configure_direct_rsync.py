#!/usr/bin/env python3
"""Migrate only the metadata pull configuration; never touch receiver services."""
import argparse
import json
from datetime import datetime,timezone
from pathlib import Path
import shutil


def configure(config,unit_dir,bwlimit_kb=128,metadata_mount=None,mirror_mount=None,service_repo=None):
    text=config.read_text()
    required=('JUMP_HOST=','PANSY_SSH_HOST=','PANSY_SSH_PORT=')
    if not all(key in text for key in required): raise ValueError('Missing direct SSH configuration')
    if bwlimit_kb < 0: raise ValueError('Bandwidth budget must be nonnegative')
    settings={'RSYNC_DIRECT_SSH':'1','RSYNC_BWLIMIT_KB':str(bwlimit_kb),'RSYNC_TOTAL_BWLIMIT_KB':str(bwlimit_kb)}
    for key,value in (('LOCAL_METADATA_MOUNT',metadata_mount),('LOCAL_MIRROR_METADATA_MOUNT',mirror_mount)):
        if value is not None:
            if not Path(value).is_absolute() or str(Path(value)) == '/': raise ValueError('Archive mount must be an absolute non-root path')
            settings[key]=str(value)
    scripts={'pansy-backup-rsync':'pansy_backup_rsync.sh',
             'pansy-local-metadata-mirror':'pansy_local_metadata_mirror.sh',
             'pansy-backup-web':'pansy_backup_status_loop.sh'}
    repo=None
    if service_repo is not None:
        repo=Path(service_repo).resolve()
        if '%' in str(repo): raise ValueError('Service checkout path must not contain systemd specifiers')
        for script in scripts.values():
            path=repo/'backup/scripts'/script
            if not path.is_file(): raise ValueError(f'Missing service script: {path}')
    obsolete={'RSYNC_RECENT_MINUTES','RSYNC_HISTORY_SLICE_SECONDS'}
    lines=[]
    for line in text.splitlines():
        key=line.partition('=')[0].strip()
        if key not in settings and key not in obsolete and not line.startswith('# Direct metadata pull:'): lines.append(line)
    while lines and not lines[-1].strip(): lines.pop()
    lines+=['','# Direct metadata pull: ordinary rsync with a shared bandwidth budget.']
    lines += [key+'='+value for key,value in settings.items()]
    backup=config.with_name(config.name+'.pre-direct-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')+'.bak')
    shutil.copy2(config,backup)
    temporary=config.with_suffix(config.suffix+'.partial')
    temporary.write_text('\n'.join(lines)+'\n');shutil.copymode(config,temporary);temporary.replace(config)
    dropin=unit_dir/'pansy-backup-rsync.service.d'/'direct-outbound.conf'
    dropin.parent.mkdir(parents=True,exist_ok=True)
    dropin.write_text('[Unit]\nAfter=\nAfter=network-online.target\nWants=\nWants=network-online.target\n')
    if repo is not None:
        for unit,script in scripts.items():
            directory=unit_dir/(unit+'.service.d')
            directory.mkdir(parents=True,exist_ok=True)
            (directory/'transfer-checkout.conf').write_text(
                '[Service]\nEnvironment='+json.dumps('PANSY_RECEIVER_REPO='+str(repo))+
                '\nExecStart=\nExecStart=/usr/bin/bash '+json.dumps(str(repo/'backup/scripts'/script))+'\n')
    print('Backed up metadata transfer configuration to',backup)
    print(f'Enabled ordinary rsync over direct outbound SSH, metadata budget {bwlimit_kb} KiB/s; no wall-clock cutoff. No services restarted.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True,type=Path)
    parser.add_argument('--unit-dir',required=True,type=Path)
    parser.add_argument('--bwlimit-kb',type=int,default=128,help='Total metadata budget in KiB/s (default: 128; 0: unlimited)')
    parser.add_argument('--metadata-mount',help='Required mount for the primary archive')
    parser.add_argument('--mirror-mount',help='Required mount for the local mirror')
    parser.add_argument('--service-repo',type=Path,help='Use this checkout for pull, local mirror and status publisher')
    args=parser.parse_args();configure(args.config,args.unit_dir,args.bwlimit_kb,args.metadata_mount,args.mirror_mount,args.service_repo)
