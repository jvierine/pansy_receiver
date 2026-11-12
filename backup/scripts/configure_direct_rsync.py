#!/usr/bin/env python3
"""Migrate only the metadata pull configuration; never touch receiver services."""
import argparse
from datetime import datetime,timezone
from pathlib import Path
import shutil


def configure(config,unit_dir,bwlimit_kb=128):
    text=config.read_text()
    required=('JUMP_HOST=','PANSY_SSH_HOST=','PANSY_SSH_PORT=')
    if not all(key in text for key in required): raise ValueError('Missing direct SSH configuration')
    if bwlimit_kb < 0: raise ValueError('Bandwidth budget must be nonnegative')
    settings={'RSYNC_DIRECT_SSH':'1','RSYNC_BWLIMIT_KB':str(bwlimit_kb),'RSYNC_TOTAL_BWLIMIT_KB':str(bwlimit_kb),
              'RSYNC_RECENT_MINUTES':'20','RSYNC_HISTORY_SLICE_SECONDS':'180'}
    lines=[]
    for line in text.splitlines():
        key=line.partition('=')[0].strip()
        if key not in settings: lines.append(line)
    lines+=['','# Direct metadata pull: recent files first, bounded historical catch-up.']
    lines += [key+'='+value for key,value in settings.items()]
    backup=config.with_name(config.name+'.pre-direct-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')+'.bak')
    shutil.copy2(config,backup)
    temporary=config.with_suffix(config.suffix+'.partial')
    temporary.write_text('\n'.join(lines)+'\n');shutil.copymode(config,temporary);temporary.replace(config)
    dropin=unit_dir/'pansy-backup-rsync.service.d'/'direct-outbound.conf'
    dropin.parent.mkdir(parents=True,exist_ok=True)
    dropin.write_text('[Unit]\nAfter=\nAfter=network-online.target\nWants=\nWants=network-online.target\n')
    print('Backed up metadata transfer configuration to',backup)
    print(f'Enabled direct outbound SSH, metadata budget {bwlimit_kb} KiB/s and recent-first catch-up. No services restarted.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True,type=Path)
    parser.add_argument('--unit-dir',required=True,type=Path)
    parser.add_argument('--bwlimit-kb',type=int,default=128,help='Total metadata budget in KiB/s (default: 128; 0: unlimited)')
    args=parser.parse_args();configure(args.config,args.unit_dir,args.bwlimit_kb)
