"""Outbound-only native HDF5 publication and archive pulls for PANSY RVD.

Syowa pushes stable snapshots to juha.no. Revontuli pulls from juha.no.
No SSH port forwarding or inbound connection into Revontuli is used.
"""
import argparse
import json
from pathlib import Path
import shlex
import subprocess
import time

import h5py

from .replicate_sparse import rsync, save_metrics, historical_turn
from .web_retention import sha256


SSH = ['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=20',
       '-o', 'StrictHostKeyChecking=yes', '-o', 'ControlMaster=auto',
       '-o', 'ControlPersist=60', '-o', 'ControlPath=/tmp/pansy-rvd-outbound-%C']


def check_archive(host, root):
    command = 'mountpoint -q /mnt/shovel && test -d ' + shlex.quote(root)
    subprocess.run(SSH + [host, command], check=True, timeout=30)


def source_files(root):
    return sorted(root.glob('20*T*/pmse@*.h5'), reverse=True)


def pending_files(root, receipt):
    """Compare stable snapshot versions with completed publication receipts."""
    files = source_files(root)
    if not receipt.exists():
        return files
    pending = []
    with h5py.File(receipt, 'r') as h:
        for path in files:
            group = h.get('files/' + str(path.relative_to(root)))
            st = path.stat()
            if (group is None or int(group.attrs.get('mtime_ns', -1)) != st.st_mtime_ns
                    or int(group.attrs.get('bytes', -1)) != st.st_size):
                pending.append(path)
    return pending


def update_receipt(root, receipt, files):
    receipt.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(receipt, 'a') as h:
        groups = h.require_group('files')
        for path in files:
            st = path.stat()
            group = groups.require_group(str(path.relative_to(root)))
            if (int(group.attrs.get('mtime_ns', -1)) != st.st_mtime_ns
                    or int(group.attrs.get('bytes', -1)) != st.st_size):
                group.attrs.update(sha256=sha256(path), bytes=st.st_size,
                                   mtime_ns=st.st_mtime_ns)
        h.attrs.update(source=str(root), checked_unix_s=time.time(),
                       generator='pmse/archive_sync.py')


def publication_selection(pending, now, turn, frontier=None):
    """Keep live chunks fresh while filling earlier holes before old backfill."""
    # Do not let a missing chunk age out of the live queue during catch-up.
    today=int(now)//86400*86400
    cutoff=min(today,now-3*3600)
    recent=[p for p in pending if int(p.stem.split('@')[1])>=cutoff]
    historical=[p for p in pending if p not in recent]
    selected=recent[:1]
    catchup=recent[-1:] if len(recent)>1 else []
    selected+=catchup
    # Historical work keeps a share of the link, but cannot strand live holes.
    if historical and (not catchup or turn%3==0):
        names=[str(p) for p in historical]
        choice=historical_turn(names,turn,processed_frontier_us=frontier)[0]
        selected.append(Path(choice))
    return selected


def publish_once(a, turn):
    from .snapshot_archive import snapshot
    snapshot(a.source, a.root)
    try:
        from .processing_status import status
        temporary = a.root / 'processing-status.json.partial'
        temporary.write_text(json.dumps(status()))
        temporary.replace(a.root / 'processing-status.json')
    except Exception as error:
        print(f'Optional processing status unavailable: {error}', flush=True)
    check_archive(a.host, a.archive_root)
    destination = f'{a.host}:{a.archive_root}/'
    metadata = [name for name in ['dmd_properties.h5', 'processing-status.json']
                if (a.root / name).is_file()]
    if metadata:
        rsync(str(a.root) + '/', destination, a.bandwidth, shlex.join(SSH),
              files_from=metadata, preserve_times=False, force=True,
              timeout=60, partial_ok=True)
    pending = pending_files(a.root, a.receipt)
    if not pending:
        return
    frontier = None
    try:
        frontier = json.loads((a.root / 'processing-status.json').read_text()).get(
            'backfill', {}).get('current_start_us')
    except (OSError, ValueError):
        pass
    selected=publication_selection(pending,time.time(),turn,frontier)
    for path in selected:
        relative = str(path.relative_to(a.root))
        print(f'Publishing stable snapshot: {relative}', flush=True)
        values = rsync(str(a.root) + '/', destination, a.bandwidth,
                       shlex.join(SSH), files_from=[relative],
                       preserve_times=False, force=True, timeout=120,
                       partial_ok=True)
        save_metrics(a.metrics, 'syowa_to_juha_no', values)
        if not values['partial']:
            update_receipt(a.root, a.receipt, [path])


def pull_once(a):
    check_archive(a.host, a.archive_root)
    values = rsync(f'{a.host}:{a.archive_root}/', str(a.root) + '/',
                   a.bandwidth, shlex.join(SSH), checksum=True,
                   timeout=600, partial_ok=True)
    save_metrics(a.metrics, 'juha_no_to_revontuli', values)
    if values['partial']:
        return
    update_receipt(a.root, a.receipt, source_files(a.root))
    rsync(str(a.receipt), f'{a.host}:/mnt/shovel/rvd/revontuli-receipt.h5',
          remote_shell=shlex.join(SSH), preserve_times=False, force=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['publish', 'pull'])
    parser.add_argument('--host', default='j@juha.no')
    parser.add_argument('--archive-root', default='/mnt/shovel/rvd/pmse_sparse16')
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--source', type=Path,
                        default=Path('/media/analysis/metadata/pmse_sparse16'))
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--metrics', type=Path, required=True)
    parser.add_argument('--bandwidth', type=int, required=True,
                        help='Transfer cap in KiB/s; source publication remains capped at 10')
    parser.add_argument('--loop', action='store_true')
    a = parser.parse_args()
    if a.bandwidth <= 0 or (a.mode == 'publish' and a.bandwidth > 10):
        parser.error('Require positive bandwidth; Syowa publication may not exceed 10 KiB/s')
    a.root.mkdir(parents=True, exist_ok=True)
    turn = 1
    while True:
        try:
            if a.mode == 'publish':
                publish_once(a, turn)
                turn += 1
            else:
                pull_once(a)
        except Exception:
            import traceback
            traceback.print_exc()
            if not a.loop:
                raise
        if not a.loop:
            break
        time.sleep(30)


if __name__ == '__main__':
    main()
