"""Oldest-first rectangular reprocessing; recoverable atomic Hann replacement.

Run only in the dedicated single-core service. Acquisition is never controlled
here. Existing rectangular records are immutable; historical Hann groups are
backed up outside the live archive before an atomically replaced file is exposed.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import resource
import shutil
import tempfile
import time
from dataclasses import replace

import h5py
import numpy as np

from .sparse_record import Config, closed_modes, process_mode, read_complete, write_record


def native_path(root, key):
    sec = int(key) // 1_000_000
    hour = datetime.fromtimestamp(sec // 3600 * 3600, timezone.utc).strftime('%Y-%m-%dT%H-%M-%S')
    return Path(root) / hour / f'pmse@{sec // 600 * 600}.h5'


def is_rectangular(config, weights):
    cfg = json.loads(config)
    w = np.asarray(weights)
    if cfg.get('window_kind') != 'rectangular' or w.ndim != 2 or not w.shape[0]:
        return False
    a = w.max(axis=1)
    return bool(np.isfinite(w).all() and np.all(a > 0) and
                np.all((w == 0) | (w == a[:, None])) and
                np.allclose(np.sum(w.astype(float)**2, axis=1), 1., atol=2e-6))


def rectangular_record(root, key):
    path = native_path(root, key)
    if not path.exists():
        return False
    with (Path(root) / '.writer.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_SH)
        with h5py.File(path, 'r') as h:
            if str(key) not in h:
                return False
            g = h[str(key)]
            return (int(g['complete'][()]) == 1 and
                    is_rectangular(g['config_json'][()], g['look_window'][()]))


def sync_file(path):
    with open(path, 'rb') as stream:
        os.fsync(stream.fileno())


def sync_directory(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def promote(root, key, record, backup_root, reserve_gib=200, max_archive_gib=500):
    """Validate a real native write, preserve old group, atomically replace it.

    The same advisory lock is used by the live producer and snapshot process.
    Other groups in the chunk, including incomplete ones, are copied unchanged.
    A crash before os.replace leaves the old file; after it, the new complete
    file is authoritative. A durable backup precedes replacement in either case.
    """
    root = Path(root).resolve()
    backup_root = Path(backup_root).resolve()
    if backup_root.is_relative_to(root) or root.is_relative_to(backup_root):
        raise ValueError('Recovery backups must be separate from the live archive')
    if not is_rectangular(record['config_json'], record['look_window']):
        raise ValueError('Replacement must have validated rectangular weights')
    if record.get('quality') == 'insufficient_looks' or not np.isfinite(record['spectra']).all():
        raise ValueError('Do not replace historical data with an invalid/insufficient record')
    if int(record['mode_start_sample']) != int(key):
        raise ValueError('Replacement start/key mismatch')
    root.mkdir(parents=True, exist_ok=True)
    path = native_path(root, key)
    if shutil.disk_usage(root).free < reserve_gib * 2**30 + 4 * (path.stat().st_size if path.exists() else 0):
        raise RuntimeError('Backfill free-space reserve reached; original retained')
    if sum(p.stat().st_size for p in root.glob('*T*/pmse@*.h5')) >= max_archive_gib * 2**30:
        raise RuntimeError('Backfill archive-size budget reached; original retained')
    with tempfile.TemporaryDirectory(prefix='pmse-rectangular-stage-', dir=root.parent) as stage:
        write_record(stage, key, record)
        saved = read_complete(stage, key)
        for name in ('spectra', 'sparse_index', 'look_window'):
            np.testing.assert_array_equal(saved[name], record[name])
        candidate = native_path(stage, key)
        temporary = None
        with (root / '.writer.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            exists = False
            if path.exists():
                with h5py.File(path, 'r') as src:
                    exists = str(key) in src
                    if exists:
                        old = src[str(key)]
                        if int(old['complete'][()]) != 1:
                            raise ValueError('Incomplete live group is owned by its writer; retry later')
                        if is_rectangular(old['config_json'][()], old['look_window'][()]):
                            return 'already_rectangular'
                        if set(old) != set(saved):
                            raise ValueError('Native field layout differs; original retained')
                        if int(saved['mode_stop_sample_exclusive']) < int(old['mode_stop_sample_exclusive'][()]):
                            raise ValueError('Replacement truncates original mode; original retained')
                        if float(saved['coverage_seconds']) < .95 * float(old['coverage_seconds'][()]):
                            raise ValueError('Replacement loses coverage; original retained')
                        for name in ('range_sample', 'doppler_frequency_hz'):
                            np.testing.assert_array_equal(old[name][()], saved[name])
                        backup = backup_root / path.parent.name / f'{key}.h5'
                        backup.parent.mkdir(parents=True, exist_ok=True)
                        if not backup.exists():
                            backup_tmp = backup.with_suffix('.h5.partial')
                            with h5py.File(backup_tmp, 'w') as dst:
                                src.copy(str(key), dst)
                                dst.attrs['original_native_file'] = str(path)
                            sync_file(backup_tmp)
                            os.replace(backup_tmp, backup)
                            sync_directory(backup.parent)
                        with h5py.File(backup, 'r') as recovery:
                            np.testing.assert_array_equal(recovery[str(key)]['spectra'][()], old['spectra'][()])
                            if int(recovery[str(key)]['complete'][()]) != 1:
                                raise ValueError('Invalid recovery copy; original retained')
                fd, name = tempfile.mkstemp(prefix=path.name + '.rectangular-', suffix='.partial', dir=path.parent)
                os.close(fd)
                temporary = Path(name)
                try:
                    with h5py.File(path, 'r') as src, h5py.File(candidate, 'r') as new, h5py.File(temporary, 'w') as dst:
                        for name, value in src.attrs.items():
                            dst.attrs[name] = value
                        for name in src:
                            if name != str(key):
                                src.copy(name, dst)
                        new.copy(str(key), dst)
                    sync_file(temporary)
                    os.replace(temporary, path)
                    sync_directory(path.parent)
                finally:
                    if temporary.exists():
                        temporary.unlink()
                return 'replaced_hann' if exists else 'created'
        # No bucket existed: the ordinary writer creates its channel/properties
        # under its own lock. A racing live writer wins idempotently.
        return 'created' if write_record(root, key, record) else 'already_rectangular'


def pending(root, choices):
    """Strict oldest-first ordering; skip only already-rectangular records."""
    return [(start, last) for start, last in sorted(choices) if not rectangular_record(root, start)]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw-root', type=Path, default=Path('/media/archive'))
    p.add_argument('--tx-root', type=Path, default=Path('/media/analysis/metadata/tx'))
    p.add_argument('--mode-root', type=Path, default=Path('/media/analysis/metadata/mesomode'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--backup-root', type=Path, required=True)
    p.add_argument('--reserve-free-gib', type=float, default=200)
    p.add_argument('--max-archive-gib', type=float, default=500)
    p.add_argument('--replace-hann', action='store_true', help='Explicitly authorize recoverable historical replacement')
    a = p.parse_args()
    if not a.replace_hann or a.reserve_free_gib < 0 or a.max_archive_gib <= 0:
        p.error('Require --replace-hann and a nonnegative storage reserve')
    a.output.mkdir(parents=True, exist_ok=True)
    cfg = replace(Config(), nfft=128, complementary_pairs=True, velocity_limit_mps=None)
    with (a.output / '.rectangular-backfill.lock').open('a') as worker_lock:
        fcntl.flock(worker_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        print('Indexing closed modes from the oldest currently retained raw samples', flush=True)
        choices = pending(a.output, closed_modes(a.mode_root, a.tx_root, a.raw_root, 0))
        print(json.dumps(dict(pending=len(choices), oldest=choices[0][0] if choices else None)), flush=True)
        failures = 0
        for number, (start, last) in enumerate(choices, 1):
            if rectangular_record(a.output, start):
                continue
            print(f'Rectangular backfill {number}/{len(choices)}: {start}..{last}', flush=True)
            begun = time.monotonic()
            try:
                record = process_mode(a.raw_root, a.tx_root, start, last, cfg, backend='rust')
                record['processing_wall_seconds'] = time.monotonic() - begun
                record['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
                result = promote(a.output, start, record, a.backup_root, a.reserve_free_gib, a.max_archive_gib)
                print(json.dumps(dict(start=start, result=result, selected=len(record['sparse_index']),
                                      seconds=time.monotonic() - begun)), flush=True)
            except RuntimeError:
                # Resource guards pause the whole worker; do not process every
                # remaining raw mode only to fail its publication again.
                raise
            except Exception:
                import traceback
                traceback.print_exc()
                failures += 1
        if failures:
            raise RuntimeError(f'{failures} modes failed; their original data remain intact')
        print('Rectangular backfill complete for indexed retained raw buffer', flush=True)


if __name__ == '__main__':
    main()
