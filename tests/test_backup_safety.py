"""Failure-injection checks for the ordinary rsync pull and status publisher."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import random
import signal
import subprocess
import time

import h5py
import pytest

REPO = Path(__file__).resolve().parents[1]
PULL = REPO / "backup/scripts/pansy_backup_rsync.sh"
GUARD = REPO / "backup/scripts/pansy_storage_guard.sh"
MIRROR = REPO / "backup/scripts/pansy_local_metadata_mirror.sh"
spec = importlib.util.spec_from_file_location("backup_status", REPO / "backup/scripts/pansy_backup_status.py")
status = importlib.util.module_from_spec(spec)
spec.loader.exec_module(status)


def wait_until(predicate, timeout=10):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(0.05)
    raise AssertionError("condition did not become true")


def stop(process):
    if process.poll() is None:
        os.killpg(process.pid, signal.SIGTERM)
    process.communicate(timeout=6)


def pull_config(tmp_path, extra=""):
    source = tmp_path / "source"
    source.mkdir(exist_ok=True)
    config = tmp_path / "config.env"
    config.write_text(
        f'LOCAL_METADATA_ROOT="{tmp_path}/data"\nSTATE_DIR="{tmp_path}/state"\n'
        f'REMOTE_METADATA_ROOT="{source}"\nMETADATA_CHANNELS="cut"\n'
        'RSYNC_TOTAL_BWLIMIT_KB=128\nRSYNC_SLEEP_SECONDS=0.1\nLOCAL_TUNNEL_PORT=2222\n'
        # Obsolete settings must never revive the old wall-clock cutoff.
        'RSYNC_RECENT_MINUTES=20\nRSYNC_HISTORY_SLICE_SECONDS=1\n' + extra
    )
    return dict(os.environ, PANSY_BACKUP_CONFIG=str(config)), source


def test_empty_directories_and_mtimes_do_not_fake_fresh_data(tmp_path):
    hour = tmp_path / "2026-10-31T01-00-00"
    hour.mkdir()
    timestamp = 1793411460.25
    with h5py.File(hour / "cut@1793411460.h5", "w") as handle:
        handle.create_group(str(int(timestamp * 1e6)))
    (tmp_path / "2026-11-13T04-00-00").mkdir()
    newest, day, count, size = status.channel_summary(tmp_path)
    assert newest == timestamp
    assert day == "2026-10-31"
    assert count is None and size is None  # No invented file counts from directories.


def test_recent_only_success_is_not_full_archive_verification(tmp_path):
    path = tmp_path / "status.tsv"
    path.write_text("cut\tstart\tfinish\t0\ntx\tstart\tfinish\t0\tfull_channel\n")
    records = status.read_rsync_status(path)
    assert records["cut"]["last_rsync_exit_code"] is None
    assert records["cut"]["archive_verified_at_pass_start"] is None
    assert records["tx"]["last_rsync_exit_code"] == 0
    assert records["tx"]["archive_verified_at_pass_start"] == "start"


def test_dead_worker_is_not_reported_active(tmp_path):
    (tmp_path / "cut.tsv").write_text("cut\t2026-11-13T00:00:00Z\t999999999\n")
    assert status.read_active_rsync(tmp_path, time.time()) == {}


def test_status_document_exposes_unknown_inventory_and_missing_mount(tmp_path, monkeypatch):
    config = tmp_path / 'status.env'
    config.write_text(f'LOCAL_METADATA_ROOT={tmp_path}/data\nLOCAL_MIRROR_METADATA_ROOT={tmp_path}/mirror\n'
                      f'STATE_DIR={tmp_path}/state\nWEB_BUILD_DIR={tmp_path}/web\n'
                      f'LOCAL_METADATA_MOUNT={tmp_path}/missing-mount\nMETADATA_CHANNELS=cut\n')
    monkeypatch.setenv('PANSY_BACKUP_CONFIG', str(config))
    status.main()
    document = json.loads((tmp_path / 'web/status.json').read_text())
    assert document['channels'][0]['age_seconds'] is None
    assert document['channels'][0]['entry_count'] is None
    assert document['local_mirror']['channels'][0]['entry_count_delta'] is None
    assert document['storage'][0]['mount_available'] is False
    assert document['storage'][0]['critical'] is True
    assert not list((tmp_path / 'web').glob('*.tmp'))


def test_mount_missing_fails_before_creating_archive_or_starting_rsync(tmp_path):
    env, _ = pull_config(tmp_path, f'LOCAL_METADATA_MOUNT="{tmp_path}/missing-mount"\n')
    result = subprocess.run(["bash", str(PULL)], env=env, capture_output=True, text=True, timeout=3)
    assert result.returncode != 0
    assert "Required archive mount is unavailable" in result.stderr
    assert not (tmp_path / "data").exists()


def test_mount_guard_rejects_destination_outside_mount(tmp_path):
    mount = tmp_path / "mount"
    mount.mkdir()
    result = subprocess.run(
        ["bash", "-c", 'source "$1"; mountpoint() { return 0; }; pansy_require_mount "$2" "$3"',
         "_", str(GUARD), str(tmp_path / "outside"), str(mount)], capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "not under required mount" in result.stderr


def test_duplicate_instance_cannot_clear_status_or_start_another_transfer(tmp_path):
    env, _ = pull_config(tmp_path)
    binaries = tmp_path / "bin"
    binaries.mkdir()
    mock = binaries / "rsync"
    mock.write_text("#!/bin/sh\nsleep 30\n")
    mock.chmod(0o755)
    env["PATH"] = str(binaries) + os.pathsep + env["PATH"]
    first = subprocess.Popen(["bash", str(PULL)], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
    try:
        active = tmp_path / "state/rsync_active/cut.tsv"
        wait_until(active.exists, 3)
        original = active.read_text()
        second = subprocess.run(["bash", str(PULL)], env=env, capture_output=True, text=True, timeout=3)
        assert second.returncode != 0
        assert "already running" in second.stderr
        assert active.read_text() == original
        assert first.poll() is None
    finally:
        stop(first)


def test_real_rsync_interruption_resume_and_no_source_or_destination_deletion(tmp_path):
    """Run the real rsync protocol over a local stand-in for SSH, including restart."""
    env, source = pull_config(tmp_path)
    binaries = tmp_path / "bin"
    binaries.mkdir()
    ssh = binaries / "ssh"
    ssh.write_text(
        "#!/usr/bin/env python3\nimport os,sys\n"
        "start=sys.argv.index('rsync')\nos.execvp('rsync',sys.argv[start:])\n"
    )
    ssh.chmod(0o755)
    env["PATH"] = str(binaries) + os.pathsep + env["PATH"]
    source_hour = source / "cut/2026-10-31T01-00-00"
    source_hour.mkdir(parents=True)
    payload = random.Random(12).randbytes(1024 * 1024)
    original = source_hour / "cut@1793411460.h5"
    original.write_bytes(payload)
    expected = hashlib.sha256(payload).hexdigest()
    destination = tmp_path / "data/cut/2026-10-31T01-00-00/cut@1793411460.h5"
    keep = tmp_path / "data/cut/local-only.h5"
    keep.parent.mkdir(parents=True)
    keep.write_bytes(b"never delete local-only data")
    first = subprocess.Popen(["bash", str(PULL)], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
    try:
        wait_until(lambda: (tmp_path / "state/rsync_active/cut.tsv").exists(), 3)
        time.sleep(1.5)
        assert not destination.exists()
    finally:
        stop(first)
    partial = destination.parent / ".rsync-partial" / destination.name
    assert partial.exists() and 0 < partial.stat().st_size < len(payload)
    second = subprocess.Popen(["bash", str(PULL)], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
    try:
        wait_until(destination.exists, 15)
        assert hashlib.sha256(destination.read_bytes()).hexdigest() == expected
        assert hashlib.sha256(original.read_bytes()).hexdigest() == expected
        assert keep.read_bytes() == b"never delete local-only data"
        assert not partial.exists()
        wait_until(lambda: "full_channel" in (tmp_path / "state/rsync_status.tsv").read_text(), 3)
    finally:
        stop(second)


def test_mirror_keeps_partial_files_hidden_and_does_not_copy_source_partials(tmp_path):
    env, _ = pull_config(tmp_path,
        f'LOCAL_MIRROR_METADATA_ROOT="{tmp_path}/mirror"\n'
        'LOCAL_MIRROR_BWLIMIT_KB=128\nLOCAL_MIRROR_SLEEP_SECONDS=0.1\n')
    source = tmp_path / "data/cut"
    source.mkdir(parents=True)
    original = source / "example.h5"
    original.write_bytes(random.Random(13).randbytes(1024 * 1024))
    incomplete = source / ".rsync-partial/incomplete.h5"
    incomplete.parent.mkdir()
    incomplete.write_bytes(b"must not be mirrored")
    destination = tmp_path / "mirror/cut/example.h5"
    first = subprocess.Popen(["bash", str(MIRROR)], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
    try:
        wait_until(lambda: (tmp_path / "state/local_mirror_active.tsv").exists(), 3)
        time.sleep(1.5)
        assert not destination.exists()
    finally:
        stop(first)
    assert (destination.parent / ".rsync-partial/example.h5").exists()
    second = subprocess.Popen(["bash", str(MIRROR)], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
    try:
        wait_until(destination.exists, 15)
        assert destination.read_bytes() == original.read_bytes()
        assert not (tmp_path / "mirror/cut/.rsync-partial/incomplete.h5").exists()
    finally:
        stop(second)


def test_mirror_rejects_nested_destination(tmp_path):
    env, _ = pull_config(tmp_path, f'LOCAL_MIRROR_METADATA_ROOT="{tmp_path}/data/recursive-mirror"\n')
    result = subprocess.run(["bash", str(MIRROR)], env=env, capture_output=True, text=True, timeout=3)
    assert result.returncode != 0
    assert "overlap" in result.stderr


@pytest.mark.parametrize("channels", ["cut cut", "../escape"])
def test_invalid_channel_configuration_is_rejected(tmp_path, channels):
    env, _ = pull_config(tmp_path, f'METADATA_CHANNELS="{channels}"\n')
    result = subprocess.run(["bash", str(PULL)], env=env, capture_output=True, text=True, timeout=3)
    assert result.returncode != 0
    assert "Invalid or duplicate" in result.stderr
