import importlib.util
from pathlib import Path
import pytest

spec=importlib.util.spec_from_file_location('configure_direct_rsync',Path(__file__).resolve().parents[1]/'backup/scripts/configure_direct_rsync.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_migration_preserves_paths_backs_up_config_and_sets_shared_cap(tmp_path):
    config=tmp_path/'pansy.env'
    original='JUMP_HOST=j@example.org\nPANSY_SSH_HOST=radar@localhost\nPANSY_SSH_PORT=3131\nRSYNC_BWLIMIT_KB=15\nRSYNC_TOTAL_BWLIMIT_KB=50\nRSYNC_RECENT_MINUTES=20\nRSYNC_HISTORY_SLICE_SECONDS=180\nLOCAL_METADATA_ROOT=/keep/my/data\n'
    config.write_text(original)
    module.configure(config,tmp_path/'units')
    text=config.read_text()
    assert 'RSYNC_BWLIMIT_KB=128\n' in text
    assert 'RSYNC_TOTAL_BWLIMIT_KB=128\n' in text
    assert 'RSYNC_DIRECT_SSH=1\n' in text
    assert 'RSYNC_RECENT_MINUTES=' not in text
    assert 'RSYNC_HISTORY_SLICE_SECONDS=' not in text
    assert 'LOCAL_METADATA_ROOT=/keep/my/data\n' in text
    assert next(tmp_path.glob('*.bak')).read_text()==original
    dropin=(tmp_path/'units/pansy-backup-rsync.service.d/direct-outbound.conf').read_text()
    assert 'Wants=\nWants=network-online.target' in dropin
    module.configure(config,tmp_path/'units')
    assert config.read_text().count('RSYNC_BWLIMIT_KB=')==1
    assert config.read_text()==text


def test_explicit_bandwidth_override(tmp_path):
    config=tmp_path/'pansy.env'
    config.write_text('JUMP_HOST=j@example.org\nPANSY_SSH_HOST=radar@localhost\nPANSY_SSH_PORT=3131\n')
    module.configure(config,tmp_path/'units',256)
    assert 'RSYNC_TOTAL_BWLIMIT_KB=256\n' in config.read_text()


def test_mount_guards_and_checkout_are_explicit_and_persistent(tmp_path):
    config=tmp_path/'pansy.env'
    config.write_text('JUMP_HOST=j@example.org\nPANSY_SSH_HOST=radar@localhost\nPANSY_SSH_PORT=3131\n')
    repo=Path(__file__).resolve().parents[1]
    module.configure(config,tmp_path/'units',128,'/mnt/data','/urdr/data',repo)
    assert 'LOCAL_METADATA_MOUNT=/mnt/data\n' in config.read_text()
    assert 'LOCAL_MIRROR_METADATA_MOUNT=/urdr/data\n' in config.read_text()
    for unit in ('pansy-backup-rsync','pansy-local-metadata-mirror','pansy-backup-web'):
        dropin=(tmp_path/'units'/(unit+'.service.d')/'transfer-checkout.conf').read_text()
        assert str(repo) in dropin
        assert 'ExecStart=\nExecStart=' in dropin


def test_missing_checkout_rejected_without_changing_configuration(tmp_path):
    config=tmp_path/'pansy.env'
    original='JUMP_HOST=j@example.org\nPANSY_SSH_HOST=radar@localhost\nPANSY_SSH_PORT=3131\n'
    config.write_text(original)
    with pytest.raises(ValueError,match='Missing service script'):
        module.configure(config,tmp_path/'units',service_repo=tmp_path/'missing-checkout')
    assert config.read_text()==original
    assert not list(tmp_path.glob('*.bak'))
