import importlib.util
from pathlib import Path

spec=importlib.util.spec_from_file_location('configure_direct_rsync',Path(__file__).resolve().parents[1]/'backup/scripts/configure_direct_rsync.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_migration_preserves_paths_backs_up_config_and_sets_shared_cap(tmp_path):
    config=tmp_path/'pansy.env'
    original='JUMP_HOST=j@example.org\nPANSY_SSH_HOST=radar@localhost\nPANSY_SSH_PORT=3131\nRSYNC_BWLIMIT_KB=15\nRSYNC_TOTAL_BWLIMIT_KB=50\nLOCAL_METADATA_ROOT=/keep/my/data\n'
    config.write_text(original)
    module.configure(config,tmp_path/'units')
    text=config.read_text()
    assert 'RSYNC_BWLIMIT_KB=128\n' in text
    assert 'RSYNC_TOTAL_BWLIMIT_KB=128\n' in text
    assert 'RSYNC_DIRECT_SSH=1\n' in text
    assert 'LOCAL_METADATA_ROOT=/keep/my/data\n' in text
    assert next(tmp_path.glob('*.bak')).read_text()==original
    dropin=(tmp_path/'units/pansy-backup-rsync.service.d/direct-outbound.conf').read_text()
    assert 'Wants=\nWants=network-online.target' in dropin
    module.configure(config,tmp_path/'units')
    assert config.read_text().count('RSYNC_BWLIMIT_KB=')==1


def test_explicit_bandwidth_override(tmp_path):
    config=tmp_path/'pansy.env'
    config.write_text('JUMP_HOST=j@example.org\nPANSY_SSH_HOST=radar@localhost\nPANSY_SSH_PORT=3131\n')
    module.configure(config,tmp_path/'units',256)
    assert 'RSYNC_TOTAL_BWLIMIT_KB=256\n' in config.read_text()
