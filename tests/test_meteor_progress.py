import glob

import h5py

from meteor_progress import read_completed_frontier, read_progress, write_progress


def test_progress_checkpoint_is_replaced_atomically(tmp_path):
    checkpoint = tmp_path / "meteor_mf_complete.h5"

    write_progress(checkpoint, 123)
    assert read_progress(checkpoint) == 123

    write_progress(checkpoint, 456)
    assert read_progress(checkpoint) == 456
    assert glob.glob(str(checkpoint) + ".*.tmp") == []


def test_progress_checkpoint_has_expected_hdf5_schema(tmp_path):
    checkpoint = tmp_path / "meteor_mf_complete.h5"
    write_progress(checkpoint, 789)

    with h5py.File(checkpoint, "r") as progress_file:
        assert list(progress_file.keys()) == ["latest"]
        assert int(progress_file["latest"][()]) == 789


def test_committed_frontier_overrides_stale_idle_rank(tmp_path):
    rank_pattern = str(tmp_path / "meteor_mf_[0-9]*.h5")
    complete = tmp_path / "meteor_mf_complete.h5"
    write_progress(tmp_path / "meteor_mf_0.h5", 500)
    write_progress(tmp_path / "meteor_mf_5.h5", 100)

    assert read_completed_frontier(complete, rank_pattern) == 100

    write_progress(complete, 500)
    assert read_completed_frontier(complete, rank_pattern) == 500
