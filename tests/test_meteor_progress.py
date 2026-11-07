import glob

import h5py

from meteor_progress import read_completed_frontier, read_progress, write_progress
from meteor_progress import inclusive_chunks
from meteor_progress import tx_scan_frontier
from meteor_progress import search_checkpoint


def test_chunks_include_tail_once_and_allow_empty_input():
    assert inclusive_chunks(5, 5, 60) == [(5, 5)]
    assert inclusive_chunks(5, 4, 60) == []
    chunks=inclusive_chunks(5, 130, 60)
    assert chunks == [(5,64),(65,124),(125,130)]
    assert [v for a,b in chunks for v in range(a,b+1)] == list(range(5,131))


def test_mpi_chunks_cover_every_sample_once():
    chunks=inclusive_chunks(1_000_000,4_032_000,1_000_000)
    assigned=[chunks[r::6] for r in range(6)]
    assert sum(map(len,assigned)) == len(chunks)
    assert chunks[-1][1] == 4_032_000


def test_scanned_empty_interval_advances_mode_closure(tmp_path):
    progress=tmp_path/"scan.json"
    assert tx_scan_frontier(100,progress)==100
    progress.write_text('{"next_sample": 500}')
    assert tx_scan_frontier(100,progress)==500


def test_search_prefers_reboot_persistent_checkpoint(tmp_path,monkeypatch):
    monkeypatch.setenv("HOME",str(tmp_path))
    assert search_checkpoint()=="/tmp/meteor_mf_complete.h5"
    durable=tmp_path/".local/state/pansy_receiver/mf_complete_v1.h5"
    durable.parent.mkdir(parents=True)
    write_progress(durable,123)
    assert search_checkpoint()==str(durable)


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
