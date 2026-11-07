from types import SimpleNamespace
from unittest.mock import patch
import pytest
import process_cut_meteor as fit
from meteor_progress import read_progress, write_progress


@pytest.mark.parametrize("failure", [None, "read", "fit"])
def test_rejected_fit_advances_but_io_failures_do_not(tmp_path, monkeypatch, failure):
    monkeypatch.setenv("HOME",str(tmp_path))
    state=tmp_path/".local/state/pansy_receiver"
    state.mkdir(parents=True)
    write_progress(state/"cut_complete_v1.h5",2_000_000)
    class Reader:
        def read(self,lo,hi):
            assert (lo,hi)==(1_000_000,1_999_999)
            if failure=="read": raise OSError("read failed")
            return {1_000_000:{}}
    def bounds(path,label):
        return (Reader(),(1_000_000,1_000_000)) if label=="cut" else (None,(-1,-1))
    comm=SimpleNamespace(bcast=lambda value,root:value,allreduce=lambda value,op:value)
    with patch.object(fit,"comm",comm), patch.object(fit,"rank",0), patch.object(fit,"size",1), \
         patch.object(fit,"metadata_bounds",side_effect=bounds), \
         patch.object(fit.pc,"simple_fit_metadata_dir",str(tmp_path/"fit")), \
         patch.object(fit,"simple_fit_writer",return_value=object()), \
         patch.object(fit,"process_cut",side_effect=OSError("fit failed") if failure=="fit" else None,return_value=False):
        assert fit.process_latest()==(1 if failure is None else 0)
    assert read_progress(state/"fit_complete_v1.h5")== (2_000_000 if failure is None else 1_000_000)
