from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import pytest
import cut_meteors as cut
from meteor_progress import read_progress, write_progress


@pytest.mark.parametrize("fail", [False, True, "empty"])
def test_short_final_window_and_failed_write_checkpoint(tmp_path, monkeypatch, fail):
    monkeypatch.setenv("HOME",str(tmp_path))
    state=tmp_path/".local/state/pansy_receiver"
    state.mkdir(parents=True)
    first=1_000_000
    end=2_000_000
    write_progress(state/"cluster_phase_v1.h5",end)
    samples=first+np.arange(6)*8000
    record=dict(tx_idx=samples,xhat=np.zeros(3),range=np.ones(6)*100,
                doppler=np.ones(6)*-20000,snr=np.ones(6)*20,beam=np.zeros(6))
    reads=[]
    class Reader:
        def read(self,lo,hi,*args):
            reads.append((lo,hi))
            return {first:record}
    class Writer:
        def write(self,*args):
            if fail: raise OSError("disk failure")
    def bounds(path,label):
        return (Reader(),(first,first)) if label=="detections" else (None,(-1,-1))
    with patch.object(cut.pc,"cut_metadata_dir",str(tmp_path/"cut")), \
         patch.object(cut,"metadata_bounds",side_effect=bounds), \
         patch.object(cut.drf,"DigitalRFReader",return_value=SimpleNamespace(get_bounds=lambda ch:(0,10_000_000))), \
         patch.object(cut.drf,"DigitalMetadataWriter",return_value=Writer()), \
         patch.object(cut.cmf,"fit_obs",return_value=(0,0,np.zeros(3),0,lambda t:100,lambda t:-20)), \
         patch.object(cut,"cut_raw_voltage",return_value={} if fail=="empty" else {"tx_idx":[samples]}):
        if fail=="empty":
            with pytest.raises(RuntimeError,match="no mesosphere TX batches"): cut.cut_block()
        elif fail:
            with pytest.raises(OSError,match="disk failure"): cut.cut_block()
        else:
            assert cut.cut_block()==1
    assert reads==[(first,end-1)]
    assert read_progress(state/"cut_complete_v1.h5")== (first if fail else end)
