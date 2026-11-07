import os
from unittest.mock import patch
import digital_rf
import h5py
import numpy as np
import pytest
from meteor_metadata import StrictDigitalMetadataReader
from meteor_metadata import matched_filter_record_exists, MF_FIELDS
from meteor_metadata import write_complete_record


def make_reader(tmp_path):
    writer=digital_rf.DigitalMetadataWriter(str(tmp_path),3600,60,1,1,"test")
    writer.write([100,101],{"value":np.array([3,4],dtype=np.int32)})
    return StrictDigitalMetadataReader(str(tmp_path))


def test_strict_reader_matches_digital_metadata_columns(tmp_path):
    reader=make_reader(tmp_path)
    assert list(reader.read(100,101,"value").values())==[3,4]
    assert reader.read(101,101,["value"])[101]["value"]==4
    assert reader.read(100,100)[100]["value"]==3
    assert reader.read(102,102,"value",method="pad")[101]==4


def test_failed_open_raises_and_never_deletes_existing_file(tmp_path):
    reader=make_reader(tmp_path)
    path=next(tmp_path.glob("*/test@*.h5"))
    os.utime(path,(1,1))
    with patch("meteor_metadata.h5py.File",side_effect=OSError("read failure")):
        with pytest.raises(OSError,match="read failure"):
            reader.read(100,101,"value")
    assert path.exists()


def test_partial_record_raises_instead_of_returning_incomplete_rows(tmp_path):
    reader=make_reader(tmp_path)
    path=next(tmp_path.glob("*/test@*.h5"))
    with h5py.File(path,"a") as h:
        del h["101/value"]
    with pytest.raises(KeyError): reader.read(100,101,["value"])


def test_interrupted_mf_row_replays_without_removing_valid_neighbour(tmp_path):
    writer=digital_rf.DigitalMetadataWriter(str(tmp_path),3600,60,1,1,"test")
    def payload(key):
        data={name:np.ones(20,dtype=np.float32) for name in MF_FIELDS}
        data["tx_idxs"]=key+np.arange(20)*1600
        return data
    writer.write(100,payload(100)); writer.write(101,payload(101))
    reader=StrictDigitalMetadataReader(str(tmp_path))
    assert matched_filter_record_exists(reader,100)
    path=next(tmp_path.glob("*/test@*.h5"))
    with h5py.File(path,"a") as h: del h["100/max_snr"]
    assert not matched_filter_record_exists(reader,100)
    assert matched_filter_record_exists(reader,101)
    writer.write(100,payload(100))
    assert matched_filter_record_exists(reader,100)


def test_duplicate_and_partial_output_record_recovery(tmp_path):
    writer=digital_rf.DigitalMetadataWriter(str(tmp_path),3600,60,1,1,"test")
    payload={"a":np.arange(3),"b":np.arange(4)}
    assert write_complete_record(writer,tmp_path,100,payload)
    assert write_complete_record(writer,tmp_path,101,payload)
    assert not write_complete_record(writer,tmp_path,100,payload)
    path=next(tmp_path.glob("*/test@*.h5"))
    with h5py.File(path,"a") as h: del h["100/b"]
    assert write_complete_record(writer,tmp_path,100,payload)
    rows=StrictDigitalMetadataReader(str(tmp_path)).read(100,101)
    assert list(rows)==[100,101]
    np.testing.assert_array_equal(rows[100]["b"],payload["b"])
