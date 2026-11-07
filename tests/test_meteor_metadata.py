import os
from unittest.mock import patch
import digital_rf
import h5py
import numpy as np
import pytest
from meteor_metadata import StrictDigitalMetadataReader


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
