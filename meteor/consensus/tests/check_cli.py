"""Exercise native MF metadata reading, boundary pulses, TX gate and HDF5 output."""
from pathlib import Path
import subprocess
import tempfile

import digital_rf
import h5py
import numpy as np

crate = Path(__file__).resolve().parents[1]
subprocess.run(["cargo", "build", "--release", "--manifest-path", str(crate/"Cargo.toml")], check=True)
with tempfile.TemporaryDirectory(prefix="pansy-consensus-") as tmp:
    root = Path(tmp)
    metadata = root/"mf"
    metadata.mkdir()
    writer = digital_rf.DigitalMetadataWriter(str(metadata), 3600, 60, 1_000_000, 1, "mf")
    epoch = 1_700_000_000_000_000
    for batch in range(4):
        t = epoch + batch*32000 + np.arange(20, dtype=np.int64)*1600
        power = np.full(20, 2e9)
        if batch == 3:
            power[0] = 0
        writer.write([int(t[0])], {
            "tx_idxs": [t], "max_range": [(100+(t-epoch)/1e6*20).astype(np.float32)],
            "max_dopvel": [np.full(20, 20000, dtype=np.float32)],
            "max_snr": [np.full(20, 10, dtype=np.float32)], "tx_pwr": [power],
        })
    del writer
    output = root/"result.h5"
    subprocess.run([str(crate/"target/release/pansy-meteor-consensus"), str(metadata),
                    str(epoch+1600), str(epoch+128000), "--output", str(output)], check=True)
    with h5py.File(output) as f:
        assert len(f["sample_idx"]) == 59  # first batch tails included; bad TX batch rejected
        assert f["sample_idx"][0] == epoch+1600
        assert np.all(f["track_id"][:] == 0)
        assert f["sample_idx"].dtype == np.dtype("int64")
        assert f["track_id"].dtype == np.dtype("int32")
        for name in ("range_km", "radial_velocity_m_s", "snr_db", "track_range_at_epoch_km"):
            assert f[name].dtype == np.dtype("float32"), name
    print("PASS: boundary pulse retention, TX gate, membership and HDF5 dtypes")
