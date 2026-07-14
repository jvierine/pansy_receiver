import sys
from pathlib import Path

import h5py
import numpy as np


PANSY_RECEIVER_ROOT = Path(__file__).resolve().parents[1]
if str(PANSY_RECEIVER_ROOT) not in sys.path:
    sys.path.insert(0, str(PANSY_RECEIVER_ROOT))


def test_orbit_path_rows_preserve_selection_keep():
    import orbit_metadata_table as omt

    _, _, paths = omt.payload_to_rows(
        123,
        {
            "path_t_rel_s": [0.0, 1.0],
            "path_position_enu_km": [[0.0, 0.0, 100.0], [0.0, 0.0, 40.0]],
            "path_doppler_mps": [1000.0, 1000.0],
            "path_snr": [10.0, 5.0],
            "path_beam_id": [0, 0],
            "path_selection_keep": [True, False],
        },
    )

    np.testing.assert_array_equal(paths["selection_keep"], [True, False])


def test_orbit_path_rows_reject_missing_selection_mask():
    import orbit_metadata_table as omt

    _, _, paths = omt.payload_to_rows(
        123,
        {
            "path_t_rel_s": [0.0],
            "path_position_enu_km": [[0.0, 0.0, 100.0]],
            "path_doppler_mps": [1000.0],
            "path_snr": [10.0],
            "path_beam_id": [0],
        },
    )

    np.testing.assert_array_equal(paths["selection_keep"], [False])


def test_backfill_reads_winning_robust_fit_mask(tmp_path):
    from backfill_orbit_path_selection_keep import load_selection_keep

    path = tmp_path / "state.h5"
    with h5py.File(path, "w") as handle:
        losing = handle.create_group("H01")
        losing.attrs["combined_rank"] = 1
        losing.create_dataset("path_selection_keep", data=[True, True])
        winning = handle.create_group("H02")
        winning.attrs["combined_rank"] = 0
        winning.create_dataset("path_selection_keep", data=[True, False])

    np.testing.assert_array_equal(load_selection_keep(path, b"H02", 2), [True, False])
