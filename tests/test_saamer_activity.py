import numpy as np

from saamer_activity import sliding_counts
from shower_selection_windows import WINDOWS


def test_sliding_counts_uses_requested_window_width():
    solar = np.asarray([109.6, 110.0, 110.4, 111.1])
    centers = np.asarray([110.0, 111.0])

    counts = sliding_counts(solar, centers, window_deg=1.0)

    np.testing.assert_array_equal(counts, [3.0, 1.0])


def test_common_oes_dcs_window_contains_both_individual_windows():
    common = WINDOWS["OES_DCS"]
    for code in ("OES", "DCS"):
        window = WINDOWS[code]
        assert common.vg_range_km_s[0] <= window.vg_range_km_s[0]
        assert common.vg_range_km_s[1] >= window.vg_range_km_s[1]
        assert common.e_range[0] <= window.e_range[0]
        assert common.e_range[1] >= window.e_range[1]
        assert common.a_range_au[0] <= window.a_range_au[0]
        assert common.a_range_au[1] >= window.a_range_au[1]
        assert common.i_range_deg[0] <= window.i_range_deg[0]
        assert common.i_range_deg[1] >= window.i_range_deg[1]
        assert common.q_range_au[0] <= window.q_range_au[0]
        assert common.q_range_au[1] >= window.q_range_au[1]
        assert abs(common.ra_center_deg - window.ra_center_deg) + window.ra_half_width_deg <= common.ra_half_width_deg
        assert common.dec_range_deg[0] <= window.dec_range_deg[0]
        assert common.dec_range_deg[1] >= window.dec_range_deg[1]
