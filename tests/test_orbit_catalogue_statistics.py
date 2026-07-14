import sys
from pathlib import Path

import numpy as np


PANSY_RECEIVER_ROOT = Path(__file__).resolve().parents[1]
if str(PANSY_RECEIVER_ROOT) not in sys.path:
    sys.path.insert(0, str(PANSY_RECEIVER_ROOT))


def test_solar_longitude_count_density_uses_bin_width():
    from plot_orbit_catalogue_statistics import count_density_from_counts

    counts_by_year = np.asarray([[10, 20], [4, 8]], dtype=np.int32)
    edges_deg = np.asarray([0.0, 2.0, 6.0], dtype=np.float32)

    density_by_year, all_density = count_density_from_counts(counts_by_year, edges_deg)

    np.testing.assert_allclose(density_by_year, [[5.0, 5.0], [2.0, 2.0]])
    np.testing.assert_allclose(all_density, [7.0, 7.0])

