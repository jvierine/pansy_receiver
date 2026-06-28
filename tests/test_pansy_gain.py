import sys
from pathlib import Path

import numpy as np


PANSY_RECEIVER_ROOT = Path(__file__).resolve().parents[1]
if str(PANSY_RECEIVER_ROOT) not in sys.path:
    sys.path.insert(0, str(PANSY_RECEIVER_ROOT))


def test_tx_gain_is_module_gain_times_sparse_array_gain():
    import pansy_gain as pgain

    beam_vecs = pgain.tx_beam_unit_vectors()
    uvw = np.asarray(
        [
            beam_vecs[0],
            [0.02, -0.03, -np.sqrt(1.0 - 0.02**2 - 0.03**2)],
            beam_vecs[1],
        ],
        dtype=np.float64,
    )
    beam_id = np.asarray([0, 0, 1], dtype=np.int64)
    modules = pgain.tx_module_positions()
    module_centers = pgain.tx_module_center_positions(modules)
    rel_module = modules[0] - np.mean(modules[0], axis=0, keepdims=True)

    expected = np.empty(len(uvw), dtype=np.float64)
    for beam in np.unique(beam_id):
        idx = beam_id == beam
        module_gain = pgain.module_power_gain(uvw[idx], steer=beam_vecs[beam], module_positions=rel_module)
        sparse_gain = pgain.sparse_tx_array_power_gain(
            uvw[idx],
            int(beam),
            beam_vecs=beam_vecs,
            module_centers=module_centers,
        )
        expected[idx] = module_gain * sparse_gain

    np.testing.assert_allclose(pgain.tx_power_gain(uvw, beam_id, beam_vecs=beam_vecs, modules=modules), expected)
    np.testing.assert_allclose(pgain.tx_gain_db(uvw, beam_id, beam_vecs=beam_vecs), pgain.power_to_db(expected))


def test_rx_gain_uses_single_module_pattern():
    import pansy_gain as pgain

    zenith = np.asarray([[0.0, 0.0, -1.0]], dtype=np.float64)
    off_zenith = np.asarray([[0.03, 0.04, -np.sqrt(1.0 - 0.03**2 - 0.04**2)]], dtype=np.float64)

    np.testing.assert_allclose(pgain.rx_power_gain(zenith, channel=0), [1.0], rtol=1e-12, atol=1e-12)
    off_gain = pgain.rx_power_gain(off_zenith, channel=0)
    assert np.isfinite(off_gain[0])
    assert 0.0 < off_gain[0] <= 1.0


def test_two_way_gain_is_tx_gain_times_rx_module_gain():
    import pansy_gain as pgain

    uvw = np.asarray(
        [
            [0.0, 0.0, -1.0],
            [0.04, -0.02, -np.sqrt(1.0 - 0.04**2 - 0.02**2)],
        ],
        dtype=np.float64,
    )
    beam_id = np.asarray([0, 1], dtype=np.int64)
    beam_vecs = pgain.tx_beam_unit_vectors()
    expected = np.empty(len(uvw), dtype=np.float64)
    tx_gain = pgain.tx_power_gain(uvw, beam_id, beam_vecs=beam_vecs)
    for beam in np.unique(beam_id):
        idx = beam_id == beam
        expected[idx] = tx_gain[idx] * pgain.rx_power_gain(uvw[idx], steer=beam_vecs[beam])
    np.testing.assert_allclose(pgain.two_way_power_gain(uvw, beam_id), expected)
    np.testing.assert_allclose(pgain.two_way_gain_db(uvw, beam_id), pgain.power_to_db(expected))
