import numpy as np

from mesomode_boundary import closed_mode_blocks


def test_missed_cycle_does_not_split_mesosphere_mode():
    cycle = 32_000
    first = np.arange(0, 86_000_000, cycle, dtype=np.int64)
    first = np.delete(first, 1000)  # one missed detection makes a 64 ms gap
    second_start = int(first[-1] + 130_000_000)
    second = second_start + np.arange(0, 86_000_000, cycle, dtype=np.int64)

    blocks = closed_mode_blocks(
        np.concatenate((first, second)),
        process_end=int(second[-1] + 1_000_000),
        max_gap=500_000,
    )

    assert blocks == [(int(first[0]), int(first[-1])), (int(second[0]), int(second[-1]))]


def test_trailing_open_block_is_not_written():
    starts = np.arange(0, 10 * 32_000, 32_000, dtype=np.int64)

    blocks = closed_mode_blocks(
        starts,
        process_end=int(starts[-1] + 100_000),
        max_gap=500_000,
    )

    assert blocks == []
