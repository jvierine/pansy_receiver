"""Atomic progress checkpoints shared by the realtime meteor pipeline."""

import glob
import os
import json

import h5py


def inclusive_chunks(start, end, width):
    """Disjoint bounded intervals for DigitalMetadata's inclusive reads."""
    if width <= 0:
        raise ValueError("chunk width must be positive")
    return [(i, min(i + width - 1, end)) for i in range(int(start), int(end) + 1, int(width))]


def tx_scan_frontier(latest_detection, path=None):
    """Use completed scanner input, including intervals with no TX detections."""
    if path is None:
        path=os.environ.get("PANSY_FIND_MODE_PROGRESS_PATH",
            os.path.expanduser("~/.local/state/pansy-receiver/find_mode_starts.json"))
    if not os.path.exists(path):
        return int(latest_detection)
    with open(path) as stream:
        frontier=json.load(stream)["next_sample"]
    if not isinstance(frontier,int) or frontier < 0:
        raise ValueError("invalid TX scanner frontier")
    return frontier


def write_progress(path, latest):
    """Atomically replace an HDF5 progress checkpoint."""
    temporary_path = "%s.%d.tmp" % (path, os.getpid())
    try:
        with h5py.File(temporary_path, "w") as progress_file:
            progress_file["latest"] = int(latest)
        os.replace(temporary_path, path)
    finally:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)


def read_progress(path):
    """Read a progress checkpoint and return its integer sample index."""
    with h5py.File(path, "r") as progress_file:
        return int(progress_file["latest"][()])


def read_completed_frontier(complete_path, rank_pattern):
    """Return the committed frontier, with legacy rank files as a fallback."""
    if os.path.exists(complete_path):
        return read_progress(complete_path)

    rank_paths = glob.glob(rank_pattern)
    if not rank_paths:
        return None
    return min(read_progress(path) for path in rank_paths)
