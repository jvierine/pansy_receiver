"""Lightweight sparse-archive monitoring, without reading complex spectra.

Missing cells are not zero power, but use the zero-level display color. Only
retained zenith detections are drawn; complete noise-only records still advance
the processing frontier.
"""
from datetime import datetime, timezone
from pathlib import Path
import warnings

import h5py
import numpy as np


MONITOR_FIELDS = (
    "complete", "mode_stop_sample_exclusive", "sparse_index",
    "detector_snr_db", "nominal_altitude_km",
)


def read_sparse_status(root, duration_us=24 * 3600 * 1_000_000):
    # Read HDF5 explicitly in mode 'r'. DigitalMetadataReader's recovery path
    # can delete an unreadable old file; a monitor must never repair its input.
    files = []
    for path in Path(root).glob("*/pmse@*.h5"):
        try:
            files.append((int(path.stem.split("@", 1)[1]), path))
        except ValueError:
            continue
    if not files:
        return [-1, -1], []
    files.sort()
    cutoff = files[-1][0] - duration_us / 1e6
    records = []
    bounds = [-1, -1]
    for file_second, path in files:
        if file_second + 600 < cutoff:
            continue
        try:
            with h5py.File(path, "r") as handle:
                for key in sorted(handle, key=int):
                    group = handle[key]
                    if "complete" not in group or int(group["complete"][()]) != 1:
                        continue
                    row = {name: group[name][()] for name in MONITOR_FIELDS}
                    start, stop = int(key), int(row["mode_stop_sample_exclusive"])
                    if stop <= start:
                        continue
                    # Freshness is independent of selection and plot availability.
                    bounds[0] = start if bounds[0] == -1 else min(bounds[0], start)
                    bounds[1] = max(bounds[1], stop)
                    records.append((start, stop, row))
        except (OSError, KeyError, ValueError) as exc:
            warnings.warn("Sparse monitor skipped unreadable metadata %s: %s" % (path.name, exc))
    records.sort(key=lambda record: record[0])
    return bounds, records


def zenith_cells(row):
    """Altitude edges and Doppler-maximum SNR for retained zenith cells only."""
    index = np.asarray(row["sparse_index"], dtype=np.int64).reshape(-1, 3)
    snr = np.asarray(row["detector_snr_db"], dtype=float).reshape(-1)
    altitude = np.asarray(row["nominal_altitude_km"], dtype=float)[0]
    if len(index) != len(snr) or altitude.ndim != 1 or len(altitude) < 2:
        raise ValueError("Invalid sparse monitor axes")
    if not np.isfinite(altitude).all() or not np.all(np.diff(altitude) > 0):
        raise ValueError("Invalid sparse altitude axis")
    edges = np.r_[altitude[0] - (altitude[1] - altitude[0]) / 2,
                  (altitude[:-1] + altitude[1:]) / 2,
                  altitude[-1] + (altitude[-1] - altitude[-2]) / 2]
    selected = (index[:, 0] == 0) & np.isfinite(snr)
    ranges = index[selected, 2]
    if np.any((ranges < 0) | (ranges >= len(altitude))):
        raise ValueError("Sparse range index outside altitude axis")
    peak = np.full(len(altitude), -np.inf)
    np.maximum.at(peak, ranges, snr[selected])
    retained = np.flatnonzero(np.isfinite(peak))
    return edges[retained], edges[retained + 1], peak[retained]


def plot_sparse_status(fig, ax, records, now_seconds):
    from matplotlib import colormaps
    import matplotlib.dates as dates
    from matplotlib.collections import PolyCollection

    # Sparse gaps have no polygons. Match their background and explicit NaNs to
    # the zero-level color without inserting fabricated zero-valued samples.
    cmap = colormaps["viridis"].copy()
    cmap.set_bad(cmap(0.0))
    ax.set_facecolor(cmap(0.0))

    def date_number(sample):
        return dates.date2num(datetime.fromtimestamp(sample / 1e6, timezone.utc))

    start = (now_seconds - 24 * 3600) * 1e6
    end = now_seconds * 1e6
    polygons, colors = [], []
    for left, right, row in records:
        if right <= start or left >= end:
            continue
        low, high, peak = zenith_cells(row)
        x0, x1 = date_number(max(left, start)), date_number(min(right, end))
        polygons.extend([(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
                        for y0, y1 in zip(low, high))
        colors.extend(peak)
    if polygons:
        collection = PolyCollection(polygons, array=np.asarray(colors), cmap=cmap,
                                    edgecolors="none")
        collection.set_clim(0, max(10, float(np.max(colors))))
        ax.add_collection(collection)
        ax.autoscale_view(scalex=False)
        fig.colorbar(collection, ax=ax, label="Peak retained SNR (dB)")
    else:
        ax.text(0.5, 0.5, "No retained zenith detections in the last 24 hours",
                ha="center", va="center", color="white", transform=ax.transAxes)
    ax.set_xlim(date_number(start), date_number(end))
    ax.xaxis_date(tz=timezone.utc)
    ax.set_xlabel("Date (UTC)")
    ax.set_ylabel("Nominal altitude (km)")
    ax.set_title("Compressed spectra: selected zenith echoes (not dense power)")
    fig.autofmt_xdate()
