"""Seven-parameter shower windows inferred from the Zenodo CSV selections.

Each interval is the full selected-event envelope in current Revontuli orbit
metadata, rounded outward to a human-readable boundary.  The intervals
therefore retain every event in the corresponding authoritative CSV while
providing a reproducible joint cut for activity-rate estimates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ShowerSelectionWindow:
    code: str
    csv_name: str
    csv_count: int
    vg_range_km_s: tuple[float, float]
    e_range: tuple[float, float]
    a_range_au: tuple[float, float]
    i_range_deg: tuple[float, float]
    q_range_au: tuple[float, float]
    ra_center_deg: float
    ra_half_width_deg: float
    dec_range_deg: tuple[float, float]


WINDOWS = {
    "JUE": ShowerSelectionWindow(
        code="JUE",
        csv_name="-1_JUE_pansy01.csv",
        csv_count=83,
        vg_range_km_s=(45.5, 52.5),
        e_range=(0.625, 0.830),
        a_range_au=(2.2, 4.8),
        i_range_deg=(84.0, 99.0),
        q_range_au=(0.79, 0.96),
        ra_center_deg=54.0,
        ra_half_width_deg=9.0,
        dec_range_deg=(-34.5, -25.5),
    ),
    "CAP": ShowerSelectionWindow(
        code="CAP",
        csv_name="001_CAP_pansy01.csv",
        csv_count=283,
        vg_range_km_s=(19.5, 27.0),
        e_range=(0.690, 0.830),
        a_range_au=(1.6, 3.0),
        i_range_deg=(4.0, 12.0),
        q_range_au=(0.40, 0.69),
        ra_center_deg=301.0,
        ra_half_width_deg=17.0,
        dec_range_deg=(-15.5, -5.0),
    ),
    "DCS": ShowerSelectionWindow(
        code="DCS",
        csv_name="115_DCS_pansy01.csv",
        csv_count=249,
        vg_range_km_s=(22.0, 28.5),
        e_range=(0.725, 0.845),
        a_range_au=(1.7, 3.0),
        i_range_deg=(4.0, 10.5),
        q_range_au=(0.40, 0.56),
        ra_center_deg=312.0,
        ra_half_width_deg=11.0,
        dec_range_deg=(-30.5, -22.0),
    ),
    "OES": ShowerSelectionWindow(
        code="OES",
        csv_name="1208_OES_pansy01.csv",
        csv_count=186,
        vg_range_km_s=(18.5, 25.5),
        e_range=(0.660, 0.825),
        a_range_au=(1.6, 3.6),
        i_range_deg=(5.0, 10.5),
        q_range_au=(0.50, 0.68),
        ra_center_deg=301.0,
        ra_half_width_deg=14.0,
        dec_range_deg=(-36.5, -27.5),
    ),
}


def wrap180(values: np.ndarray) -> np.ndarray:
    return (np.asarray(values, dtype=np.float64) + 180.0) % 360.0 - 180.0


def selection_mask(
    window: ShowerSelectionWindow,
    *,
    vg_km_s: np.ndarray,
    kepler: np.ndarray,
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
) -> np.ndarray:
    """Return the joint seven-parameter selection mask."""
    vg = np.asarray(vg_km_s, dtype=np.float64)
    kep = np.asarray(kepler, dtype=np.float64)
    ra = np.asarray(ra_deg, dtype=np.float64)
    dec = np.asarray(dec_deg, dtype=np.float64)
    a, e, inc, q = kep[:, 0], kep[:, 1], kep[:, 2], kep[:, 6]
    finite = (
        np.isfinite(vg)
        & np.isfinite(a)
        & np.isfinite(e)
        & np.isfinite(inc)
        & np.isfinite(q)
        & np.isfinite(ra)
        & np.isfinite(dec)
    )
    keep = finite.copy()
    for values, bounds in (
        (vg, window.vg_range_km_s),
        (e, window.e_range),
        (a, window.a_range_au),
        (inc, window.i_range_deg),
        (q, window.q_range_au),
        (dec, window.dec_range_deg),
    ):
        keep &= (values >= bounds[0]) & (values <= bounds[1])
    keep &= np.abs(wrap180(ra - window.ra_center_deg)) <= window.ra_half_width_deg
    return keep


def format_window(window: ShowerSelectionWindow) -> str:
    return (
        f"{window.code}: vg={window.vg_range_km_s}, e={window.e_range}, "
        f"a={window.a_range_au}, i={window.i_range_deg}, q={window.q_range_au}, "
        f"RA={window.ra_center_deg} +/- {window.ra_half_width_deg}, "
        f"Dec={window.dec_range_deg}"
    )
