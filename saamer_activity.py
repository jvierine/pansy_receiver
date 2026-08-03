"""Reusable SAAMER shower-selection cache for paper activity profiles."""

from __future__ import annotations

import zipfile
from pathlib import Path

import h5py
import numpy as np

from shower_selection_windows import WINDOWS, selection_mask


DEFAULT_ARCHIVES = sorted(Path.home().joinpath("Desktop").glob("iaumdcSAAMER20*.zip"))
DEFAULT_CACHE = Path(__file__).resolve().parent / "figs" / "saamer_shower_activity.h5"
SHOWER_CODES = ("JUE", "CAP", "DCS", "OES", "OES_DCS")


def _source_signature(paths: list[Path]) -> str:
    return "|".join(
        f"{path.resolve()}:{path.stat().st_size}:{path.stat().st_mtime_ns}" for path in paths
    )


def build_cache(paths: list[Path], cache_path: Path = DEFAULT_CACHE) -> None:
    """Read full-precision SAAMER archives once and retain selected shower rows."""
    selected: dict[str, dict[str, list[np.ndarray]]] = {
        code: {name: [] for name in ("solar", "ra", "dec", "vg", "a", "e", "i", "q")}
        for code in SHOWER_CODES
    }
    for path in paths:
        with zipfile.ZipFile(path) as archive:
            for member in sorted(archive.namelist()):
                if not member.lower().endswith(".dat"):
                    continue
                with archive.open(member) as stream:
                    values = np.loadtxt(stream, usecols=(4, 6, 7, 8, 10, 11, 12, 13))
                if values.ndim == 1:
                    values = values[None, :]
                solar, ra, dec, vg, q, eccentricity, semimajor_axis, inclination = values.T
                kepler = np.full((solar.size, 7), np.nan, dtype=np.float64)
                kepler[:, 0] = semimajor_axis
                kepler[:, 1] = eccentricity
                kepler[:, 2] = inclination
                kepler[:, 6] = q
                columns = (solar, ra, dec, vg, semimajor_axis, eccentricity, inclination, q)
                for code in SHOWER_CODES:
                    keep = selection_mask(
                        WINDOWS[code],
                        vg_km_s=vg,
                        kepler=kepler,
                        ra_deg=ra,
                        dec_deg=dec,
                    )
                    for name, column in zip(selected[code], columns, strict=True):
                        selected[code][name].append(column[keep])
        print(f"read SAAMER archive {path.name}", flush=True)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = cache_path.with_suffix(cache_path.suffix + ".tmp")
    with h5py.File(temporary, "w") as h5:
        h5.attrs["source_signature"] = _source_signature(paths)
        h5.attrs["selection"] = "same seven-parameter vg,e,a,i,q,RA,Dec windows used for PANSY"
        for code in SHOWER_CODES:
            group = h5.create_group(code)
            for name, parts in selected[code].items():
                values = np.concatenate(parts) if parts else np.asarray([], dtype=np.float64)
                group.create_dataset(name, data=values, compression="gzip", shuffle=True)
    temporary.replace(cache_path)


def load_selected(
    code: str,
    paths: list[Path] | None = None,
    cache_path: Path = DEFAULT_CACHE,
) -> dict[str, np.ndarray]:
    """Return full-precision SAAMER rows passing a paper shower selection."""
    archives = DEFAULT_ARCHIVES if paths is None else paths
    if not archives:
        raise FileNotFoundError("no SAAMER IAU MDC archives found")
    signature = _source_signature(archives)
    rebuild = not cache_path.exists()
    if not rebuild:
        with h5py.File(cache_path, "r") as h5:
            rebuild = h5.attrs.get("source_signature", "") != signature or code not in h5
    if rebuild:
        build_cache(archives, cache_path)
    with h5py.File(cache_path, "r") as h5:
        return {name: np.asarray(values) for name, values in h5[code].items()}


def sliding_counts(solar_longitude_deg, centers_deg, window_deg: float = 1.0) -> np.ndarray:
    solar = np.asarray(solar_longitude_deg, dtype=np.float64)
    centers = np.asarray(centers_deg, dtype=np.float64)
    half_window = 0.5 * float(window_deg)
    difference = (solar[:, None] - centers[None, :] + 180.0) % 360.0 - 180.0
    return np.sum(np.abs(difference) <= half_window, axis=0).astype(np.float64)
