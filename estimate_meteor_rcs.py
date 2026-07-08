#!/usr/bin/env python3
"""Estimate meteor radar cross section from PANSY orbit metadata."""

from __future__ import annotations

import argparse
from collections import deque
from pathlib import Path

import h5py
import numpy as np
import scipy.constants as sc
import scipy.special as ss

import fit_cut_noise_pygdsm as noise_fit
import orbit_metadata_table as omt
import pansy_config as pc
import pansy_gain as pg

DEFAULT_TREC_K = 535.613
DEFAULT_TSYS_FRACTIONAL_ERROR = 0.057
DEFAULT_ELEMENT_PATTERN_BLEND = 0.1
DEFAULT_TX_POWER_W = 500e3
DEFAULT_YAGI_GAIN_DBI = 7.2
DEFAULT_TX_MODULES = 54
DEFAULT_YAGIS_PER_MODULE = 19
DEFAULT_CODE_BITS = 16
DEFAULT_BAUD_S = 8e-6
DEFAULT_SUMMED_RX_CHANNELS = 7


def db_to_linear(db: float | np.ndarray) -> float | np.ndarray:
    return 10.0 ** (np.asarray(db) / 10.0)


def linear_to_db(value: float | np.ndarray) -> float | np.ndarray:
    return 10.0 * np.log10(np.maximum(np.asarray(value), 1e-300))


def decode_attr(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def absolute_tx_peak_gain_dbi(element_gain_dbi: float, modules: int, yagis_per_module: int) -> float:
    return float(element_gain_dbi + 10.0 * np.log10(float(modules * yagis_per_module)))


def absolute_rx_module_gain_dbi(element_gain_dbi: float, yagis_per_module: int) -> float:
    return float(element_gain_dbi + 10.0 * np.log10(float(yagis_per_module)))


def gamma_power_median_to_mean(shape: float) -> float:
    """Median/mean factor for a gamma-distributed summed power estimate."""
    shape = float(shape)
    if shape <= 0.0:
        raise ValueError("gamma shape must be positive")
    return float(ss.gammaincinv(shape, 0.5) / shape)


def median_referenced_snr_to_mean_referenced_snr(snr: np.ndarray, median_to_mean: float) -> np.ndarray:
    """Convert ``(P - median_noise) / median_noise`` to ``(P - mean_noise) / mean_noise``."""
    return (np.asarray(snr, dtype=np.float64) + 1.0) * float(median_to_mean) - 1.0


class MainLobeMask:
    """Connected -3 dB transmit main-lobe masks in direction-cosine space."""

    def __init__(
        self,
        beam_vecs: np.ndarray,
        threshold_db: float,
        grid_half_width_dc: float,
        grid_n: int,
    ) -> None:
        self.beam_vecs = np.asarray(beam_vecs, dtype=np.float64)
        self.threshold_linear = float(db_to_linear(threshold_db))
        self.grid_half_width_dc = float(grid_half_width_dc)
        self.grid_n = int(grid_n)
        if self.grid_n < 5 or self.grid_n % 2 == 0:
            raise ValueError("main-lobe grid size must be an odd integer >= 5")
        self.du = 2.0 * self.grid_half_width_dc / float(self.grid_n - 1)
        self.component_masks: list[np.ndarray] = []
        self.u0: list[float] = []
        self.v0: list[float] = []
        self.peak_gain: list[float] = []
        self._build()

    def _build(self) -> None:
        offset = np.linspace(-self.grid_half_width_dc, self.grid_half_width_dc, self.grid_n)
        du, dv = np.meshgrid(offset, offset, indexing="xy")
        center_i = self.grid_n // 2
        for beam, steer in enumerate(self.beam_vecs):
            u = steer[0] + du
            v = steer[1] + dv
            valid = u * u + v * v < 1.0
            w = np.full_like(u, np.nan)
            w[valid] = -np.sqrt(np.maximum(0.0, 1.0 - u[valid] * u[valid] - v[valid] * v[valid]))
            dirs = np.stack([u, v, w], axis=-1)
            peak = float(pg.tx_power_gain(steer[None, :], beam, beam_vecs=self.beam_vecs)[0])
            rel = np.full(u.shape, np.nan, dtype=np.float64)
            rel[valid] = pg.tx_power_gain(dirs[valid], beam, beam_vecs=self.beam_vecs) / max(peak, 1e-300)
            above = np.isfinite(rel) & (rel >= self.threshold_linear)
            component = np.zeros_like(above, dtype=bool)
            if above[center_i, center_i]:
                q: deque[tuple[int, int]] = deque([(center_i, center_i)])
                component[center_i, center_i] = True
                while q:
                    y, x = q.popleft()
                    for yy, xx in ((y - 1, x), (y + 1, x), (y, x - 1), (y, x + 1)):
                        if 0 <= yy < self.grid_n and 0 <= xx < self.grid_n and above[yy, xx] and not component[yy, xx]:
                            component[yy, xx] = True
                            q.append((yy, xx))
            self.component_masks.append(component)
            self.u0.append(float(steer[0] - self.grid_half_width_dc))
            self.v0.append(float(steer[1] - self.grid_half_width_dc))
            self.peak_gain.append(peak)

    def contains(self, uvw: np.ndarray, beam_id: np.ndarray) -> np.ndarray:
        dirs = np.asarray(uvw, dtype=np.float64)
        beams = np.asarray(beam_id, dtype=np.int64)
        out = np.zeros(len(dirs), dtype=bool)
        for beam in np.unique(beams[(beams >= 0) & (beams < len(self.beam_vecs))]):
            idx = np.where(beams == beam)[0]
            ix = np.rint((dirs[idx, 0] - self.u0[beam]) / self.du).astype(np.int64)
            iy = np.rint((dirs[idx, 1] - self.v0[beam]) / self.du).astype(np.int64)
            inside = (ix >= 0) & (ix < self.grid_n) & (iy >= 0) & (iy < self.grid_n)
            if np.any(inside):
                out[idx[inside]] = self.component_masks[beam][iy[inside], ix[inside]]
        return out


def selected_hypothesis_group(handle: h5py.File) -> h5py.Group:
    selected = decode_attr(handle.attrs.get("selected_hypothesis", ""))
    hypotheses = handle["hypotheses"]
    if selected and selected in hypotheses:
        return hypotheses[selected]
    best_name = None
    best_rank = np.inf
    for name, group in hypotheses.items():
        rank = float(group.attrs.get("combined_rank", np.inf))
        if rank < best_rank:
            best_rank = rank
            best_name = name
    if best_name is None:
        raise RuntimeError("no selected hypothesis and no ranked hypotheses")
    return hypotheses[best_name]


def read_event(path: Path) -> dict[str, np.ndarray | float | str | int]:
    with h5py.File(path, "r") as handle:
        group = selected_hypothesis_group(handle)
        epoch = float(handle.attrs.get("sample_epoch_unix", handle.attrs.get("event_epoch_unix", np.nan)))
        sample_idx = int(handle.attrs.get("sample_idx", round(epoch * 1e6) if np.isfinite(epoch) else -1))
        selected = group.name.rsplit("/", 1)[-1]
        return {
            "path": str(path),
            "sample_idx": sample_idx,
            "epoch_unix": epoch,
            "selected_hypothesis": selected,
            "t_rel_s": np.asarray(group["t_rel_s"][()], dtype=np.float64),
            "range_km": np.asarray(group["range_km"][()], dtype=np.float64),
            "snr": np.asarray(group["snr"][()], dtype=np.float64),
            "beam_id": np.asarray(group["beam_id"][()], dtype=np.int64),
            "uvw": np.asarray(group["direction_cosines_uvw"][()], dtype=np.float64),
            "selection_keep": np.asarray(group.get("selection_keep", np.ones(len(group["snr"]), dtype=bool)), dtype=bool),
        }


def orbit_files_for_day(root: Path, day: str) -> list[Path]:
    return sorted(Path(root).glob(f"{day}T*/orbit@*.h5"))


def read_orbit_metadata_file(path: Path) -> list[dict[str, np.ndarray | float | str | int]]:
    out: list[dict[str, np.ndarray | float | str | int]] = []
    with h5py.File(path, "r") as handle:
        events = handle["events"][()] if "events" in handle else np.zeros(0, dtype=omt.EVENT_DTYPE)
        paths = handle["paths"][()] if "paths" in handle else np.zeros(0, dtype=omt.PATH_DTYPE)
    events = np.asarray(events)
    paths = np.asarray(paths)
    for event in events:
        sample_idx = int(event["sample_idx"])
        p = paths[paths["sample_idx"] == sample_idx] if len(paths) else np.zeros(0, dtype=omt.PATH_DTYPE)
        if len(p) == 0:
            continue
        position = np.asarray(p["position_enu_km"], dtype=np.float64)
        range_km = np.linalg.norm(position, axis=1)
        uvw = np.full_like(position, np.nan, dtype=np.float64)
        good_range = range_km > 0.0
        uvw[good_range, 0] = position[good_range, 0] / range_km[good_range]
        uvw[good_range, 1] = position[good_range, 1] / range_km[good_range]
        uvw[good_range, 2] = -position[good_range, 2] / range_km[good_range]
        selected = event["selected_hypothesis"]
        out.append(
            {
                "path": str(path),
                "sample_idx": sample_idx,
                "epoch_unix": sample_idx / 1e6,
                "selected_hypothesis": decode_attr(selected).strip("\x00"),
                "t_rel_s": np.asarray(p["t_rel_s"], dtype=np.float64),
                "range_km": range_km,
                "snr": np.asarray(p["snr"], dtype=np.float64),
                "beam_id": np.asarray(p["beam_id"], dtype=np.int64),
                "uvw": uvw,
                "selection_keep": np.ones(len(p), dtype=bool),
            }
        )
    return out


def sky_temperature_at_events(
    times_s: np.ndarray,
    beam_id: np.ndarray,
    n_az: int,
    n_el: int,
    min_elevation_deg: float,
    element_pattern_blend: float,
    sky_step_seconds: float,
) -> np.ndarray:
    if len(times_s) == 0:
        return np.asarray([], dtype=np.float64)
    step = float(sky_step_seconds)
    t0 = np.floor(np.nanmin(times_s) / step) * step
    t1 = np.ceil(np.nanmax(times_s) / step) * step
    model_times = np.arange(t0, t1 + 0.5 * step, step, dtype=np.float64)
    if len(model_times) == 1:
        model_times = np.array([model_times[0], model_times[0] + step], dtype=np.float64)
    model_tsky = noise_fit.receive_sky_temperature_matrix(
        model_times,
        n_az=n_az,
        n_el=n_el,
        min_elevation_deg=min_elevation_deg,
        include_element_pattern=True,
        element_pattern_blend=element_pattern_blend,
    )
    out = np.full(len(times_s), np.nan, dtype=np.float64)
    for beam in np.unique(beam_id[(beam_id >= 0) & (beam_id < model_tsky.shape[1])]):
        idx = beam_id == beam
        out[idx] = np.interp(times_s[idx], model_times, model_tsky[:, beam])
    return out


def estimate_event_rcs(
    event: dict[str, np.ndarray | float | str | int],
    lobe: MainLobeMask,
    args: argparse.Namespace,
) -> dict[str, np.ndarray | float | str | int]:
    beam_id = np.asarray(event["beam_id"], dtype=np.int64)
    uvw = np.asarray(event["uvw"], dtype=np.float64)
    snr_median_noise = np.asarray(event["snr"], dtype=np.float64)
    snr = median_referenced_snr_to_mean_referenced_snr(snr_median_noise, float(args.snr_median_to_mean))
    range_m = np.asarray(event["range_km"], dtype=np.float64) * 1e3
    times_s = float(event["epoch_unix"]) + np.asarray(event["t_rel_s"], dtype=np.float64)
    selection_keep = np.asarray(event["selection_keep"], dtype=bool)

    rel_tx_gain = np.full(len(snr), np.nan, dtype=np.float64)
    for beam in np.unique(beam_id[(beam_id >= 0) & (beam_id < len(lobe.beam_vecs))]):
        idx = beam_id == beam
        rel_tx_gain[idx] = pg.tx_power_gain(uvw[idx], int(beam), beam_vecs=lobe.beam_vecs) / max(lobe.peak_gain[beam], 1e-300)
    main_lobe = lobe.contains(uvw, beam_id)
    finite = (
        selection_keep
        & main_lobe
        & np.isfinite(snr)
        & (snr > 0.0)
        & np.isfinite(range_m)
        & (range_m > 0.0)
        & np.isfinite(rel_tx_gain)
        & (rel_tx_gain > 0.0)
    )

    tsky_k = sky_temperature_at_events(
        times_s,
        beam_id,
        n_az=args.sky_n_az,
        n_el=args.sky_n_el,
        min_elevation_deg=args.sky_min_elevation_deg,
        element_pattern_blend=args.element_pattern_blend,
        sky_step_seconds=args.sky_step_seconds,
    )
    tsys_k = tsky_k + float(args.trec_k)
    noise_power_w = sc.Boltzmann * tsys_k * float(args.analysis_bandwidth_hz)
    received_echo_power_w = snr * noise_power_w

    tx_gain_linear = db_to_linear(float(args.tx_peak_gain_dbi)) * rel_tx_gain
    rx_gain_linear = db_to_linear(float(args.rx_module_gain_dbi))
    sigma_m2 = (
        received_echo_power_w
        * (4.0 * np.pi) ** 3
        * range_m**4
        / (
            float(args.tx_power_w)
            * tx_gain_linear
            * rx_gain_linear
            * pc.wavelength**2
        )
    )
    sigma_m2[~finite] = np.nan
    sigma_dbsm = linear_to_db(sigma_m2)
    sigma_fractional_error = np.full(len(snr), np.nan, dtype=np.float64)
    sigma_fractional_error[finite] = float(args.system_noise_fractional_error)
    return {
        **event,
        "times_s": times_s,
        "tx_relative_gain": rel_tx_gain,
        "tx_relative_gain_db": linear_to_db(rel_tx_gain),
        "in_tx_main_lobe_3db": main_lobe,
        "used_for_rcs": finite,
        "tsky_k": tsky_k,
        "tsys_k": tsys_k,
        "snr_median_noise": snr_median_noise,
        "snr_mean_noise": snr,
        "noise_power_w": noise_power_w,
        "received_echo_power_w": received_echo_power_w,
        "sigma_m2": sigma_m2,
        "sigma_dbsm": sigma_dbsm,
        "sigma_fractional_error": sigma_fractional_error,
    }


def write_event_group(parent: h5py.Group, name: str, result: dict[str, np.ndarray | float | str | int]) -> None:
    group = parent.create_group(name)
    for key, value in result.items():
        if isinstance(value, np.ndarray):
            group.create_dataset(key, data=value, compression="gzip", shuffle=True)
        elif isinstance(value, (float, int, np.floating, np.integer)):
            group.attrs[key] = value
        else:
            group.attrs[key] = str(value)
    used = np.asarray(result["used_for_rcs"], dtype=bool)
    sigma = np.asarray(result["sigma_m2"], dtype=np.float64)
    sigma_db = np.asarray(result["sigma_dbsm"], dtype=np.float64)
    group.attrs["n_points"] = int(len(used))
    group.attrs["n_used_for_rcs"] = int(np.count_nonzero(used))
    group.attrs["median_sigma_m2"] = float(np.nanmedian(sigma[used])) if np.any(used) else np.nan
    group.attrs["median_sigma_dbsm"] = float(np.nanmedian(sigma_db[used])) if np.any(used) else np.nan
    group.attrs["mean_sigma_m2"] = float(np.nanmean(sigma[used])) if np.any(used) else np.nan


def find_event_files(args: argparse.Namespace) -> list[Path]:
    files: list[Path] = []
    for item in args.event_file:
        files.append(Path(item))
    if args.events_dir and args.day:
        day_dir = Path(args.events_dir) / args.day
        files.extend(sorted(day_dir.glob("pansy_disambiguation_diagnostics_*.h5")))
    if args.limit:
        files = files[: int(args.limit)]
    if not files:
        raise FileNotFoundError("no event HDF5 files selected")
    return files


def iter_input_events(args: argparse.Namespace):
    if args.orbit_metadata_dir is not None:
        orbit_files: list[Path] = []
        for item in args.orbit_file:
            orbit_files.append(Path(item))
        if args.day:
            orbit_files.extend(orbit_files_for_day(Path(args.orbit_metadata_dir), args.day))
        if args.limit_files:
            orbit_files = orbit_files[: int(args.limit_files)]
        n_events = 0
        for path in orbit_files:
            for event in read_orbit_metadata_file(path):
                yield event
                n_events += 1
                if args.limit and n_events >= args.limit:
                    return
        if n_events:
            return
    for path in find_event_files(args):
        yield read_event(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-metadata-dir", type=Path, default=Path("/mnt/data/juha/pansy/metadata/orbit"))
    parser.add_argument("--orbit-file", action="append", default=[])
    parser.add_argument("--events-dir", type=Path, default=Path("/mnt/data/juha/pansy/events"))
    parser.add_argument("--day", default=None)
    parser.add_argument("--event-file", action="append", default=[])
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--limit-files", type=int, default=0)
    parser.add_argument("--output", type=Path, default=Path("figs/meteor_rcs_estimates.h5"))
    parser.add_argument("--tx-power-w", type=float, default=DEFAULT_TX_POWER_W)
    parser.add_argument("--tx-modules", type=int, default=DEFAULT_TX_MODULES)
    parser.add_argument("--yagis-per-module", type=int, default=DEFAULT_YAGIS_PER_MODULE)
    parser.add_argument("--yagi-gain-dbi", type=float, default=DEFAULT_YAGI_GAIN_DBI)
    parser.add_argument("--tx-peak-gain-dbi", type=float, default=None)
    parser.add_argument("--rx-module-gain-dbi", type=float, default=None)
    parser.add_argument("--baud-s", type=float, default=DEFAULT_BAUD_S)
    parser.add_argument("--code-bits", type=int, default=DEFAULT_CODE_BITS)
    parser.add_argument("--analysis-bandwidth-hz", type=float, default=None)
    parser.add_argument("--summed-rx-channels", type=float, default=DEFAULT_SUMMED_RX_CHANNELS)
    parser.add_argument("--snr-median-to-mean", type=float, default=None)
    parser.add_argument("--trec-k", type=float, default=DEFAULT_TREC_K)
    parser.add_argument("--system-noise-fractional-error", type=float, default=DEFAULT_TSYS_FRACTIONAL_ERROR)
    parser.add_argument("--element-pattern-blend", type=float, default=DEFAULT_ELEMENT_PATTERN_BLEND)
    parser.add_argument("--sky-n-az", type=int, default=180)
    parser.add_argument("--sky-n-el", type=int, default=60)
    parser.add_argument("--sky-min-elevation-deg", type=float, default=0.5)
    parser.add_argument("--sky-step-seconds", type=float, default=60.0)
    parser.add_argument("--main-lobe-threshold-db", type=float, default=-3.0)
    parser.add_argument("--main-lobe-grid-half-width-dc", type=float, default=0.05)
    parser.add_argument("--main-lobe-grid-n", type=int, default=401)
    args = parser.parse_args()

    args.tx_peak_gain_dbi = (
        absolute_tx_peak_gain_dbi(args.yagi_gain_dbi, args.tx_modules, args.yagis_per_module)
        if args.tx_peak_gain_dbi is None
        else float(args.tx_peak_gain_dbi)
    )
    args.rx_module_gain_dbi = (
        absolute_rx_module_gain_dbi(args.yagi_gain_dbi, args.yagis_per_module)
        if args.rx_module_gain_dbi is None
        else float(args.rx_module_gain_dbi)
    )
    args.analysis_bandwidth_hz = (
        1.0 / (float(args.baud_s) * float(args.code_bits))
        if args.analysis_bandwidth_hz is None
        else float(args.analysis_bandwidth_hz)
    )
    args.snr_median_to_mean = (
        gamma_power_median_to_mean(float(args.summed_rx_channels))
        if args.snr_median_to_mean is None
        else float(args.snr_median_to_mean)
    )

    if args.orbit_metadata_dir is None and not args.event_file:
        raise FileNotFoundError("select --day/--orbit-file orbit metadata, or pass --event-file")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    beam_vecs = pg.tx_beam_unit_vectors()
    lobe = MainLobeMask(
        beam_vecs=beam_vecs,
        threshold_db=args.main_lobe_threshold_db,
        grid_half_width_dc=args.main_lobe_grid_half_width_dc,
        grid_n=args.main_lobe_grid_n,
    )
    with h5py.File(args.output, "w") as out:
        out.attrs["description"] = "Meteor RCS estimates from selected disambiguation hypotheses."
        out.attrs["tx_peak_gain_dbi"] = float(args.tx_peak_gain_dbi)
        out.attrs["rx_module_gain_dbi"] = float(args.rx_module_gain_dbi)
        out.attrs["tx_power_w"] = float(args.tx_power_w)
        out.attrs["analysis_bandwidth_hz"] = float(args.analysis_bandwidth_hz)
        out.attrs["snr_source"] = "stored path_snr = (range-Doppler power - median noise floor) / median noise floor"
        out.attrs["summed_rx_channels"] = float(args.summed_rx_channels)
        out.attrs["snr_median_to_mean"] = float(args.snr_median_to_mean)
        out.attrs["trec_k"] = float(args.trec_k)
        out.attrs["system_noise_fractional_error"] = float(args.system_noise_fractional_error)
        out.attrs["wavelength_m"] = float(pc.wavelength)
        out.attrs["main_lobe_threshold_db"] = float(args.main_lobe_threshold_db)
        out.attrs["formula"] = "sigma = SNR_mean k_B Tsys B (4pi)^3 R^4 / (Pt Gt Gr lambda^2)"
        events = out.create_group("events")
        n_written = 0
        for idx, event in enumerate(iter_input_events(args)):
            try:
                result = estimate_event_rcs(event, lobe, args)
                write_event_group(events, f"event_{idx:06d}", result)
                used = int(np.count_nonzero(result["used_for_rcs"]))
                print(f"{result['path']} sample {result['sample_idx']} used {used}/{len(result['snr'])}")
                n_written += 1
            except Exception as exc:
                failed = events.create_group(f"event_{idx:06d}")
                failed.attrs["path"] = str(event.get("path", ""))
                failed.attrs["sample_idx"] = int(event.get("sample_idx", -1))
                failed.attrs["error"] = repr(exc)
                print(f"{event.get('path', '')} failed {exc!r}")
                n_written += 1
        if n_written == 0:
            raise FileNotFoundError("no orbit/event records with path data were selected")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
