#!/usr/bin/env python3
"""Interactively select upper/lower Figure 6 height bands and compare subsets."""

from __future__ import annotations

import argparse
import datetime as dt
from dataclasses import dataclass
from pathlib import Path

import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.path import Path as MplPath
from matplotlib.widgets import Button

import orbit_metadata_table as omt
from healpix_hammer import render_healpix_hammer
from plot_orbit_catalogue_statistics import coerce_structured, iter_orbit_tables
from plot_paper_radiant_results import HEALPIX_NSIDE, PLOT_CENTER_LONGITUDE_DEG, style_hammer


DEFAULT_STATISTICS_H5 = Path("figs/paper_refresh_20260721_current/paper_orbit_catalogue_statistics.h5")
DEFAULT_RADIANTS_H5 = Path("figs/paper_refresh_20260721_current/paper_radiant_results.h5")
DEFAULT_SELECTION_H5 = Path("figs/height_band_selection.h5")
DEFAULT_COMPARISON_PNG = Path("figs/height_band_comparison.png")


@dataclass
class BandData:
    rows: np.ndarray
    sample_idx: np.ndarray
    height_km: np.ndarray
    speed_km_s: np.ndarray
    mask_upper: np.ndarray
    mask_lower: np.ndarray


def read_polygon(h5: h5py.File, name: str) -> np.ndarray:
    if f"polygons/{name}" not in h5:
        return np.empty((0, 2), dtype=np.float64)
    return np.asarray(h5[f"polygons/{name}"], dtype=np.float64).reshape(-1, 2)


def load_polygons(path: Path) -> dict[str, np.ndarray]:
    if not path.exists():
        return {"upper": np.empty((0, 2)), "lower": np.empty((0, 2))}
    with h5py.File(path, "r") as h5:
        return {"upper": read_polygon(h5, "upper"), "lower": read_polygon(h5, "lower")}


def mask_from_polygon(x: np.ndarray, y: np.ndarray, polygon: np.ndarray) -> np.ndarray:
    if len(polygon) < 3:
        return np.zeros(len(x), dtype=bool)
    points = np.column_stack((np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)))
    return MplPath(np.asarray(polygon, dtype=np.float64)).contains_points(points)


def load_band_data(radiants_h5: Path, polygons: dict[str, np.ndarray]) -> BandData:
    with h5py.File(radiants_h5, "r") as h5:
        rows = h5["radiants"][()]
    height = np.asarray(rows["first_alt_km"], dtype=np.float64)
    speed = np.asarray(rows["speed_km_s"], dtype=np.float64)
    good = np.isfinite(height) & np.isfinite(speed) & (height > 0.0) & (speed > 0.0)
    rows = rows[good]
    height = height[good]
    speed = speed[good]
    return BandData(
        rows=rows,
        sample_idx=np.asarray(rows["sample_idx"], dtype=np.int64),
        height_km=height,
        speed_km_s=speed,
        mask_upper=mask_from_polygon(speed, height, polygons["upper"]),
        mask_lower=mask_from_polygon(speed, height, polygons["lower"]),
    )


def load_statistics(statistics_h5: Path):
    with h5py.File(statistics_h5, "r") as h5:
        return (
            np.asarray(h5["height_velocity_count"], dtype=np.float64),
            np.asarray(h5["height_edges_km"], dtype=np.float64),
            np.asarray(h5["speed_edges_km_s"], dtype=np.float64),
        )


def write_band_selection_group(h5: h5py.File, name: str, rows: np.ndarray, mask: np.ndarray) -> np.ndarray:
    selected = rows[mask]
    sample_idx = np.asarray(selected["sample_idx"], dtype=np.int64)
    group = h5.create_group(name)
    group.attrs["count"] = int(len(selected))
    group.attrs["event_id_definition"] = "event_id is the orbit/radiant event sample_idx at 1 MHz sample rate"
    group.create_dataset("event_id", data=sample_idx, compression="gzip", shuffle=True)
    group.create_dataset("sample_idx", data=sample_idx, compression="gzip", shuffle=True)
    group.create_dataset("radiants", data=selected, compression="gzip", shuffle=True)
    scalar_fields = {
        "initial_detection_height_km": "first_alt_km",
        "v_g_km_s": "speed_km_s",
        "radiant_ra_gcrs_deg": "radiant_ra_gcrs_deg",
        "radiant_dec_gcrs_deg": "radiant_dec_gcrs_deg",
        "radiant_lambda_ecliptic_deg": "radiant_lambda_ecliptic_deg",
        "radiant_beta_ecliptic_deg": "radiant_beta_ecliptic_deg",
        "lambda_minus_sun_deg": "lambda_minus_sun_deg",
        "plot_longitude_deg": "plot_longitude_deg",
        "solar_longitude_deg": "sun_lambda_ecliptic_deg",
        "combined_score": "combined_score",
    }
    for out_name, row_name in scalar_fields.items():
        if row_name in selected.dtype.names:
            group.create_dataset(out_name, data=np.asarray(selected[row_name]), compression="gzip", shuffle=True)
    if "hypothesis" in selected.dtype.names:
        group.create_dataset("hypothesis", data=np.asarray(selected["hypothesis"]), compression="gzip", shuffle=True)
    return sample_idx


def save_selection(path: Path, polygons: dict[str, np.ndarray], data: BandData, radiants_h5: Path, statistics_h5: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h5:
        h5.attrs["script"] = Path(__file__).name
        h5.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        h5.attrs["radiants_h5"] = str(radiants_h5)
        h5.attrs["statistics_h5"] = str(statistics_h5)
        h5.attrs["polygon_columns"] = "speed_km_s, initial_detection_height_km"
        h5.attrs["schema"] = "height_band_selection_v2"
        h5.attrs["event_id_definition"] = "event_id is the orbit/radiant event sample_idx at 1 MHz sample rate"
        h5.attrs["selection_source"] = "closed polygons in geocentric speed versus initial detection height"
        h5.attrs["upper_count"] = int(np.sum(data.mask_upper))
        h5.attrs["lower_count"] = int(np.sum(data.mask_lower))
        grp = h5.create_group("polygons")
        for name, polygon in polygons.items():
            grp.create_dataset(name, data=np.asarray(polygon, dtype=np.float32))
        upper_sample_idx = write_band_selection_group(h5, "upper", data.rows, data.mask_upper)
        lower_sample_idx = write_band_selection_group(h5, "lower", data.rows, data.mask_lower)
        h5.create_dataset("upper_sample_idx", data=upper_sample_idx, compression="gzip", shuffle=True)
        h5.create_dataset("lower_sample_idx", data=lower_sample_idx, compression="gzip", shuffle=True)
        h5.create_dataset("upper_event_id", data=upper_sample_idx, compression="gzip", shuffle=True)
        h5.create_dataset("lower_event_id", data=lower_sample_idx, compression="gzip", shuffle=True)
        selected_event_id = np.concatenate((upper_sample_idx, lower_sample_idx))
        selected_band = np.concatenate(
            (
                np.full(len(upper_sample_idx), b"upper", dtype="S8"),
                np.full(len(lower_sample_idx), b"lower", dtype="S8"),
            )
        )
        order = np.argsort(selected_event_id)
        h5.create_dataset("selected_event_id", data=selected_event_id[order], compression="gzip", shuffle=True)
        h5.create_dataset("selected_sample_idx", data=selected_event_id[order], compression="gzip", shuffle=True)
        h5.create_dataset("selected_band", data=selected_band[order], compression="gzip", shuffle=True)


def healpix_counts(rows: np.ndarray, mask: np.ndarray) -> np.ndarray:
    lon = np.asarray(rows["lambda_minus_sun_deg"], dtype=np.float64)[mask]
    lat = np.asarray(rows["radiant_beta_ecliptic_deg"], dtype=np.float64)[mask]
    good = np.isfinite(lon) & np.isfinite(lat)
    pix = hp.ang2pix(HEALPIX_NSIDE, lon[good], lat[good], lonlat=True)
    return np.bincount(pix, minlength=hp.nside2npix(HEALPIX_NSIDE)).astype(np.float64)


def selected_event_table(orbit_metadata_dir: Path | None, sample_idx: np.ndarray) -> np.ndarray:
    if orbit_metadata_dir is None or len(sample_idx) == 0:
        return np.zeros(0, dtype=omt.EVENT_DTYPE)
    wanted = set(int(s) for s in np.asarray(sample_idx, dtype=np.int64))
    chunks = []
    for _path, events, _paths, _aliases in iter_orbit_tables(orbit_metadata_dir, include_paths=False, include_aliases=False):
        if len(events) == 0:
            continue
        events = coerce_structured(events, omt.EVENT_DTYPE)
        sample = np.asarray(events["sample_idx"], dtype=np.int64)
        keep = np.fromiter((int(s) in wanted for s in sample), dtype=bool, count=len(sample))
        if np.any(keep):
            chunks.append(events[keep])
    if not chunks:
        return np.zeros(0, dtype=omt.EVENT_DTYPE)
    out = np.concatenate(chunks)
    order = np.argsort(out["sample_idx"])
    return out[order]


def summarize_band(events: np.ndarray) -> dict[str, float]:
    out: dict[str, float] = {"n_orbit": float(len(events))}
    if len(events) == 0:
        return out
    for name, scale in (("ceplecha_initial_radius_m", 1e6), ("v_g_km_s", 1.0), ("initial_detection_height_km", 1.0)):
        x = np.asarray(events[name], dtype=np.float64) * scale
        x = x[np.isfinite(x) & (x > 0.0)]
        if len(x):
            out[f"{name}_median"] = float(np.nanmedian(x))
    kep = np.asarray(events["kepler"], dtype=np.float64)
    labels = ("a_au", "e", "i_deg", "raan_deg", "argp_deg", "nu_deg", "q_au")
    for i, label in enumerate(labels):
        x = kep[:, i]
        x = x[np.isfinite(x)]
        if len(x):
            out[f"{label}_median"] = float(np.nanmedian(x))
    return out


def write_summary_text(path: Path, upper_events: np.ndarray, lower_events: np.ndarray, data: BandData) -> None:
    upper = summarize_band(upper_events)
    lower = summarize_band(lower_events)
    lines = [
        "band,n_radiants,n_orbit,median_height_km,median_vg_km_s,median_r0_um,median_q_au,median_i_deg,median_a_au,median_e",
    ]
    for name, events, mask, summary in (
        ("upper", upper_events, data.mask_upper, upper),
        ("lower", lower_events, data.mask_lower, lower),
    ):
        lines.append(
            ",".join(
                [
                    name,
                    str(int(np.sum(mask))),
                    str(int(summary.get("n_orbit", 0.0))),
                    f"{summary.get('initial_detection_height_km_median', np.nan):.6g}",
                    f"{summary.get('v_g_km_s_median', np.nan):.6g}",
                    f"{summary.get('ceplecha_initial_radius_m_median', np.nan):.6g}",
                    f"{summary.get('q_au_median', np.nan):.6g}",
                    f"{summary.get('i_deg_median', np.nan):.6g}",
                    f"{summary.get('a_au_median', np.nan):.6g}",
                    f"{summary.get('e_median', np.nan):.6g}",
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def step_hist(ax, values: np.ndarray, bins: np.ndarray, label: str, color: str) -> None:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return
    hist, edges = np.histogram(values, bins=bins)
    if np.sum(hist) > 0:
        hist = hist / np.sum(hist)
    ax.step(0.5 * (edges[:-1] + edges[1:]), hist, where="mid", label=label, color=color, lw=1.5)


def make_comparison_plot(
    output: Path,
    statistics_h5: Path,
    polygons: dict[str, np.ndarray],
    data: BandData,
    upper_events: np.ndarray,
    lower_events: np.ndarray,
) -> None:
    hv, height_edges, speed_edges = load_statistics(statistics_h5)
    plot_hv = hv.copy()
    plot_hv[plot_hv <= 0.0] = np.nan
    upper_hp = healpix_counts(data.rows, data.mask_upper)
    lower_hp = healpix_counts(data.rows, data.mask_lower)
    positive = np.r_[upper_hp[upper_hp > 0.0], lower_hp[lower_hp > 0.0]]
    radiant_norm = LogNorm(vmin=1.0, vmax=max(1.0, np.nanpercentile(positive, 99.5) if len(positive) else 1.0))

    fig = plt.figure(figsize=(12.0, 10.0), constrained_layout=True)
    gs = fig.add_gridspec(3, 3)
    ax_sel = fig.add_subplot(gs[0, 0])
    ax_up = fig.add_subplot(gs[0, 1], projection="hammer")
    ax_lo = fig.add_subplot(gs[0, 2], projection="hammer")
    ax_qi = fig.add_subplot(gs[1, 0])
    ax_ae = fig.add_subplot(gs[1, 1])
    ax_r0 = fig.add_subplot(gs[1, 2])
    ax_season = fig.add_subplot(gs[2, :])

    ax_sel.pcolormesh(speed_edges, height_edges, plot_hv, shading="auto", cmap="magma", norm=LogNorm(vmin=1.0, vmax=np.nanmax(plot_hv)))
    for name, color in (("upper", "#54a0ff"), ("lower", "#ff9f43")):
        poly = polygons[name]
        if len(poly):
            closed = np.vstack((poly, poly[0]))
            ax_sel.plot(closed[:, 0], closed[:, 1], color=color, lw=1.8, label=name)
    ax_sel.set_xlabel(r"$v_g$ (km s$^{-1}$)")
    ax_sel.set_ylabel("Initial detection height (km)")
    ax_sel.set_title("Selected bands")
    ax_sel.legend(frameon=False)

    render_healpix_hammer(ax_up, upper_hp, HEALPIX_NSIDE, cmap="magma", norm=radiant_norm)
    render_healpix_hammer(ax_lo, lower_hp, HEALPIX_NSIDE, cmap="magma", norm=radiant_norm)
    for ax, title in ((ax_up, f"Upper band, N={np.sum(data.mask_upper):d}"), (ax_lo, f"Lower band, N={np.sum(data.mask_lower):d}")):
        style_hammer(ax)
        ax.set_title(title)

    def kep(events: np.ndarray) -> np.ndarray:
        return np.asarray(events["kepler"], dtype=np.float64) if len(events) else np.empty((0, 7))

    ku = kep(upper_events)
    kl = kep(lower_events)
    ax_qi.scatter(ku[:, 6], ku[:, 2], s=2, alpha=0.16, color="#54a0ff", label="upper")
    ax_qi.scatter(kl[:, 6], kl[:, 2], s=2, alpha=0.16, color="#ff9f43", label="lower")
    ax_qi.set_xlabel("Perihelion distance q (au)")
    ax_qi.set_ylabel("Inclination i (deg)")
    ax_qi.set_xlim(0.0, 1.5)
    ax_qi.set_ylim(0.0, 180.0)
    ax_qi.legend(frameon=False, markerscale=4)

    au = ku[:, 0]
    al = kl[:, 0]
    good_u = np.isfinite(au) & (au > 0.0) & (au < 20.0)
    good_l = np.isfinite(al) & (al > 0.0) & (al < 20.0)
    ax_ae.scatter(au[good_u], ku[:, 1][good_u], s=2, alpha=0.16, color="#54a0ff")
    ax_ae.scatter(al[good_l], kl[:, 1][good_l], s=2, alpha=0.16, color="#ff9f43")
    ax_ae.set_xlabel("Semi-major axis a (au)")
    ax_ae.set_ylabel("Eccentricity e")
    ax_ae.set_xlim(0.0, 10.0)
    ax_ae.set_ylim(0.0, 1.1)

    r0u = np.asarray(upper_events["ceplecha_initial_radius_m"], dtype=np.float64) * 1e6 if len(upper_events) else np.asarray([])
    r0l = np.asarray(lower_events["ceplecha_initial_radius_m"], dtype=np.float64) * 1e6 if len(lower_events) else np.asarray([])
    bins_r0 = np.logspace(0, 4, 80)
    step_hist(ax_r0, r0u[(r0u > 0) & np.isfinite(r0u)], bins_r0, "upper", "#54a0ff")
    step_hist(ax_r0, r0l[(r0l > 0) & np.isfinite(r0l)], bins_r0, "lower", "#ff9f43")
    ax_r0.set_xscale("log")
    ax_r0.set_xlabel(r"$r_0$ ($\mu$m)")
    ax_r0.set_ylabel("Fraction")
    ax_r0.legend(frameon=False)

    season_bins = np.arange(0.0, 361.0, 5.0)
    for mask, label, color in ((data.mask_upper, "upper", "#54a0ff"), (data.mask_lower, "lower", "#ff9f43")):
        step_hist(ax_season, np.asarray(data.rows["sun_lambda_ecliptic_deg"], dtype=np.float64)[mask], season_bins, label, color)
    ax_season.set_xlim(0.0, 360.0)
    ax_season.set_xlabel(r"Solar longitude, $\lambda_\odot$ (deg)")
    ax_season.set_ylabel("Fraction per 5 deg")
    ax_season.legend(frameon=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


class BandSelector:
    def __init__(self, statistics_h5: Path, radiants_h5: Path, selection_h5: Path, comparison_png: Path, orbit_metadata_dir: Path | None):
        self.statistics_h5 = statistics_h5
        self.radiants_h5 = radiants_h5
        self.selection_h5 = selection_h5
        self.comparison_png = comparison_png
        self.orbit_metadata_dir = orbit_metadata_dir
        self.polygons = load_polygons(selection_h5)
        self.active = "upper"
        self.current_vertices: list[tuple[float, float]] = []
        self.current_artists = []
        self.fig, self.ax = plt.subplots(figsize=(8.5, 6.0))
        plt.subplots_adjust(bottom=0.18)
        self._draw_base()
        self._make_buttons()
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)

    def _draw_base(self) -> None:
        hv, height_edges, speed_edges = load_statistics(self.statistics_h5)
        hv[hv <= 0.0] = np.nan
        self.ax.clear()
        self.ax.pcolormesh(speed_edges, height_edges, hv, shading="auto", cmap="magma", norm=LogNorm(vmin=1.0, vmax=np.nanmax(hv)))
        self.ax.set_xlabel(r"Geocentric velocity, $v_g$ (km s$^{-1}$)")
        self.ax.set_ylabel("Initial detection height (km)")
        self.ax.set_title("Draw polygon; choose upper/lower, then save")
        self.ax.set_xlim(float(speed_edges[0]), float(speed_edges[-1]))
        self.ax.set_ylim(float(height_edges[0]), float(height_edges[-1]))
        self._draw_polygons()

    def _draw_polygons(self) -> None:
        drew_any = False
        for name, color in (("upper", "#54a0ff"), ("lower", "#ff9f43")):
            poly = self.polygons[name]
            if len(poly):
                closed = np.vstack((poly, poly[0]))
                self.ax.plot(closed[:, 0], closed[:, 1], color=color, lw=2.0, label=name)
                drew_any = True
        if drew_any:
            self.ax.legend(frameon=False, loc="upper right")

    def _make_buttons(self) -> None:
        specs = [
            ("Upper", 0.12, self.set_upper),
            ("Lower", 0.24, self.set_lower),
            ("Clear", 0.36, self.clear_active),
            ("Save", 0.50, self.save),
            ("Compare", 0.62, self.compare),
            ("Close", 0.76, lambda _event: plt.close(self.fig)),
        ]
        self.buttons = []
        for label, x0, callback in specs:
            axb = self.fig.add_axes([x0, 0.04, 0.10, 0.06])
            button = Button(axb, label)
            button.on_clicked(callback)
            self.buttons.append(button)

    def reset_current(self) -> None:
        self.current_vertices = []
        for artist in self.current_artists:
            try:
                artist.remove()
            except ValueError:
                pass
        self.current_artists = []
        self.fig.canvas.draw_idle()

    def redraw_current(self) -> None:
        for artist in self.current_artists:
            try:
                artist.remove()
            except ValueError:
                pass
        self.current_artists = []
        if not self.current_vertices:
            self.fig.canvas.draw_idle()
            return
        color = "#54a0ff" if self.active == "upper" else "#ff9f43"
        vertices = np.asarray(self.current_vertices, dtype=np.float64)
        line, = self.ax.plot(vertices[:, 0], vertices[:, 1], color=color, lw=2.2, marker="o", ms=4.5, zorder=20)
        self.current_artists.append(line)
        self.fig.canvas.draw_idle()

    def finish_current(self) -> None:
        if len(self.current_vertices) >= 3:
            self.polygons[self.active] = np.asarray(self.current_vertices, dtype=np.float64)
        self.reset_current()
        self._draw_base()
        self.fig.canvas.draw_idle()

    def on_click(self, event) -> None:
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return
        if event.button == 3 or event.dblclick:
            self.finish_current()
            return
        if event.button != 1:
            return
        self.current_vertices.append((float(event.xdata), float(event.ydata)))
        if len(self.current_vertices) >= 3:
            dx = self.current_vertices[-1][0] - self.current_vertices[0][0]
            dy = self.current_vertices[-1][1] - self.current_vertices[0][1]
            if np.hypot(dx, dy) < 1.5:
                self.current_vertices.pop()
                self.finish_current()
                return
        self.redraw_current()

    def set_upper(self, _event) -> None:
        self.reset_current()
        self.active = "upper"
        self.ax.set_title("Drawing upper-band polygon: left-click points, right-click/double-click to finish")
        self.fig.canvas.draw_idle()

    def set_lower(self, _event) -> None:
        self.reset_current()
        self.active = "lower"
        self.ax.set_title("Drawing lower-band polygon: left-click points, right-click/double-click to finish")
        self.fig.canvas.draw_idle()

    def clear_active(self, _event) -> None:
        self.reset_current()
        self.polygons[self.active] = np.empty((0, 2), dtype=np.float64)
        self._draw_base()
        self.fig.canvas.draw_idle()

    def band_data(self) -> BandData:
        return load_band_data(self.radiants_h5, self.polygons)

    def save(self, _event=None) -> None:
        data = self.band_data()
        save_selection(self.selection_h5, self.polygons, data, self.radiants_h5, self.statistics_h5)
        print(f"saved {self.selection_h5} upper={np.sum(data.mask_upper)} lower={np.sum(data.mask_lower)}")

    def compare(self, _event=None) -> None:
        self.save()
        data = self.band_data()
        upper_events = selected_event_table(self.orbit_metadata_dir, data.sample_idx[data.mask_upper])
        lower_events = selected_event_table(self.orbit_metadata_dir, data.sample_idx[data.mask_lower])
        make_comparison_plot(self.comparison_png, self.statistics_h5, self.polygons, data, upper_events, lower_events)
        write_summary_text(self.comparison_png.with_suffix(".csv"), upper_events, lower_events, data)
        print(f"wrote {self.comparison_png}")
        print(f"wrote {self.comparison_png.with_suffix('.csv')}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--statistics-h5", type=Path, default=DEFAULT_STATISTICS_H5)
    parser.add_argument("--radiants-h5", type=Path, default=DEFAULT_RADIANTS_H5)
    parser.add_argument("--selection-h5", type=Path, default=DEFAULT_SELECTION_H5)
    parser.add_argument("--comparison-output", type=Path, default=DEFAULT_COMPARISON_PNG)
    parser.add_argument("--orbit-metadata-dir", type=Path, default=None)
    parser.add_argument("--no-gui", action="store_true", help="Use an existing selection HDF5 and only regenerate comparison outputs.")
    args = parser.parse_args()

    if args.no_gui:
        polygons = load_polygons(args.selection_h5)
        data = load_band_data(args.radiants_h5, polygons)
        save_selection(args.selection_h5, polygons, data, args.radiants_h5, args.statistics_h5)
        upper_events = selected_event_table(args.orbit_metadata_dir, data.sample_idx[data.mask_upper])
        lower_events = selected_event_table(args.orbit_metadata_dir, data.sample_idx[data.mask_lower])
        make_comparison_plot(args.comparison_output, args.statistics_h5, polygons, data, upper_events, lower_events)
        write_summary_text(args.comparison_output.with_suffix(".csv"), upper_events, lower_events, data)
        print(f"upper={np.sum(data.mask_upper)} lower={np.sum(data.mask_lower)}")
        print(args.selection_h5)
        print(args.comparison_output)
        print(args.comparison_output.with_suffix(".csv"))
        return

    BandSelector(args.statistics_h5, args.radiants_h5, args.selection_h5, args.comparison_output, args.orbit_metadata_dir)
    plt.show()


if __name__ == "__main__":
    main()
