#!/usr/bin/env python3
"""Interactive editor for symmetric sporadic-source flux regions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.widgets import Button, Slider


DEFAULT_SIDECAR = Path("figs/paper_radiant_results_current/paper_radiant_results.h5")
DEFAULT_OUTPUT = Path("figs/sporadic_source_regions_manual.json")
PLOT_CENTER_LONGITUDE_DEG = -90.0


DEFAULT = {
    "ha_half_lon_deg": 30.0,
    "ha_beta_abs_deg": 20.0,
    "apex_center_lon_deg": 270.0,
    "apex_half_lon_deg": 30.0,
    "apex_beta_abs_deg": 25.0,
    "narrow_apex_center_lon_deg": 270.0,
    "narrow_apex_half_lon_deg": 15.0,
    "narrow_apex_beta_abs_deg": 10.0,
    "toroidal_center_lon_deg": 270.0,
    "toroidal_half_lon_deg": 45.0,
    "toroidal_beta_center_abs_deg": 50.0,
    "toroidal_beta_half_width_deg": 15.0,
}

COLORS = {
    "helion": "#00a6d6",
    "antihelion": "#00a6d6",
    "apex": "#009e73",
    "narrow_apex": "#d55e00",
    "northern_toroidal": "#cc79a7",
    "southern_toroidal": "#cc79a7",
}


def wrap180(deg):
    return (np.asarray(deg, dtype=np.float64) + 180.0) % 360.0 - 180.0


def wrap360(deg):
    return np.asarray(deg, dtype=np.float64) % 360.0


def physical_lon_from_plot_x(plot_x_deg):
    return wrap360(PLOT_CENTER_LONGITUDE_DEG - plot_x_deg)


def centered_plot_longitude_deg(lambda_minus_sun_deg):
    signed = wrap180(lambda_minus_sun_deg)
    return -wrap180(signed - PLOT_CENTER_LONGITUDE_DEG)


def load_params(path):
    params = dict(DEFAULT)
    if not path.exists():
        return params
    with path.open("r", encoding="utf-8") as fh:
        payload = json.load(fh)
    if "parameters" in payload:
        loaded = dict(payload["parameters"])
        if "toroidal_beta_abs_min_deg" in loaded and "toroidal_beta_abs_max_deg" in loaded:
            tor_min = loaded.pop("toroidal_beta_abs_min_deg")
            tor_max = loaded.pop("toroidal_beta_abs_max_deg")
            loaded["toroidal_beta_center_abs_deg"] = 0.5 * (tor_min + tor_max)
            loaded["toroidal_beta_half_width_deg"] = 0.5 * (tor_max - tor_min)
        params.update(loaded)
    elif "symmetric_groups" in payload:
        groups = payload["symmetric_groups"]
        if "helion_antihelion" in groups:
            params["ha_half_lon_deg"] = groups["helion_antihelion"]["half_lon_deg"]
            params["ha_beta_abs_deg"] = groups["helion_antihelion"]["beta_abs_deg"]
        if "apex" in groups:
            params["apex_center_lon_deg"] = groups["apex"]["center_lon_deg"]
            params["apex_half_lon_deg"] = groups["apex"]["half_lon_deg"]
            params["apex_beta_abs_deg"] = groups["apex"]["beta_abs_deg"]
        if "narrow_apex" in groups:
            params["narrow_apex_center_lon_deg"] = groups["narrow_apex"]["center_lon_deg"]
            params["narrow_apex_half_lon_deg"] = groups["narrow_apex"]["half_lon_deg"]
            params["narrow_apex_beta_abs_deg"] = groups["narrow_apex"]["beta_abs_deg"]
        if "toroidal" in groups:
            params["toroidal_center_lon_deg"] = groups["toroidal"]["center_lon_deg"]
            params["toroidal_half_lon_deg"] = groups["toroidal"]["half_lon_deg"]
            params["toroidal_beta_center_abs_deg"] = 0.5 * (groups["toroidal"]["beta_abs_min_deg"] + groups["toroidal"]["beta_abs_max_deg"])
            params["toroidal_beta_half_width_deg"] = 0.5 * (groups["toroidal"]["beta_abs_max_deg"] - groups["toroidal"]["beta_abs_min_deg"])
    return params


def load_sidecar(path):
    with h5py.File(path, "r") as h5:
        raw = np.asarray(h5["raw_count"], dtype=np.float64)
        rate = np.asarray(h5["debiased_rate_h_inv"], dtype=np.float64)
        exposure = np.asarray(h5["radiant_exposure_hours"], dtype=np.float64)
        xedges = np.asarray(h5["plot_longitude_edges_deg"], dtype=np.float64)
        yedges = np.asarray(h5["ecliptic_latitude_edges_deg"], dtype=np.float64)
    return raw, rate, exposure, xedges, yedges


def regions_from_params(p):
    tor_min = max(0.0, p["toroidal_beta_center_abs_deg"] - p["toroidal_beta_half_width_deg"])
    tor_max = min(80.0, p["toroidal_beta_center_abs_deg"] + p["toroidal_beta_half_width_deg"])
    return {
        "helion": {
            "center_lon_deg": 0.0,
            "half_lon_deg": p["ha_half_lon_deg"],
            "beta_min_deg": -p["ha_beta_abs_deg"],
            "beta_max_deg": p["ha_beta_abs_deg"],
        },
        "antihelion": {
            "center_lon_deg": 180.0,
            "half_lon_deg": p["ha_half_lon_deg"],
            "beta_min_deg": -p["ha_beta_abs_deg"],
            "beta_max_deg": p["ha_beta_abs_deg"],
        },
        "apex": {
            "center_lon_deg": p["apex_center_lon_deg"],
            "half_lon_deg": p["apex_half_lon_deg"],
            "beta_min_deg": -p["apex_beta_abs_deg"],
            "beta_max_deg": p["apex_beta_abs_deg"],
        },
        "narrow_apex": {
            "center_lon_deg": p["narrow_apex_center_lon_deg"],
            "half_lon_deg": p["narrow_apex_half_lon_deg"],
            "beta_min_deg": -p["narrow_apex_beta_abs_deg"],
            "beta_max_deg": p["narrow_apex_beta_abs_deg"],
        },
        "northern_toroidal": {
            "center_lon_deg": p["toroidal_center_lon_deg"],
            "half_lon_deg": p["toroidal_half_lon_deg"],
            "beta_min_deg": tor_min,
            "beta_max_deg": tor_max,
        },
        "southern_toroidal": {
            "center_lon_deg": p["toroidal_center_lon_deg"],
            "half_lon_deg": p["toroidal_half_lon_deg"],
            "beta_min_deg": -tor_max,
            "beta_max_deg": -tor_min,
        },
    }


def style_hammer(ax):
    ticks = np.asarray([-90.0, 0.0, 90.0])
    labels = [f"{int(wrap360(PLOT_CENTER_LONGITUDE_DEG - tick))}$^\\circ$" for tick in ticks]
    ax.set_xticks(np.deg2rad(ticks))
    ax.set_xticklabels(labels)
    ax.set_yticks(np.deg2rad([-60.0, -30.0, 0.0, 30.0, 60.0]))
    ax.set_yticklabels([r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$30^\circ$", r"$60^\circ$"])
    ax.grid(True, alpha=0.24, lw=0.45)


def add_source_markers(ax):
    for lon_deg, marker, color, size in ((0.0, "o", "#ffd21f", 20), (270.0, r"$\otimes$", "black", 42), (180.0, "o", "black", 20)):
        ax.scatter(np.deg2rad(centered_plot_longitude_deg(lon_deg)), 0.0, marker=marker, s=size, color=color, edgecolor="black" if marker == "o" else None, linewidth=0.3 if marker == "o" else 0.0, zorder=10)


def region_mask(region, lon, beta, valid):
    return valid & (np.abs(wrap180(lon - region["center_lon_deg"])) <= region["half_lon_deg"]) & (beta >= region["beta_min_deg"]) & (beta <= region["beta_max_deg"])


def plot_box(ax, region, color, lw=1.35):
    center = region["center_lon_deg"]
    half = region["half_lon_deg"]
    beta_min = region["beta_min_deg"]
    beta_max = region["beta_max_deg"]
    lon_edge = np.linspace(center - half, center + half, 180)
    beta_edge = np.linspace(beta_min, beta_max, 80)
    artists = []
    artists.extend(ax.plot(np.deg2rad(centered_plot_longitude_deg(lon_edge)), np.deg2rad(np.full_like(lon_edge, beta_min)), color=color, lw=lw, alpha=0.95))
    artists.extend(ax.plot(np.deg2rad(centered_plot_longitude_deg(lon_edge)), np.deg2rad(np.full_like(lon_edge, beta_max)), color=color, lw=lw, alpha=0.95))
    artists.extend(ax.plot(np.deg2rad(centered_plot_longitude_deg(np.full_like(beta_edge, center - half))), np.deg2rad(beta_edge), color=color, lw=lw, alpha=0.95))
    artists.extend(ax.plot(np.deg2rad(centered_plot_longitude_deg(np.full_like(beta_edge, center + half))), np.deg2rad(beta_edge), color=color, lw=lw, alpha=0.95))
    return artists


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sidecar", type=Path, default=DEFAULT_SIDECAR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    raw, rate, exposure, xedges, yedges = load_sidecar(args.sidecar)
    params = load_params(args.output)

    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])
    xgrid, beta = np.meshgrid(xcenters, ycenters)
    lon = physical_lon_from_plot_x(xgrid)
    valid = exposure > 0.0

    raw_norm = LogNorm(vmin=1.0, vmax=max(2.0, float(np.nanmax(raw[raw > 0.0]))))
    rate_positive = rate[np.isfinite(rate) & (rate > 0.0)]
    rate_norm = LogNorm(vmin=float(np.nanpercentile(rate_positive, 5.0)), vmax=float(np.nanpercentile(rate_positive, 99.5)))

    fig = plt.figure(figsize=(14.2, 8.2))
    ax_raw = fig.add_axes([0.04, 0.43, 0.39, 0.52], projection="hammer")
    ax_rate = fig.add_axes([0.50, 0.43, 0.39, 0.52], projection="hammer")
    ax_raw.pcolormesh(np.deg2rad(xedges), np.deg2rad(yedges), np.where(raw > 0.0, raw, raw_norm.vmin), shading="auto", cmap="magma", norm=raw_norm)
    ax_rate.pcolormesh(np.deg2rad(xedges), np.deg2rad(yedges), np.where(rate > 0.0, rate, rate_norm.vmin), shading="auto", cmap="magma", norm=rate_norm)
    for ax, title in ((ax_raw, "Raw counts"), (ax_rate, "Speed-debiased rate")):
        style_hammer(ax)
        add_source_markers(ax)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(r"Sun-centered ecliptic longitude, $\lambda-\lambda_\odot$")
        ax.set_ylabel(r"Ecliptic latitude, $\beta$")

    slider_specs = [
        ("ha_half_lon_deg", "H/A half lon", 1.0, 120.0, 0.29, 0.05),
        ("ha_beta_abs_deg", "H/A beta abs", 1.0, 75.0, 0.25, 0.05),
        ("apex_center_lon_deg", "APEX center lon", 220.0, 320.0, 0.36, 0.05),
        ("apex_half_lon_deg", "APEX half lon", 1.0, 120.0, 0.32, 0.05),
        ("apex_beta_abs_deg", "APEX beta abs", 1.0, 75.0, 0.28, 0.05),
        ("narrow_apex_center_lon_deg", "NAR center lon", 220.0, 320.0, 0.36, 0.52),
        ("narrow_apex_half_lon_deg", "NAR half lon", 1.0, 80.0, 0.32, 0.52),
        ("narrow_apex_beta_abs_deg", "NAR beta abs", 1.0, 50.0, 0.28, 0.52),
        ("toroidal_center_lon_deg", "TOR center lon", 220.0, 320.0, 0.20, 0.52),
        ("toroidal_half_lon_deg", "TOR half lon", 1.0, 120.0, 0.16, 0.52),
        ("toroidal_beta_center_abs_deg", "TOR beta center", 0.0, 80.0, 0.12, 0.52),
        ("toroidal_beta_half_width_deg", "TOR beta half width", 0.5, 45.0, 0.08, 0.52),
    ]
    sliders = {}
    for key, label, lo, hi, y, x in slider_specs:
        sliders[key] = Slider(fig.add_axes([x, y, 0.34, 0.024]), label, lo, hi, valinit=params[key], valstep=0.5)

    stats_text = fig.text(0.05, 0.06, "", fontsize=9, family="monospace", va="top")
    artists = []

    def current_params():
        return {key: float(slider.val) for key, slider in sliders.items()}

    def redraw(_=None):
        nonlocal artists
        for artist in artists:
            artist.remove()
        artists = []
        p = current_params()
        regions = regions_from_params(p)
        for name, region in regions.items():
            artists.extend(plot_box(ax_raw, region, COLORS[name], lw=1.7 if name == "apex" else 1.2))
            artists.extend(plot_box(ax_rate, region, COLORS[name], lw=1.7 if name == "apex" else 1.2))
        masks = {name: region_mask(region, lon, beta, valid) for name, region in regions.items()}
        masks["apex"] = masks["apex"] & ~masks["narrow_apex"]
        source_masks = {
            "helion": masks["helion"],
            "antihelion": masks["antihelion"],
            "apex_remainder": masks["apex"],
            "narrow_apex": masks["narrow_apex"],
            "toroidal_total": masks["northern_toroidal"] | masks["southern_toroidal"],
        }
        total_rate = sum(float(np.nansum(rate[m])) for m in source_masks.values())
        total_raw = sum(int(np.nansum(raw[m])) for m in source_masks.values())
        lines = ["region                  rate      raw"]
        for name, mask in source_masks.items():
            value = float(np.nansum(rate[mask]))
            count = int(np.nansum(raw[mask]))
            rp = 100.0 * value / total_rate if total_rate > 0.0 else 0.0
            cp = 100.0 * count / total_raw if total_raw > 0 else 0.0
            lines.append(f"{name:20s} {rp:6.2f}% {cp:6.2f}%")
        stats_text.set_text("\n".join(lines))
        fig.canvas.draw_idle()

    def save(_=None):
        p = current_params()
        payload = {
            "description": "Symmetric manual sporadic-source selection regions from interactive_sporadic_source_regions.py",
            "sidecar": str(args.sidecar),
            "parameters": p,
            "regions": regions_from_params(p),
            "apex_region_note": "Flux estimator should use apex minus narrow_apex for non-overlapping diagnostic regions.",
            "toroidal_region_note": "Toroidal is mirrored north/south and reported as total unless separated explicitly.",
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"saved {args.output}")

    for slider in sliders.values():
        slider.on_changed(redraw)
    save_button = Button(fig.add_axes([0.80, 0.015, 0.07, 0.04]), "Save")
    close_button = Button(fig.add_axes([0.89, 0.015, 0.07, 0.04]), "Close")
    save_button.on_clicked(save)
    close_button.on_clicked(lambda _: plt.close(fig))
    fig._pansy_region_editor_widgets = {
        "sliders": sliders,
        "save_button": save_button,
        "close_button": close_button,
    }

    redraw()
    plt.show()


if __name__ == "__main__":
    main()
