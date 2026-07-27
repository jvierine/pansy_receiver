#!/usr/bin/env python3
"""Interactive tuner for radiant spatial-frequency filtering."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.widgets import Button, Slider

from radiant_spatial_frequency_filter import tapered_spatial_frequency_split


SIDECAR = Path("figs/paper_radiant_results_current/paper_radiant_results.h5")


def add_source_markers(ax) -> None:
    for lon_deg, marker, color, size in (
        (0.0, "o", "#ffd21f", 22),
        (-90.0, r"$\otimes$", "black", 45),
        (180.0, "o", "black", 22),
    ):
        ax.scatter(
            np.deg2rad(lon_deg),
            0.0,
            marker=marker,
            s=size,
            color=color,
            edgecolor="black" if marker == "o" else None,
            linewidth=0.3 if marker == "o" else 0.0,
            zorder=12,
        )


def style_hammer(ax):
    ax.set_xticks(np.deg2rad([-90.0, 0.0, 90.0]))
    ax.set_xticklabels([r"$270^\circ$", r"$0^\circ$", r"$90^\circ$"])
    ax.set_yticks(np.deg2rad([-60.0, -30.0, 0.0, 30.0, 60.0]))
    ax.set_yticklabels([r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$30^\circ$", r"$60^\circ$"])
    ax.grid(True, alpha=0.25, lw=0.45)


def plot_grid(ax, values, xedges_deg, yedges_deg, *, cmap, norm):
    x, y = np.meshgrid(np.deg2rad(xedges_deg), np.deg2rad(yedges_deg))
    mesh = ax.pcolormesh(x, y, values, shading="auto", cmap=cmap, norm=norm, rasterized=True)
    style_hammer(ax)
    add_source_markers(ax)
    return mesh


def positive_log_count_norm(values: np.ndarray, percentile: float = 99.8) -> LogNorm:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite) & (finite > 0.0)]
    vmax = max(2.0, float(np.nanpercentile(finite, percentile))) if finite.size else 2.0
    return LogNorm(vmin=1.0, vmax=vmax)


def signed_log_count_norm(values: np.ndarray, percentile: float = 99.5) -> SymLogNorm:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    vmax = max(2.0, float(np.nanpercentile(np.abs(finite), percentile))) if finite.size else 2.0
    return SymLogNorm(linthresh=1.0, linscale=0.7, vmin=-vmax, vmax=vmax, base=10.0)


class FilterGui:
    def __init__(self, raw, xedges, yedges):
        self.raw = raw
        self.xedges = xedges
        self.yedges = yedges
        self.fig = plt.figure(figsize=(13.0, 7.6))
        self.fig.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.23, wspace=0.16)
        self.ax_low = self.fig.add_subplot(1, 2, 1, projection="hammer")
        self.ax_high = self.fig.add_subplot(1, 2, 2, projection="hammer")
        self.sliders = self.make_sliders()
        self.save_button = Button(self.fig.add_axes([0.82, 0.045, 0.11, 0.045]), "Print args")
        self.save_button.on_clicked(self.print_args)
        self.mesh_low = None
        self.mesh_high = None
        for slider in self.sliders.values():
            slider.on_changed(lambda _value: self.update())
        self.update()

    def make_sliders(self):
        specs = [
            ("zonal_pass", "Zonal pass", 8.0, 1.0, 40.0, 0.05),
            ("zonal_taper", "Zonal taper", 10.0, 1.0, 40.0, 0.10),
            ("meridional_pass", "Meridional pass", 5.0, 1.0, 24.0, 0.05),
            ("meridional_taper", "Meridional taper", 7.0, 1.0, 24.0, 0.10),
        ]
        sliders = {}
        for i, (name, label, value, vmin, vmax, step) in enumerate(specs):
            ax = self.fig.add_axes([0.10, 0.145 - i * 0.035, 0.62, 0.02])
            sliders[name] = Slider(ax, label, vmin, vmax, valinit=value, valstep=step)
        return sliders

    def values(self):
        return {name: float(slider.val) for name, slider in self.sliders.items()}

    def update(self):
        values = self.values()
        lowpass, highpass = tapered_spatial_frequency_split(
            self.raw,
            values["zonal_pass"],
            values["zonal_taper"],
            values["meridional_pass"],
            values["meridional_taper"],
        )
        self.ax_low.clear()
        self.ax_high.clear()
        self.mesh_low = plot_grid(
            self.ax_low,
            np.clip(lowpass, 0.0, None),
            self.xedges,
            self.yedges,
            cmap="magma",
            norm=positive_log_count_norm(lowpass),
        )
        self.mesh_high = plot_grid(
            self.ax_high,
            highpass,
            self.xedges,
            self.yedges,
            cmap="coolwarm",
            norm=signed_log_count_norm(highpass),
        )
        self.ax_low.set_title(
            f"Low spatial frequencies: zonal {values['zonal_pass']:.1f}+{values['zonal_taper']:.1f}, "
            f"meridional {values['meridional_pass']:.1f}+{values['meridional_taper']:.1f}"
        )
        self.ax_high.set_title("High spatial frequencies, signed count residual")
        self.fig.canvas.draw_idle()

    def print_args(self, _event=None):
        values = self.values()
        print(
            "--zonal-pass {zonal_pass:.2f} --zonal-taper {zonal_taper:.2f} "
            "--meridional-pass {meridional_pass:.2f} --meridional-taper {meridional_taper:.2f}".format(**values)
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sidecar", type=Path, default=SIDECAR)
    args = parser.parse_args()
    with h5py.File(args.sidecar, "r") as h5:
        raw = np.asarray(h5["raw_count"], dtype=np.float64)
        xedges = np.asarray(h5["plot_longitude_edges_deg"], dtype=np.float64)
        yedges = np.asarray(h5["ecliptic_latitude_edges_deg"], dtype=np.float64)
    raw = np.where(np.isfinite(raw) & (raw > 0.0), raw, 0.0)
    FilterGui(raw, xedges, yedges)
    plt.show()


if __name__ == "__main__":
    main()
