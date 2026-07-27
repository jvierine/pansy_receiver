#!/usr/bin/env python3
"""Qt tuner for iterative spherical-harmonic radiant filtering."""

from __future__ import annotations

import argparse
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Qt5Agg")

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PyQt5 import QtCore, QtWidgets

from healpix_hammer import centered_plot_longitude_deg, render_healpix_hammer
from radiant_spatial_frequency_filter import (
    iterative_positive_spherical_harmonic_split,
    positive_log_norm,
)


SIDECAR = Path("figs/paper_radiant_results_current/paper_radiant_results.h5")
PLOT_CENTER_LONGITUDE_DEG = -90.0


def add_source_markers(ax) -> None:
    for lon_deg, marker, color, edgecolor, linewidth, size in (
        (0.0, "o", "#ffd21f", "black", 0.3, 22),
        (270.0, r"$\otimes$", "black", "none", 0.0, 45),
        (180.0, "o", "black", "white", 0.8, 22),
    ):
        ax.scatter(
            np.deg2rad(centered_plot_longitude_deg(lon_deg, PLOT_CENTER_LONGITUDE_DEG)),
            0.0,
            marker=marker,
            s=size,
            color=color,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=12,
        )


def style_hammer(ax) -> None:
    ax.set_xticks(np.deg2rad([-90.0, 0.0, 90.0]))
    ax.set_xticklabels([])
    ax.set_yticks(np.deg2rad([-60.0, -30.0, 0.0, 30.0, 60.0]))
    ax.set_yticklabels(
        [r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$30^\circ$", r"$60^\circ$"]
    )
    ax.grid(True, alpha=0.25, lw=0.45)
    add_source_markers(ax)


def plot_healpix(ax, values, nside) -> None:
    render_healpix_hammer(
        ax,
        np.where(values > 0.0, values, np.nan),
        nside,
        cmap="magma",
        norm=positive_log_norm(values),
        center_longitude_deg=PLOT_CENTER_LONGITUDE_DEG,
    )
    style_hammer(ax)


CONTROL_SPECS = (
    ("zonal_pass", "Zonal pass |m|", 1.0, 0.0, 40.0, 0.05, 2),
    ("zonal_taper", "Zonal taper", 19.1, 0.1, 40.0, 0.10, 1),
    ("meridional_pass", "Meridional pass l-|m|", 6.4, 0.0, 40.0, 0.05, 2),
    ("meridional_taper", "Meridional taper", 23.9, 0.1, 40.0, 0.10, 1),
    ("positive_iterations", "Recursive iterations", 12.0, 1.0, 30.0, 1.0, 0),
    ("positive_fraction", "Positive fraction", 0.5, 0.05, 1.0, 0.05, 2),
)


class FilterGui(QtWidgets.QMainWindow):
    def __init__(self, raw, nside):
        super().__init__()
        self.raw = raw
        self.nside = int(nside)
        self.executor = ThreadPoolExecutor(max_workers=1)
        self.future: Future | None = None
        self.computed_values = None
        self.spinboxes = {}
        self.sliders = {}

        self.setWindowTitle("Radiant spherical-harmonic filter")
        self.resize(1400, 900)

        central = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(central)
        self.setCentralWidget(central)

        self.figure = Figure(figsize=(13.0, 6.2))
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas, stretch=1)
        self.ax_low = self.figure.add_subplot(1, 2, 1, projection="hammer")
        self.ax_high = self.figure.add_subplot(1, 2, 2, projection="hammer")
        self.figure.subplots_adjust(
            left=0.06, right=0.98, top=0.90, bottom=0.08, wspace=0.16
        )
        for ax, title in (
            (self.ax_low, "Recursive low-frequency background"),
            (self.ax_high, "Final positive high-frequency excess"),
        ):
            ax.set_title(title)
            ax.text(
                0.5,
                0.5,
                "Adjust controls, then press Recompute",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )

        controls = QtWidgets.QGroupBox("Filter settings")
        grid = QtWidgets.QGridLayout(controls)
        layout.addWidget(controls)
        self.make_controls(grid)

        buttons = QtWidgets.QHBoxLayout()
        buttons.addStretch(1)
        self.recompute_button = QtWidgets.QPushButton("Recompute")
        self.print_button = QtWidgets.QPushButton("Print args")
        self.recompute_button.clicked.connect(self.recompute)
        self.print_button.clicked.connect(self.print_args)
        buttons.addWidget(self.recompute_button)
        buttons.addWidget(self.print_button)
        grid.addLayout(buttons, len(CONTROL_SPECS), 0, 1, 3)

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.poll_result)
        self.timer.start(100)

    def make_controls(self, layout) -> None:
        for row, (name, label, value, vmin, vmax, step, decimals) in enumerate(
            CONTROL_SPECS
        ):
            text = QtWidgets.QLabel(label)
            slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
            slider.setRange(0, int(round((vmax - vmin) / step)))
            slider.setSingleStep(1)
            slider.setPageStep(max(1, int(round(1.0 / step))))
            slider.setTracking(True)
            spinbox = QtWidgets.QDoubleSpinBox()
            spinbox.setRange(vmin, vmax)
            spinbox.setSingleStep(step)
            spinbox.setDecimals(decimals)
            spinbox.setValue(value)
            slider.setValue(int(round((value - vmin) / step)))

            slider.valueChanged.connect(
                lambda index, box=spinbox, lo=vmin, delta=step: box.setValue(
                    lo + index * delta
                )
            )
            spinbox.valueChanged.connect(
                lambda current, control=slider, lo=vmin, delta=step: control.setValue(
                    int(round((current - lo) / delta))
                )
            )

            layout.addWidget(text, row, 0)
            layout.addWidget(slider, row, 1)
            layout.addWidget(spinbox, row, 2)
            self.sliders[name] = slider
            self.spinboxes[name] = spinbox

    def values(self):
        values = {
            name: float(spinbox.value()) for name, spinbox in self.spinboxes.items()
        }
        values["positive_iterations"] = int(round(values["positive_iterations"]))
        return values

    def recompute(self) -> None:
        if self.future is not None and not self.future.done():
            return
        values = self.values()
        self.computed_values = values
        self.recompute_button.setText("Working...")
        self.recompute_button.setEnabled(False)
        self.future = self.executor.submit(
            iterative_positive_spherical_harmonic_split,
            self.raw,
            self.nside,
            values["zonal_pass"],
            values["zonal_taper"],
            values["meridional_pass"],
            values["meridional_taper"],
            iterations=values["positive_iterations"],
            positive_fraction=values["positive_fraction"],
        )

    def poll_result(self) -> None:
        if self.future is None or not self.future.done():
            return
        future = self.future
        self.future = None
        try:
            lowpass, highpass, _, history = future.result()
        except Exception as exc:
            self.recompute_button.setText("Recompute")
            self.recompute_button.setEnabled(True)
            QtWidgets.QMessageBox.critical(self, "Filter failed", str(exc))
            return

        values = self.computed_values
        self.ax_low.clear()
        self.ax_high.clear()
        plot_healpix(self.ax_low, lowpass, self.nside)
        plot_healpix(self.ax_high, highpass, self.nside)
        self.ax_low.set_title(
            f"Recursive low-frequency background\n"
            f"|m| {values['zonal_pass']:.1f}+{values['zonal_taper']:.1f}, "
            f"l-|m| {values['meridional_pass']:.1f}+{values['meridional_taper']:.1f}"
        )
        self.ax_high.set_title(
            f"Final positive high-frequency excess\n"
            f"{values['positive_iterations']} iterations, "
            f"fraction {values['positive_fraction']:.2f}; "
            f"last residual sum {history[-1, 0]:.0f}"
        )
        self.canvas.draw_idle()
        self.recompute_button.setText(
            "Recompute" if self.values() == values else "Recompute *"
        )
        self.recompute_button.setEnabled(True)

    def print_args(self) -> None:
        values = self.values()
        print(
            "--zonal-pass {zonal_pass:.2f} --zonal-taper {zonal_taper:.2f} "
            "--meridional-pass {meridional_pass:.2f} "
            "--meridional-taper {meridional_taper:.2f} "
            "--positive-iterations {positive_iterations:d} "
            "--positive-fraction {positive_fraction:.2f}".format(**values),
            flush=True,
        )

    def closeEvent(self, event) -> None:
        self.timer.stop()
        self.executor.shutdown(wait=False, cancel_futures=True)
        super().closeEvent(event)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sidecar", type=Path, default=SIDECAR)
    args = parser.parse_args()
    with h5py.File(args.sidecar, "r") as h5:
        raw = np.asarray(h5["healpix_raw_count"], dtype=np.float64)
        nside = int(h5.attrs["healpix_nside"])
    raw = np.where(np.isfinite(raw) & (raw > 0.0), raw, 0.0)

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    gui = FilterGui(raw, nside)
    gui.show()
    gui.raise_()
    gui.activateWindow()
    app.exec_()


if __name__ == "__main__":
    main()
