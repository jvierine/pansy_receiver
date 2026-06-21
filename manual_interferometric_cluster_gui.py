#!/usr/bin/env python
"""Manual labeling tool for interferometric direction-cosine clusters."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.widgets import Button, LassoSelector
import numpy as np

import pansy_interferometry as pint
import plot_interferometric_disambiguation as disamb


def build_candidates(sample_idx, cut_dir, grid_n, coherence_threshold, snr_threshold, recompute):
    cut = disamb.load_cut(cut_dir, sample_idx)
    obs_all = disamb.recompute_cut_observables(cut) if recompute else disamb.cached_cut_observables(cut)
    good = obs_all["snr"] > snr_threshold
    obs = {
        key: val[good] if isinstance(val, np.ndarray) and len(val) == len(good) else val
        for key, val in obs_all.items()
    }

    antpos = pint.get_antpos()
    ch_pairs = np.asarray(list(itertools.combinations(np.arange(7), 2)))
    dmat = pint.pair_mat(ch_pairs, antpos)
    phasecal = pint.get_phasecal()
    u, v, w, valid = disamb.horizon_grid(grid_n)
    steering, _uvw = disamb.steering_matrix(dmat, u, v, w, valid)

    candidates = []
    if len(obs["tx_idx"]) == 0:
        return obs, candidates
    t_rel = obs["tx_idx"] / 1e6 - obs["tx_idx"][0] / 1e6
    for pulse_i in range(len(obs["snr"])):
        coh = disamb.coherence_map(
            obs["xc"][pulse_i],
            int(obs["beam_id"][pulse_i]),
            phasecal,
            ch_pairs,
            steering,
            valid,
            u.shape,
        )
        ii, jj = disamb.local_coherence_peaks(coh, coherence_threshold)
        for row, col in zip(ii, jj):
            candidates.append(
                {
                    "u": float(u[row, col]),
                    "v": float(v[row, col]),
                    "w": float(w[row, col]),
                    "grid_row": int(row),
                    "grid_col": int(col),
                    "coherence": float(coh[row, col]),
                    "t_rel": float(t_rel[pulse_i]),
                    "pulse": int(pulse_i),
                    "range_km": float(obs["range_km"][pulse_i]),
                    "doppler_mps": float(obs["doppler_mps"][pulse_i]),
                    "snr": float(obs["snr"][pulse_i]),
                    "beam_id": int(obs["beam_id"][pulse_i]),
                }
            )
    return obs, candidates


def load_labels(path, n_candidates):
    if not path.exists():
        return np.zeros(n_candidates, dtype=np.int16)
    with h5py.File(path, "r") as h:
        if "label" not in h or len(h["label"]) != n_candidates:
            return np.zeros(n_candidates, dtype=np.int16)
        return h["label"][:].astype(np.int16)


def save_labels(path, sample_idx, grid_n, coherence_threshold, snr_threshold, candidates, labels):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = {
        "u": np.asarray([c["u"] for c in candidates], dtype=np.float64),
        "v": np.asarray([c["v"] for c in candidates], dtype=np.float64),
        "w": np.asarray([c["w"] for c in candidates], dtype=np.float64),
        "t_rel_s": np.asarray([c["t_rel"] for c in candidates], dtype=np.float64),
        "pulse": np.asarray([c["pulse"] for c in candidates], dtype=np.int64),
        "range_km": np.asarray([c["range_km"] for c in candidates], dtype=np.float64),
        "doppler_mps": np.asarray([c["doppler_mps"] for c in candidates], dtype=np.float64),
        "snr": np.asarray([c["snr"] for c in candidates], dtype=np.float64),
        "beam_id": np.asarray([c["beam_id"] for c in candidates], dtype=np.int64),
        "coherence": np.asarray([c["coherence"] for c in candidates], dtype=np.float64),
        "grid_row": np.asarray([c["grid_row"] for c in candidates], dtype=np.int64),
        "grid_col": np.asarray([c["grid_col"] for c in candidates], dtype=np.int64),
        "label": np.asarray(labels, dtype=np.int16),
    }
    with h5py.File(path, "w") as h:
        h.attrs["sample_idx"] = int(sample_idx)
        h.attrs["grid_n"] = int(grid_n)
        h.attrs["coherence_threshold"] = float(coherence_threshold)
        h.attrs["snr_threshold"] = float(snr_threshold)
        h.attrs["format"] = "manual_interferometric_cluster_labels_v1"
        for name, value in fields.items():
            h.create_dataset(name, data=value)
        labels_unique = sorted(int(x) for x in np.unique(labels) if int(x) > 0)
        clusters = h.create_group("clusters")
        for label in labels_unique:
            clusters.create_dataset(f"C{label:02d}", data=np.flatnonzero(labels == label).astype(np.int64))


class ClusterGui:
    def __init__(self, args, candidates, labels):
        self.args = args
        self.candidates = candidates
        self.labels = labels
        self.current_label = 1
        self.undo_stack = []

        self.u = np.asarray([c["u"] for c in candidates], dtype=np.float64)
        self.v = np.asarray([c["v"] for c in candidates], dtype=np.float64)
        self.t = np.asarray([c["t_rel"] for c in candidates], dtype=np.float64)
        self.pulse = np.asarray([c["pulse"] for c in candidates], dtype=np.int64)
        self.snr = np.asarray([c["snr"] for c in candidates], dtype=np.float64)
        self.points_uv = np.column_stack([self.u, self.v])

        self.fig, self.ax = plt.subplots(figsize=(9.0, 8.4))
        self.fig.subplots_adjust(bottom=0.16)
        self.scatter = None
        self.text = self.fig.text(0.02, 0.035, "", ha="left", va="bottom", family="monospace", fontsize=9)
        self._add_buttons()
        self.lasso = LassoSelector(self.ax, onselect=self.on_lasso)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.redraw()

    def _add_buttons(self):
        buttons = [
            ("Save", 0.58, self.save),
            ("New", 0.68, self.new_cluster),
            ("Clear", 0.78, self.clear_current),
            ("Undo", 0.88, self.undo),
        ]
        for label, left, callback in buttons:
            ax_button = self.fig.add_axes([left, 0.035, 0.08, 0.05])
            button = Button(ax_button, label)
            button.on_clicked(callback)
            setattr(self, f"_button_{label.lower()}", button)

    def colors(self):
        cmap = plt.get_cmap("tab20")
        colors = np.full((len(self.labels), 4), (0.55, 0.55, 0.55, 0.28), dtype=np.float64)
        for label in sorted(int(x) for x in np.unique(self.labels) if int(x) > 0):
            colors[self.labels == label] = cmap((label - 1) % 20)
        return colors

    def redraw(self):
        self.ax.clear()
        self.ax.add_patch(plt.Circle((0.0, 0.0), 1.0, color="black", fill=False, linewidth=1.0))
        if len(self.candidates):
            sizes = np.clip(8.0 + 1.8 * self.snr, 14.0, 70.0)
            self.scatter = self.ax.scatter(
                self.u,
                self.v,
                c=self.colors(),
                s=sizes,
                edgecolors="none",
            )
        beam_vecs = disamb.tx_beam_unit_vectors()
        for bid, (u_b, v_b, _w_b) in enumerate(beam_vecs):
            self.ax.scatter([u_b], [v_b], marker="x", s=55, color=f"C{bid}", linewidths=1.4, zorder=5)
            self.ax.text(u_b, v_b, f" {disamb.TX_BEAM_SHORT_NAMES[bid]}", fontsize=8, color=f"C{bid}")
        self._draw_cluster_lines()
        self.ax.set_xlim(-1.03, 1.03)
        self.ax.set_ylim(-1.03, 1.03)
        self.ax.set_aspect("equal")
        self.ax.grid(True, alpha=0.2)
        self.ax.set_xlabel("u/east direction cosine")
        self.ax.set_ylabel("v/north direction cosine")
        self.ax.set_title(f"Manual interferometric clusters: {self.args.sample_idx}")
        counts = {int(label): int(np.count_nonzero(self.labels == label)) for label in np.unique(self.labels) if label > 0}
        self.text.set_text(
            f"Current cluster: C{self.current_label:02d} | "
            "lasso=assign, right-click lasso unavailable | keys: 1-9 select, n new, 0 erase lasso, c clear current, u undo, s save, q quit\n"
            f"Candidates: {len(self.candidates)} | labeled: {int(np.count_nonzero(self.labels))} | counts: {counts}"
        )
        self.fig.canvas.draw_idle()

    def _draw_cluster_lines(self):
        for label in sorted(int(x) for x in np.unique(self.labels) if int(x) > 0):
            idx = np.flatnonzero(self.labels == label)
            if len(idx) < 2:
                continue
            order = idx[np.argsort(self.t[idx], kind="mergesort")]
            color = plt.get_cmap("tab20")((label - 1) % 20)
            self.ax.plot(self.u[order], self.v[order], "-", lw=1.4, color=color, alpha=0.85)
            mid = order[len(order) // 2]
            self.ax.text(
                self.u[mid],
                self.v[mid],
                f"C{label:02d}",
                fontsize=8,
                ha="center",
                va="center",
                bbox={"facecolor": "white", "edgecolor": color, "alpha": 0.85, "pad": 1.2},
            )

    def snapshot(self):
        self.undo_stack.append(self.labels.copy())
        if len(self.undo_stack) > 50:
            self.undo_stack.pop(0)

    def on_lasso(self, verts):
        if len(self.candidates) == 0:
            return
        path = MplPath(verts)
        selected = path.contains_points(self.points_uv)
        if not np.any(selected):
            return
        self.snapshot()
        self.labels[selected] = int(self.current_label)
        self.redraw()

    def on_key(self, event):
        if event.key is None:
            return
        key = event.key.lower()
        if key in {str(i) for i in range(1, 10)}:
            self.current_label = int(key)
            self.redraw()
        elif key == "0":
            self.current_label = 0
            self.redraw()
        elif key == "n":
            self.new_cluster(None)
        elif key == "c":
            self.clear_current(None)
        elif key == "u":
            self.undo(None)
        elif key == "s":
            self.save(None)
        elif key == "q":
            plt.close(self.fig)

    def new_cluster(self, _event):
        self.current_label = max([1, *[int(x) for x in np.unique(self.labels)]]) + 1
        self.redraw()

    def clear_current(self, _event):
        self.snapshot()
        self.labels[self.labels == self.current_label] = 0
        self.redraw()

    def undo(self, _event):
        if not self.undo_stack:
            return
        self.labels[:] = self.undo_stack.pop()
        self.redraw()

    def save(self, _event):
        save_labels(
            self.args.output_h5,
            self.args.sample_idx,
            self.args.grid_n,
            self.args.coherence_threshold,
            self.args.snr_threshold,
            self.candidates,
            self.labels,
        )
        png = self.args.output_h5.with_suffix(".png")
        self.fig.savefig(png, dpi=180)
        print(f"wrote {self.args.output_h5}")
        print(f"wrote {png}")


def main():
    parser = argparse.ArgumentParser(description="Manually label interferometric alias clusters in direction-cosine space.")
    parser.add_argument("--sample-idx", type=int, required=True)
    parser.add_argument("--cut-dir", type=Path, default=Path("data/metadata/cut"))
    parser.add_argument("--output-h5", type=Path, default=None)
    parser.add_argument("--grid-n", type=int, default=501)
    parser.add_argument("--coherence-threshold", type=float, default=0.9)
    parser.add_argument("--snr-threshold", type=float, default=9.0)
    parser.add_argument("--recompute-cut-observables", action="store_true")
    args = parser.parse_args()
    if args.output_h5 is None:
        args.output_h5 = Path("test_plots") / f"manual_interferometric_clusters_{args.sample_idx}.h5"

    _obs, candidates = build_candidates(
        args.sample_idx,
        args.cut_dir,
        args.grid_n,
        args.coherence_threshold,
        args.snr_threshold,
        args.recompute_cut_observables,
    )
    labels = load_labels(args.output_h5, len(candidates))
    gui = ClusterGui(args, candidates, labels)
    plt.show()
    return gui


if __name__ == "__main__":
    main()
