#!/usr/bin/env python3
"""Fast GUI for annotating issues in PANSY interferometric event plots."""

from __future__ import annotations

import argparse
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import numpy as np


ISSUE_LABELS = [
    "interferometric hypothesis clustering",
    "too few RTI points selected despite clear trace",
    "multiple simultaneous meteors not separated",
    "wrong hypothesis selected",
    "bad range fit",
    "bad doppler fit",
]


def sample_idx_from_path(path: Path) -> int:
    return int(path.stem.split("_")[-1])


def discover_plots(plot_dir: Path, pattern: str):
    paths = sorted(plot_dir.glob(pattern), key=sample_idx_from_path)
    if not paths:
        raise FileNotFoundError(f"No event plot PNGs matching {plot_dir / pattern}")
    return paths


def load_annotations(path: Path, plots: list[Path]):
    sample_idx = np.asarray([sample_idx_from_path(p) for p in plots], dtype=np.int64)
    rel_paths = np.asarray([str(p) for p in plots], dtype=object)
    issue_flags = np.zeros((len(plots), len(ISSUE_LABELS)), dtype=np.bool_)
    reviewed = np.zeros(len(plots), dtype=np.bool_)
    current_index = 0
    if not path.exists():
        return sample_idx, rel_paths, issue_flags, reviewed, current_index
    with h5py.File(path, "r") as h:
        old_sample = h["sample_idx"][()].astype(np.int64) if "sample_idx" in h else np.empty(0, dtype=np.int64)
        old_flags = h["issue_flags"][()].astype(np.bool_) if "issue_flags" in h else np.zeros((0, len(ISSUE_LABELS)), dtype=np.bool_)
        old_reviewed = h["reviewed"][()].astype(np.bool_) if "reviewed" in h else np.zeros(len(old_sample), dtype=np.bool_)
        current_index = int(h.attrs.get("current_index", 0))
    lookup = {int(s): i for i, s in enumerate(old_sample)}
    for i, sample in enumerate(sample_idx):
        old_i = lookup.get(int(sample))
        if old_i is None:
            continue
        n_issue = min(issue_flags.shape[1], old_flags.shape[1])
        issue_flags[i, :n_issue] = old_flags[old_i, :n_issue]
        if old_i < len(old_reviewed):
            reviewed[i] = old_reviewed[old_i]
    current_index = int(np.clip(current_index, 0, len(sample_idx) - 1))
    return sample_idx, rel_paths, issue_flags, reviewed, current_index


def write_annotations(path: Path, sample_idx, plot_paths, issue_flags, reviewed, current_index):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(tmp, "w") as h:
        h.attrs["schema"] = "pansy_event_plot_issue_annotations_v1"
        h.attrs["updated_utc"] = datetime.now(timezone.utc).isoformat()
        h.attrs["current_index"] = int(current_index)
        h.create_dataset("sample_idx", data=np.asarray(sample_idx, dtype=np.int64))
        h.create_dataset("plot_path", data=np.asarray(plot_paths, dtype=object), dtype=string_dtype)
        h.create_dataset("issue_flags", data=np.asarray(issue_flags, dtype=np.bool_))
        h.create_dataset("reviewed", data=np.asarray(reviewed, dtype=np.bool_))
        h.create_dataset("issue_ids", data=np.arange(1, len(ISSUE_LABELS) + 1, dtype=np.int16))
        h.create_dataset("issue_labels", data=np.asarray(ISSUE_LABELS, dtype=object), dtype=string_dtype)
    tmp.replace(path)


class ImageCache:
    def __init__(self, max_images: int):
        self.max_images = max(1, int(max_images))
        self.cache: OrderedDict[Path, np.ndarray] = OrderedDict()

    def get(self, path: Path):
        if path in self.cache:
            image = self.cache.pop(path)
            self.cache[path] = image
            return image
        image = mpimg.imread(path)
        self.cache[path] = image
        while len(self.cache) > self.max_images:
            self.cache.popitem(last=False)
        return image


class EventIssueGui:
    def __init__(self, args):
        self.args = args
        self.plots = discover_plots(args.plot_dir, args.pattern)
        self.sample_idx, self.plot_paths, self.issue_flags, self.reviewed, self.index = load_annotations(args.output_h5, self.plots)
        self.cache = ImageCache(args.cache_images)

        self.fig, self.ax = plt.subplots(figsize=(14.0, 9.0))
        self.fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.12)
        self.ax.set_axis_off()
        self.image_artist = None
        self.status_text = self.fig.text(0.01, 0.055, "", ha="left", va="bottom", family="monospace", fontsize=9)
        self.issue_text = self.fig.text(0.99, 0.055, "", ha="right", va="bottom", family="monospace", fontsize=9)
        self._add_buttons()
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.save()
        self.show_current()

    def _add_buttons(self):
        buttons = [
            ("Prev", 0.36, self.prev),
            ("Next", 0.43, self.next),
            ("Clear", 0.50, self.clear_current),
            ("Reviewed", 0.58, self.toggle_reviewed),
            ("Unreviewed", 0.68, self.next_unreviewed),
            ("Save", 0.80, self.save),
        ]
        for label, left, callback in buttons:
            ax_button = self.fig.add_axes([left, 0.02, 0.065, 0.04])
            button = Button(ax_button, label)
            button.on_clicked(callback)
            setattr(self, f"_button_{label.lower()}", button)

    def current_path(self) -> Path:
        return self.plots[self.index]

    def current_issue_summary(self) -> str:
        active = [f"{i + 1}" for i, flag in enumerate(self.issue_flags[self.index]) if flag]
        return ",".join(active) if active else "none"

    def current_issue_text(self) -> str:
        lines = []
        for i, label in enumerate(ISSUE_LABELS):
            marker = "x" if self.issue_flags[self.index, i] else " "
            lines.append(f"[{marker}] {i + 1}: {label}")
        return "\n".join(lines)

    def show_current(self):
        path = self.current_path()
        image = self.cache.get(path)
        if self.image_artist is None:
            self.image_artist = self.ax.imshow(image)
        else:
            self.image_artist.set_data(image)
        self.ax.set_axis_off()
        n_reviewed = int(np.count_nonzero(self.reviewed))
        n_flagged = int(np.count_nonzero(np.any(self.issue_flags, axis=1)))
        sample = int(self.sample_idx[self.index])
        self.fig.suptitle(
            f"{self.index + 1}/{len(self.plots)}  sample {sample}  reviewed={bool(self.reviewed[self.index])}  issues={self.current_issue_summary()}",
            fontsize=12,
        )
        self.status_text.set_text(
            "keys: left/right browse | 1-6 toggle issue | 0 clear | r reviewed | u next unreviewed | q quit\n"
            f"reviewed {n_reviewed}/{len(self.plots)} | flagged {n_flagged} | autosave: {self.args.output_h5}"
        )
        self.issue_text.set_text(self.current_issue_text())
        self.fig.canvas.draw_idle()
        self.prefetch_neighbors()

    def prefetch_neighbors(self):
        for offset in (-1, 1):
            j = self.index + offset
            if 0 <= j < len(self.plots):
                try:
                    self.cache.get(self.plots[j])
                except OSError:
                    pass

    def save(self, _event=None):
        write_annotations(
            self.args.output_h5,
            self.sample_idx,
            self.plot_paths,
            self.issue_flags,
            self.reviewed,
            self.index,
        )

    def mark_reviewed(self):
        self.reviewed[self.index] = True

    def goto(self, index: int):
        self.index = int(np.clip(index, 0, len(self.plots) - 1))
        self.save()
        self.show_current()

    def next(self, _event=None):
        self.goto(self.index + 1)

    def prev(self, _event=None):
        self.goto(self.index - 1)

    def next_unreviewed(self, _event=None):
        if not np.any(~self.reviewed):
            return
        for step in range(1, len(self.reviewed) + 1):
            j = (self.index + step) % len(self.reviewed)
            if not self.reviewed[j]:
                self.goto(j)
                return

    def toggle_issue(self, issue_id: int):
        col = issue_id - 1
        if not (0 <= col < self.issue_flags.shape[1]):
            return
        self.issue_flags[self.index, col] = ~self.issue_flags[self.index, col]
        self.mark_reviewed()
        self.save()
        self.show_current()

    def clear_current(self, _event=None):
        self.issue_flags[self.index, :] = False
        self.mark_reviewed()
        self.save()
        self.show_current()

    def toggle_reviewed(self, _event=None):
        self.reviewed[self.index] = ~self.reviewed[self.index]
        self.save()
        self.show_current()

    def on_key(self, event):
        key = (event.key or "").lower()
        if key in {"right", "d", "n", " "}:
            self.next()
        elif key in {"left", "a", "p", "backspace"}:
            self.prev()
        elif key in {str(i) for i in range(1, len(ISSUE_LABELS) + 1)}:
            self.toggle_issue(int(key))
        elif key in {"0", "delete", "c"}:
            self.clear_current()
        elif key == "r":
            self.toggle_reviewed()
        elif key == "u":
            self.next_unreviewed()
        elif key == "home":
            self.goto(0)
        elif key == "end":
            self.goto(len(self.plots) - 1)
        elif key in {"q", "escape"}:
            self.save()
            plt.close(self.fig)


def main():
    parser = argparse.ArgumentParser(description="Annotate PANSY event plot issues and autosave to HDF5.")
    parser.add_argument("--plot-dir", type=Path, default=Path("test_plots"))
    parser.add_argument("--pattern", default="pansy_interferometer_disambiguation_summary_*.png")
    parser.add_argument("--output-h5", type=Path, default=Path("test_plots/event_plot_issue_annotations.h5"))
    parser.add_argument("--cache-images", type=int, default=5)
    args = parser.parse_args()

    gui = EventIssueGui(args)
    plt.show()
    gui.save()


if __name__ == "__main__":
    main()
