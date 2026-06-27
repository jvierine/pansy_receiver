#!/usr/bin/env python3
"""Create memo figures explaining PANSY interferometric hypothesis formation."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from scipy import ndimage as ndi

import pansy_interferometry as pint
from interferometer_alias_diagnostics import load_cut, recompute_cut_observables
import plot_interferometric_disambiguation as disamb


def filtered_observations(obs_all: dict, snr_threshold: float) -> dict:
    good = np.asarray(obs_all["snr"]) > float(snr_threshold)
    return {
        key: val[good] if isinstance(val, np.ndarray) and len(val) == len(good) else val
        for key, val in obs_all.items()
    }


def make_interferometer_context(grid_n: int) -> dict:
    antpos = pint.get_antpos()
    ch_pairs = np.asarray(list(itertools.combinations(np.arange(7), 2)))
    dmat = pint.pair_mat(ch_pairs, antpos)
    phasecal = disamb.pint.get_phasecal()
    u, v, w, valid = disamb.horizon_grid(grid_n)
    steering, _uvw = disamb.steering_matrix(dmat, u, v, w, valid)
    return {
        "ch_pairs": ch_pairs,
        "phasecal": phasecal,
        "u": u,
        "v": v,
        "w": w,
        "valid": valid,
        "steering": steering,
    }


def load_diagnostics(path: Path) -> tuple[dict[str, np.ndarray], list[dict], dict]:
    with h5py.File(path, "r") as handle:
        cand_group = handle["candidates"]
        candidates = {key: np.asarray(cand_group[key]) for key in cand_group.keys()}
        tracks = []
        for label in sorted(handle["hypotheses"].keys()):
            group = handle["hypotheses"][label]
            if "candidate_indices" not in group:
                continue
            idx = np.asarray(group["candidate_indices"], dtype=np.int64)
            seed_idx = np.asarray(group["seed_candidate_indices"], dtype=np.int64) if "seed_candidate_indices" in group else idx
            tracks.append(
                {
                    "label": label,
                    "idx": idx,
                    "seed_idx": seed_idx,
                    "attrs": dict(group.attrs),
                }
            )
        attrs = dict(handle.attrs)
    return candidates, tracks, attrs


def candidate_dicts(candidates: dict[str, np.ndarray]) -> list[dict]:
    n = len(candidates["u"])
    rows = []
    for i in range(n):
        rows.append(
            {
                "u": float(candidates["u"][i]),
                "v": float(candidates["v"][i]),
                "pulse": int(candidates["pulse"][i]),
                "t_rel": float(candidates["t_rel_s"][i]),
                "coherence": float(candidates["coherence"][i]),
            }
        )
    return rows


def local_maxima_2d(image: np.ndarray, threshold: float, max_peaks: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    filled = np.nan_to_num(image, nan=-np.inf)
    max_img = ndi.maximum_filter(filled, size=7, mode="constant", cval=-np.inf)
    mask = (filled == max_img) & (filled >= threshold)
    ii, jj = np.where(mask)
    if max_peaks is not None and len(ii) > max_peaks:
        order = np.argsort(filled[ii, jj])[::-1][:max_peaks]
        ii = ii[order]
        jj = jj[order]
    return ii, jj


def peak_support(
    peak_u: np.ndarray,
    peak_v: np.ndarray,
    candidates: dict[str, np.ndarray],
    radius: float,
    used_candidate_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cand_u = np.asarray(candidates["u"], dtype=np.float64)
    cand_v = np.asarray(candidates["v"], dtype=np.float64)
    cand_pulse = np.asarray(candidates["pulse"], dtype=np.int64)
    support = np.zeros(len(peak_u), dtype=np.int64)
    used_support = np.zeros(len(peak_u), dtype=np.int64)
    nearest_used = np.full(len(peak_u), np.nan, dtype=np.float64)
    used_u = cand_u[used_candidate_mask]
    used_v = cand_v[used_candidate_mask]
    for i, (uu, vv) in enumerate(zip(peak_u, peak_v)):
        dist = np.hypot(cand_u - uu, cand_v - vv)
        near = dist <= radius
        support[i] = len(np.unique(cand_pulse[near]))
        used_near = near & used_candidate_mask
        used_support[i] = len(np.unique(cand_pulse[used_near]))
        if len(used_u):
            nearest_used[i] = float(np.min(np.hypot(used_u - uu, used_v - vv)))
    return support, used_support, nearest_used


def draw_horizon(ax):
    ax.add_patch(plt.Circle((0.0, 0.0), 1.0, fill=False, color="0.15", linewidth=0.8))
    ax.set_aspect("equal")
    ax.set_xlim(-1.03, 1.03)
    ax.set_ylim(-1.03, 1.03)
    ax.set_xlabel("u/east direction cosine")
    ax.set_ylabel("v/north direction cosine")
    ax.grid(True, color="0.0", alpha=0.08, linewidth=0.5)


def plot_hypothesis_formation(
    output_png: Path,
    audit_h5: Path,
    sample_idx: int,
    u_grid: np.ndarray,
    v_grid: np.ndarray,
    summed: np.ndarray,
    candidates: dict[str, np.ndarray],
    tracks: list[dict],
    peak_i: np.ndarray,
    peak_j: np.ndarray,
    support: np.ndarray,
    used_support: np.ndarray,
    nearest_used: np.ndarray,
    min_unique_pulses: int,
    support_radius: float,
    coherence_threshold: float,
):
    output_png.parent.mkdir(parents=True, exist_ok=True)
    audit_h5.parent.mkdir(parents=True, exist_ok=True)

    cand_u = np.asarray(candidates["u"], dtype=np.float64)
    cand_v = np.asarray(candidates["v"], dtype=np.float64)
    cand_coh = np.asarray(candidates["coherence"], dtype=np.float64)
    cand_pulse = np.asarray(candidates["pulse"], dtype=np.int64)
    peak_u = u_grid[peak_i, peak_j]
    peak_v = v_grid[peak_i, peak_j]
    peak_value = summed[peak_i, peak_j]

    used_mask = np.zeros(len(cand_u), dtype=bool)
    seed_mask = np.zeros(len(cand_u), dtype=bool)
    for track in tracks:
        used_mask[np.asarray(track["idx"], dtype=np.int64)] = True
        seed_mask[np.asarray(track["seed_idx"], dtype=np.int64)] = True

    with h5py.File(audit_h5, "w") as handle:
        handle.attrs["schema_version"] = "pansy_hypothesis_formation_audit_v1"
        handle.attrs["source_program"] = "plot_hypothesis_formation_memo.py"
        handle.attrs["sample_idx"] = int(sample_idx)
        handle.attrs["coherence_threshold"] = float(coherence_threshold)
        handle.attrs["support_radius_direction_cosine"] = float(support_radius)
        handle.attrs["min_unique_pulses"] = int(min_unique_pulses)
        handle.create_dataset("candidate_u", data=cand_u)
        handle.create_dataset("candidate_v", data=cand_v)
        handle.create_dataset("candidate_pulse", data=cand_pulse)
        handle.create_dataset("candidate_coherence", data=cand_coh)
        handle.create_dataset("candidate_used_by_completed_track", data=used_mask.astype(np.uint8))
        handle.create_dataset("candidate_used_by_seed_track", data=seed_mask.astype(np.uint8))
        peaks = handle.create_group("summed_response_peaks")
        peaks.create_dataset("u", data=peak_u)
        peaks.create_dataset("v", data=peak_v)
        peaks.create_dataset("normalized_response", data=peak_value)
        peaks.create_dataset("support_unique_pulses", data=support)
        peaks.create_dataset("used_support_unique_pulses", data=used_support)
        peaks.create_dataset("nearest_used_candidate_distance", data=nearest_used)

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 6.0), constrained_layout=True)
    cmap = plt.get_cmap("Blues").copy()
    cmap.set_bad("white")

    ax = axes[0]
    ax.pcolormesh(u_grid, v_grid, summed, shading="auto", cmap=cmap, vmin=0.0, vmax=1.0)
    ax.scatter(cand_u, cand_v, s=2.0, c="black", alpha=0.34, linewidths=0, label="per-pulse maxima")
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, max(1, len(tracks))))
    for color, track in zip(colors, tracks):
        idx = np.asarray(track["idx"], dtype=np.int64)
        if len(idx) == 0:
            continue
        x = cand_u[idx]
        y = cand_v[idx]
        pad = 0.018
        x0 = float(np.nanmin(x) - pad)
        y0 = float(np.nanmin(y) - pad)
        width = float(np.nanmax(x) - np.nanmin(x) + 2 * pad)
        height = float(np.nanmax(y) - np.nanmin(y) + 2 * pad)
        ax.add_patch(Rectangle((x0, y0), width, height, fill=False, edgecolor=color, linewidth=1.1))
        ax.text(
            x0,
            y0 + height,
            track["label"],
            ha="left",
            va="bottom",
            fontsize=7,
            color="black",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.3, "pad": 0.8},
        )
    draw_horizon(ax)
    ax.set_title("A. Per-pulse maxima grouped into track hypotheses")
    ax.text(
        0.02,
        0.02,
        f"{len(cand_u)} candidates, {len(np.unique(cand_pulse))} pulses, {len(tracks)} tracks",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72},
    )

    ax = axes[1]
    ax.pcolormesh(u_grid, v_grid, summed, shading="auto", cmap=cmap, vmin=0.0, vmax=1.0)
    ok = support >= min_unique_pulses
    sc = ax.scatter(
        peak_u[~ok],
        peak_v[~ok],
        c=np.maximum(support[~ok], 1),
        s=26,
        cmap="plasma",
        vmin=1,
        vmax=max(min_unique_pulses, int(np.nanmax(support)) if len(support) else 1),
        marker="x",
        linewidths=0.9,
        label="below pulse support",
    )
    ax.scatter(
        peak_u[ok],
        peak_v[ok],
        c=np.maximum(support[ok], 1),
        s=46,
        cmap="plasma",
        vmin=1,
        vmax=max(min_unique_pulses, int(np.nanmax(support)) if len(support) else 1),
        edgecolors="black",
        linewidths=0.5,
        label="enough pulse support",
    )
    draw_horizon(ax)
    ax.set_title("B. Summed-response peaks need pulse support")
    cb = fig.colorbar(sc, ax=ax, shrink=0.84)
    cb.set_label("nearby unique pulses")
    ax.legend(loc="lower left", fontsize=8, frameon=True, framealpha=0.75)
    ax.text(
        0.02,
        0.98,
        f"support radius={support_radius:.2f}, required pulses={min_unique_pulses}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72},
    )

    fig.suptitle(f"Interferometric hypothesis formation audit, sample {sample_idx}", fontsize=12)
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-idx", type=int, default=1746490036904845)
    parser.add_argument("--cut-dir", type=Path, default=Path("data/metadata/cut"))
    parser.add_argument("--diagnostics-h5", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=Path("../pansy_paper/memos/figures"))
    parser.add_argument("--grid-n", type=int, default=501)
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    parser.add_argument("--coherence-threshold", type=float, default=0.80)
    parser.add_argument("--summed-peak-threshold", type=float, default=0.85)
    parser.add_argument("--summed-max-peaks", type=int, default=80)
    parser.add_argument("--support-radius", type=float, default=0.04)
    args = parser.parse_args()

    diagnostics_h5 = args.diagnostics_h5
    if diagnostics_h5 is None:
        diagnostics_h5 = Path("test_plots") / f"pansy_disambiguation_diagnostics_{args.sample_idx}.h5"

    candidates, tracks, _attrs = load_diagnostics(diagnostics_h5)
    obs = filtered_observations(recompute_cut_observables(load_cut(args.cut_dir, args.sample_idx)), args.snr_threshold)
    context = make_interferometer_context(args.grid_n)
    rows = candidate_dicts(candidates)
    summed = disamb.selected_pulse_interferometer_response(
        rows,
        obs,
        context["phasecal"],
        context["ch_pairs"],
        context["steering"],
        context["valid"],
        context["u"].shape,
    )
    if summed is None:
        raise RuntimeError("could not form summed receiver interferometer response")

    peak_i, peak_j = local_maxima_2d(summed, args.summed_peak_threshold, max_peaks=args.summed_max_peaks)
    used_mask = np.zeros(len(candidates["u"]), dtype=bool)
    for track in tracks:
        used_mask[np.asarray(track["idx"], dtype=np.int64)] = True
    support, used_support, nearest_used = peak_support(
        context["u"][peak_i, peak_j],
        context["v"][peak_i, peak_j],
        candidates,
        args.support_radius,
        used_mask,
    )

    total_unique_pulses = len(np.unique(candidates["pulse"]))
    min_unique_pulses = int(min(total_unique_pulses, 10, max(5, np.ceil(0.75 * total_unique_pulses))))
    output_png = args.output_dir / f"hypothesis_formation_audit_{args.sample_idx}.png"
    audit_h5 = args.output_dir / f"hypothesis_formation_audit_{args.sample_idx}.h5"
    plot_hypothesis_formation(
        output_png,
        audit_h5,
        args.sample_idx,
        context["u"],
        context["v"],
        summed,
        candidates,
        tracks,
        peak_i,
        peak_j,
        support,
        used_support,
        nearest_used,
        min_unique_pulses,
        args.support_radius,
        args.coherence_threshold,
    )
    print(output_png)
    print(audit_h5)
    print(f"candidates {len(candidates['u'])}")
    print(f"tracks {len(tracks)}")
    print(f"summed_response_peaks {len(peak_i)}")
    print(f"summed_response_peaks_with_support_ge_{min_unique_pulses} {int(np.sum(support >= min_unique_pulses))}")


if __name__ == "__main__":
    main()
