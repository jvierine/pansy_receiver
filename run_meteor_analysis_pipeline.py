#!/usr/bin/env python3
"""Run the standardized PANSY meteor alias/orbit pipeline on multiple cuts."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import itertools
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

import pansy_interferometry as pint
from interferometer_alias_diagnostics import load_cut, recompute_cut_observables
import plot_interferometric_disambiguation as disamb


@dataclass
class PipelineConfig:
    cut_dir: Path
    output_dir: Path
    table_dir: Path
    grid_n: int
    coherence_threshold: float
    snr_threshold: float
    linearity_angular_sigma_deg: float
    linearity_range_sigma_km: float
    linearity_p_threshold: float
    ballistic_p_threshold: float
    target_alt_km: float


STANDARD_PIPELINE_NAME = "interferometric_disambiguation_v1"


def sample_utc(sample_idx: int) -> str:
    return dt.datetime.fromtimestamp(sample_idx / 1e6, tz=dt.timezone.utc).isoformat(timespec="milliseconds")


def discover_cut_samples(cut_dir: Path) -> list[dict]:
    """Discover candidate meteor cut keys without scanning the whole metadata channel."""
    rows = []
    for path in sorted(cut_dir.glob("*/*.h5")):
        with h5py.File(path, "r") as handle:
            for key in handle.keys():
                group = handle[key]
                c_snr = np.asarray(group.get("c_snr", []), dtype=np.float64)
                tx_idx = np.asarray(group.get("tx_idx", []), dtype=np.float64)
                if len(c_snr) == 0 or len(tx_idx) == 0:
                    continue
                rows.append(
                    {
                        "sample_idx": int(key),
                        "file": str(path),
                        "n_cut_pulses": int(len(tx_idx)),
                        "n_cluster_points": int(len(c_snr)),
                        "max_cluster_snr_db": float(10.0 * np.log10(np.nanmax(np.maximum(c_snr, 1e-6)))),
                    }
                )
    rows.sort(key=lambda r: (r["max_cluster_snr_db"], r["n_cut_pulses"], r["n_cluster_points"]), reverse=True)
    return rows


def make_context(grid_n: int):
    antpos = pint.get_antpos()
    ch_pairs = np.asarray(list(itertools.combinations(np.arange(7), 2)))
    dmat = pint.pair_mat(ch_pairs, antpos)
    phasecal = disamb.pint.get_phasecal()
    u, v, w, valid = disamb.horizon_grid(grid_n)
    steering, _uvw = disamb.steering_matrix(dmat, u, v, w, valid)
    tx_gain_maps = disamb.precompute_tx_array_gain_maps(u, v, w, valid)
    return {
        "ch_pairs": ch_pairs,
        "phasecal": phasecal,
        "u": u,
        "v": v,
        "w": w,
        "valid": valid,
        "steering": steering,
        "tx_gain_maps": tx_gain_maps,
    }


def build_candidates(obs: dict, context: dict, coherence_threshold: float) -> tuple[list[dict], tuple[np.ndarray, np.ndarray], np.ndarray]:
    """Find all high-coherence interferometer aliases for one event."""
    u = context["u"]
    v = context["v"]
    valid = context["valid"]
    steering = context["steering"]
    phasecal = context["phasecal"]
    ch_pairs = context["ch_pairs"]

    strongest = int(np.argmax(obs["snr"]))
    single = disamb.coherence_map(
        obs["xc"][strongest], int(obs["beam_id"][strongest]), phasecal, ch_pairs, steering, valid, u.shape
    )
    peak_ij = disamb.local_coherence_peaks(single, coherence_threshold)

    candidates = []
    t_rel = obs["tx_idx"] / 1e6 - obs["tx_idx"][0] / 1e6
    for i in range(len(obs["snr"])):
        coh = disamb.coherence_map(obs["xc"][i], int(obs["beam_id"][i]), phasecal, ch_pairs, steering, valid, u.shape)
        ii, jj = disamb.local_coherence_peaks(coh, coherence_threshold)
        for row, col in zip(ii, jj):
            candidates.append(
                {
                    "u": float(u[row, col]),
                    "v": float(v[row, col]),
                    "w": float(context["w"][row, col]),
                    "grid_row": int(row),
                    "grid_col": int(col),
                    "coherence": float(coh[row, col]),
                    "t_rel": float(t_rel[i]),
                    "pulse": int(i),
                    "range_km": float(obs["range_km"][i]),
                    "doppler_mps": float(obs["doppler_mps"][i]),
                    "snr": float(obs["snr"][i]),
                    "beam_id": int(obs["beam_id"][i]),
                }
            )
    return candidates, peak_ij, single


def filter_observations_by_snr(obs_all: dict, snr_threshold: float) -> dict:
    """Keep pulse-wise cut observables above the standard SNR threshold."""
    good = obs_all["snr"] > snr_threshold
    return {
        key: val[good] if isinstance(val, np.ndarray) and len(val) == len(good) else val
        for key, val in obs_all.items()
    }


def run_standard_disambiguation(obs_all: dict, sample_idx: int, config: PipelineConfig, context: dict) -> dict:
    """Run the standard interferometric meteor-cut disambiguation pipeline."""
    obs = filter_observations_by_snr(obs_all, config.snr_threshold)
    if len(obs["snr"]) < 10:
        raise RuntimeError(f"only {len(obs['snr'])} pulses above SNR threshold")

    candidates, peak_ij, single = build_candidates(obs, context, config.coherence_threshold)
    tracks = disamb.fit_candidate_tracks(candidates)
    t_rel = obs["tx_idx"] / 1e6 - obs["tx_idx"][0] / 1e6
    t_start = float(np.min(t_rel))
    t_end = float(np.max(t_rel))
    tracks = [disamb.classify_track_visibility(track, candidates, t_start, t_end) for track in tracks]
    tracks = [
        disamb.classify_track_linearity(
            track,
            candidates,
            angular_sigma_deg=config.linearity_angular_sigma_deg,
            range_sigma_km=config.linearity_range_sigma_km,
            p_threshold=config.linearity_p_threshold,
        )
        if track["reason"] == "kept"
        else track
        for track in tracks
    ]
    tracks = [
        disamb.classify_track_descent(track)
        if track["reason"] == "kept"
        and not track.get("linearity_reject", False)
        and not track.get("low_detection_altitude_reject", False)
        else track
        for track in tracks
    ]
    sigma_pos, sigma_dop = disamb.fit_ballistic_survivors(
        tracks,
        candidates,
        event_epoch_unix=sample_idx / 1e6,
        p_threshold=config.ballistic_p_threshold,
    )
    disamb.score_tx_beam_consistency(tracks, candidates)
    disamb.score_candidate_diagnostics(tracks, candidates, tx_gain_maps=context["tx_gain_maps"])
    tx_sigma = disamb.score_combined_hypotheses(tracks)
    ranked = sorted([t for t in tracks if "combined_score" in t], key=lambda t: t["combined_rank"])
    if not ranked:
        raise RuntimeError("no combined-ranked hypotheses survived")

    for track in ranked[:3]:
        track["candidate_orbit"] = disamb.orbit_for_candidate_track(track, sample_idx / 1e6)

    return {
        "pipeline": STANDARD_PIPELINE_NAME,
        "obs": obs,
        "candidates": candidates,
        "peak_ij": peak_ij,
        "single": single,
        "tracks": tracks,
        "ranked": ranked,
        "sigma_pos": sigma_pos,
        "sigma_dop": sigma_dop,
        "tx_sigma": tx_sigma,
    }


def analyze_event(sample_idx: int, config: PipelineConfig, context: dict, make_event_plots: bool = True) -> dict:
    cut = load_cut(config.cut_dir, sample_idx)
    obs_all = recompute_cut_observables(cut)
    result = run_standard_disambiguation(obs_all, sample_idx, config, context)
    obs = result["obs"]
    candidates = result["candidates"]
    peak_ij = result["peak_ij"]
    single = result["single"]
    tracks = result["tracks"]
    ranked = result["ranked"]
    sigma_pos = result["sigma_pos"]
    sigma_dop = result["sigma_dop"]
    tx_sigma = result["tx_sigma"]

    event_dir = config.output_dir / str(sample_idx)
    if make_event_plots:
        event_dir.mkdir(parents=True, exist_ok=True)
        disamb.plot_single_echo(context["u"], context["v"], single, peak_ij, obs, int(np.argmax(obs["snr"])), event_dir / "single_echo_coherence.png")
        disamb.plot_all_candidates(context["u"], context["v"], obs, candidates, event_dir / "alias_candidates.png")
        disamb.plot_visibility_rejections(candidates, tracks, event_dir / "visibility_rejections.png")
        disamb.plot_linearity_rejections(candidates, tracks, event_dir / "linearity_rejections.png")
        disamb.plot_ballistic_ranking(tracks, event_dir / "ballistic_tx_ranking.png")
        disamb.plot_tx_beam_consistency(tracks, event_dir / "tx_beam_consistency.png")
        disamb.plot_tx_array_snr_consistency(tracks, event_dir / "tx_array_snr_consistency.png")

    best = ranked[0]
    second = ranked[1] if len(ranked) > 1 else None
    kep = best["candidate_orbit"]["kepler"]
    v_gcrs = np.linalg.norm(best["candidate_orbit"]["state_gcrs_m_mps"][3:]) / 1e3
    return {
        "sample_idx": int(sample_idx),
        "utc": sample_utc(sample_idx),
        "n_pulses": int(len(obs["snr"])),
        "max_snr_db": float(10.0 * np.log10(np.nanmax(np.maximum(obs["snr"], 1e-6)))),
        "single_echo_peaks": int(len(peak_ij[0])),
        "coherence_candidates": int(len(candidates)),
        "path_hypotheses": int(len(tracks)),
        "after_horizon": int(sum(t["reason"] == "kept" for t in tracks)),
        "after_linearity": int(sum(t["reason"] == "kept" and not t.get("linearity_reject", False) for t in tracks)),
        "after_descent": int(
            sum(
                t["reason"] == "kept" and not t.get("linearity_reject", False) and not t.get("descent_reject", False)
                for t in tracks
            )
        ),
        "ranked_hypotheses": int(len(ranked)),
        "best_score": float(best["combined_score"]),
        "second_score": float(second["combined_score"]) if second is not None else np.nan,
        "delta_score": float(second["combined_score"] - best["combined_score"]) if second is not None else np.inf,
        "odds_best_second": float(best.get("combined_odds_best_vs_second", np.inf)),
        "best_ballistic_chi2_nu": float(best["ballistic_reduced_chi2"]),
        "best_tx_rms_deg": float(best["tx_beam_weighted_rms_deg"]),
        "best_coherence_mean": float(best.get("coherence_weighted_mean", np.nan)),
        "best_coherence_p10": float(best.get("coherence_p10", np.nan)),
        "best_track_pulse_fraction": float(best.get("track_unique_pulse_fraction", np.nan)),
        "best_tx_array_snr_rms_db": float(best.get("tx_array_snr_rms_db", np.nan)),
        "best_tx_array_snr_corr": float(best.get("tx_array_snr_corr", np.nan)),
        "best_tx_array_gain_span_db": float(best.get("tx_array_gain_span_db", np.nan)),
        "sigma_pos_km": float(sigma_pos),
        "sigma_dop_km_s": float(sigma_dop),
        "tx_sigma_deg": float(tx_sigma),
        "v_gcrs_km_s": float(v_gcrs),
        "a_au": float(kep[0]),
        "e": float(kep[1]),
        "i_deg": float(kep[2]),
        "raan_deg": float(kep[3]),
        "argp_deg": float(kep[4]),
        "nu_deg": float(kep[5]),
        "q_au": float(kep[6]),
        "hyperbolic": bool(kep[1] > 1.0 or kep[0] < 0.0),
        "event_plot_dir": str(event_dir),
    }


def write_csv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_latex_table(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\\begin{tabular}{rrrrrrrrrrr}",
        "\\hline",
        "Event & Pulses & Paths & Survivors & $\\chi^2_{\\nu,0}$ & $\\Delta\\chi^2_\\nu$ & $\\bar{C}$ & $v_{\\mathrm{GCRS}}$ & $a$ & $e$ & $i$ \\\\",
        " & & & & & & & (\\si{km.s^{-1}}) & (AU) & & (deg) \\\\",
        "\\hline",
    ]
    for i, row in enumerate(rows, start=1):
        lines.append(
            f"{i} & {row['n_pulses']} & {row['path_hypotheses']} & {row['after_descent']} & "
            f"{row['best_score']:.2f} & {row['delta_score']:.2f} & {row['best_coherence_mean']:.3f} & {row['v_gcrs_km_s']:.2f} & "
            f"{row['a_au']:.3g} & {row['e']:.3f} & {row['i_deg']:.1f} \\\\"
        )
    lines.extend(["\\hline", "\\end{tabular}", ""])
    path.write_text("\n".join(lines))


def plot_summary(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    idx = np.arange(1, len(rows) + 1)
    delta = np.asarray([r["delta_score"] for r in rows], dtype=np.float64)
    speed = np.asarray([r["v_gcrs_km_s"] for r in rows], dtype=np.float64)
    ecc = np.asarray([r["e"] for r in rows], dtype=np.float64)
    after = np.asarray([r["after_descent"] for r in rows], dtype=np.float64)
    ranked = np.asarray([r["ranked_hypotheses"] for r in rows], dtype=np.float64)

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2), constrained_layout=True)
    ax = axes[0, 0]
    ax.bar(idx, ranked, color="0.75", label="ranked")
    ax.bar(idx, after, color="tab:blue", label="after descent")
    ax.set_xlabel("Event")
    ax.set_ylabel("Number of hypotheses")
    ax.set_title("Hypotheses reaching final ranking")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.bar(idx, delta, color="tab:green")
    ax.set_xlabel("Event")
    ax.set_ylabel("$\\Delta\\chi^2_\\nu$ to second-best")
    ax.set_title("Ballistic score separation")
    ax.grid(True, axis="y", alpha=0.25)

    ax = axes[1, 0]
    colors = ["tab:red" if r["hyperbolic"] else "tab:blue" for r in rows]
    ax.scatter(idx, speed, c=colors, s=60)
    ax.set_xlabel("Event")
    ax.set_ylabel("GCRS speed at 140 km (km/s)")
    ax.set_title("Selected state speed")
    ax.grid(True, alpha=0.25)

    ax = axes[1, 1]
    ax.scatter([r["q_au"] for r in rows], ecc, c=colors, s=70)
    ax.axhline(1.0, color="black", lw=1.0, ls="--")
    ax.set_xlabel("Perihelion distance q (AU)")
    ax.set_ylabel("Eccentricity")
    ax.set_title("Preliminary selected-orbit elements")
    ax.grid(True, alpha=0.25)

    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run PANSY meteor analysis pipeline on multiple cut events.")
    parser.add_argument("--cut-dir", type=Path, default=Path("data/metadata/cut"))
    parser.add_argument("--output-dir", type=Path, default=Path("../pansy_paper/memos/figures/pansy_meteor_pipeline"))
    parser.add_argument("--table-dir", type=Path, default=Path("../pansy_paper/memos/tables"))
    parser.add_argument("--n-events", type=int, default=10)
    parser.add_argument("--grid-n", type=int, default=401)
    parser.add_argument("--coherence-threshold", type=float, default=0.9)
    parser.add_argument("--snr-threshold", type=float, default=9.0)
    parser.add_argument("--linearity-angular-sigma-deg", type=float, default=0.25)
    parser.add_argument("--linearity-range-sigma-km", type=float, default=0.15)
    parser.add_argument("--linearity-p-threshold", type=float, default=0.01)
    parser.add_argument("--ballistic-p-threshold", type=float, default=0.01)
    parser.add_argument("--target-alt-km", type=float, default=140.0)
    parser.add_argument("--sample-idx", type=int, nargs="*", help="Explicit sample indices to process.")
    parser.add_argument("--no-event-plots", action="store_true")
    args = parser.parse_args()

    config = PipelineConfig(
        cut_dir=args.cut_dir,
        output_dir=args.output_dir,
        table_dir=args.table_dir,
        grid_n=args.grid_n,
        coherence_threshold=args.coherence_threshold,
        snr_threshold=args.snr_threshold,
        linearity_angular_sigma_deg=args.linearity_angular_sigma_deg,
        linearity_range_sigma_km=args.linearity_range_sigma_km,
        linearity_p_threshold=args.linearity_p_threshold,
        ballistic_p_threshold=args.ballistic_p_threshold,
        target_alt_km=args.target_alt_km,
    )

    samples = args.sample_idx
    if not samples:
        discovered = discover_cut_samples(args.cut_dir)
        samples = [row["sample_idx"] for row in discovered[: args.n_events]]

    context = make_context(args.grid_n)
    rows = []
    for sample_idx in samples:
        print(f"processing {sample_idx} {sample_utc(sample_idx)}")
        try:
            row = analyze_event(sample_idx, config, context, make_event_plots=not args.no_event_plots)
            rows.append(row)
            print(
                f"  score={row['best_score']:.2f} delta={row['delta_score']:.2f} "
                f"e={row['e']:.3f} a={row['a_au']:.3g} speed={row['v_gcrs_km_s']:.2f} km/s"
            )
        except Exception as exc:
            print(f"  failed: {exc}")

    if not rows:
        raise RuntimeError("no events completed")

    csv_path = args.table_dir / "pansy_meteor_pipeline_10.csv"
    tex_path = args.table_dir / "pansy_meteor_pipeline_10.tex"
    summary_path = args.output_dir / "pansy_meteor_pipeline_10_summary.png"
    write_csv(rows, csv_path)
    write_latex_table(rows, tex_path)
    plot_summary(rows, summary_path)
    print(f"wrote {csv_path}")
    print(f"wrote {tex_path}")
    print(f"wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
