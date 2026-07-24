#!/usr/bin/env python3
"""Plot radiant uncertainty distributions from compact PANSY orbit metadata."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_ORBIT_METADATA_DIR = Path("/mnt/data/juha/pansy/metadata/orbit")
DEFAULT_OUTPUT_DIR = Path("test_plots/radiant_uncertainty")
PERCENTILES = np.asarray([1, 5, 10, 25, 50, 75, 90, 95, 99], dtype=np.float64)
ANGLE_THRESHOLDS_DEG = np.asarray([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0], dtype=np.float64)
UNCERTAINTY_DTYPE = np.dtype(
    [
        ("combined_score", "f8"),
        ("initial_detection_speed_km_s", "f8"),
        ("n_uncertainty_samples", "i8"),
        ("initial_state_position_sigma_m", "f8"),
        ("initial_state_velocity_sigma_mps", "f8"),
        ("initial_state_radiant_angle_sigma_deg", "f8"),
    ]
)


def finite_quality_rows(rows: np.ndarray, max_combined_score: float, min_uncertainty_samples: int) -> np.ndarray:
    good = np.isfinite(rows["initial_state_radiant_angle_sigma_deg"])
    good &= np.isfinite(rows["initial_state_velocity_sigma_mps"])
    good &= np.isfinite(rows["initial_state_position_sigma_m"])
    good &= np.isfinite(rows["initial_detection_speed_km_s"]) & (rows["initial_detection_speed_km_s"] > 0.0)
    good &= np.isfinite(rows["combined_score"]) & (rows["combined_score"] <= max_combined_score)
    good &= rows["n_uncertainty_samples"] >= int(min_uncertainty_samples)
    return rows[good]


def threshold_counts(values: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    return np.asarray([np.sum(values <= threshold) for threshold in thresholds], dtype=np.int64)


def rows_from_event_table(events: np.ndarray) -> np.ndarray:
    rows = np.zeros(len(events), dtype=UNCERTAINTY_DTYPE)
    rows["combined_score"] = np.asarray(events["combined_score"], dtype=np.float64)
    rows["initial_detection_speed_km_s"] = np.linalg.norm(
        np.asarray(events["fit_parameters"][:, 3:6], dtype=np.float64),
        axis=1,
    ) / 1e3
    rows["n_uncertainty_samples"] = np.asarray(events["n_uncertainty_samples"], dtype=np.int64)
    cov = np.asarray(events["fit_parameter_covariance"][:, :6, :6], dtype=np.float64)
    pos_var = np.trace(cov[:, :3, :3], axis1=1, axis2=2)
    vel_var = np.trace(cov[:, 3:6, 3:6], axis1=1, axis2=2)
    rows["initial_state_position_sigma_m"] = np.sqrt(np.where(pos_var >= 0.0, pos_var, np.nan))
    rows["initial_state_velocity_sigma_mps"] = np.sqrt(np.where(vel_var >= 0.0, vel_var, np.nan))
    rows["initial_state_radiant_angle_sigma_deg"] = np.rad2deg(
        np.arctan2(
            rows["initial_state_velocity_sigma_mps"],
            np.maximum(rows["initial_detection_speed_km_s"] * 1e3, 1.0),
        )
    )
    return rows


def rows_from_radiant_table(radiants: np.ndarray) -> np.ndarray:
    """Coerce a radiant sidecar that explicitly stores the local initial speed."""
    required = {
        "combined_score",
        "initial_detection_speed_km_s",
        "n_uncertainty_samples",
        "initial_state_position_sigma_m",
        "initial_state_velocity_sigma_mps",
        "initial_state_radiant_angle_sigma_deg",
    }
    names = set(radiants.dtype.names or ())
    missing = sorted(required - names)
    if missing:
        raise ValueError(
            "radiant sidecar predates the v0 uncertainty definition; rebuild it. "
            f"Missing fields: {', '.join(missing)}"
        )
    rows = np.zeros(len(radiants), dtype=UNCERTAINTY_DTYPE)
    for name in UNCERTAINTY_DTYPE.names:
        rows[name] = radiants[name]
    return rows


def collect_uncertainty_rows(orbit_metadata_dir: Path) -> tuple[np.ndarray, int, int]:
    chunks = []
    files_read = 0
    files_skipped = 0
    for path in sorted(orbit_metadata_dir.glob("**/orbit@*.h5")):
        try:
            with h5py.File(path, "r") as h:
                if "events" not in h:
                    files_skipped += 1
                    continue
                chunks.append(rows_from_event_table(h["events"][()]))
                files_read += 1
        except OSError:
            files_skipped += 1
    if chunks:
        return np.concatenate(chunks), files_read, files_skipped
    return np.zeros(0, dtype=UNCERTAINTY_DTYPE), files_read, files_skipped


def write_summary_h5(path: Path, rows: np.ndarray, files_read: int, files_skipped: int, source: str) -> None:
    angle = rows["initial_state_radiant_angle_sigma_deg"]
    vel = rows["initial_state_velocity_sigma_mps"]
    pos = rows["initial_state_position_sigma_m"]
    speed = rows["initial_detection_speed_km_s"]
    fractional_velocity_sigma = vel / np.maximum(speed * 1e3, 1.0)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h:
        h.attrs["script"] = Path(__file__).name
        h.attrs["source"] = source
        h.attrs["files_read"] = int(files_read)
        h.attrs["files_skipped"] = int(files_skipped)
        h.attrs["n_radiants"] = int(len(rows))
        h.attrs["radiant_angle_sigma_definition"] = (
            "atan(sqrt(trace(local initial-state velocity covariance)) / local initial-detection speed), degrees"
        )
        h.attrs["velocity_sigma_definition"] = "sqrt(trace(local initial-state velocity covariance)), m/s"
        h.attrs["fractional_velocity_sigma_definition"] = (
            "sqrt(trace(local initial-state velocity covariance)) / local initial-detection speed"
        )
        h.create_dataset("percentiles", data=PERCENTILES)
        h.create_dataset("radiant_angle_sigma_deg_percentiles", data=np.nanpercentile(angle, PERCENTILES))
        h.create_dataset("initial_state_velocity_sigma_mps_percentiles", data=np.nanpercentile(vel, PERCENTILES))
        h.create_dataset("initial_state_position_sigma_m_percentiles", data=np.nanpercentile(pos, PERCENTILES))
        h.create_dataset("initial_detection_speed_km_s_percentiles", data=np.nanpercentile(speed, PERCENTILES))
        h.create_dataset("angle_thresholds_deg", data=ANGLE_THRESHOLDS_DEG)
        h.create_dataset("angle_threshold_counts", data=threshold_counts(angle, ANGLE_THRESHOLDS_DEG))
        h.create_dataset("initial_state_velocity_sigma_mps", data=vel, compression="gzip", shuffle=True)
        h.create_dataset("initial_detection_speed_km_s", data=speed, compression="gzip", shuffle=True)
        h.create_dataset("fractional_initial_state_velocity_sigma", data=fractional_velocity_sigma, compression="gzip", shuffle=True)
        h.create_dataset("initial_state_radiant_angle_sigma_deg", data=angle, compression="gzip", shuffle=True)


def histogram_density(values: np.ndarray, bins: np.ndarray, logarithmic: bool) -> tuple[np.ndarray, np.ndarray]:
    """Return a unit-area density in linear units or probability per decade."""
    counts, edges = np.histogram(values, bins=bins)
    coordinate_edges = np.log10(edges) if logarithmic else edges
    widths = np.diff(coordinate_edges)
    density = counts / np.maximum(np.sum(counts) * widths, np.finfo(float).tiny)
    return density, edges


def plot_distribution(
    rows: np.ndarray,
    output_png: Path,
    output_pdf: Path | None,
    chosen_angle_deg: float,
    single_panel: bool,
    linear_xaxis: bool,
) -> None:
    angle = rows["initial_state_radiant_angle_sigma_deg"]
    velocity_sigma_mps = rows["initial_state_velocity_sigma_mps"]
    speed_mps = rows["initial_detection_speed_km_s"] * 1e3
    fractional_velocity_sigma = velocity_sigma_mps / np.maximum(speed_mps, 1.0)
    total = max(1, len(rows))
    angle_counts = threshold_counts(angle, ANGLE_THRESHOLDS_DEG)

    if single_panel:
        fig, ax0 = plt.subplots(1, 1, figsize=(4.8, 3.4), constrained_layout=True)
        axes = np.asarray([ax0])
    else:
        fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6), constrained_layout=True)

    if linear_xaxis:
        bins = np.linspace(0.0, max(3.05, np.nanpercentile(angle, 99.8)), 55)
    else:
        bins = np.geomspace(0.15, 60.0, 70)
    density, hist_edges = histogram_density(angle, bins, logarithmic=not linear_xaxis)
    peak_index = int(np.nanargmax(density)) if len(density) else 0
    if linear_xaxis:
        peak_angle_deg = 0.5 * (hist_edges[peak_index] + hist_edges[peak_index + 1])
    else:
        peak_angle_deg = np.sqrt(hist_edges[peak_index] * hist_edges[peak_index + 1])
    axes[0].stairs(
        density,
        hist_edges,
        color="#2f5f9f",
        linewidth=1.8,
    )
    axes[0].axvline(peak_angle_deg, color="black", lw=1.4, ls="--")
    if not linear_xaxis:
        axes[0].set_xscale("log")
    axes[0].set_xlabel("Velocity-vector angular uncertainty (deg)")
    axes[0].set_ylabel("Probability density per decade" if not linear_xaxis else r"Probability density (deg$^{-1}$)")
    axes[0].grid(True, which="both", alpha=0.25)

    good_fractional_sigma = fractional_velocity_sigma[
        np.isfinite(fractional_velocity_sigma) & (fractional_velocity_sigma > 0.0)
    ]
    top_axis = axes[0].twiny()
    if linear_xaxis:
        fractional_bins = np.linspace(0.0, max(1e-3, np.nanpercentile(good_fractional_sigma, 99.8)), 55)
    else:
        fractional_bins = np.geomspace(
            max(1e-5, np.nanpercentile(good_fractional_sigma, 0.1)),
            max(2e-5, np.nanpercentile(good_fractional_sigma, 99.8)),
            70,
        )
        top_axis.set_xscale("log")
    fractional_density, fractional_edges = histogram_density(
        good_fractional_sigma,
        fractional_bins,
        logarithmic=not linear_xaxis,
    )
    fractional_peak_index = int(np.nanargmax(fractional_density)) if len(fractional_density) else 0
    if linear_xaxis:
        peak_fractional_sigma = 0.5 * (
            fractional_edges[fractional_peak_index] + fractional_edges[fractional_peak_index + 1]
        )
    else:
        peak_fractional_sigma = np.sqrt(
            fractional_edges[fractional_peak_index] * fractional_edges[fractional_peak_index + 1]
        )
    top_axis.stairs(
        fractional_density,
        fractional_edges,
        color="#c45a1a",
        linewidth=1.5,
    )
    top_axis.set_xlabel(r"Fractional velocity uncertainty $\sigma_v/v_0$")
    top_axis.tick_params(axis="x", colors="#9f4613")
    top_axis.xaxis.label.set_color("#9f4613")
    axes[0].text(
        0.97,
        0.95,
        f"N = {len(rows):,}\nangle peak = {peak_angle_deg:.2f} deg\n$\\sigma_v/v_0$ peak = {peak_fractional_sigma:.3f}\n{chosen_angle_deg:.1f} deg ($\\sigma_v/v_0={np.tan(np.deg2rad(chosen_angle_deg)):.3f}$) keeps {100.0*np.sum(angle <= chosen_angle_deg)/total:.1f}%",
        transform=axes[0].transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
    )

    if not single_panel:
        axes[1].plot(ANGLE_THRESHOLDS_DEG, 100.0 * angle_counts / total, "o-", color="#3b6fb6", label="angular")
        axes[1].axvline(chosen_angle_deg, color="black", lw=1.4, ls="--")
        axes[1].set_xlabel("Radiant angular uncertainty threshold (deg)")
        axes[1].set_ylabel("Cumulative retained fraction (%)")
        axes[1].set_ylim(0, 100)
        axes[1].grid(True, alpha=0.25)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220)
    if output_pdf is not None:
        fig.savefig(output_pdf)
    plt.close(fig)


def write_tex_snippet(path: Path, figure_name: str, data_name: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        r"""\newif\ifshowscriptprovenance
\showscriptprovenancetrue

\begin{figure}
  \centering
  \includegraphics[width=\linewidth]{"""
        + figure_name
        + r"""}
  \caption{Distribution of the initial-state velocity uncertainty expressed as an upper bound
  on the radiant angular uncertainty, $\sigma_\theta \simeq \tan^{-1}(\sigma_v/v_g)$.
  The vertical dashed line marks the default high-quality threshold used for the aggregate
  radiant histogram.
  \ifshowscriptprovenance\emph{Script: \texttt{plot_radiant_uncertainty_distribution.py};
  data product: \texttt{"""
        + data_name
        + r"""}.}\fi}
  \label{fig:pansy-radiant-uncertainty}
\end{figure}
"""
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-metadata-dir", type=Path, default=DEFAULT_ORBIT_METADATA_DIR)
    parser.add_argument(
        "--radiants-h5",
        type=Path,
        default=None,
        help="Use a collected radiant sidecar containing a /radiants table instead of scanning orbit metadata.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--max-combined-score", type=float, default=1.5)
    parser.add_argument("--min-uncertainty-samples", type=int, default=3)
    parser.add_argument("--chosen-angle-deg", type=float, default=3.0)
    parser.add_argument("--single-panel", action="store_true", help="Only plot the histogram panel.")
    parser.add_argument("--linear-xaxis", action="store_true", help="Use a linear x-axis for the histogram.")
    args = parser.parse_args()

    if args.radiants_h5 is not None:
        with h5py.File(args.radiants_h5, "r") as h:
            rows = rows_from_radiant_table(h["radiants"][()])
        files_read = 1
        files_skipped = 0
        source = f"collected radiant sidecar: {args.radiants_h5}"
    else:
        rows, files_read, files_skipped = collect_uncertainty_rows(args.orbit_metadata_dir)
        source = "compact hourly orbit metadata events tables, direct covariance read"
    rows = finite_quality_rows(rows, args.max_combined_score, args.min_uncertainty_samples)
    out_h5 = args.output_dir / "radiant_uncertainty_distribution.h5"
    out_png = args.output_dir / "radiant_uncertainty_distribution.png"
    out_pdf = args.output_dir / "radiant_uncertainty_distribution.pdf"
    out_tex = args.output_dir / "radiant_uncertainty_distribution_snippet.tex"
    write_summary_h5(out_h5, rows, files_read, files_skipped, source)
    plot_distribution(rows, out_png, out_pdf, args.chosen_angle_deg, args.single_panel, args.linear_xaxis)
    write_tex_snippet(out_tex, out_png.name, out_h5.name)

    angle = rows["initial_state_radiant_angle_sigma_deg"]
    print(f"radiant_uncertainty_rows {len(rows)} files_read {files_read} files_skipped {files_skipped}")
    print(f"median_angle_deg {np.nanmedian(angle):.4f}")
    for threshold in ANGLE_THRESHOLDS_DEG:
        print(f"angle_threshold_deg {threshold:.2f} retained {int(np.sum(angle <= threshold))} fraction {np.sum(angle <= threshold) / max(1, len(angle)):.4f}")
    print(out_png)
    print(out_h5)


if __name__ == "__main__":
    main()
