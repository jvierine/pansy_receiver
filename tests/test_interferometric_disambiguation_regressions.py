import math
import subprocess
import sys
import itertools
from pathlib import Path

import h5py
import numpy as np
import pytest


PANSY_RECEIVER_ROOT = Path(__file__).resolve().parents[1]
if str(PANSY_RECEIVER_ROOT) not in sys.path:
    sys.path.insert(0, str(PANSY_RECEIVER_ROOT))
MULTI_METEOR_SAMPLE_581983 = 1746492577581983
TRAIL_ECHO_SAMPLE_742042 = 1746493024742042
TRAIL_ECHO_SAMPLE_606335 = 1746495226606335
RANGE_FIT_SAMPLE_450693 = 1746575731450693
RANGE_FIT_SAMPLE_895405 = 1746575072895405
REGRESSION_EVENTS = [
    pytest.param(1746489745288806, "88806", 101, id="event_88806"),
    pytest.param(1746489819007216, "07216", 501, id="event_07216"),
    pytest.param(1746490224199270, "99270", 701, id="event_99270"),
    pytest.param(1746492411469961, "469961", 501, id="event_469961"),
    pytest.param(1746492839697218, "97218", 501, id="event_97218"),
]

TRAIL_ECHO_EVENTS = [
    pytest.param(TRAIL_ECHO_SAMPLE_742042, "742042", True, id="event_742042"),
    pytest.param(TRAIL_ECHO_SAMPLE_606335, "606335", False, id="event_606335"),
]


def _subset_observations(obs, idx):
    return {
        key: val[idx] if isinstance(val, np.ndarray) and len(val) == len(obs["snr"]) else val
        for key, val in obs.items()
    }


def _fit_best_component(obs, sample_idx, grid_n=501):
    import pansy_interferometry as pint
    import plot_interferometric_disambiguation as disamb

    antpos = pint.get_antpos()
    ch_pairs = np.asarray(list(itertools.combinations(np.arange(7), 2)))
    dmat = pint.pair_mat(ch_pairs, antpos)
    phasecal = pint.get_phasecal()
    u, v, w, valid = disamb.horizon_grid(grid_n)
    steering, _uvw = disamb.steering_matrix(dmat, u, v, w, valid)
    tx_gain_maps = disamb.precompute_tx_array_gain_maps(u, v, w, valid)
    tx_lobe_centroids, _tx_lobe_maps = disamb.tx_grating_lobe_centroids_from_gain_maps(tx_gain_maps, u, v, valid)

    candidates = []
    t_rel = obs["tx_idx"] / 1e6 - obs["tx_idx"][0] / 1e6
    for i in range(len(obs["snr"])):
        coh = disamb.coherence_map(obs["xc"][i], int(obs["beam_id"][i]), phasecal, ch_pairs, steering, valid, u.shape)
        ii, jj = disamb.local_coherence_peaks(coh, 0.9)
        for row, col in zip(ii, jj):
            candidates.append(
                {
                    "u": float(u[row, col]),
                    "v": float(v[row, col]),
                    "w": float(w[row, col]),
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

    tracks = disamb.complete_tracks_to_common_pulses(disamb.fit_candidate_tracks(candidates), candidates)
    t_start = float(np.min(t_rel))
    t_end = float(np.max(t_rel))
    tracks = [disamb.classify_track_visibility(track, candidates, t_start, t_end) for track in tracks]
    tracks = [
        disamb.classify_track_descent(track) if track["reason"] == "kept" else track
        for track in tracks
    ]
    tracks = [
        disamb.classify_track_linearity(track, candidates)
        if track["reason"] == "kept" and not track.get("descent_reject", False)
        else track
        for track in tracks
    ]
    sigma_pos, sigma_dop = disamb.fit_ballistic_survivors(tracks, candidates, event_epoch_unix=sample_idx / 1e6)
    rho_of_alt_m, _msis_meta = disamb.pbal.density_interpolator(sample_idx / 1e6)
    disamb.fit_fixed_velocity_survivors(tracks, rho_of_alt_m, sigma_pos, sigma_dop)
    disamb.score_tx_beam_consistency(tracks, candidates)
    disamb.score_candidate_diagnostics(
        tracks,
        candidates,
        tx_gain_maps=tx_gain_maps,
        tx_lobe_centroids=tx_lobe_centroids,
    )
    disamb.score_combined_hypotheses(tracks)
    ranked = sorted([track for track in tracks if "combined_score" in track], key=lambda track: track["combined_rank"])
    return ranked[0] if ranked else None


@pytest.mark.slow
@pytest.mark.parametrize(("sample_idx", "suffix", "grid_n"), REGRESSION_EVENTS)
def test_event_finds_valid_trajectory(sample_idx, suffix, grid_n, tmp_path):
    cut_dir = PANSY_RECEIVER_ROOT / "data" / "metadata" / "cut"
    if not cut_dir.exists():
        pytest.skip(f"local cut metadata not available: {cut_dir}")

    output_dir = tmp_path / f"event_{suffix}"
    cmd = [
        sys.executable,
        str(PANSY_RECEIVER_ROOT / "plot_interferometric_disambiguation.py"),
        "--sample-idx",
        str(sample_idx),
        "--cut-dir",
        str(cut_dir),
        "--output-dir",
        str(output_dir),
        "--grid-n",
        str(grid_n),
        "--overview-only",
        "--orbit-samples",
        "0",
    ]
    result = subprocess.run(
        cmd,
        cwd=PANSY_RECEIVER_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=300,
        check=False,
    )
    assert result.returncode == 0, result.stdout

    diagnostics_h5 = output_dir / f"pansy_disambiguation_diagnostics_{sample_idx}.h5"
    state_h5 = output_dir / f"pansy_candidate_orbit_states_{sample_idx}.h5"
    summary_png = output_dir / f"pansy_interferometer_disambiguation_summary_{sample_idx}.png"
    assert diagnostics_h5.exists(), result.stdout
    assert state_h5.exists(), result.stdout
    assert summary_png.exists(), result.stdout

    with h5py.File(diagnostics_h5, "r") as h:
        assert h.attrs["event_status"] == "processed"
        selected = str(h.attrs["selected_hypothesis"])
        selected_score = float(h.attrs["selected_combined_score"])
        assert selected.startswith("H")
        assert math.isfinite(selected_score)

        hypothesis = h["hypotheses"][selected]
        assert int(hypothesis.attrs["combined_rank"]) == 0
        assert not bool(hypothesis.attrs["combined_reject"])
        assert bool(hypothesis.attrs["combined_good_fit"])
        assert str(hypothesis.attrs["selection_model_type"]) in {"fixed_velocity", "msis_drag"}
        assert math.isfinite(float(hypothesis.attrs["selection_reduced_chi2"]))
        assert float(hypothesis.attrs["selection_reduced_chi2"]) < 1.5
        assert math.isfinite(float(hypothesis.attrs["combined_score"]))

    with h5py.File(state_h5, "r") as h:
        assert selected in h
        assert int(h[selected].attrs["combined_rank"]) == 0


@pytest.mark.slow
def test_event_581983_identifies_two_meteors_separately():
    import plot_interferometric_disambiguation as disamb

    cut_dir = PANSY_RECEIVER_ROOT / "data" / "metadata" / "cut"
    if not cut_dir.exists():
        pytest.skip(f"local cut metadata not available: {cut_dir}")

    cut = disamb.load_cut(cut_dir, MULTI_METEOR_SAMPLE_581983)
    obs_all = disamb.recompute_cut_observables(cut)
    good = obs_all["snr"] > 9.0
    obs = {
        key: val[good] if isinstance(val, np.ndarray) and len(val) == len(good) else val
        for key, val in obs_all.items()
    }

    segments = disamb.split_observations_by_range_time(obs, min_points=3)
    assert [len(segment) for segment in segments] == [6, 3]

    best_tracks = []
    for segment in segments:
        segment_obs = _subset_observations(obs, segment)
        assert float(segment_obs["range_km"][-1]) < float(segment_obs["range_km"][0])
        best = _fit_best_component(segment_obs, MULTI_METEOR_SAMPLE_581983, grid_n=501)
        assert best is not None
        assert not bool(best.get("combined_reject", True))
        assert str(best["selection_model_type"]) in {"fixed_velocity", "msis_drag"}
        assert math.isfinite(float(best["selection_reduced_chi2"]))
        assert int(best["unique_pulses"]) == len(segment)
        best_tracks.append(best)

    assert float(best_tracks[0]["selection_reduced_chi2"]) < 1.5
    assert float(best_tracks[1]["selection_reduced_chi2"]) < 3.0


@pytest.mark.slow
@pytest.mark.parametrize(("sample_idx", "suffix", "expect_speed_reject"), TRAIL_ECHO_EVENTS)
def test_trail_echo_is_not_accepted_as_head_echo(sample_idx, suffix, expect_speed_reject, tmp_path):
    cut_dir = PANSY_RECEIVER_ROOT / "data" / "metadata" / "cut"
    if not cut_dir.exists():
        pytest.skip(f"local cut metadata not available: {cut_dir}")

    output_dir = tmp_path / f"event_{suffix}"
    cmd = [
        sys.executable,
        str(PANSY_RECEIVER_ROOT / "plot_interferometric_disambiguation.py"),
        "--sample-idx",
        str(sample_idx),
        "--cut-dir",
        str(cut_dir),
        "--output-dir",
        str(output_dir),
        "--grid-n",
        "501",
        "--overview-only",
        "--orbit-samples",
        "0",
    ]
    result = subprocess.run(
        cmd,
        cwd=PANSY_RECEIVER_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=300,
        check=False,
    )
    assert result.returncode == 0, result.stdout

    diagnostics_h5 = output_dir / f"pansy_disambiguation_diagnostics_{sample_idx}.h5"
    assert diagnostics_h5.exists(), result.stdout

    with h5py.File(diagnostics_h5, "r") as h:
        assert h.attrs["event_status"] == "processed"
        ranked = [
            hypothesis
            for hypothesis in h["hypotheses"].values()
            if int(hypothesis.attrs.get("combined_rank", -1)) >= 0
        ]
        assert not any(
            (not bool(hypothesis.attrs["combined_reject"])) and bool(hypothesis.attrs["combined_good_fit"])
            for hypothesis in ranked
        )
        if expect_speed_reject:
            assert ranked, result.stdout
            assert all(bool(hypothesis.attrs["head_echo_speed_reject"]) for hypothesis in ranked)
            assert max(float(hypothesis.attrs["selection_speed_km_s"]) for hypothesis in ranked) < 5.0
        else:
            assert not ranked


@pytest.mark.slow
def test_event_450693_correct_hypothesis_has_good_range_fit(tmp_path):
    cut_dir = PANSY_RECEIVER_ROOT / "data" / "metadata" / "cut"
    if not cut_dir.exists():
        pytest.skip(f"local cut metadata not available: {cut_dir}")

    output_dir = tmp_path / "event_450693"
    cmd = [
        sys.executable,
        str(PANSY_RECEIVER_ROOT / "plot_interferometric_disambiguation.py"),
        "--sample-idx",
        str(RANGE_FIT_SAMPLE_450693),
        "--cut-dir",
        str(cut_dir),
        "--output-dir",
        str(output_dir),
        "--grid-n",
        "501",
        "--overview-only",
        "--orbit-samples",
        "0",
    ]
    result = subprocess.run(
        cmd,
        cwd=PANSY_RECEIVER_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=300,
        check=False,
    )
    assert result.returncode == 0, result.stdout

    diagnostics_h5 = output_dir / f"pansy_disambiguation_diagnostics_{RANGE_FIT_SAMPLE_450693}.h5"
    assert diagnostics_h5.exists(), result.stdout

    with h5py.File(diagnostics_h5, "r") as h:
        assert h.attrs["event_status"] == "processed"
        selected = str(h.attrs["selected_hypothesis"])
        assert selected == "H04"

        hypothesis = h["hypotheses"][selected]
        assert int(hypothesis.attrs["combined_rank"]) == 0
        assert not bool(hypothesis.attrs["combined_reject"])
        assert bool(hypothesis.attrs["combined_good_fit"])
        assert str(hypothesis.attrs["selection_model_type"]) == "fixed_velocity"
        assert float(hypothesis.attrs["selection_reduced_chi2"]) < 1.0

        points = hypothesis["position_enu_km"][:]
        model = hypothesis["selection_model"][:]
        keep = hypothesis["selection_keep"][:]
        range_res_m = (np.linalg.norm(points, axis=1) - np.linalg.norm(model, axis=1)) * 1e3
        kept_res = range_res_m[keep]
        assert len(kept_res) >= 30
        assert abs(float(np.mean(kept_res))) < 25.0
        assert float(np.sqrt(np.mean(kept_res**2))) < 100.0


@pytest.mark.slow
def test_event_895405_correct_hypothesis_has_good_range_fit(tmp_path):
    cut_dir = PANSY_RECEIVER_ROOT / "data" / "metadata" / "cut"
    if not cut_dir.exists():
        pytest.skip(f"local cut metadata not available: {cut_dir}")

    output_dir = tmp_path / "event_895405"
    cmd = [
        sys.executable,
        str(PANSY_RECEIVER_ROOT / "plot_interferometric_disambiguation.py"),
        "--sample-idx",
        str(RANGE_FIT_SAMPLE_895405),
        "--cut-dir",
        str(cut_dir),
        "--output-dir",
        str(output_dir),
        "--grid-n",
        "501",
        "--overview-only",
        "--orbit-samples",
        "0",
    ]
    result = subprocess.run(
        cmd,
        cwd=PANSY_RECEIVER_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=300,
        check=False,
    )
    assert result.returncode == 0, result.stdout

    diagnostics_h5 = output_dir / f"pansy_disambiguation_diagnostics_{RANGE_FIT_SAMPLE_895405}.h5"
    assert diagnostics_h5.exists(), result.stdout

    with h5py.File(diagnostics_h5, "r") as h:
        assert h.attrs["event_status"] == "processed"
        selected = str(h.attrs["selected_hypothesis"])
        assert selected == "H01"

        hypothesis = h["hypotheses"][selected]
        assert int(hypothesis.attrs["combined_rank"]) == 0
        assert not bool(hypothesis.attrs["combined_reject"])
        assert bool(hypothesis.attrs["combined_good_fit"])
        assert str(hypothesis.attrs["selection_model_type"]) == "fixed_velocity"
        assert float(hypothesis.attrs["selection_reduced_chi2"]) < 1.0

        points = hypothesis["position_enu_km"][:]
        model = hypothesis["selection_model"][:]
        keep = hypothesis["selection_keep"][:]
        range_res_m = (np.linalg.norm(points, axis=1) - np.linalg.norm(model, axis=1)) * 1e3
        kept_res = range_res_m[keep]
        assert len(kept_res) >= 65
        assert abs(float(np.mean(kept_res))) < 25.0
        assert float(np.sqrt(np.mean(kept_res**2))) < 100.0
