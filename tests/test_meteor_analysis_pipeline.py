import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import run_meteor_analysis_pipeline as pipeline


def config(tmp_path):
    return pipeline.PipelineConfig(
        cut_dir=tmp_path / "cut",
        output_dir=tmp_path / "out",
        table_dir=tmp_path / "tables",
        grid_n=11,
        coherence_threshold=0.9,
        max_peaks_per_pulse=32,
        snr_threshold=9.0,
        linearity_angular_sigma_deg=0.25,
        linearity_range_sigma_km=0.15,
        linearity_p_threshold=0.01,
        ballistic_p_threshold=0.01,
        target_alt_km=140.0,
    )


def observations(snr):
    snr = np.asarray(snr, dtype=np.float64)
    n = len(snr)
    return {
        "snr": snr,
        "tx_idx": 1_700_000_000_000_000 + np.arange(n, dtype=np.float64) * 10_000,
        "range_km": np.linspace(90.0, 100.0, n),
        "doppler_mps": np.linspace(-12000.0, -11800.0, n),
        "beam_id": np.zeros(n, dtype=np.int64),
        "xc": np.ones((n, 2), dtype=np.complex64),
        "non_pulse_metadata": "kept",
    }


def test_filter_observations_by_snr_keeps_pulse_aligned_arrays():
    obs = observations([5.0, 10.0, 12.0, 4.0])

    filtered = pipeline.filter_observations_by_snr(obs, 9.0)

    assert filtered["snr"].tolist() == [10.0, 12.0]
    assert filtered["range_km"].tolist() == pytest.approx([93.3333333333, 96.6666666667])
    assert filtered["beam_id"].tolist() == [0, 0]
    assert filtered["non_pulse_metadata"] == "kept"


def test_standard_disambiguation_rejects_short_events(tmp_path):
    with pytest.raises(RuntimeError, match="only 9 pulses above SNR threshold"):
        pipeline.run_standard_disambiguation(
            observations(np.full(9, 12.0)),
            sample_idx=1_700_000_000_000_000,
            config=config(tmp_path),
            context={"tx_gain_maps": np.zeros((5, 1, 1))},
        )


def test_standard_disambiguation_runs_interferometric_ranking(monkeypatch, tmp_path):
    obs = observations(np.full(12, 12.0))
    candidates = [{"pulse": i, "snr": 12.0} for i in range(12)]
    tracks = [
        {"hypothesis_id": 1, "idx": np.arange(12), "reason": "kept"},
        {"hypothesis_id": 2, "idx": np.arange(12), "reason": "kept"},
    ]
    calls = []

    def fake_build_candidates(filtered, context, coherence_threshold, max_peaks_per_pulse=32):
        calls.append(("build", len(filtered["snr"]), coherence_threshold, max_peaks_per_pulse))
        return candidates, (np.array([0]), np.array([0])), np.ones((1, 1), dtype=np.float32)

    def fake_fit_candidate_tracks(got_candidates):
        calls.append(("fit_tracks", got_candidates is candidates))
        return tracks

    def fake_visibility(track, got_candidates, t_start, t_end):
        calls.append(("visibility", track["hypothesis_id"], t_start, t_end))
        return track

    def fake_linearity(track, got_candidates, **kwargs):
        calls.append(("linearity", track["hypothesis_id"], kwargs["p_threshold"]))
        track["linearity_reject"] = False
        track["low_detection_altitude_reject"] = False
        return track

    def fake_descent(track):
        calls.append(("descent", track["hypothesis_id"]))
        track["descent_reject"] = False
        return track

    def fake_ballistic(got_tracks, got_candidates, event_epoch_unix, p_threshold):
        calls.append(("ballistic", event_epoch_unix, p_threshold))
        got_tracks[0].update({"ballistic_reduced_chi2": 4.0, "combined_rank": 1, "combined_score": 4.0})
        got_tracks[1].update({"ballistic_reduced_chi2": 1.0, "combined_rank": 0, "combined_score": 1.0})
        return 0.1, 0.2

    def fake_score_combined(got_tracks):
        calls.append(("combined", [t["hypothesis_id"] for t in got_tracks]))
        return 0.3

    def fake_orbit(track, sample_epoch_unix):
        calls.append(("orbit", track["hypothesis_id"], sample_epoch_unix))
        return {
            "kepler": np.array([1.0, 0.5, 30.0, 0.0, 0.0, 0.0, 0.8]),
            "state_gcrs_m_mps": np.array([0.0, 0.0, 0.0, 1000.0, 2000.0, 2000.0]),
        }

    monkeypatch.setattr(pipeline, "build_candidates", fake_build_candidates)
    monkeypatch.setattr(pipeline.disamb, "fit_candidate_tracks", fake_fit_candidate_tracks)
    monkeypatch.setattr(pipeline.disamb, "classify_track_visibility", fake_visibility)
    monkeypatch.setattr(pipeline.disamb, "classify_track_linearity", fake_linearity)
    monkeypatch.setattr(pipeline.disamb, "classify_track_descent", fake_descent)
    monkeypatch.setattr(pipeline.disamb, "fit_ballistic_survivors", fake_ballistic)
    monkeypatch.setattr(pipeline.disamb, "score_tx_beam_consistency", lambda *args, **kwargs: calls.append(("tx_beam",)))
    monkeypatch.setattr(pipeline.disamb, "score_candidate_diagnostics", lambda *args, **kwargs: calls.append(("diagnostics",)))
    monkeypatch.setattr(pipeline.disamb, "score_combined_hypotheses", fake_score_combined)
    monkeypatch.setattr(pipeline.disamb, "orbit_for_candidate_track", fake_orbit)

    result = pipeline.run_standard_disambiguation(
        obs,
        sample_idx=1_700_000_000_000_000,
        config=config(tmp_path),
        context={"tx_gain_maps": np.zeros((5, 1, 1))},
    )

    assert result["pipeline"] == pipeline.STANDARD_PIPELINE_NAME
    assert result["ranked"][0]["hypothesis_id"] == 2
    assert result["ranked"][1]["hypothesis_id"] == 1
    assert result["sigma_pos"] == pytest.approx(0.1)
    assert result["sigma_dop"] == pytest.approx(0.2)
    assert result["tx_sigma"] == pytest.approx(0.3)
    assert result["ranked"][0]["candidate_orbit"]["kepler"][1] == pytest.approx(0.5)
    assert ("build", 12, 0.9, 32) in calls
    assert ("fit_tracks", True) in calls
    assert ("tx_beam",) in calls
    assert ("diagnostics",) in calls
    assert ("combined", [1, 2]) in calls
    assert ("orbit", 2, 1_700_000_000.0) in calls
