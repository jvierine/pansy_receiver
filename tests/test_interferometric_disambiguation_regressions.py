import math
import subprocess
import sys
from pathlib import Path

import h5py
import pytest


PANSY_RECEIVER_ROOT = Path(__file__).resolve().parents[1]
REGRESSION_EVENTS = [
    pytest.param(1746489745288806, "88806", 101, id="event_88806"),
    pytest.param(1746489819007216, "07216", 501, id="event_07216"),
]


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
