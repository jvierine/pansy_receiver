from pathlib import Path

import h5py
import numpy as np
import pytest

import tx_phase_quality as txpq


def write_phase_file(phase_dir: Path, hour: str, samples: list[int]) -> None:
    directory = phase_dir / hour
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"txphase@{samples[0] // 1_000_000}.h5"
    with h5py.File(path, "w") as h5:
        for sample in samples:
            group = h5.create_group(str(sample))
            phase = np.linspace(0.0, 0.07, 8, dtype=np.float32)
            group.create_dataset("xphase", data=np.exp(1j * phase).astype(np.complex64))


def test_preflight_refreshes_stale_quality_table(tmp_path: Path) -> None:
    phase_dir = tmp_path / "phase"
    quality_h5 = tmp_path / "tx_phase_quality.h5"
    old = [1_780_000_000_000_000, 1_780_000_060_000_000]
    new = [1_780_010_000_000_000, 1_780_010_060_000_000]
    write_phase_file(phase_dir, "2026-05-27T00-00-00", old)
    txpq.write_quality_h5(phase_dir, quality_h5)
    assert txpq.quality_table_coverage(quality_h5)["last_sample_idx"] == old[-1]

    write_phase_file(phase_dir, "2026-05-27T02-00-00", new)
    result = txpq.ensure_quality_table(
        quality_h5,
        old[0],
        new[-1],
        phase_dir=phase_dir,
    )
    assert result["action"] == "refreshed"
    assert result["last_sample_idx"] == new[-1]


def test_preflight_keeps_current_quality_table(tmp_path: Path) -> None:
    phase_dir = tmp_path / "phase"
    quality_h5 = tmp_path / "tx_phase_quality.h5"
    samples = [1_780_000_000_000_000, 1_780_000_060_000_000]
    write_phase_file(phase_dir, "2026-05-27T00-00-00", samples)
    txpq.write_quality_h5(phase_dir, quality_h5)

    result = txpq.ensure_quality_table(
        quality_h5,
        samples[0],
        samples[-1],
        phase_dir=phase_dir,
    )
    assert result["action"] == "current"


def test_preflight_fails_when_phase_metadata_cannot_cover_request(tmp_path: Path) -> None:
    phase_dir = tmp_path / "phase"
    quality_h5 = tmp_path / "tx_phase_quality.h5"
    samples = [1_780_000_000_000_000, 1_780_000_060_000_000]
    write_phase_file(phase_dir, "2026-05-27T00-00-00", samples)
    txpq.write_quality_h5(phase_dir, quality_h5, max_nearest_age_s=10.0)

    with pytest.raises(RuntimeError, match="did not cover"):
        txpq.ensure_quality_table(
            quality_h5,
            samples[0],
            samples[-1] + 60_000_000,
            phase_dir=phase_dir,
            max_age_s=10.0,
        )
