#!/usr/bin/env python3
"""Build the mass-analysis manifest for the CAP/DCS orbit-panel events."""

from __future__ import annotations

import argparse
import datetime as dt
import os
from pathlib import Path

import h5py
import numpy as np

import plot_capricornid_conjugate_stream as capplot


def diagnostics_path(events_dir: Path, sample_idx: int) -> Path:
    day = dt.datetime.fromtimestamp(
        sample_idx / 1e6, tz=dt.timezone.utc
    ).date().isoformat()
    return events_dir / day / f"pansy_disambiguation_diagnostics_{sample_idx}.h5"


def selected_rows(catalogue: Path) -> list[tuple[int, capplot.Passage, np.ndarray]]:
    selections = []
    for passage_index, passage in enumerate(capplot.PASSAGES):
        rows = capplot.load_passage_rows(
            capplot.RADVIEW_DATA,
            passage,
            capplot.CLUSTER_SOLAR_WINDOW_DEG / 2.0,
            catalogue_path=catalogue,
            source="PANSY",
        )
        rows = capplot.cluster_filter(
            rows, capplot.CLUSTER_VG_RANGE, capplot.CLUSTER_E_RANGE
        )
        rows = capplot.select_associated(rows, passage, 5.0, 10.0)
        selections.append((passage_index, passage, rows))
    return selections


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, default=capplot.DEFAULT_CATALOGUE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--events-dir",
        type=Path,
        help="If supplied, require a diagnostics file for every selected event.",
    )
    args = parser.parse_args()

    parts = selected_rows(args.catalogue)
    passage_index = np.concatenate(
        [np.full(len(rows), index, dtype=np.int8) for index, _passage, rows in parts]
    )
    rows = np.concatenate([rows for _index, _passage, rows in parts])
    sample_idx = np.rint(rows["epoch"] * 1e6).astype(np.int64)
    if len(np.unique(sample_idx)) != len(sample_idx):
        values, counts = np.unique(sample_idx, return_counts=True)
        duplicate = values[counts > 1]
        raise RuntimeError(f"duplicate event identifiers in shower selections: {duplicate}")

    diagnostics = None
    if args.events_dir is not None:
        diagnostics = np.asarray(
            [str(diagnostics_path(args.events_dir, int(sample))) for sample in sample_idx],
            dtype=object,
        )
        missing = [path for path in diagnostics if not Path(path).is_file()]
        if missing:
            preview = "\n".join(missing[:10])
            raise FileNotFoundError(
                f"missing diagnostics for {len(missing)} of {len(sample_idx)} selected events:\n{preview}"
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["schema"] = "pansy.cap_dcs_three_pulse_mass_manifest.v1"
        handle.attrs["created_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        handle.attrs["catalogue"] = str(args.catalogue)
        handle.attrs["selection_source"] = (
            "the unchanged PANSY CAP/DCS orbit-panel selection in "
            "plot_capricornid_conjugate_stream.py"
        )
        handle.attrs["mass_blinding"] = (
            "selection is fixed before mass fitting; no fitted mass or orbit-contraction "
            "criterion is used for membership or acceptance"
        )
        handle.attrs["sample_idx_conversion"] = "round(epoch_unix * 1e6)"
        handle.create_dataset("sample_idx", data=sample_idx)
        handle.create_dataset("passage_index", data=passage_index)
        handle.create_dataset("epoch_unix", data=rows["epoch"])
        handle.create_dataset("solar_longitude_deg", data=rows["solar_lon"])
        handle.create_dataset("sun_centered_lon_deg", data=rows["sun_centered_lon"])
        handle.create_dataset("ecliptic_latitude_deg", data=rows["beta"])
        handle.create_dataset("geocentric_speed_km_s", data=rows["vg"])
        handle.create_dataset("kepler", data=rows["kepler"])
        handle["kepler"].attrs["columns"] = (
            "a_AU,e,i_deg,Omega_deg,omega_deg,nu_deg,q_AU"
        )
        names = np.asarray([passage.name for passage in capplot.PASSAGES], dtype=object)
        handle.create_dataset("passage_names", data=names, dtype=string_dtype)
        if diagnostics is not None:
            handle.create_dataset("diagnostics_h5", data=diagnostics, dtype=string_dtype)
        for index, passage, selected in parts:
            handle.attrs[f"passage_{index}_name"] = passage.name
            handle.attrs[f"passage_{index}_count"] = len(selected)
    os.replace(temporary, args.output)

    print(f"manifest={args.output} events={len(sample_idx)}")
    for index, passage, selected in parts:
        print(f"passage={index} name={passage.name} events={len(selected)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
