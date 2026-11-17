"""Validate standalone blind detection and new cuts against frozen native results."""

import argparse
from pathlib import Path
import h5py
import numpy as np
from scipy.io import loadmat
from meteor_analysis.voltage import load_mat_cut, read_iq


def write_group(group, values):
    for key, value in values.items():
        if isinstance(value, dict):
            write_group(group.create_group(key), value)
        else:
            value = np.asarray(value)
            filters = (
                dict(compression="gzip", compression_opts=4, shuffle=True)
                if value.ndim and value.size
                else {}
            )
            group.create_dataset(key, data=value, **filters)


def compare(actual, reference, original_cuts, output):
    assert not output.exists(), "Use a new verification output"
    records = []
    file_rows = []
    with h5py.File(reference / "catalogue.h5") as catalogue:
        for source_id, source in sorted(catalogue["sources"].items()):
            source_index = int(source_id)
            stem = Path(bytes(source["source_path"][()]).decode()).stem
            directory = actual / stem
            result = loadmat(directory / "catalogue.mat", simplify_cells=True)["result"]
            assert result["complete"] and result["processed_pulses"] == int(
                source["processed_pulses"][()]
            )
            assert (
                result["candidate_count"] == source["candidate_points/pulse_index"].size
            ), (stem, "candidate count")
            assert (
                result["code_alias_veto_count"]
                == source["code_alias_vetoes/pulse_index"].size
            ), (stem, "code veto count")
            expected = [
                name
                for name, event in sorted(catalogue["events"].items())
                if int(event["source_file_index"][()]) == source_index
            ]
            assert result["event_count"] == len(expected)
            assert len(list(directory.glob("event_*.mat"))) == len(expected)
            file_rows.append(
                [
                    source_index,
                    result["processed_pulses"],
                    result["candidate_count"],
                    result["code_alias_veto_count"],
                    len(expected),
                    result["processing_seconds"],
                ]
            )
            for local_id, name in enumerate(expected):
                path = directory / f"event_{local_id:05}.mat"
                record = load_mat_cut(path)
                measured = record["detections"]
                event = catalogue["events/" + name]
                for key in [
                    "pulse_index",
                    "delay_samples",
                    "doppler_bin",
                    "strongest_full_array_rx",
                    "tx_time_us",
                ]:
                    np.testing.assert_array_equal(
                        measured[key], event[key][()], err_msg=f"{stem}/{name}/{key}"
                    )
                fd_error = float(
                    np.max(abs(measured["doppler_hz"] - event["doppler_hz"][()]))
                )
                statistic_error = float(
                    np.max(
                        abs(
                            measured["detector_excess_over_median_db"]
                            - event["detector_excess_over_median_db"][()]
                        )
                    )
                )
                assert fd_error < 0.01, (stem, name, "Doppler", fd_error)
                assert statistic_error < 0.00001, (
                    stem,
                    name,
                    "statistic",
                    statistic_error,
                )
                amplitude = measured["matched_amplitude_re_im"]
                target = event["matched_amplitude_re_im"][:, 5:]
                delta = np.linalg.norm(amplitude - target, axis=-1)
                relative = delta / np.maximum(np.linalg.norm(target, axis=-1), 1.0)
                amplitude_error = float(np.max(relative))
                assert amplitude_error < 0.000003, (
                    stem,
                    name,
                    "amplitude",
                    amplitude_error,
                )
                components = 0
                with h5py.File(original_cuts / "cuts" / f"{name}.h5") as cut:
                    for key in [
                        "pulse_index",
                        "ipp_start_time_us",
                        "tx_history_valid",
                        "tx_history_decoded_baud_sign",
                        "raw_header",
                        "source_channel_indices",
                        "tx_alpha_per_pulse",
                        "tx_beta_per_pulse",
                        "tx_beam_id",
                    ]:
                        np.testing.assert_array_equal(record[key], cut[key][()])
                    for key in ["zrx_echoes", "ztx_history"]:
                        a = read_iq(record, key)
                        b = read_iq(cut, key)
                        np.testing.assert_array_equal(
                            a, b, err_msg=f"{stem}/{name}/{key}"
                        )
                        components += a.size
                records.append(
                    (
                        name,
                        record,
                        [
                            len(measured["pulse_index"]),
                            components,
                            fd_error,
                            statistic_error,
                            amplitude_error,
                        ],
                    )
                )
                print(
                    name,
                    "PASS",
                    "points",
                    len(measured["pulse_index"]),
                    "max fd difference",
                    fd_error,
                    flush=True,
                )
    output.mkdir(parents=True)
    replay = output / "native_replay"
    replay.mkdir()
    for name, record, _ in records:
        with h5py.File(replay / f"{name}.h5", "w") as h:
            write_group(h, record)
    with h5py.File(output / "verification.h5", "w") as h:
        h["file_metrics"] = file_rows
        h["event_metrics"] = [row for _, _, row in records]
        h.attrs["file_columns"] = (
            "source, IPPs, candidates, code vetoes, events, processing seconds"
        )
        h.attrs["event_columns"] = (
            "detections, exact IQ components, max fd difference Hz, max detector-statistic difference dB, max relative amplitude difference"
        )
        h.attrs["scope"] = (
            "Standalone Octave blind search of all original raw records; exact detection-index/bin/voltage comparison with frozen Rust results. Windows MATLAB execution remains untested."
        )
        h["complete"] = np.uint8(1)
    metrics = np.asarray([row for _, _, row in records])
    print(
        "PASS",
        len(file_rows),
        "files;",
        len(records),
        "events;",
        int(metrics[:, 0].sum()),
        "detections;",
        int(metrics[:, 1].sum()),
        "exact IQ components",
    )
    print(
        "Maximum fd difference",
        metrics[:, 2].max(),
        "Hz; maximum amplitude difference",
        metrics[:, 4].max(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ["actual", "reference", "original_cuts", "output"]:
        parser.add_argument(key, type=Path)
    args = parser.parse_args()
    compare(args.actual, args.reference, args.original_cuts, args.output)
