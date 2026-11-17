#!/usr/bin/env python
"""Publish native ten-channel cuts using the standard Digital RF metadata writer."""
import argparse
from pathlib import Path
import hashlib
import math
import digital_rf
import h5py
import numpy as np

NAMES = np.array(
    ["S0", "S1", "S2", "L0", "L1", "L2", "L3", "L4", "L5", "L6"], dtype="S2"
)


def unpack(group):
    return {k: v[()] for k, v in group.items() if isinstance(v, h5py.Dataset)}


def publish(root, run, geometry):
    with h5py.File(root / "index.h5") as index:
        assert index["complete"][()] == 1
        rates = {int(g["sample_rate_hz"][()]) for g in index["sources"].values()}
        if len(rates) != 1:
            raise ValueError("Mixed sample rates require separate metadata channels")
        fs = rates.pop()
        divisor = 1_000_000 // math.gcd(fs, 1_000_000)
        multiplier = fs // math.gcd(fs, 1_000_000)

        def sample_index(us):
            if np.any(us % divisor):
                raise ValueError("Nonintegral metadata sample timestamps")
            return (us // divisor) * multiplier

        directory = root / "metadata.partial"
        if directory.exists() or (root / "metadata").exists():
            raise FileExistsError("Metadata destination exists")
        directory.mkdir()

        def writer(name):
            path = directory / name
            path.mkdir()
            return digital_rf.DigitalMetadataWriter(str(path), 3600, 60, fs, 1, name)

        cut_writer = writer("cut")
        detection_writer = writer("detections")
        used = set()
        first_key = None
        for name, event in index["events"].items():
            cut_path = root / bytes(event["cut_file"][()]).decode()
            with h5py.File(cut_path) as cut:
                record = unpack(cut)
                data = unpack(cut["detections"])
                for field in [
                    "rx_alpha_code",
                    "rx_beta_code",
                    "rx_direction_uv_inferred",
                    "phase_weight_index",
                    "amplitude_weight_index",
                ]:
                    record["source_" + field] = record[field]
                    record[field] = record[field][5:]
                compact = "zrx_echoes_iq_i16" in record
                if compact:
                    tx_count = record["ztx_history_iq_i16"].shape[0]
                else:
                    iq = record.pop("zrx_echoes_iq_i32")
                    tx = record.pop("ztx_history_iq_i32")
                    tx_count = len(tx)
                ipp_us = record["ipp_start_time_us"]
                ipp_idx = sample_index(ipp_us)
                tx_idx = ipp_idx + int(record["tx_start_gate"])
                actual_key = int(data["tx_time_us"][0]) * fs // 1_000_000 + int(
                    record["tx_start_gate"]
                )
                key = actual_key
                while key in used:
                    key += 1
                used.add(key)
                first_key = (
                    int(tx_idx[0])
                    if first_key is None
                    else min(first_key, int(tx_idx[0]))
                )
                record.update(
                    tx_idx=tx_idx,
                    ipp_idx=ipp_idx,
                    beam_id=record["tx_beam_id"],
                    channels=NAMES,
                    rx_source_channels=np.arange(5, 15, dtype=np.uint16),
                    tx_channel=np.uint16(2),
                    delays=np.full(
                        tx_count,
                        int(record["range_start_gate"]) - int(record["tx_start_gate"]),
                        np.int32,
                    ),
                    metadata_key_offset_samples=np.int64(key - actual_key),
                    c_tx_idx=sample_index(data["tx_time_us"])
                    + int(record["tx_start_gate"]),
                    c_range_km_nominal=data["range_km_nominal"],
                    c_doppler_hz=data["doppler_hz"],
                    c_detection_statistic_db=data["detector_excess_over_median_db"],
                    c_alias_matched_power=data[
                        "alias_matched_power_at_observed_delay_doppler"
                    ],
                    cut_sha256=np.bytes_(
                        hashlib.file_digest(cut_path.open("rb"), "sha256").hexdigest()
                    ),
                )
                if not compact:
                    record.update(
                        zrx_echoes_re=iq[..., 0],
                        zrx_echoes_im=iq[..., 1],
                        ztx_pulses_re=tx[:, 0, :, 0],
                        ztx_pulses_im=tx[:, 0, :, 1],
                        ztx_history_re=tx[..., 0],
                        ztx_history_im=tx[..., 1],
                    )
                for field in [
                    "alias_rd_peak_power",
                    "alias_rd_peak_delay_samples_mod_ipp",
                    "alias_rd_peak_doppler_hz",
                ]:
                    if field in data:
                        record["c_" + field] = data[field]
                data.update(
                    event_id=np.bytes_(name),
                    metadata_key_offset_samples=np.int64(key - actual_key),
                )
                cut_writer.write(key, {k: [v] for k, v in record.items()})
                detection_writer.write(key, {k: [v] for k, v in data.items()})
        if first_key is None:
            raise ValueError("No events to publish")
        with h5py.File(geometry) as g:
            record = unpack(g)
            record.update(
                channels=NAMES,
                source_channel_indices=np.arange(5, 15, dtype=np.uint16),
                active_mask_confirmed=np.uint8(0),
                channel_map_hardware_confirmed=np.uint8(0),
            )
            writer("geometry").write(first_key, {k: [v] for k, v in record.items()})
        cal_path = run / "phasecal_boresight_provisional.h5"
        if cal_path.exists():
            with h5py.File(cal_path) as cal:
                record = unpack(cal)
                # Preserve a ten-channel calibration, never raw phase-corrected voltages.
                for k, v in list(record.items()):
                    if isinstance(v, np.ndarray) and v.shape and v.shape[0] == 15:
                        record[k] = v[5:]
                record.update(
                    channels=NAMES,
                    reference_source_channel=np.uint16(8),
                    reference_cut_channel=np.uint16(3),
                )
                record["reference_channel"] = np.uint32(3)
                writer("phasecal_provisional").write(
                    first_key, {k: [v] for k, v in record.items()}
                )
        for subdir in directory.iterdir():
            for path in subdir.rglob("*.h5"):
                if path.name != "dmd_properties.h5":
                    compress_metadata(path)
            reader = digital_rf.DigitalMetadataReader(str(subdir))
            lo, hi = reader.get_bounds()
            count = len(reader.read(lo, hi))
            print(subdir.name, count, "records")
            expected = (
                len(index["events"]) if subdir.name in ["cut", "detections"] else 1
            )
            assert count == expected
        directory.rename(root / "metadata")


def compress_metadata(path):
    temporary = path.with_suffix(".h5.compressed")

    def copy(source, destination):
        destination.attrs.update(source.attrs)
        for key, value in source.items():
            if isinstance(value, h5py.Group):
                copy(value, destination.create_group(key))
            else:
                data = value[()]
                filters = (
                    dict(compression="gzip", compression_opts=4, shuffle=True)
                    if value.ndim and value.size and value.dtype.kind != "O"
                    else {}
                )
                dataset = destination.create_dataset(key, data=data, **filters)
                dataset.attrs.update(value.attrs)

    with h5py.File(path) as source, h5py.File(temporary, "w") as destination:
        copy(source, destination)
    temporary.replace(path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cuts", type=Path)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument(
        "--geometry", type=Path, default=Path("outputs/interferometer/design.h5")
    )
    args = p.parse_args()
    publish(args.cuts, args.run, args.geometry)


if __name__ == "__main__":
    main()
