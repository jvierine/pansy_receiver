"""Optional standard Digital Metadata publication of completed event fits."""

import argparse
from pathlib import Path
import h5py
import numpy as np
import digital_rf


def publish(run, sample_rate=500000):
    destination = run / "metadata" / "meteor_trajectory"
    if destination.exists():
        raise FileExistsError(destination)
    destination.mkdir(parents=True)
    writer = digital_rf.DigitalMetadataWriter(
        str(destination), 3600, 60, sample_rate, 1, "meteor_trajectory"
    )
    used = set()
    for path in sorted((run / "events").glob("event_*.h5")):
        with h5py.File(path) as h:
            m = h["measurement"]
            r = h["result"]
            g = r["preferred_geocentric_conditional"]
            epoch = m["epoch_unix"][()]
            actual = int(round(float(epoch[0]) * sample_rate))
            key = actual
            while key in used:
                key += 1
            used.add(key)
            record = dict(
                event_id=np.bytes_(path.stem),
                epoch_unix=epoch,
                time_s=m["time_s"][()],
                position_array_m=h["model/position_m"][()],
                velocity_array_mps=h["model/velocity_mps"][()],
                measured_range_m=m["range_m"][()],
                single_pulse_radial_velocity_mps=m["radial_velocity_mps"][()],
                geocentric_speed_km_s=g["geocentric_speed_km_s"][()],
                ra_deg_illustrative=g["ra_deg"][()],
                dec_deg_illustrative=g["dec_deg"][()],
                array_radiant_azimuth_deg=r["array_radiant_azimuth_deg"][()],
                array_radiant_zenith_angle_deg=r["array_radiant_zenith_angle_deg"][()],
                initial_radius_m=(
                    np.nan if bool(r["upper_open"][()]) else r["best_radius_m"][()]
                ),
                profile_best_grid_radius_m=r["best_radius_m"][()],
                profile_lower_radius_m=r["profile_lower_radius_m"][()],
                profile_upper_radius_m=r["profile_upper_radius_m"][()],
                upper_open=r["upper_open"][()],
                phase_velocity_time_s=h["pulse_phase/time_s"][()],
                phase_velocity_mps=h["pulse_phase/velocity_mps"][()],
                phase_pair_keep=h["pulse_phase/keep"][()],
                celestial_orientation_confirmed=np.uint8(0),
                range_alias_order_confirmed=np.uint8(0),
                metadata_key_offset_samples=np.int64(key - actual),
                status=np.bytes_(
                    "conditional: nominal range branch, provisional phase calibration, unconfirmed geographic axes; working radius profiles"
                ),
            )
            writer.write(key, {k: [v] for k, v in record.items()})
    reader = digital_rf.DigitalMetadataReader(str(destination))
    lo, hi = reader.get_bounds()
    assert len(reader.read(lo, hi)) == len(used)
    print("Verified", len(used), "derived Digital Metadata records")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    a = p.parse_args()
    publish(a.run)


if __name__ == "__main__":
    main()
