"""Check complete event products, energy/radius invariants, metadata, and PDF ink."""

import argparse
from pathlib import Path
import h5py
import numpy as np
import digital_rf
import pymupdf
from .coordinates import MU_EARTH


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    args = p.parse_args()
    events = sorted((args.run / "events").glob("event_*.h5"))
    counts = dict(
        events=len(events),
        detections=0,
        valid_local_covariance=0,
        phase_pairs_retained=0,
        upper_open_radius_profiles=0,
        same_code_pairs_retained=0,
    )
    for path in events:
        with h5py.File(path) as h:
            m = h["measurement"]
            r = h["result"]
            t = m["time_s"][()]
            assert np.all(np.diff(t) > 0)
            for key in [
                "range_m",
                "radial_velocity_mps",
                "amplitude",
                "noise_covariance",
            ]:
                assert np.all(np.isfinite(m[key][()])), (path, key)
            covariance = m["range_velocity_covariance"][()]
            assert np.all(np.linalg.eigvalsh(covariance) > 0)
            position = h["model/position_m"][()]
            model_range = np.linalg.norm(position, axis=1)
            keep = h["model/pulse_keep"][()].astype(bool)
            rms = lambda x: np.sqrt(np.mean(x**2))
            np.testing.assert_allclose(
                r["range_rms_m"][()], rms((model_range - m["range_m"][()])[keep])
            )
            uv = h["original_conditional_uv"][()]
            observed_xy = uv * m["range_m"][()][:, None]
            for j, key in enumerate(["horizontal_x_rms_m", "horizontal_y_rms_m"]):
                np.testing.assert_allclose(
                    r[key][()], rms((observed_xy - position[:, :2])[keep, j])
                )
            np.testing.assert_allclose(
                r["radial_velocity_rms_mps"][()],
                rms(
                    (h["model/radial_velocity_mps"][()] - m["radial_velocity_mps"][()])[
                        keep
                    ]
                ),
            )
            pp = h["pulse_phase"]
            selected = pp["keep"][()].astype(bool)
            if np.any(selected):
                previous, current = pp["previous"][()], pp["current"][()]
                prediction = (model_range[current] - model_range[previous]) / (
                    t[current] - t[previous]
                )
                np.testing.assert_allclose(
                    r["pulse_phase_velocity_rms_mps"][()],
                    rms((prediction - pp["velocity_mps"][()])[selected]),
                )
            assert np.all(np.isfinite(m["matched_filter_snr_db"][()]))
            radius = h["model/radius_m"][()]
            assert np.all(radius > 0) and np.all(np.diff(radius) <= 1e-14)
            assert np.isclose(radius[0], r["best_radius_m"][()])
            assert r["model_lookup_max_phase_error_rad"][()] < 0.002
            geo = r["preferred_geocentric_conditional"]
            state = geo["state_gcrs_m_mps"][()]
            speed = (
                np.sqrt(
                    np.sum(state[3:] ** 2) - 2 * MU_EARTH / np.linalg.norm(state[:3])
                )
                / 1000
            )
            np.testing.assert_allclose(
                speed, geo["geocentric_speed_km_s"][()], rtol=1e-12
            )
            n = len(t)
            counts["detections"] += n
            counts["valid_local_covariance"] += int(np.sum(m["refinement_valid"][()]))
            counts["phase_pairs_retained"] += int(r["selected_ppp_pairs"][()])
            counts["same_code_pairs_retained"] += int(
                r["selected_same_code_two_ipp_pairs"][()]
            )
            counts["upper_open_radius_profiles"] += int(r["upper_open"][()])
            for group in ["profile_single_pulse", "profile_with_pulse_phase"]:
                assert np.isfinite(np.min(h[group + "/chi2"][()]))
            with h5py.File(h.attrs["source"]) as cut:
                assert len(cut["detections/pulse_index"]) == n
            for suffix in ["", "_channels"]:
                assert (args.run / "figures" / f"{path.stem}{suffix}.pdf").exists()
    reader = digital_rf.DigitalMetadataReader(
        str(args.run / "metadata/meteor_trajectory")
    )
    lo, hi = reader.get_bounds()
    records = reader.read(lo, hi)
    assert len(records) == len(events)
    for record in records.values():
        assert record["celestial_orientation_confirmed"] == 0
        assert record["range_alias_order_confirmed"] == 0
        if record["upper_open"]:
            assert np.isnan(record["initial_radius_m"])
    atlas = pymupdf.open(args.run / "tibet_meteor_event_atlas.pdf")
    assert len(atlas) in (2 + len(events), 2 + 2 * len(events))
    ink = []
    for page in atlas:
        pix = page.get_pixmap(matrix=pymupdf.Matrix(0.2, 0.2), alpha=False)
        values = np.frombuffer(pix.samples, np.uint8)
        ink.append(float(np.mean(values < 245)))
        assert ink[-1] > 0.003
        text = page.get_text()
        assert "/Users/" not in text and ".py" not in text and "Rustplotlib" not in text
    with h5py.File(args.run / "verification.h5", "w") as h:
        for k, v in counts.items():
            h[k] = v
        h["atlas_pages"] = len(atlas)
        h["atlas_ink_fraction"] = ink
        h.attrs["scope"] = (
            "consistency and correctness checks; not calibrated uncertainty coverage or absolute celestial/range validation"
        )
    print("PASS", counts, "atlas pages", len(atlas))


if __name__ == "__main__":
    main()
