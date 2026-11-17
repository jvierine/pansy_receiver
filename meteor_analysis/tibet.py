"""Offline Tibet cut analysis using reusable site/DSP/trajectory modules."""

import argparse
from pathlib import Path
import datetime as dt
import hashlib
import time
import numpy as np
import h5py
from scipy.integrate import solve_ivp
from .stations import TIBET_MST
from .atmosphere import build_profile
from .measurements import load_tibet_cut, refine, C
from .geometry import ArrayModel
from .trajectory import EventFit, profile_bounds, radial
from .pulse_phase import measure_pairs
from .coordinates import radiant, configure_eop, load_eop
from .plot_data import range_time_power


def write(group, values):
    for k, v in values.items():
        if isinstance(v, dict):
            write(group.create_group(k), v)
        elif isinstance(v, str):
            group[k] = np.bytes_(v)
        elif v is not None:
            group[k] = v


def above_atmosphere(state, radius, atmosphere, station, target=120000.0):
    alt, rho = atmosphere

    def height(y):
        return station.altitude_m + y[2] + (y[0] ** 2 + y[1] ** 2) / (2 * 6371000.0)

    def rhs(t, y):
        h = height(y)
        air = np.exp(np.interp(h, alt, np.log(rho)))
        speed = np.linalg.norm(y[3:6])
        drag = -0.75 * air / (3000 * y[6]) * speed * y[3:6]
        drag[2] -= 9.80665 * (6371000 / (6371000 + h)) ** 2
        return np.r_[y[3:6], drag, -air * 1e-8 * speed**3 / (4 * 3000)]

    def top(t, y):
        return height(y) - target

    top.terminal = True
    solution = solve_ivp(
        rhs,
        (0.0, -15.0),
        np.r_[state, radius],
        events=top,
        rtol=1e-8,
        atol=np.r_[np.full(6, 1e-6), 1e-12],
        max_step=0.01,
    )
    if not solution.success or not len(solution.t_events[0]):
        raise ValueError("120 km backward extrapolation failed")
    if np.linalg.norm(solution.y[3:6, -1]) > 100000:
        raise ValueError("back-extrapolated speed exceeds 100 km/s")
    return solution.t[-1], solution.y[:, -1]


def analyze(name, args, array, atmosphere, phase, calstd):
    cut = load_tibet_cut(args.cuts / "cuts" / f"{name}.h5")
    m = refine(cut)
    with h5py.File(args.run / "interferometry_provisional.h5") as h:
        uv = h["events"][name]["conditional_uv_within_4deg_box"][()]
    fit = EventFit(m, array, TIBET_MST, atmosphere, phase, calstd)
    cv = fit.initialize(uv)
    pos, vel, _ = fit.model(cv.x)
    cv_cov = np.linalg.pinv(cv.jac.T @ cv.jac) * max(
        1.0, np.sum(cv.fun**2) / max(len(cv.fun) - len(cv.x), 1)
    )
    beam_signature = np.column_stack(
        [cut["tx_alpha"][m["raw_rows"]], cut["tx_beta"][m["raw_rows"]]]
    )
    ppp = measure_pairs(
        m,
        array,
        pos,
        phase,
        cut["detections"]["pulse_index"],
        C / cut["frequency"],
        cv_cov[:6, :6],
        beam_signature,
    )
    ppp_two = measure_pairs(
        m,
        array,
        pos,
        phase,
        cut["detections"]["pulse_index"],
        C / cut["frequency"],
        cv_cov[:6, :6],
        beam_signature,
        pulse_stride=2,
    )
    fit.ppp = ppp
    radii = np.geomspace(1e-6, 0.1, args.radius_points)
    base = fit.profile(cv.x, radii, use_ppp=False)
    combined = fit.profile(cv.x, radii, use_ppp=True)
    # Refine around both minima without changing the global radius domain.
    extra = []
    for profile in [base, combined]:
        j = int(np.argmin(profile["chi2"]))
        lo = max(0, j - 1)
        hi = min(len(radii) - 1, j + 1)
        extra.extend(np.geomspace(radii[lo], radii[hi], 13))
    refined = np.unique(np.r_[radii, extra])
    base = fit.profile(
        base["parameters"][np.argmin(base["chi2"])], refined, use_ppp=False
    )
    combined = fit.profile(
        combined["parameters"][np.argmin(combined["chi2"])], refined, use_ppp=True
    )
    radii = refined
    k = int(np.argmin(combined["chi2"]))
    best = combined["parameters"][k]
    radius = float(radii[k])
    position, velocity, rr = fit.model(best, radius)
    geometry_error = array.check(
        position[:: max(1, len(position) // 10), :2]
        / np.linalg.norm(position[:: max(1, len(position) // 10)], axis=1)[:, None]
    )
    result = profile_bounds(combined)
    keep = fit.keep
    observed_position = (
        np.column_stack([uv, np.sqrt(np.maximum(1 - np.sum(uv**2, axis=1), 0))])
        * m["range_m"][:, None]
    )
    position_residual = observed_position - position
    model_range = np.linalg.norm(position, axis=1)
    pair_model = (model_range[ppp["current"]] - model_range[ppp["previous"]]) / (
        m["time_s"][ppp["current"]] - m["time_s"][ppp["previous"]]
    )
    rms = lambda values: float(np.sqrt(np.mean(np.asarray(values) ** 2)))
    phase_rms = (
        rms((pair_model - ppp["velocity_mps"])[ppp["keep"]])
        if np.any(ppp["keep"])
        else np.nan
    )
    result.update(
        initial_radius_epoch="first retained scattering measurement, not pre-atmospheric radius",
        material_density_kg_m3=3000.0,
        ablation_sigma_kg_j=1e-8,
        speed_first_km_s=np.linalg.norm(velocity[0]) / 1000,
        initial_mass_kg=4 * np.pi / 3 * 3000 * radius**3,
        model_lookup_max_phase_error_rad=geometry_error,
        range_rms_m=rms((model_range - m["range_m"])[keep]),
        horizontal_x_rms_m=rms(position_residual[keep, 0]),
        horizontal_y_rms_m=rms(position_residual[keep, 1]),
        horizontal_position_rms_m=rms(
            np.linalg.norm(position_residual[keep, :2], axis=1)
        ),
        height_rms_m=rms(position_residual[keep, 2]),
        radial_velocity_rms_mps=rms(
            (radial(position, velocity) - m["radial_velocity_mps"])[keep]
        ),
        pulse_phase_velocity_rms_mps=phase_rms,
        rms_definition="unweighted sqrt(mean(residual squared)); single-pulse and accepted nonoverlapping phase pairs reported separately; horizontal residuals use the plotted conditional AoA positions",
        selected_ppp_pairs=int(np.sum(ppp["keep"])),
        calibration_and_orientation_confirmed=False,
        reference_frame="array x/y/up; geographic ENU orientation unconfirmed",
        nominal_range_alias_order=0,
        even_ipp_range_alias_resolved=False,
    )
    if len(ppp_two["time_s"]) and len(ppp["time_s"]):
        difference = ppp_two["velocity_mps"] - np.interp(
            ppp_two["time_s"], ppp["time_s"], ppp["velocity_mps"]
        )
        result["one_vs_two_ipp_median_velocity_difference_mps"] = float(
            np.median(difference)
        )
        result["one_vs_two_ipp_velocity_difference_rms_mps"] = float(np.std(difference))
    result["selected_same_code_two_ipp_pairs"] = int(np.sum(ppp_two["keep"]))
    result["radius_identification_status"] = (
        "open upper profile limit"
        if result["upper_open"]
        else (
            "bounded conditional profile; same-code phase check available"
            if np.sum(ppp_two["keep"]) >= 3
            else "bounded conditional profile; insufficient same-code phase validation"
        )
    )
    plotting = range_time_power(cut, m, velocity)
    response = array(position[:, :2] / np.linalg.norm(position, axis=1)[:, None])
    observed = np.angle(m["amplitude"] * np.exp(-1j * phase[None, :]))
    common = np.angle(
        m["amplitude"][:, 3] * np.exp(-1j * phase[3]) * response[:, 3].conj()
    )
    plotting.update(
        channel_phase_rad=np.angle(np.exp(1j * (observed - common[:, None]))),
        model_channel_phase_rad=np.angle(
            response * np.exp(1j * np.r_[best[6:9], 0.0, best[9:]][None, :])
        ),
    )
    result["conditional_speed_profile_range_km_s"] = [
        np.min(
            np.linalg.norm(
                combined["parameters"][
                    combined["chi2"] <= np.min(combined["chi2"]) + 3.84, 3:6
                ],
                axis=1,
            )
        ),
        np.max(
            np.linalg.norm(
                combined["parameters"][
                    combined["chi2"] <= np.min(combined["chi2"]) + 3.84, 3:6
                ],
                axis=1,
            )
        ),
    ]
    # This is an illustrative ENU conversion, not a verified astronomical radiant.
    first_rad = radiant(TIBET_MST, np.r_[position[0], velocity[0]], m["epoch_unix"][0])
    result["geocentric_first_detected_conditional"] = first_rad
    orientation_speeds = []
    orientation_radiants = []
    for angle in np.arange(0, 360, 10):
        a = np.radians(angle)
        rotation = np.array(
            [[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]]
        )
        r = radiant(
            TIBET_MST,
            np.r_[rotation @ position[0], rotation @ velocity[0]],
            m["epoch_unix"][0],
        )
        orientation_speeds.append(r["geocentric_speed_km_s"])
        orientation_radiants.append(
            [r["sun_centered_ecliptic_longitude_deg"], r["ecliptic_latitude_deg"]]
        )
    result["unknown_axis_rotation_geocentric_speed_range_km_s"] = [
        min(orientation_speeds),
        max(orientation_speeds),
    ]
    result["array_radiant_azimuth_deg"] = (
        np.degrees(np.arctan2(-velocity[0, 1], -velocity[0, 0])) % 360
    )
    result["array_radiant_zenith_angle_deg"] = np.degrees(
        np.arccos(-velocity[0, 2] / np.linalg.norm(velocity[0]))
    )
    result["physical_speed_valid"] = bool(11000 < np.linalg.norm(velocity[0]) < 100000)
    try:
        top_time, top = above_atmosphere(
            np.r_[position[0], velocity[0]], radius, atmosphere, TIBET_MST
        )
        result["at_120km_conditional"] = radiant(
            TIBET_MST, top[:6], m["epoch_unix"][0] + top_time
        )
        result["preferred_geocentric_conditional"] = result["at_120km_conditional"]
        result["preferred_geocentric_epoch"] = (
            "back-extrapolated 120 km state; fixed MSIS/material/ablation model"
        )
        result["at_120km_radius_m"] = top[6]
        result["at_120km_time_offset_s"] = top_time
    except ValueError as exc:
        result["at_120km_status"] = str(exc)
        result["preferred_geocentric_conditional"] = first_rad
        result["preferred_geocentric_epoch"] = (
            "first observed state; atmospheric deceleration not removed"
        )
    with h5py.File(args.output / "events" / f"{name}.h5", "w") as h:
        write(
            h,
            dict(
                measurement=m,
                calibration=dict(phase_rad=phase, ensemble_phase_std_rad=calstd),
                constant_velocity=dict(
                    parameters=cv.x,
                    position_m=pos,
                    velocity_mps=vel,
                    conditional_parameter_covariance=cv_cov,
                ),
                model=dict(
                    position_m=position,
                    velocity_mps=velocity,
                    radius_m=rr,
                    radial_velocity_mps=radial(position, velocity),
                    parameters=best,
                    pulse_keep=fit.keep,
                    range_velocity_working_scale=fit.rv_scale,
                    phase_working_scale=fit.phase_scale,
                ),
                pulse_phase=ppp,
                pulse_phase_same_code_two_ipp=ppp_two,
                profile_single_pulse=base,
                profile_with_pulse_phase=combined,
                plotting=plotting,
                result=result,
                unknown_orientation_radiant_locus=np.array(orientation_radiants),
                original_conditional_uv=uv,
                cut_source_sha256=hashlib.sha256(
                    cut["cut_path"].read_bytes()
                ).hexdigest(),
            ),
        )
        h.attrs.update(
            status="PROVISIONAL: unknown array orientation/sign, boresight ensemble calibration, nominal alias order; conditional radius profile",
            source=str(cut["cut_path"]),
            station_name=TIBET_MST.name,
            latitude_deg=TIBET_MST.latitude_deg,
            longitude_deg=TIBET_MST.longitude_deg,
            altitude_m=TIBET_MST.altitude_m,
        )
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cuts", type=Path)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--geometry", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--max-events", type=int, default=100000)
    p.add_argument("--radius-points", type=int, default=25)
    p.add_argument(
        "--environment-from",
        type=Path,
        help="reuse saved observation-day MSIS/EOP HDF5 inputs",
    )
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    (args.output / "events").mkdir(parents=True)
    with h5py.File(args.cuts / "index.h5") as h:
        names = sorted(h["events"])[: args.max_events]
    sample = load_tibet_cut(args.cuts / "cuts" / f"{names[0]}.h5")
    if args.environment_from:
        import shutil

        with h5py.File(args.environment_from / "atmosphere.h5") as h:
            when = dt.datetime.fromisoformat(h.attrs["utc"]).timestamp()
            if (
                abs(when - float(sample["ipp_time_us"][0]) / 1e6) > 3600
                or h.attrs["latitude_deg"] != TIBET_MST.latitude_deg
                or h.attrs["longitude_deg"] != TIBET_MST.longitude_deg
            ):
                raise ValueError(
                    "cached atmosphere does not match the observation time/site"
                )
            atmosphere = (h["altitude_m"][()], h["density_kg_m3"][()])
        for name in ["earth_orientation.h5", "atmosphere.h5"]:
            shutil.copy2(args.environment_from / name, args.output / name)
        load_eop(args.output / "earth_orientation.h5")
    else:
        configure_eop(
            float(sample["ipp_time_us"][0]) / 1e6, args.output / "earth_orientation.h5"
        )
        atmosphere = build_profile(
            TIBET_MST,
            float(sample["ipp_time_us"][0]) / 1e6,
            args.output / "atmosphere.h5",
        )
    array = ArrayModel(args.geometry, C / sample["frequency"])
    with h5py.File(args.run / "phasecal_boresight_provisional.h5") as h:
        phase = h["phasecal_rad"][()][5:]
        concentration = h["circular_concentration"][()][5:]
        calstd = np.sqrt(-2 * np.log(np.clip(concentration, 1e-6, 1))) / np.sqrt(
            float(h["selected_events"][()])
        )
    with h5py.File(args.output / "catalogue.h5", "w") as catalog:
        events = catalog.create_group("events")
        for name in names:
            started = time.monotonic()
            try:
                result = analyze(name, args, array, atmosphere, phase, calstd)
                write(events.create_group(name), result)
                print(
                    name,
                    "v=",
                    round(result["speed_first_km_s"], 2),
                    "km/s, r0=",
                    round(result["best_radius_m"] * 1e6, 1),
                    "um, upper open=",
                    result["upper_open"],
                    "PPP=",
                    result["selected_ppp_pairs"],
                    "elapsed",
                    round(time.monotonic() - started, 1),
                    flush=True,
                )
            except Exception as exc:
                g = events.create_group(name)
                g["failure"] = np.bytes_(f"{type(exc).__name__}: {exc}")
                import traceback

                traceback.print_exc()
            catalog.flush()
        catalog.attrs.update(
            complete=True,
            coordinate_status="array frame is unconfirmed relative to ENU; celestial conversion is illustrative only",
        )


if __name__ == "__main__":
    main()
