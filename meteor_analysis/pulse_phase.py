"""Nonoverlapping inter-pulse Doppler pairs with explicit ambiguity branches."""

import numpy as np


def measure_pairs(
    measurement,
    array,
    position,
    phasecal,
    pulse_indices,
    wavelength,
    state_covariance=None,
    beam_signature=None,
    pulse_stride=1,
):
    m = measurement
    distance = np.linalg.norm(position, axis=1)
    pattern = array(position[:, :2] / distance[:, None])
    covariance = m["noise_covariance"][3:, 3:]
    inverse = np.linalg.inv(covariance)
    common = []
    variance = []
    for i, a in enumerate(m["amplitude"]):
        # Seven geometrically disjoint large sections; the overlapping small
        # groups remain in angular fitting, not counted as independent power here.
        steering = abs(a[3:]) * np.exp(1j * (np.angle(pattern[i, 3:]) + phasecal[3:]))
        w = inverse @ steering
        gain = float(np.vdot(steering, w).real)
        common.append(np.vdot(w, a[3:]) / gain)
        variance.append(1 / gain)
    common = np.array(common)
    variance = np.array(variance)
    result = []
    last_used = -1
    for i in range(len(common) - 1):
        j = int(np.searchsorted(pulse_indices, pulse_indices[i] + pulse_stride))
        if (
            i <= last_used
            or j >= len(pulse_indices)
            or pulse_indices[j] - pulse_indices[i] != pulse_stride
        ):
            continue
        if beam_signature is not None and not np.array_equal(
            beam_signature[i], beam_signature[j]
        ):
            continue
        dt = m["time_s"][j] - m["time_s"][i]
        beat = common[j] * common[i].conj()
        angle = np.angle(beat)
        wrapped = -wavelength * angle / (4 * np.pi * dt)
        spacing = wavelength / (2 * dt)
        anchor = (distance[j] - distance[i]) / dt
        branch = int(np.rint((anchor - wrapped) / spacing))
        velocity = wrapped + branch * spacing
        phasevar = variance[i] / (2 * abs(common[i]) ** 2) + variance[j] / (
            2 * abs(common[j]) ** 2
        )
        sigma = wavelength * np.sqrt(phasevar) / (4 * np.pi * dt)
        coarse = np.mean(m["radial_velocity_mps"][[i, j]])
        coarse_sigma = np.sqrt(np.sum(m["range_velocity_covariance"][[i, j], 1, 1])) / 2
        if state_covariance is None:
            anchor_sigma = coarse_sigma
        else:
            ui = position[i] / distance[i]
            uj = position[j] / distance[j]
            gradient = (
                1000
                * np.r_[
                    (uj - ui) / dt, (m["time_s"][j] * uj - m["time_s"][i] * ui) / dt
                ]
            )
            anchor_sigma = np.sqrt(
                max(float(gradient @ state_covariance @ gradient), 0.0)
            )
        margin = spacing / 2 - abs(velocity - anchor)
        # These are quality/branch flags, not a calibrated false-alarm probability.
        unique = margin > 3 * anchor_sigma
        result.append(
            (
                i,
                j,
                (m["time_s"][i] + m["time_s"][j]) / 2,
                angle,
                wrapped,
                spacing,
                branch,
                velocity,
                sigma,
                anchor,
                coarse,
                unique,
                anchor_sigma,
            )
        )
        last_used = j
    if not result:
        return {
            k: np.array(
                [],
                dtype=(
                    bool
                    if k in ["keep", "branch_unique"]
                    else (int if k in ["previous", "current", "branch"] else float)
                ),
            )
            for k in [
                "previous",
                "current",
                "time_s",
                "phase_rad",
                "wrapped_velocity_mps",
                "branch_spacing_mps",
                "branch",
                "velocity_mps",
                "formal_sigma_mps",
                "model_anchor_mps",
                "single_pulse_anchor_mps",
                "branch_unique",
                "sigma_mps",
                "keep",
                "trajectory_anchor_sigma_mps",
                "empirical_model_discrepancy_mps",
                "adjacent_branch_minus_mps",
                "adjacent_branch_plus_mps",
            ]
        }
    result = np.asarray(result)
    keys = [
        "previous",
        "current",
        "time_s",
        "phase_rad",
        "wrapped_velocity_mps",
        "branch_spacing_mps",
        "branch",
        "velocity_mps",
        "formal_sigma_mps",
        "model_anchor_mps",
        "single_pulse_anchor_mps",
        "branch_unique",
        "trajectory_anchor_sigma_mps",
    ]
    pairs = {k: result[:, i] for i, k in enumerate(keys)}
    for k in ["previous", "current", "branch"]:
        pairs[k] = pairs[k].astype(int)
    residual = pairs["velocity_mps"] - pairs["model_anchor_mps"]
    centre = np.median(residual)
    scatter = 1.4826 * np.median(abs(residual - centre))
    pairs["sigma_mps"] = np.sqrt(pairs["formal_sigma_mps"] ** 2 + scatter**2)
    pairs["empirical_model_discrepancy_mps"] = np.full(len(result), scatter)
    pairs["keep"] = (pairs["branch_unique"] > 0) & (
        abs(residual) < np.maximum(200.0, 4 * pairs["sigma_mps"])
    )
    pairs["adjacent_branch_minus_mps"] = (
        pairs["velocity_mps"] - pairs["branch_spacing_mps"]
    )
    pairs["adjacent_branch_plus_mps"] = (
        pairs["velocity_mps"] + pairs["branch_spacing_mps"]
    )
    return pairs
