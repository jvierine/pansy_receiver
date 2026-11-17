"""Nonoverlapping inter-pulse Doppler pairs with explicit ambiguity branches."""

import numpy as np


def _doppler_reference(measurement, pair_times):
    """Independent within-pulse Doppler anchors; curvature is not forced to zero."""
    t = measurement["time_s"]
    centre = np.mean(t)
    span = max(np.ptp(t), 1e-6)
    degree = min(3, max(0, len(t) // 6 - 1))
    x = np.polynomial.polynomial.polyvander((t - centre) / span, degree)
    target = np.polynomial.polynomial.polyvander((pair_times - centre) / span, degree)
    sigma = np.sqrt(measurement["range_velocity_covariance"][:, 1, 1])
    robust = np.ones(len(t))
    for _ in range(5):
        weight = np.sqrt(robust) / sigma
        coefficients = np.linalg.lstsq(
            x * weight[:, None], measurement["radial_velocity_mps"] * weight, rcond=None
        )[0]
        residual = (measurement["radial_velocity_mps"] - x @ coefficients) / sigma
        scale = max(1.0, 1.4826 * np.median(abs(residual - np.median(residual))))
        robust = np.minimum(1.0, 2.5 * scale / np.maximum(abs(residual), 1e-12))
    covariance = (
        np.linalg.pinv((x * weight[:, None]).T @ (x * weight[:, None])) * scale**2
    )
    uncertainty = np.sqrt(
        np.maximum(np.einsum("ni,ij,nj->n", target, covariance, target), 0.0)
    )
    return target @ coefficients, uncertainty


def _resolve_event_branches(pairs, measurement):
    """Minimum-cost continuous branch paths, anchored to independent Doppler data."""
    t = pairs["time_s"]
    wrapped = pairs["wrapped_velocity_mps"]
    spacing = pairs["branch_spacing_mps"]
    reference, sigma_reference = _doppler_reference(measurement, t)
    centre = np.rint((reference - wrapped) / spacing).astype(int)
    candidates = centre[:, None] + np.arange(-2, 3)[None, :]
    velocity = wrapped[:, None] + candidates * spacing[:, None]
    variance = sigma_reference**2 + pairs["formal_sigma_mps"] ** 2
    local = (velocity - reference[:, None]) ** 2 / np.maximum(variance[:, None], 1.0)
    forward = local.copy()
    parents = np.zeros(local.shape, int)
    transitions = []
    gap_limit = max(0.02, 5 * np.median(np.diff(t))) if len(t) > 1 else 0.02
    for j in range(1, len(t)):
        if t[j] - t[j - 1] > gap_limit:
            transition = np.zeros((5, 5))
        else:
            change = (velocity[j, None, :] - velocity[j - 1, :, None]) - (
                reference[j] - reference[j - 1]
            )
            variance_change = (
                pairs["formal_sigma_mps"][j] ** 2
                + pairs["formal_sigma_mps"][j - 1] ** 2
                + 50.0**2
            )
            # Bounded influence rejects a bad phase pair rather than moving the whole path.
            transition = 2 * (np.sqrt(1 + change**2 / variance_change) - 1)
        transitions.append(transition)
        cost = forward[j - 1, :, None] + transition
        parents[j] = np.argmin(cost, axis=0)
        forward[j] += np.min(cost, axis=0)
    selected = np.zeros(len(t), int)
    selected[-1] = np.argmin(forward[-1])
    for j in range(len(t) - 1, 0, -1):
        selected[j - 1] = parents[j, selected[j]]
    backward = np.zeros_like(local)
    for j in range(len(t) - 2, -1, -1):
        backward[j] = np.min(
            transitions[j] + local[j + 1, None, :] + backward[j + 1, None, :], axis=1
        )
    total = forward + backward
    chosen = total[np.arange(len(t)), selected]
    total[np.arange(len(t)), selected] = np.inf
    margin = np.min(total, axis=1) - chosen
    pairs["branch"] = candidates[np.arange(len(t)), selected]
    pairs["velocity_mps"] = velocity[np.arange(len(t)), selected]
    pairs["doppler_reference_mps"] = reference
    pairs["doppler_reference_sigma_mps"] = sigma_reference
    pairs["branch_alternative_cost_gap"] = margin
    pairs["branch_unique"] = margin > 9.0
    residual = pairs["velocity_mps"] - reference
    differences = np.diff(residual)
    scatter = (
        1.4826 * np.median(abs(differences - np.median(differences))) / np.sqrt(2.0)
        if len(differences)
        else 0.0
    )
    pairs["sigma_mps"] = np.sqrt(pairs["formal_sigma_mps"] ** 2 + scatter**2)
    pairs["empirical_model_discrepancy_mps"] = np.full(len(t), scatter)
    pairs["keep"] = pairs["branch_unique"] & (
        abs(residual)
        < np.maximum(200.0, 4 * np.sqrt(pairs["sigma_mps"] ** 2 + sigma_reference**2))
    )


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
                "doppler_reference_mps",
                "doppler_reference_sigma_mps",
                "branch_alternative_cost_gap",
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
    pairs["independent_nearest_cv_velocity_mps"] = pairs["velocity_mps"].copy()
    _resolve_event_branches(pairs, m)
    pairs["adjacent_branch_minus_mps"] = (
        pairs["velocity_mps"] - pairs["branch_spacing_mps"]
    )
    pairs["adjacent_branch_plus_mps"] = (
        pairs["velocity_mps"] + pairs["branch_spacing_mps"]
    )
    return pairs
