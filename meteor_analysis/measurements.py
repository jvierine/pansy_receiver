"""Lossless Tibet cut adapter and common-reference range/Doppler refinement."""

from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import h5py
from .native import matched

C = 299792458.0


def load_tibet_cut(path):
    with h5py.File(path) as h:
        iq = h["zrx_echoes_iq_i32"][()]
        tx = h["ztx_history_iq_i32"][:, 0]
        d = {k: v[()] for k, v in h["detections"].items()}
        return dict(
            rx=iq[..., 0].astype(float) + 1j * iq[..., 1],
            tx=tx[..., 0].astype(float) + 1j * tx[..., 1],
            pulse=h["pulse_index"][()],
            ipp_time_us=h["ipp_start_time_us"][()],
            tx_alpha=h["tx_alpha_per_pulse"][()],
            tx_beta=h["tx_beta_per_pulse"][()],
            gate0=int(h["range_start_gate"][()]),
            tx_gate=int(h["tx_start_gate"][()]),
            fs=float(h["sample_rate_hz"][()]),
            frequency=float(h["frequency_hz"][()]),
            detections=d,
            cut_path=Path(path),
        )


def noise_covariance(cut):
    # Only range margins, outside the echo code and its nearby sidelobes.
    samples = (
        np.concatenate([cut["rx"][:, :, :24], cut["rx"][:, :, -24:]], axis=-1)
        .transpose(0, 2, 1)
        .reshape(-1, 10)
    )
    energy = np.sum(abs(samples) ** 2, axis=1)
    samples = samples[energy <= np.quantile(energy, 0.95)]
    samples = samples - np.mean(samples, axis=0)
    cov = samples.T @ samples.conj() / len(samples)
    eigen, vec = np.linalg.eigh(cov)
    floor = np.max(eigen) * 1e-4
    cov = (vec * np.maximum(eigen, floor)) @ vec.conj().T
    return cov, len(samples)


def refine(cut):
    d = cut["detections"]
    rows = np.searchsorted(cut["pulse"], d["pulse_index"])
    cov, nsamples = noise_covariance(cut)
    inverse = np.linalg.inv(cov[3:, 3:])
    output = []
    amplitudes = []
    joint = []
    success = []
    reference_ipp = int(cut["ipp_time_us"][rows[0]])
    for j, row in enumerate(rows):
        gate = float(d["delay_samples"][j]) + cut["tx_gate"] - cut["gate0"]
        fd = float(d["doppler_hz"][j])
        rx = cut["rx"][row]
        tx = cut["tx"][row]

        def objective(x):
            a = matched(rx[3:], tx, gate + x[0], fd + x[1] * 1000, cut["fs"])
            return -2 * float(np.vdot(a, inverse @ a).real)

        fit = minimize(
            objective,
            [0.0, 0.0],
            method="L-BFGS-B",
            bounds=[(-1.5, 1.5), (-8.0, 8.0)],
            options={"maxiter": 45, "ftol": 1e-11, "gtol": 1e-5},
        )
        x = fit.x
        a = matched(rx, tx, gate + x[0], fd + x[1] * 1000, cut["fs"])
        steps = np.array([0.01, 0.025])
        hes = np.empty((2, 2))
        f0 = objective(x)
        for k in range(2):
            e = np.eye(2)[k] * steps[k]
            hes[k, k] = (objective(x + e) - 2 * f0 + objective(x - e)) / steps[k] ** 2
        e = np.array([steps[0], 0.0])
        v = np.array([0.0, steps[1]])
        hes[0, 1] = hes[1, 0] = (
            objective(x + e + v)
            - objective(x + e - v)
            - objective(x - e + v)
            + objective(x - e - v)
        ) / (4 * np.prod(steps))
        valid = (
            np.all(np.linalg.eigvalsh(hes) > 0)
            and abs(x[0]) < 1.49
            and abs(x[1]) < 7.99
        )
        covariance = 2 * np.linalg.inv(hes) if valid else np.diag([0.5**2, 2.0**2])
        scale = np.diag([C / (2 * cut["fs"]), -C / cut["frequency"] / 2 * 1000])
        joint.append(scale @ covariance @ scale.T)
        refined_delay = float(d["delay_samples"][j]) + x[0]
        energy = abs(tx) ** 2
        centre = float(np.sum(np.arange(len(tx)) * energy) / sum(energy))
        # Scattering epoch: receive sample centre minus the one-way light time.
        light = refined_delay / (2 * cut["fs"])
        relative = (
            (int(cut["ipp_time_us"][row]) - reference_ipp) / 1e6
            + (gate + x[0] + cut["gate0"] + centre) / cut["fs"]
            - light
        )
        output.append(
            [
                relative,
                refined_delay * C / (2 * cut["fs"]),
                -(fd + x[1] * 1000) * C / cut["frequency"] / 2,
                fd + x[1] * 1000,
            ]
        )
        amplitudes.append(a)
        success.append(valid)
    output = np.asarray(output)
    return dict(
        epoch_unix=reference_ipp / 1e6 + output[:, 0],
        time_s=output[:, 0] - output[0, 0],
        range_m=output[:, 1],
        radial_velocity_mps=output[:, 2],
        doppler_hz=output[:, 3],
        amplitude=np.array(amplitudes),
        range_velocity_covariance=np.array(joint),
        refinement_valid=np.array(success),
        noise_covariance=cov,
        noise_samples=nsamples,
        raw_rows=rows,
        detector_statistic_db=d["detector_excess_over_median_db"],
        original_delay_samples=d["delay_samples"],
        original_doppler_hz=d["doppler_hz"],
    )
