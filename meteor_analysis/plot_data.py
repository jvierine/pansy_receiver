"""Scientific plotting arrays derived from the retained measured voltages."""

import numpy as np


def range_time_power(cut, measurement, model_velocity):
    width = cut["rx"].shape[-1]
    length = cut["tx"].shape[-1]
    t = (
        cut["ipp_time_us"].astype(np.int64)
        - int(cut["ipp_time_us"][measurement["raw_rows"][0]])
    ) / 1e6
    fd = np.interp(t, measurement["time_s"], measurement["doppler_hz"])
    energy = abs(cut["tx"]) ** 2
    centre = np.sum(np.arange(length)[None, :] * energy, axis=1) / np.sum(
        energy, axis=1
    )
    template = cut["tx"].conj() / np.sqrt(np.sum(energy, axis=1))[:, None]
    template *= np.exp(
        -2j
        * np.pi
        * fd[:, None]
        * (np.arange(length)[None, :] - centre[:, None])
        / cut["fs"]
    )
    windows = np.lib.stride_tricks.sliding_window_view(
        cut["rx"][:, 3:, :], length, axis=-1
    )
    amp = np.einsum("ncgk,nk->ncg", windows, template, optimize=True)
    inverse = np.linalg.inv(measurement["noise_covariance"][3:, 3:])
    power = (
        np.einsum("nig,ij,njg->ng", amp.conj(), inverse, amp, optimize=True).real / 7
    )
    # Vertices describe full cells; no plot interpolation creates extra samples.
    time_edges = np.r_[
        t[0] - (t[1] - t[0]) / 2, (t[:-1] + t[1:]) / 2, t[-1] + (t[-1] - t[-2]) / 2
    ]
    dr = 299792458 / (2 * cut["fs"]) / 1000
    range_edges = (
        cut["gate0"] - cut["tx_gate"] + np.arange(width - length + 2) - 0.5
    ) * dr
    return dict(
        time_edges_s=time_edges,
        range_edges_km=range_edges,
        power_over_noise_db=10 * np.log10(np.maximum(power.T, 1e-3)),
    )
