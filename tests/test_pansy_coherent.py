from pathlib import Path

import h5py
import numpy as np
import pytest
import digital_rf as drf

import pansy_config as pc
from pansy_coherent import (
    aoa_steering_vector,
    coherent_add_modules,
    enu_line_of_sight_to_arrival_uvw,
    estimate_event_module_voltage_gain,
)


SAMPLE_IDX = 1739754984380607
REPO = Path(__file__).resolve().parents[1]


def test_synthetic_coherent_addition_and_enu_conversion():
    position = np.array([3.0, 4.0, 12.0])
    arrival = enu_line_of_sight_to_arrival_uvw(position)
    np.testing.assert_allclose(arrival, [3.0 / 13.0, 4.0 / 13.0, -12.0 / 13.0])

    antpos = np.array([[0.0, 0.0, 0.0], [2.0, -1.0, 0.0], [-3.0, 4.0, 0.0]])
    phasecal = np.array([0.2, -0.4, 0.7])
    steering = aoa_steering_vector(
        arrival, 0, phasecal=phasecal, antpos_m=antpos, wavelength_m=6.0
    )
    signal = np.array([1.0 + 2.0j, -0.5 + 0.25j])
    gain = np.array([0.8, 1.1, 1.4])
    received = gain[:, None] * np.conj(steering)[:, None] * signal[None, :]
    summed = coherent_add_modules(
        received,
        arrival,
        0,
        module_voltage_gain=gain,
        phasecal=phasecal,
        antpos_m=antpos,
        wavelength_m=6.0,
    )
    np.testing.assert_allclose(summed, 3.0 * signal, atol=1e-12)


def test_local_raw_cut_recovers_near_ideal_seven_module_gain():
    metadata = REPO / "data/metadata/cut"
    diagnostics = REPO / "data/server_initial_fits" / (
        f"pansy_disambiguation_diagnostics_{SAMPLE_IDX}.h5"
    )
    initial = REPO / "data/server_initial_fits" / f"highres_fft_i2_p16_{SAMPLE_IDX}.h5"
    if not metadata.exists() or not diagnostics.exists() or not initial.exists():
        pytest.skip("local one-cut coherent-addition regression fixture is unavailable")

    cut = drf.DigitalMetadataReader(str(metadata)).read(
        SAMPLE_IDX - 1, SAMPLE_IDX + 1
    )[SAMPLE_IDX]
    with h5py.File(initial, "r") as handle:
        raw_idx = np.asarray(handle["raw_idx"], dtype=int)
        precise_range_km = np.asarray(handle["range_km"], dtype=float)
        precise_doppler_mps = np.asarray(handle["doppler_mps"], dtype=float)
    with h5py.File(diagnostics, "r") as handle:
        label = handle.attrs["selected_hypothesis"]
        label = label.decode() if isinstance(label, bytes) else str(label)
        group = handle["hypotheses"][label]
        model_km = np.asarray(group["physics_ceplecha_model"], dtype=float)

    z_rx = np.asarray(cut["zrx_echoes_re"], np.float64) + 1j * np.asarray(
        cut["zrx_echoes_im"], np.float64
    )
    z_tx = np.asarray(cut["ztx_pulses_re"], np.float64) + 1j * np.asarray(
        cut["ztx_pulses_im"], np.float64
    )
    with h5py.File(REPO / "data/mesocal.h5", "r") as handle:
        power = np.asarray(handle["pwr"])
    z_rx *= np.real(np.sqrt(power[0]) / np.sqrt(power[:7]))[None, :, None]
    delays = np.asarray(cut["delays"], dtype=float)
    range_sample_km = 299792458.0 / 2.0 / 1e6 / 1e3
    range_gate = precise_range_km / range_sample_km - delays[raw_idx]
    fast_time_s = np.arange(z_tx.shape[1]) / 1e6
    response = np.empty((len(raw_idx), 7), dtype=np.complex128)
    for observation, pulse in enumerate(raw_idx):
        start = int(np.floor(range_gate[observation]))
        fraction = range_gate[observation] - start
        frequency = np.fft.fftfreq(z_rx.shape[-1])
        shifted = np.fft.ifft(
            np.fft.fft(z_rx[pulse], axis=-1)
            * np.exp(2j * np.pi * frequency * fraction)[None, :],
            axis=-1,
        )
        echo = shifted[:, start : start + z_tx.shape[1]]
        template = z_tx[pulse]
        envelope = np.abs(template) ** 2
        use = envelope > 0.05 * np.max(envelope)
        doppler_hz = 2.0 * precise_doppler_mps[observation] / pc.wavelength
        derotation = np.exp(-2j * np.pi * doppler_hz * fast_time_s)
        response[observation] = np.mean(
            echo[:, use] * np.conj(template[use])[None, :] * derotation[use],
            axis=1,
        )

    arrival = enu_line_of_sight_to_arrival_uvw(model_km)
    module_gain = estimate_event_module_voltage_gain(response)

    selected = np.array([9, 10, 11])
    measured_gain = []
    for observation in selected:
        beam_id = int(np.asarray(cut["beam_id"])[raw_idx[observation]])
        coherent = coherent_add_modules(
            response[observation],
            arrival[observation],
            beam_id,
            module_voltage_gain=module_gain,
        )
        calibrated = response[observation] / module_gain
        measured_gain.append(
            np.abs(coherent) ** 2 / np.sum(np.abs(calibrated) ** 2)
        )
    measured_gain = np.asarray(measured_gain)
    np.testing.assert_allclose(
        measured_gain,
        [6.95796204, 6.94937748, 6.97347518],
        atol=0.02,
    )
    assert np.all(measured_gain > 6.8)
