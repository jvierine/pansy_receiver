"""Independent checks of the reusable analysis kernels and coordinate conventions."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from meteor_analysis import native
from meteor_analysis.stations import TIBET_MST
from meteor_analysis.coordinates import incoming_asymptote, hammer, MU_EARTH
from meteor_analysis.pulse_phase import measure_pairs

pytestmark = pytest.mark.skipif(
    not native.LIB.exists(), reason="build meteor/analysis first"
)


def test_site_and_enu_are_explicit():
    assert TIBET_MST.latitude_deg == 30.1 and TIBET_MST.longitude_deg == 90.5
    r = TIBET_MST.enu_to_ecef()
    np.testing.assert_allclose(r.T @ r, np.eye(3), atol=1e-14)
    assert np.linalg.det(r) > 0


def test_native_integer_matches_independent_direct_sum():
    rng = np.random.default_rng(421)
    rx = rng.normal(size=(10, 160)) + 1j * rng.normal(size=(10, 160))
    tx = rng.normal(size=32) + 1j * rng.normal(size=32)
    fd = 23567.0
    fs = 500000.0
    centre = np.sum(np.arange(32) * abs(tx) ** 2) / np.sum(abs(tx) ** 2)
    expected = np.sum(
        rx[:, 50:82]
        * tx.conj()
        * np.exp(-2j * np.pi * fd * (np.arange(32) - centre) / fs),
        axis=1,
    ) / np.linalg.norm(tx)
    np.testing.assert_allclose(
        native.matched(rx, tx, 50, fd, fs), expected, rtol=2e-13, atol=2e-13
    )


def test_rk4_against_independent_scipy_drag_ablation():
    state = np.array([1000.0, 500.0, 100000.0, 12000.0, 20000.0, -40000.0])
    radius = 0.0004
    times = np.linspace(0, 0.15, 31)
    rho = 8e-8

    def rhs(t, y):
        v = y[3:6]
        speed = np.linalg.norm(v)
        return np.r_[
            v,
            -0.75 * rho / (3000 * y[6]) * speed * v,
            -rho * 1e-8 * speed**3 / (4 * 3000),
        ]

    truth = solve_ivp(
        rhs, (0, times[-1]), np.r_[state, radius], t_eval=times, rtol=1e-11, atol=1e-13
    ).y.T
    got = native.propagate(
        state, radius, times, [0, 200000], [rho, rho], 4287, gravity=False
    )
    np.testing.assert_allclose(got, truth, rtol=1e-10, atol=1e-7)
    half = native.propagate(
        state, radius, times, [0, 200000], [rho, rho], 4287, gravity=False, dt=0.00025
    )
    np.testing.assert_allclose(got, half, rtol=1e-10, atol=1e-7)


def test_geocentric_energy_and_incoming_asymptote():
    state = np.array([6500000.0, 0.0, 0.0, -40000.0, 10000.0, 2000.0])
    mu = MU_EARTH
    got = incoming_asymptote(state)
    np.testing.assert_allclose(
        np.dot(got, got),
        np.dot(state[3:], state[3:]) - 2 * mu / np.linalg.norm(state[:3]),
        rtol=1e-14,
    )

    def rhs(t, y):
        return np.r_[y[3:], -mu * y[:3] / np.linalg.norm(y[:3]) ** 3]

    orbit = solve_ivp(rhs, (0, -100000.0), state, rtol=2e-10, atol=0.001)
    far = orbit.y[3:, -1]
    assert (
        np.degrees(
            np.arccos(
                np.clip(
                    np.dot(got, far) / np.linalg.norm(got) / np.linalg.norm(far), -1, 1
                )
            )
        )
        < 0.001
    )


def test_radiant_projection_apex_and_handedness():
    np.testing.assert_allclose(hammer(270, 0), [0, 0], atol=1e-14)
    assert hammer(0, 0)[0] < 0
    assert abs(hammer(90, 0)[0]) > 2.8


def test_nonoverlapping_pulse_phase_sign_and_alias_branch():
    n = 80
    t = np.arange(n) * 0.00128
    wavelength = 6.0
    velocity = -53210.0
    pos = np.column_stack([np.zeros(n), np.zeros(n), 100000 + velocity * t])
    phase = -4 * np.pi * np.linalg.norm(pos, axis=1) / wavelength
    amp = np.exp(1j * phase)[:, None] * np.ones((n, 10)) * 100
    m = dict(
        time_s=t,
        amplitude=amp,
        noise_covariance=np.eye(10),
        radial_velocity_mps=np.full(n, velocity),
        range_velocity_covariance=np.broadcast_to(np.diag([10.0, 10000.0]), (n, 2, 2)),
    )
    array = lambda uv: np.ones((len(uv), 10), complex)
    result = measure_pairs(
        m, array, pos, np.zeros(10), np.arange(n), wavelength, np.eye(6) * 1e-8
    )
    assert len(result["velocity_mps"]) == 40
    np.testing.assert_allclose(result["velocity_mps"], velocity, atol=1e-7)
    assert np.all(result["keep"])
    assert len(set(result["previous"]) & set(result["current"])) == 0
    switching = np.arange(n) % 2
    rejected = measure_pairs(
        m,
        array,
        pos,
        np.zeros(10),
        np.arange(n),
        wavelength,
        np.eye(6) * 1e-8,
        switching,
    )
    assert len(rejected["velocity_mps"]) == 0


def test_density_radius_scaling_is_an_exact_model_degeneracy():
    state = [1000, 500, 100000, 12000, 20000, -40000]
    times = np.linspace(0, 0.15, 31)
    first = native.propagate(state, 0.0004, times, [0, 200000], [8e-8, 8e-8], 4287)
    scaled = native.propagate(state, 0.0008, times, [0, 200000], [16e-8, 16e-8], 4287)
    np.testing.assert_allclose(first[:, :6], scaled[:, :6], rtol=1e-12, atol=1e-8)
    np.testing.assert_allclose(2 * first[:, 6], scaled[:, 6], rtol=1e-12)


def test_joint_profile_recovers_injected_shrinking_radius():
    from meteor_analysis.trajectory import EventFit

    t = np.linspace(0, 0.15, 40)
    state = np.array([1000.0, -700.0, 100000.0, 12000.0, 15000.0, -50000.0])
    truth = native.propagate(state, 0.0003, t, [0, 200000], [8e-8, 8e-8], 4287)
    position = truth[:, :3]
    velocity = truth[:, 3:6]
    distance = np.linalg.norm(position, axis=1)

    class Array:
        centres = np.array([[i * 7.3, (i * i % 7) * 8.1] for i in range(10)])

        def __call__(self, uv):
            return np.exp(2j * np.pi / 6 * (uv @ self.centres.T))

    array = Array()
    uv = position[:, :2] / distance[:, None]
    m = dict(
        time_s=t,
        range_m=distance,
        radial_velocity_mps=np.sum(position * velocity, axis=1) / distance,
        amplitude=100 * array(uv),
        noise_covariance=np.eye(10),
        range_velocity_covariance=np.broadcast_to(
            np.diag([0.5**2, 5.0**2]), (len(t), 2, 2)
        ),
    )
    fit = EventFit(
        m,
        array,
        TIBET_MST,
        ([0, 200000], [8e-8, 8e-8]),
        np.zeros(10),
        np.ones(10) * 0.05,
    )
    profile = fit.profile(
        np.r_[state / 1000, np.zeros(9)],
        np.array([0.0001, 0.0002, 0.0003, 0.0004, 0.001]),
        False,
    )
    assert profile["radius_m"][np.argmin(profile["chi2"])] == 0.0003
    assert np.min(profile["chi2"]) < 1e-4
