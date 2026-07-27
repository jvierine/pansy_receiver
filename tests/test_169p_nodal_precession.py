import numpy as np

from simulate_169p_meteoroid_nodal_precession import radiation_beta, state_elements


def test_100_um_beta_for_three_gram_per_cc_particle():
    assert np.isclose(radiation_beta(100.0, 3000.0), 0.003829, rtol=2e-3)


def test_state_elements_recovers_simple_ellipse():
    mu = 4.0 * np.pi**2
    a = 2.0
    e = 0.5
    perihelion = a * (1.0 - e)
    speed = np.sqrt(mu * (1.0 + e) / (a * (1.0 - e)))
    elements = state_elements([perihelion, 0.0, 0.0], [0.0, speed, 0.0], mu)
    assert np.isclose(elements[0], a)
    assert np.isclose(elements[1], e)
    assert np.isclose(elements[2], 0.0)
