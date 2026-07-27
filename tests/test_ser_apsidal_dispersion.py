import numpy as np

from simulate_ser_apsidal_dispersion import (
    LOW_INCLINATION_DEG,
    SER_MEAN,
    circular_mean_unwrapped_deg,
    circular_std_deg,
    integrate_ensembles,
    paired_initial_elements,
)


def test_circular_std_does_not_create_wrap_boundary_dispersion():
    values = np.asarray([359.0, 0.0, 1.0])
    assert circular_std_deg(values) < 1.0


def test_circular_mean_unwraps_across_zero():
    values = np.asarray([[358.0, 359.0, 0.0, 1.0], [0.0, 0.0, 2.0, 2.0]])
    mean = circular_mean_unwrapped_deg(values, axis=0)
    assert np.all(np.diff(mean) > 0.0)
    assert np.isclose(mean[-1] - mean[0], 2.5, atol=0.1)


def test_paired_ensembles_differ_only_in_mean_inclination():
    elements = paired_initial_elements(8, 42, 0.01, 0.002, 0.2)
    difference = elements[0] - elements[1]
    assert np.allclose(difference[:, [0, 1, 3, 4, 5]], 0.0)
    assert np.allclose(difference[:, 2], SER_MEAN[2] - LOW_INCLINATION_DEG)


def test_short_rebound_integration_preserves_shape_and_initial_state():
    initial = paired_initial_elements(2, 42, 0.001, 0.0002, 0.02)
    times, elements = integrate_ensembles(initial, duration_years=1.0, output_step_years=1.0)
    assert times.shape == (2,)
    assert elements.shape == (2, 2, 2, 7)
    assert np.allclose(elements[..., 0, 0], initial[..., 0], rtol=0.0, atol=2e-5)
    assert np.allclose(elements[..., 0, 1], initial[..., 1], rtol=0.0, atol=2e-6)
