import numpy as np

from run_catalogue_mass_profiles import density_marginalized_radius_probability


def test_density_marginalization_preserves_radius_density_degeneracy():
    radius_um = np.geomspace(1.0, 1.0e4, 241)
    conditional_delta_chi2 = (np.log(radius_um / 180.0) / 0.08) ** 2

    probability, weights, density, joint_weights, conditional = (
        density_marginalized_radius_probability(
            radius_um,
            conditional_delta_chi2,
            density_bounds_kg_m3=(100.0, 8000.0),
        )
    )

    assert np.isclose(np.sum(weights), 1.0)
    assert np.isclose(np.sum(joint_weights), 1.0)
    assert probability.shape == radius_um.shape
    assert conditional.shape == (len(radius_um), len(density))

    # Both radii can represent the same ballistic product with a permitted
    # density, so neither hypothesis may be rejected by the dynamics alone.
    index_100 = int(np.argmin(np.abs(radius_um - 100.0)))
    index_1000 = int(np.argmin(np.abs(radius_um - 1000.0)))
    assert np.min(conditional[index_100]) < 0.05
    assert np.min(conditional[index_1000]) < 0.05
    assert probability[index_1000] > 0.5 * probability[index_100]


def test_density_marginalization_rejects_too_small_radius_for_finite_density():
    radius_um = np.geomspace(1.0, 1.0e4, 241)
    conditional_delta_chi2 = (np.log(radius_um / 180.0) / 0.08) ** 2
    probability, *_ = density_marginalized_radius_probability(
        radius_um,
        conditional_delta_chi2,
        density_bounds_kg_m3=(100.0, 8000.0),
    )

    index_10 = int(np.argmin(np.abs(radius_um - 10.0)))
    index_100 = int(np.argmin(np.abs(radius_um - 100.0)))
    assert probability[index_10] < 1e-6 * probability[index_100]
