import pytest

from match_small_body_orbits import Orbit, jopek_dh


def test_jopek_dh_is_zero_for_identical_orbits():
    orbit = Orbit(
        q_au=0.899,
        e=0.712,
        i_deg=91.22,
        omega_deg=317.12,
        node_deg=290.12,
    )
    assert jopek_dh(orbit, orbit) == pytest.approx(0.0, abs=1e-12)


def test_jopek_dh_is_symmetric():
    first = Orbit(0.899, 0.712, 91.22, 317.12, 290.12)
    second = Orbit(1.0697649, 0.36915, 63.16081, 278.98083, 303.06782)
    assert jopek_dh(first, second) == pytest.approx(jopek_dh(second, first))
