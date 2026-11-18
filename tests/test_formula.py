"""formula_mass(): molecule names of LAMDA and ProDiMo files."""

from __future__ import annotations

import pytest

from pylamda.formula import formula_mass


@pytest.mark.parametrize(
    ("name", "mass"),
    [
        ("CO", 28.010),
        ("13CO", 28.999),  # mass number before the element
        ("C18O", 30.011),  # mass number of the next atom
        ("H13CN", 28.015),
        ("HN13C", 28.015),
        ("HC3N", 51.048),  # a count
        ("N2D+", 30.028),  # deuterium, charge
        ("HDO", 19.021),
        ("o-H3O+", 19.023),  # spin prefix
        ("A-CH3OH", 32.042),
        ("E-CH3OH", 32.042),
        ("c-C3H2", 38.049),
        ("Ar++", 39.95),
        ("Si+", 28.085),
    ],
)
def test_formula_mass(name, mass):
    assert formula_mass(name) == pytest.approx(mass, abs=1e-3)


@pytest.mark.parametrize("name", ["", "X", "'OH", "Xx2", "C0", "13", "C+H", "CO-pH2"])
def test_not_a_formula(name):
    assert formula_mass(name) is None
