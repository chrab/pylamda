"""The data model: rate interpolation, partner lookup, copies."""

from __future__ import annotations

import numpy as np
import pytest
from helpers import EMAA_CO, LAMDA

from pylamda import CollisionSet, read


def test_rate_interpolation():
    c = CollisionSet(
        2,
        np.array([10.0, 20.0, 40.0]),
        np.array([2, 3]),
        np.array([1, 1]),
        np.array([[1.0, 2.0, 4.0], [10.0, 10.0, 10.0]]),
    )
    assert c.partner == "p-H2" and c.text_partner is None and c.label == "p-H2"
    assert np.allclose(c.rate(20.0), [2.0, 10.0])
    assert np.allclose(c.rate(30.0), [3.0, 10.0])
    with pytest.raises(ValueError, match="outside the p-H2 grid 10-40 K"):
        c.rate(5.0)
    assert np.allclose(c.rate(5.0, clamp=True), [1.0, 10.0])
    assert np.allclose(c.rate(100.0, clamp=True), [4.0, 10.0])
    one = CollisionSet(8, np.array([10.0]), np.array([2]), np.array([1]), np.array([[3.0]]))
    assert one.partner is None and one.label == "ID 8"
    assert np.allclose(one.rate(10.0), [3.0])


def test_partner_lookup():
    d = read(LAMDA / "catom.dat")
    assert d.partner("oH2") is d.collisions[5]
    assert d.partner("electrons") is d.collisions[1]
    with pytest.raises(KeyError, match="no collision partner 'H2'"):
        d.partner("H2")
    sel = d.select_partners(["p-H2", "o-H2", "CO"])
    assert sel.partners == ["p-H2", "o-H2"] and d.n_levels == sel.n_levels
    assert read(EMAA_CO).select_partners(["H"]).collisions[0].text_partner == "CO"


def test_copy():
    d = read(LAMDA / "catom.dat")
    c = d.copy()
    c.level_energy[1] = -1.0
    c.collisions[0].rates[0, 0] = -1.0
    c.notes.append("! new")
    assert d.level_energy[1] != -1.0 and d.collisions[0].rates[0, 0] != -1.0
    assert "! new" not in d.notes
    assert c.source_map is not d.source_map and c.source_map.raw == d.source_map.raw


def test_data_without_file():
    d = read(LAMDA / "catom.dat")
    d.source_map = None
    assert d.issues == []
