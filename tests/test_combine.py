"""add_partner() and drop_partner(): ported from astrolinedata, plus the pylamda parts."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest
from helpers import EMAA_CO, LAMDA

from pylamda import CollisionSet, CombineError, LamdaData, add_partner, drop_partner, read, verify, write
from pylamda.cli import main


def _o_source() -> LamdaData:
    """An O I data set like a CHIANTI conversion: 5 levels (energies slightly different
    from LAMDA's, other quanta notation) and electron rates, also to levels 4 and 5."""
    e = CollisionSet(
        4,
        np.array([10.0, 100.0, 1000.0]),
        np.array([2, 3, 3, 4, 5, 5]),
        np.array([1, 1, 2, 1, 1, 4]),
        np.arange(18, dtype=float).reshape(6, 3) * 1e-10 + 1e-10,
        reference="CHIANTI 11.0.2: Bell et al. 1998",
        partner_text="4 O + e -- Ref: CHIANTI 11.0.2: Bell",
    )
    return LamdaData(
        "O",
        15.999,
        np.array([0.0, 158.265, 226.977, 15867.862, 33792.583]),
        np.array([5.0, 3.0, 1.0, 5.0, 1.0]),
        ["2s2.2p4_3P2", "2s2.2p4_3P1", "2s2.2p4_3P0", "2s2.2p4_1D2", "2s2.2p4_1S0"],
        np.array([2]),
        np.array([1]),
        np.array([8.5e-5]),
        np.array([4744.7]),
        np.array([227.7]),
        collisions=[e],
        source="CHIANTI 11.0.2, ion o_1",
    )


# ------------------------------------------------------------------- add_partner (from astrolinedata)
def test_add_partner():
    target = read(LAMDA / "oatom@lique.dat")  # 3 levels, no electrons (declared)
    assert "e" not in target.partners
    data, added = add_partner(target, _o_source(), "electrons")
    assert data.partners == [*target.partners, "e"] and target.partners[-1] != "e"  # target unchanged
    e = data.partner("e")
    assert (e.partner_id, e.reference) == (4, "CHIANTI 11.0.2: Bell et al. 1998")
    assert e.partner_text == "4 O + e -- Ref: CHIANTI 11.0.2: Bell et al. 1998"
    assert list(zip(e.up, e.low, strict=True)) == [(2, 1), (3, 1), (3, 2)]  # 4 -> 1, 5 -> 1, 5 -> 4 not in the target
    np.testing.assert_array_equal(e.rates, _o_source().collisions[0].rates[:3])  # unchanged
    np.testing.assert_array_equal(e.temperatures, [10.0, 100.0, 1000.0])
    assert added.level_matching == "energy" and added.level_map == {1: 1, 2: 2, 3: 3}
    assert (added.n_transitions, added.n_dropped) == (3, 3)
    assert abs(added.max_energy_difference - 0.0082492) < 1e-6
    assert added.source == "CHIANTI 11.0.2, ion o_1" and added.as_dict()["level_map"] == {"1": 1, "2": 2, "3": 3}
    assert added.warnings == []
    assert data.notes[-1] == (
        "! Collision partner e added from CHIANTI 11.0.2, ion o_1: 3 of 6 transitions, "
        "levels matched by energy (max |dE| 0.00825 cm-1)"
    )
    assert data.source_map is None  # the data no longer fit the file
    # too strict a tolerance: the levels do not match
    with pytest.raises(CombineError, match="no e transition of the source connects two levels"):
        add_partner(target, _o_source(), "e", atol=1e-3, rtol_e=0)


def test_add_partner_checks():
    target = read(LAMDA / "oatom@lique.dat")
    source = _o_source()
    with pytest.raises(CombineError, match="already has collision partner H"):
        add_partner(target, read(LAMDA / "oatom@lique.dat"), "H")
    with pytest.raises(CombineError, match="no collision partner 'H2'"):
        add_partner(target, source, "H2")
    # replacing keeps the position
    data, _ = add_partner(target, read(LAMDA / "oatom@lique.dat"), "He", replace_existing=True, source_label="copy")
    assert data.partners == target.partners and data.notes[-1].startswith("! Collision partner He added from copy")
    # matched by index, but with other weights
    bad = dataclasses.replace(source, level_g=np.array([5.0, 1.0, 3.0, 5.0, 1.0]))
    with pytest.raises(CombineError, match="matched level 2 of the target"):
        add_partner(target, bad, "e", match="index")
    # the target orders two of the levels the other way round
    flipped = dataclasses.replace(target, level_energy=np.array([0.0, 158.268741, 100.0]))
    with pytest.raises(CombineError, match="order these levels differently"):
        add_partner(flipped, source, "e", match="index")
    # H2 next to p-H2/o-H2, and another molecule name: warnings in the record
    h2 = dataclasses.replace(
        source, molecule="O2", collisions=[dataclasses.replace(source.collisions[0], partner_id=1)]
    )
    _, added = add_partner(target, h2, "H2")
    assert any("a code may count H2 collisions twice" in w for w in added.warnings)
    assert any("the molecule names differ" in w for w in added.warnings)
    assert issubclass(CombineError, ValueError)


def test_add_partner_repeated_ids():
    """EMAA's file has two blocks with ID 2 (p-H2O and p-H2): the block must be chosen."""
    emaa = read(EMAA_CO)
    target = emaa.select_partners(["o-H2"])  # without p-H2
    with pytest.raises(CombineError, match="has 2 blocks with the ID of p-H2 .* choose one with block="):
        add_partner(target, emaa, "p-H2")
    data, added = add_partner(target, emaa, "p-H2", block=3)  # block 3: the real p-H2 rates
    assert added.reference.startswith("Dagdigian") and data.partners[-1] == "p-H2"
    with pytest.raises(CombineError, match="collision block 2 of the source is o-H2, not p-H2"):
        add_partner(target, emaa, "p-H2", block=2)


def test_add_partner_writes_and_reads_back(tmp_path):
    data, _ = add_partner(read(LAMDA / "oatom@lique.dat"), _o_source(), "e")
    write(data, tmp_path / "o.dat")
    back = read(tmp_path / "o.dat")
    assert back.partners[-1] == "e" and back.issues == []
    assert not {"C02", "C03", "C10"} & {i.code for i in verify(back)}


# ------------------------------------------------------------------- drop_partner
def test_drop_partner():
    d = read(LAMDA / "catom.dat")
    out, dropped = drop_partner(d, "He")
    assert out.partners == ["H", "e", "H+", "p-H2", "o-H2"] and "He" in d.partners  # input unchanged
    assert (dropped.partner, dropped.partner_id, dropped.block, dropped.n_transitions) == ("He", 6, 4, 3)
    assert out.notes[-1] == "! Collision partner He (block 4, 3 transitions) dropped"
    assert dropped.as_dict()["block"] == 4
    with pytest.raises(CombineError, match="no collision partner 'CO'"):
        drop_partner(d, "CO")


def test_drop_partner_repeated_ids():
    emaa = read(EMAA_CO)
    with pytest.raises(CombineError, match="choose one with block="):
        drop_partner(emaa, "p-H2")
    out, dropped = drop_partner(emaa, "p-H2", block=1)  # the p-H2O block, read by its ID 2
    assert dropped.reference.startswith("Faure") and len(out.collisions) == 4
    assert [c.text_partner for c in out.collisions] == ["ortho-H2O", "para-H2", "ortho-H2", "CO"]


def test_drop_partner_outside_the_standard(tmp_path):
    d = read(LAMDA / "catom.dat")
    d.collisions[0] = dataclasses.replace(d.collisions[0], partner_id=9)
    out, dropped = drop_partner(d, "ID 9")
    assert dropped.partner == "ID 9" and len(out.collisions) == 5


# ------------------------------------------------------------------- command line
def test_combine_commands(tmp_path, capsys):
    source = tmp_path / "source.dat"
    write(_o_source(), source)
    out = tmp_path / "o.dat"
    target = str(LAMDA / "oatom@lique.dat")
    assert main(["add-partner", target, str(source), "e", "-o", str(out)]) == 0
    assert read(out).partners[-1] == "e"
    assert "Collision partner e added from CHIANTI 11.0.2, ion o_1" in capsys.readouterr().err
    assert main(["add-partner", target, str(source), "e", "-o", target]) == 2  # never the input file
    assert "must be a new file" in capsys.readouterr().err
    assert main(["add-partner", target, str(source), "H2", "-o", str(out)]) == 2
    assert "no collision partner 'H2'" in capsys.readouterr().err
    assert main(["drop-partner", str(out), "e", "-o", str(tmp_path / "back.dat")]) == 0
    assert read(tmp_path / "back.dat").partners == read(target).partners
    assert main(["drop-partner", str(EMAA_CO), "p-H2", "-o", str(tmp_path / "x.dat")]) == 2
    assert main(["drop-partner", str(EMAA_CO), "p-H2", "--block", "1", "-o", str(tmp_path / "x.dat")]) == 0
