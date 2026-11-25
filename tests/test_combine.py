"""add_partner(), drop_partner() and scale_partner(): add_partner ported from astrolinedata, plus the pylamda parts."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest
from helpers import EMAA_CO, LAMDA, collision_block, make_lamda

from pylamda import (
    CollisionSet,
    CombineError,
    LamdaData,
    __version__,
    add_partner,
    drop_partner,
    read,
    scale_partner,
    verify,
    write,
)
from pylamda.cli import main
from pylamda.combine import SCALE_DEFAULT_DIGITS, _significant_digits, patch_scale_partner
from pylamda.model import PYLAMDA_NOTE
from pylamda.reader import read_string

SCALE = f"! pylamda {__version__} scale-partner: "  # the prefix of a scale-partner note


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
        f"! pylamda {__version__} add-partner: Collision partner e added from CHIANTI 11.0.2, ion o_1: "
        "3 of 6 transitions, levels matched by energy (max |dE| 0.00825 cm-1)"
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
    assert data.partners == target.partners and data.notes[-1].startswith(
        f"! pylamda {__version__} add-partner: Collision partner He added"
    )
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
    assert (
        out.notes[-1] == f"! pylamda {__version__} drop-partner: Collision partner He (block 4, 3 transitions) dropped"
    )
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


# ------------------------------------------------------------------- scale_partner
_RATES4 = "1 2 1 1.234e-11 2.345e-11\n2 3 1 3.456e-11 4.567e-11\n3 3 2 5.678e-11 6.789e-11"
_HCN_LIKE = "1 X + H2 ! scaled (*1.37) X-He from Green (1991)"


def _scaled_file(partner: str = _HCN_LIKE, rates: str = _RATES4, after: str = "", npart: str = "1") -> LamdaData:
    """A synthetic file with rates of 4 significant digits."""
    return read_string(make_lamda(partner=partner, rates=rates, after=after, npart=npart))


def _round(x: np.ndarray, digits: int) -> np.ndarray:
    return np.array([float(f"{v:.{digits - 1}e}") for v in x.ravel()]).reshape(x.shape)


def test_scale_partner_in_place():
    d = _scaled_file()
    before = d.collisions[0].rates.copy()
    out, scaled = scale_partner(d, "H2", 1.37 / 1.39)
    c = out.collisions[0]
    assert (c.partner_id, c.partner) == (1, "H2") and out.partners == ["H2"]
    np.testing.assert_array_equal(c.rates, _round(before * (1.37 / 1.39), 4))
    np.testing.assert_array_equal(d.collisions[0].rates, before)  # input unchanged
    assert c.partner_text == "1 X + H2 -- Ref: scaled (*1.37) X-He from Green (1991)"
    assert (scaled.block, scaled.new_block, scaled.n_transitions, scaled.digits) == (1, 1, 3, 4)
    assert (scaled.kept_source, scaled.replaced_block, scaled.text_factor) == (False, None, 1.37)
    assert out.notes[-1] == f"{SCALE}Collision partner H2 (block 1, 3 transitions) scaled by 0.985612"
    assert out.source_map is None
    # the text gives 1.37, the factor is 0.9856; the reference still says "scaled (*1.37)"
    assert len(scaled.warnings) == 2
    assert "the partner text gives the factor 1.37; the factor used is 0.985612" in scaled.warnings[0]
    assert "mentions scaling" in scaled.warnings[1]
    assert scaled.as_dict()["factor"] == pytest.approx(0.985612, rel=1e-6)


def test_scale_partner_change_partner():
    """HCN-like: H2 rates that are He x 1.37 become He rates."""
    d = _scaled_file()
    out, scaled = scale_partner(d, "H2", 1 / 1.37, to="He", reference="Green (1991)")
    c = out.collisions[0]
    assert (c.partner_id, c.partner, c.reference) == (6, "He", "Green (1991)")
    assert c.partner_text == "6 X + He -- Ref: Green (1991)"
    np.testing.assert_array_equal(c.rates, _round(d.collisions[0].rates / 1.37, 4))
    assert (scaled.partner, scaled.partner_id, scaled.new_partner, scaled.new_partner_id) == ("H2", 1, "He", 6)
    assert scaled.warnings == []  # 1/1.37 fits the text; a new reference text
    assert out.notes[-1] == f"{SCALE}Collision partner H2 (block 1, 3 transitions) scaled by 0.729927 and relabelled He"
    assert not {"C02", "C03", "C12"} & {i.code for i in verify(out)}


def test_scale_partner_derive():
    """Keep the He block and add H2 = He x 1.37."""
    d = _scaled_file(partner="6 X + He ! Lique et al. (2006)")
    out, scaled = scale_partner(d, "He", 1.37, to="H2", keep=True)
    assert out.partners == ["He", "H2"]
    np.testing.assert_array_equal(out.collisions[0].rates, d.collisions[0].rates)
    np.testing.assert_array_equal(out.collisions[1].rates, _round(d.collisions[0].rates * 1.37, 4))
    assert out.collisions[1].partner_text == "1 X + H2 -- Ref: Lique et al. (2006)"
    assert (scaled.kept_source, scaled.new_block, scaled.warnings) == (True, 2, [])
    assert (
        out.notes[-1]
        == f"{SCALE}Collision partner He (block 1, 3 transitions) scaled by 1.37 and added as H2 (block 2)"
    )


def test_scale_partner_replace_existing():
    # H2 (block 1) and He (block 2)
    two = dict(partner="1 X + H2 ! Green", after=collision_block("6 X + He ! Lique"), npart="2")
    d = _scaled_file(**two)
    with pytest.raises(CombineError, match="already has collision partner H2 .use replace_existing"):
        scale_partner(d, "He", 1.37, to="H2", keep=True)
    out, scaled = scale_partner(d, "He", 1.37, to="H2", keep=True, replace_existing=True)
    assert out.partners == ["H2", "He"] and (scaled.new_block, scaled.replaced_block) == (1, 1)
    np.testing.assert_array_equal(out.collisions[0].rates, _round(d.collisions[1].rates * 1.37, 3))
    assert out.notes[-1].endswith("added as H2 (block 1); the former H2 block (block 1) replaced")
    # without keep: the scaled block stays at the place of the source, the old He block goes
    out, scaled = scale_partner(d, "H2", 1 / 1.37, to="He", replace_existing=True)
    assert out.partners == ["He"] and out.collisions[0].reference == "Green"
    assert (scaled.new_block, scaled.replaced_block) == (1, 2)
    # the source after the replaced block
    d = _scaled_file(partner="6 X + He ! Lique", after=collision_block("1 X + H2 ! Green"), npart="2")
    out, scaled = scale_partner(d, "H2", 1 / 1.37, to="He", replace_existing=True)
    assert out.partners == ["He"] and out.collisions[0].reference == "Green" and scaled.new_block == 1


def test_scale_partner_errors():
    d = _scaled_file()
    for factor in (0.0, -1.37, float("nan"), float("inf")):
        with pytest.raises(CombineError, match="finite and > 0"):
            scale_partner(d, "H2", factor)
    with pytest.raises(CombineError, match="digits must be at least 1"):
        scale_partner(d, "H2", 1.37, digits=0)
    with pytest.raises(CombineError, match="'CO' is not a standard collision partner"):
        scale_partner(d, "CO", 1.37)
    with pytest.raises(CombineError, match="'p-H2O' is not a standard collision partner"):
        scale_partner(d, "H2", 1.37, to="p-H2O")
    with pytest.raises(CombineError, match="keep=True needs another partner"):
        scale_partner(d, "H2", 1.37, to="H2", keep=True)
    with pytest.raises(CombineError, match="no collision partner 'He'"):
        scale_partner(d, "He", 1.37)
    three = collision_block("6 X + He ! a") + collision_block("6 X + He ! b")
    with pytest.raises(CombineError, match="has 2 blocks of He"):
        scale_partner(_scaled_file(after=three, npart="3"), "H2", 1 / 1.37, to="He", replace_existing=True)
    with pytest.raises(CombineError, match="choose one with block="):
        scale_partner(read(EMAA_CO), "p-H2", 1.0)
    _, scaled = scale_partner(read(EMAA_CO), "p-H2", 1.0, block=3)
    assert scaled.block == 3


def test_scale_partner_warnings():
    # C12: He in the text, ID 1, and a factor 1.39 that fits neither 1/1.37 nor 1.37
    d = _scaled_file(partner="1 X + He ! rates (Green 1995) scaled *1.39")
    _, scaled = scale_partner(d, "H2", 1 / 1.37, to="He", reference="Green 1995")
    assert len(scaled.warnings) == 2
    assert "names 'He'; the block was taken by its ID 1 as H2" in scaled.warnings[0]
    assert "gives the factor 1.39" in scaled.warnings[1]
    # unusual directions
    he = _scaled_file(partner="6 X + He ! Lique")
    _, scaled = scale_partner(he, "He", 0.73, to="H2")
    assert scaled.warnings == ["He -> H2 with a factor < 1: He rates scaled to H2 are usually larger (about 1.37)"]
    h2 = _scaled_file(partner="1 X + H2 ! Green")
    _, scaled = scale_partner(h2, "H2", 1.37, to="He")
    assert scaled.warnings[0].startswith("H2 -> He with a factor > 1")
    # a pair that is not He and H2
    e = _scaled_file(partner="4 X + e ! Faure")
    _, scaled = scale_partner(e, "e", 1.0, to="H")
    assert scaled.warnings == ["e -> H: scaling is a standard approximation only between He and H2 (p-H2, o-H2)"]
    # H2 next to p-H2 in the result (C04)
    d = _scaled_file(partner="2 X + p-H2 ! Yang", after=collision_block("6 X + He ! Lique"), npart="2")
    _, scaled = scale_partner(d, "He", 1.37, to="H2", keep=True)
    assert scaled.warnings == ["X: H2 next to p-H2; a code may count H2 collisions twice"]
    # H2 -> p-H2 (a relabelling within H2) gives no pair warning
    _, scaled = scale_partner(h2, "H2", 1.0, to="p-H2")
    assert scaled.warnings == []


def test_scale_partner_digits():
    d = _scaled_file()
    _, scaled = scale_partner(d, "H2", 1.37)
    assert scaled.digits == 4  # as written
    out, scaled = scale_partner(d, "H2", 1.37, digits=2)
    assert scaled.digits == 2 and out.collisions[0].rates[0, 0] == 1.7e-11
    out, scaled = scale_partner(dataclasses.replace(d, source_map=None), "H2", 1.37)
    assert scaled.digits == SCALE_DEFAULT_DIGITS
    np.testing.assert_array_equal(out.collisions[0].rates, _round(d.collisions[0].rates * 1.37, 6))
    # two digits as written ("1.0e-11"): at least 3
    _, scaled = scale_partner(read_string(make_lamda()), "p-H2", 1.37)
    assert scaled.digits == 3
    # rate rows that continue on the next line
    _, scaled = scale_partner(
        _scaled_file(rates="1 2 1 1.234e-11\n 2.34567e-11\n2 3 1 1e-11 1e-11\n3 3 2 1e-11 1e-11"), "H2", 2.0
    )
    assert scaled.digits == 6


@pytest.mark.parametrize(
    ("token", "digits"), [("5.47E-10", 3), ("0.00123", 3), ("10.", 2), ("1.0D-10", 2), ("7", 1), ("0.0", 0), ("abc", 0)]
)
def test_significant_digits(token, digits):
    assert _significant_digits(token) == digits


def test_scale_partner_writes_and_reads_back(tmp_path):
    out, _ = scale_partner(_scaled_file(), "H2", 1 / 1.37, to="He", reference="Green (1991)")
    write(out, tmp_path / "x.dat")
    back = read(tmp_path / "x.dat")
    assert back.partners == ["He"] and back.issues == []
    np.testing.assert_array_equal(back.collisions[0].rates, out.collisions[0].rates)
    assert not {"C02", "C03", "C12"} & {i.code for i in verify(back)}


def test_scale_partner_command(tmp_path, capsys):
    src = tmp_path / "in.dat"
    src.write_text(make_lamda(partner="2 X + p-H2 ! Yang", after=collision_block("6 X + He ! scaled *1.39"), npart="2"))
    out = tmp_path / "out.dat"
    assert main(["scale-partner", str(src), "He", "1.37", "--to", "H2", "--keep", "-o", str(out)]) == 0
    assert read(out).partners == ["p-H2", "He", "H2"]
    err = capsys.readouterr().err
    assert "warning: the partner text gives the factor 1.39" in err
    assert "warning: X: H2 next to p-H2" in err
    assert "Collision partner He (block 2, 1 transitions) scaled by 1.37 and added as H2 (block 3)" in err
    assert f"{out}:" in err and "C04 warning" in err  # the C checks of the result
    # a fraction, a new reference text and fixed digits
    assert (
        main(
            [
                "scale-partner",
                str(src),
                "p-H2",
                "1/1.37",
                "--to",
                "He",
                "--replace",
                "--reference",
                "Yang 2010",
                "--digits",
                "4",
                "-o",
                str(out),
            ]
        )
        == 0
    )
    back = read(out)
    assert back.partners == ["He"] and back.collisions[0].reference == "Yang 2010"
    assert back.collisions[0].rates[0, 0] == pytest.approx(7.299e-12, rel=0)
    capsys.readouterr()
    assert main(["scale-partner", str(src), "He", "1.37", "-o", str(src)]) == 2  # never the input file
    assert "must be a new file" in capsys.readouterr().err
    assert main(["scale-partner", str(src), "H2", "1.37", "-o", str(out)]) == 2
    assert "no collision partner 'H2'" in capsys.readouterr().err
    for bad in ("abc", "1/0", "1/2/3"):
        with pytest.raises(SystemExit) as exc:
            main(["scale-partner", str(src), "He", bad, "-o", str(out)])
        assert exc.value.code == 2
    assert "not a number or a fraction a/b" in capsys.readouterr().err


def test_notes_are_not_references():
    """The note lines of the combine commands match PYLAMDA_NOTE, as those of fix do."""
    d = read(LAMDA / "catom.dat")
    out, _ = drop_partner(d, "He")
    out, _ = scale_partner(out, "p-H2", 1.1)
    out, _ = add_partner(out, read(LAMDA / "catom.dat"), "He")
    notes = out.notes[-3:]
    assert [n.split(":")[0] for n in notes] == [
        f"! pylamda {__version__} {command}" for command in ("drop-partner", "scale-partner", "add-partner")
    ]
    assert all(PYLAMDA_NOTE.match(n) for n in notes)
    assert not PYLAMDA_NOTE.match("! pylamda is a Python package: Smith 2001")


# ------------------------------------------------------------------- scale_partner, keeping the layout
SCALED_RATES = (
    "   1    2    1  5.500E-11  0.000E+00\n   2    3    1  1.500E-11  2.0E-11\n   3    3    2  1.000E-10  4.250E-11"
)


def _he_as_h2_file(tmp_path, newline: str = "\n"):
    """A ProDiMo-style file: H2 rates that are He rates times 1.37, columns, E notation."""
    text = make_lamda(partner="1 X + H2 ! scaled (*1.37) X-He from Smith (2001)", rates=SCALED_RATES)
    path = tmp_path / "x.dat"
    path.write_bytes(text.replace("\n", newline).encode())
    return path, text.splitlines()


def test_patch_scale_partner(tmp_path):
    """Only the partner line, the changed rate values and one note line change."""
    path, before = _he_as_h2_file(tmp_path)
    text, scaled = patch_scale_partner(path, "H2", 1 / 1.37, to="He", reference="X-He from Smith (2001)")
    after = text.splitlines()
    assert len(after) == len(before) + 1 and after[-1].startswith(f"! pylamda {__version__} scale-partner: ")
    changed = [k for k, (x, y) in enumerate(zip(before, after, strict=False)) if x != y]
    partner = before.index("1 X + H2 ! scaled (*1.37) X-He from Smith (2001)")
    first = before.index("!TRANS + UP + LOW + COLLRATES(cm^3 s^-1)") + 1
    assert changed == [partner, first, first + 1, first + 2]
    assert after[partner] == "6 X + He -- Ref: X-He from Smith (2001)"
    # each number in its style and column: 4 significant digits as written; 0 stays as written
    assert after[first : first + 3] == [
        "   1    2    1  4.015E-11  0.000E+00",
        "   2    3    1  1.095E-11 1.460E-11",  # the right end stays; one blank is left
        "   3    3    2  7.299E-11  3.102E-11",
    ]
    assert scaled.digits == 4 and scaled.new_partner == "He" and scaled.warnings == []
    back = read_string(text)
    assert back.partners == ["He"] and back.issues == []
    np.testing.assert_allclose(back.collisions[0].rates[0], [4.015e-11, 0.0])
    assert "C13" not in [i.code for i in verify(back)]


def test_patch_scale_partner_keeps_crlf_and_refuses(tmp_path):
    path, _ = _he_as_h2_file(tmp_path, newline="\r\n")
    text, _ = patch_scale_partner(path, "H2", 1 / 1.37, to="He")
    assert text.count("\r\n") == text.count("\n") and text.endswith("\r\n")
    # data without a file, or a partner the data have already
    with pytest.raises(CombineError, match="only be kept for data read from a file"):
        patch_scale_partner(dataclasses.replace(read(path), source_map=None), "H2", 0.5)
    with pytest.raises(CombineError, match="already has collision partner"):
        patch_scale_partner(LAMDA / "catom.dat", "He", 1.37, to="p-H2")


def test_scale_partner_command_keeps_the_layout(tmp_path, capsys):
    path, before = _he_as_h2_file(tmp_path)
    ref = ["--reference", "X-He from Smith (2001)"]
    assert main(["scale-partner", str(path), "H2", "1/1.37", "--to", "He", *ref, "--in-place"]) == 0
    after = path.read_text().splitlines()
    partner = before.index("1 X + H2 ! scaled (*1.37) X-He from Smith (2001)")
    assert after[:partner] == before[:partner]  # the header, levels and lines as they were
    assert "   1    2    1  4.015E-11  0.000E+00" in after
    assert (
        "Collision partner H2 (block 1, 3 transitions) scaled by 0.729927 and relabelled He" in capsys.readouterr().err
    )
    # --normalise: the layout of pylamda.write
    path, _ = _he_as_h2_file(tmp_path)
    out = tmp_path / "n.dat"
    assert main(["scale-partner", str(path), "H2", "1/1.37", "--to", "He", "--normalise", "-o", str(out)]) == 0
    assert "    1     2     1 4.015e-11 0.000e+00" in out.read_text().splitlines()
    # --keep adds a block: always the normalised layout
    assert main(["scale-partner", str(path), "H2", "1/1.37", "--to", "He", "--keep", "-o", str(out)]) == 0
    assert read(out).partners == ["H2", "He"]
    with pytest.raises(SystemExit):  # -o or --in-place, not both
        main(["scale-partner", str(path), "H2", "2", "--in-place", "-o", str(out)])
