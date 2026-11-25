"""verify(): one synthetic file per check code, and the options."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from helpers import EMAA_CO, LAMDA, collision_block, make_lamda

import pylamda
from pylamda import Issue, read, verify
from pylamda.verify import HC_K, excitation_slope

LEVELS = "1 0.0 1.0 J=0\n2 3.845033 3.0 J=1\n3 11.534919 5.0 J=2"
DOI = "https://doi.org/10.1088/0004-637X/718/2/1062"


def _check(tmp_path: Path, select: str | None = None, **kw) -> list[Issue]:
    path = tmp_path / "x.dat"
    path.write_text(make_lamda(**kw))
    return verify(path, select=None if select is None else [select])


def _one(issues: list[Issue], code: str) -> Issue:
    found = [i for i in issues if i.code == code]
    assert len(found) == 1, issues
    return found[0]


def test_clean_file(tmp_path):
    assert _check(tmp_path) == []


# ------------------------------------------------------------------- L: levels
def test_l01_energy_order(tmp_path):
    levels = "1 0.0 1.0 J=0\n2 11.534919 5.0 J=2\n3 3.845033 3.0 J=1"
    issue = _one(_check(tmp_path, "L01", levels=levels), "L01")
    assert (issue.line, issue.severity, issue.fixable) == (10, "warning", True)
    assert issue.message.startswith("energies not increasing: level 3 (3.84503 cm-1) below level 2")


def test_l02_weights(tmp_path):
    issues = _check(tmp_path, "L02", levels="1 0.0 0.0 J=0\n2 3.845033 2.5 J=1\n3 11.534919 5.0 J=2")
    assert [(i.severity, i.line) for i in issues] == [("error", 8), ("info", 9)]
    assert "g <= 0" in issues[0].message and "not an integer" in issues[1].message


def test_l03_duplicate_level(tmp_path):
    levels = "1 0.0 1.0 J=0\n2 3.845033 3.0 J=1\n3 3.845033 3.0 J=1"
    issue = _one(_check(tmp_path, "L03", levels=levels), "L03")
    assert issue.line == 10 and "level 3 has the same E, g and quanta as level 2" in issue.message


def test_l04_quanta(tmp_path):
    issue = _one(_check(tmp_path, "L04", levels="1 0.0 1.0\n2 3.845033 3.0\n3 11.534919 5.0"), "L04")
    assert (issue.severity, issue.count, issue.line) == ("info", 3, 8)
    issue = _one(_check(tmp_path, "L04", levels="1 0.0 1.0 a\n2 3.845033 3.0 b\n3 11.534919 5.0  b"), "L04")
    assert issue.message == "quanta not unique: 'b' at levels 2, 3" and issue.count == 2


def test_l05_lowest_energy(tmp_path):
    levels = "1 1.0 1.0 J=0\n2 4.845033 3.0 J=1\n3 12.534919 5.0 J=2"
    issue = _one(_check(tmp_path, "L05", levels=levels), "L05")
    assert issue.message == "the lowest level (level 1) has E = 1 cm-1, not 0"


# ------------------------------------------------------------------- R: radiative transitions
def test_r01_r02_level_numbers(tmp_path):
    lines = "1 4 1 7.203e-08 115.2711894 5.532\n2 2 2 6.910e-07 230.5369826 16.596"
    issue = _one(_check(tmp_path, "R01", lines=lines), "R01")
    assert (issue.line, issue.count, issue.severity) == (14, 2, "error")
    lines = "1 1 2 7.203e-08 115.2711894 5.532\n2 3 2 6.910e-07 230.5369826 16.596"
    issue = _one(_check(tmp_path, "R02", lines=lines), "R02")
    assert issue.line == 14 and "line 1 (1 -> 2)" in issue.message


def test_r03_einstein_a(tmp_path):
    lines = "1 2 1 0.0 115.2711894 5.532\n2 3 2 NaN 230.5369826 16.596"
    issue = _one(_check(tmp_path, "R03", lines=lines), "R03")
    assert (issue.line, issue.count) == (14, 2)


def test_r04_r05_columns(tmp_path):
    lines = "1 2 1 7.203e-08 116.0 6.0\n2 3 2 6.910e-07 230.5369826 16.596"
    issues = _check(tmp_path, "R0", lines=lines)
    r04, r05 = _one(issues, "R04"), _one(issues, "R05")
    assert r04.line == 14 and r04.fixable and "line 1 (2 -> 1): 116 GHz, computed 115.2711894" in r04.message
    assert r05.line == 14 and "6 K, computed 5.532144" in r05.message
    # missing columns
    issues = _check(tmp_path, "R0", lines="1 2 1 7.203e-08\n2 3 2 6.910e-07 230.5369826 16.596")
    assert [(i.code, i.message) for i in issues] == [
        ("R04", "frequency missing: line 1 (2 -> 1)"),
        ("R05", "E_up [K] missing: line 1 (2 -> 1)"),
    ]


def test_r04_rounding(tmp_path):
    """The tolerance includes the rounding of the numbers as written."""
    path = tmp_path / "x.dat"
    path.write_text(make_lamda(lines="1 2 1 7.203e-08 115.27 5.53\n2 3 2 6.910e-07 230.54 16.6"))
    assert verify(path, select=["R04", "R05"], rtol=0) == []
    path.write_text(make_lamda(lines="1 2 1 7.203e-08 115.2800 5.532\n2 3 2 6.910e-07 230.5369826 16.596"))
    assert [i.code for i in verify(path, select=["R04", "R05"], rtol=0)] == ["R04"]
    # coarse energies widen the tolerance: 3.845 is good to 0.0005 cm-1, i.e. 0.015 GHz
    levels = "1 0.0 1.0 J=0\n2 3.845 3.0 J=1\n3 11.534919 5.0 J=2"
    path.write_text(
        make_lamda(levels=levels, lines="1 2 1 7.203e-08 115.2711894 5.532\n2 3 2 6.910e-07 230.5369826 16.596")
    )
    assert verify(path, select=["R04"], rtol=0) == []


def test_r06_duplicate_line(tmp_path):
    lines = "1 2 1 7.203e-08 115.2711894 5.532\n2 2 1 7.203e-08 115.2711894 5.532"
    issue = _one(_check(tmp_path, "R06", lines=lines), "R06")
    assert issue.line == 15 and issue.message == "line 2 (2 -> 1) repeats line 1"


# ------------------------------------------------------------------- C: collisions
def test_c01_partner_id(tmp_path):
    issue = _one(_check(tmp_path, "C", partner=f"9 X + He ! Smith 2001 {DOI}"), "C01")
    assert (issue.line, issue.severity) == (19, "error")


def test_c01_not_repeated_after_s06(tmp_path):
    """A partner line without the ID is S06; it is not reported again as C01."""
    issues = _check(tmp_path, partner=f"X + He ! Smith 2001 {DOI}")
    assert [i.code for i in issues if i.line == 19] == ["S06"]


def test_c02_name_and_id(tmp_path):
    issue = _one(_check(tmp_path, "C", partner=f"1 X + para-H2 ! Smith 2001 {DOI}"), "C02")
    assert (issue.line, issue.severity, issue.fixable) == (19, "warning", True)
    assert issue.message == 'partner "X + para-H2" (ID 1 = H2): the name says p-H2 (read as H2, ProDiMo reads p-H2)'
    assert _check(tmp_path, "C", partner=f"2 X + pH2 ! Smith 2001 {DOI}") == []


def test_c03_c04_partners(tmp_path):
    after = collision_block(f"2 X + p-H2 -- Ref: {DOI}")
    issue = _one(_check(tmp_path, "C", npart="2", after=after), "C03")
    assert issue.line == 31
    assert issue.message == 'partner "X + p-H2" (ID 2 = p-H2, block 2): the ID is also in block 1'
    after = collision_block(f"1 X + H2 -- Ref: {DOI}")
    issue = _one(_check(tmp_path, "C", npart="2", after=after), "C04")
    assert issue.line == 31
    assert issue.message.startswith('partner "X + H2" (ID 1 = H2): H2 together with p-H2')


def test_c05_temperatures(tmp_path):
    issue = _one(_check(tmp_path, "C", temps="20.0 10.0"), "C05")
    assert issue.line == 25 and issue.message.endswith("20 10")
    assert _one(_check(tmp_path, "C", temps="0.0 10.0"), "C05")


def test_c06_c07_rates(tmp_path):
    rates = "1 2 1 -1.0e-11 2.0e-11\n2 3 1 0.0 0.0\n3 3 2 NaN 6.0e-11"
    issues = _check(tmp_path, "C", rates=rates)
    c06, c07 = _one(issues, "C06"), _one(issues, "C07")
    assert (c06.line, c06.count, c06.severity) == (27, 2, "error")
    assert (c07.line, c07.count, c07.severity) == (28, 1, "info")


def test_c08_inverted_rows(tmp_path):
    rates = "1 1 2 1.0e-11 2.0e-11\n2 1 3 3.0e-11 4.0e-11\n3 3 2 5.0e-11 6.0e-11"
    issue = _one(_check(tmp_path, "C", rates=rates), "C08")
    assert (issue.line, issue.count, issue.fixable) == (27, 2, True)
    assert "row 1 (1 -> 2) and 1 more" in issue.message


def test_c09_c10_rows(tmp_path):
    rates = "1 2 1 1.0e-11 2.0e-11\n2 2 1 3.0e-11 4.0e-11\n3 4 2 5.0e-11 6.0e-11"
    issues = _check(tmp_path, "C", rates=rates)
    assert _one(issues, "C09").line == 28
    c10 = _one(issues, "C10")
    assert c10.line == 29 and "row 3 (4 -> 2)" in c10.message


def test_c11_partner_outside_standard(tmp_path):
    issue = _one(_check(tmp_path, "C", partner=f"2 X + p-H2O -- Ref: {DOI}"), "C11")
    assert (
        issue.message
        == "partner \"X + p-H2O\" (ID 2 = p-H2): 'p-H2O' is not a standard partner; read by the ID as p-H2"
    )
    issue = _one(_check(tmp_path, "C", partner=f"2 X + pH2O -- Ref: {DOI}"), "C11")  # unknown spelling
    assert "'pH2O'" in issue.message


@pytest.mark.parametrize(
    "text",
    [
        "1 SO2 + He ! rates (Green 1995 ApJS) scaled *1.39 and extrapolated to higher J",
        "1  O2 + He ! rates of Lique 2010, J. Chem. Phys., 132, 044311, scaled by 1.37'",
        "2 X + He ! He rates *1.37 (Smith 2001)",
        "3 X + He ! Smith 2001, x 1.4",
        "1 X + He -- Ref: Smith 2001, scaling of the He rates",
    ],
)
def test_c12_he_rates_scaled_to_h2(tmp_path, text):
    """The name says He, the ID is H2, and the text says scaled: C12 instead of C02."""
    issues = _check(tmp_path, "C", partner=text)
    issue = _one(issues, "C12")
    assert (issue.line, issue.severity, issue.fixable) == (19, "warning", True)
    assert "the name says He, the text says scaled: probably He rates scaled to" in issue.message
    assert issue.message.endswith("ProDiMo reads He)")
    assert [i.code for i in issues if i.code != "C13"] == ["C12"]  # C13: the hint for "scaled"


@pytest.mark.parametrize(
    "text",
    [
        "1 X + He ! Lique et al. 2010",  # not scaled: name and ID simply disagree
        "4 X + He ! scaled by 1.37",  # ID is not H2
        "1 X + p-H2 ! scaled from He",  # the name is not He
    ],
)
def test_c12_not_reported(tmp_path, text):
    assert [i.code for i in _check(tmp_path, "C", partner=text) if i.code != "C13"] == ["C02"]


def test_c12_he_with_its_id(tmp_path):
    """He with ID 6 and 'scaled' is consistent: no C12, only the hint C13 (info)."""
    assert [i.code for i in _check(tmp_path, partner=f"6 X + He ! scaled by 1.0 {DOI}")] == ["C13"]


def test_emaa_extra_partners():
    codes = [(i.code, i.line) for i in verify(EMAA_CO, select=["C"])]
    assert codes == [("C11", 41), ("C11", 105), ("C03", 169), ("C03", 269), ("C11", 369)]


# ------------------------------------------------------------------- F: references
def test_f01_no_reference(tmp_path):
    issue = _one(_check(tmp_path, "F", nlev="3"), "F01")
    assert (issue.line, issue.severity) == (6, "warning")
    # a reference in a comment line counts
    assert _check(tmp_path, "F", nlev="3", after="! levels from Smith et al. (2001), doi:10.1000/xyz\n") == []


def test_f02_f03_partner_reference(tmp_path):
    issue = _one(_check(tmp_path, "F", partner="2 X + p-H2"), "F02")
    assert (issue.line, issue.severity) == (19, "warning")
    issue = _one(_check(tmp_path, "F", partner="2 X + p-H2 ! Yang et al. (2010)"), "F03")
    assert (issue.line, issue.severity) == (19, "info") and "'Yang et al. (2010)'" in issue.message
    issue = _one(_check(tmp_path, "F", nlev="3 - Ref: Smith et al. 2001"), "F03")
    assert issue.line == 6 and "levels" in issue.message


@pytest.mark.parametrize(
    "reference",
    [
        "Yang 2010 doi:10.1088/0004-637X/718/2/1062",
        "Yang 2010 https://scixplorer.org/abs/2010ApJ...718.1062Y",
        "Yang 2010, 2010ApJ...718.1062Y",
        "Yang 2010, arXiv:1006.2722",
    ],
)
def test_reference_identifiers(tmp_path, reference):
    assert _check(tmp_path, "F", partner=f"2 X + p-H2 ! {reference}") == []


# ------------------------------------------------------------------- P: physics heuristics
def test_p01_excitation_rates(tmp_path):
    levels = "1 0.0 1.0 J=0\n2 3.845033 3.0 J=1\n3 100.0 5.0 J=2"
    temps = np.array([5.0, 10.0, 20.0, 40.0, 80.0])
    de = 100.0 * HC_K
    excitation = 1e-10 * np.exp(-de / temps)
    rates = (
        "1 2 1 " + " ".join(["1e-11"] * 5) + "\n"
        "2 3 1 " + " ".join(f"{k:.6e}" for k in excitation) + "\n"
        "3 3 2 " + " ".join(["5e-11"] * 5)
    )
    issue = _one(_check(tmp_path, "P", levels=levels, ntemp="5", temps="5 10 20 40 80", rates=rates), "P01")
    assert (issue.line, issue.count, issue.severity) == (28, 1, "info")
    assert "row 2 (3 -> 1)" in issue.message
    slope = excitation_slope(temps, excitation, de)
    assert slope == pytest.approx(-de, rel=1e-6)
    assert excitation_slope(temps, np.full(5, 1e-11), de) == pytest.approx(0.0, abs=1e-9)
    assert excitation_slope(temps, excitation, 10.0) is None  # no temperature below ΔE/2k


def test_c13_scaled(tmp_path):
    issue = _one(_check(tmp_path, "C", partner=f"1 X + H2 ! scaled from He rates, Smith 2001 {DOI}"), "C13")
    assert (issue.line, issue.severity) == (19, "info")
    assert issue.message == 'partner "X + H2" (ID 1 = H2): comment says rates are scaled "scaled"'
    issue = _one(_check(tmp_path, "C", partner=f"1 X + H2 ! derived from p-H2 rates {DOI}"), "C13")
    assert issue.message.endswith('scaled "from p-H2"')
    # the factor right after "scaled" is part of the quoted text
    for text, part in [
        ("1 X + H2 ! scaled (*1.37) X - He from Green (1991)", "scaled (*1.37)"),
        ("1 X + H2 ! scaled by 1.4 from He", "scaled by 1.4"),
        ("1 X + H2 ! He rates, scaling x1.39", "scaling x1.39"),
        ("1 X + H2 ! scaled 1991 data", "scaled"),  # a year is not a factor
    ]:
        issue = _one(_check(tmp_path, "C", partner=f"{text} {DOI}"), "C13")
        assert issue.message.endswith(f'scaled "{part}"'), issue.message
    # on by default (not a physics check)
    assert [i.code for i in _check(tmp_path, partner=f"1 X + H2 ! scaling of He rates {DOI}")] == ["C13"]
    assert _check(tmp_path, "C", partner=f"1 X + H2 ! Smith 2001 {DOI}") == []


def test_p03_large_rates(tmp_path):
    rates = "1 2 1 1.0e-11 2.0e-7\n2 3 1 3.0e-11 4.0e-11\n3 3 2 5.0e-11 6.0e-11"
    issue = _one(_check(tmp_path, "P", rates=rates), "P03")
    assert (issue.line, issue.count) == (27, 1) and "rates > 1e-08" in issue.message
    partner = f"4 X + e ! Smith 2001 {DOI}"
    assert _check(tmp_path, "P", rates=rates, partner=partner) == []  # electrons may be faster


def test_p04_mass(tmp_path):
    issue = _one(_check(tmp_path, "P", molecule="13CO", mass="28.0"), "P04")
    assert issue.line == 4 and issue.message == "mass 28 differs from the formula mass 28.999 of '13CO'"
    assert _check(tmp_path, "P", molecule="13CO", mass="29.0") == []
    assert _check(tmp_path, "P", molecule="'OH", mass="5") == []  # not a formula: no check


# ------------------------------------------------------------------- options
def test_select_ignore_error_on(tmp_path):
    path = tmp_path / "x.dat"
    # P04 (line 4): the mass of 13CO is not 28
    path.write_text(make_lamda(molecule="13CO", partner="1 X + p-H2 ! scaled from He, Yang et al. (2010)"))
    assert [i.code for i in verify(path)] == ["C02", "C13", "F03"]  # the P checks are off by default
    assert [i.code for i in verify(path, physics=True)] == ["P04", "C02", "C13", "F03"]
    assert [i.code for i in verify(path, physics=False)] == ["C02", "C13", "F03"]
    assert [i.code for i in verify(path, physics=True, ignore=["F", "P0"])] == ["C02", "C13"]
    # selecting P codes runs them, unless physics=False
    assert [i.code for i in verify(path, select=["c0", "P"])] == ["P04", "C02"]
    assert [i.code for i in verify(path, select=["P04"])] == ["P04"]
    assert [i.code for i in verify(path, select=["P"], physics=False)] == []
    raised = verify(path, select=["C02"], error_on=["C"])
    assert [(i.code, i.severity) for i in raised] == [("C02", "error")]


def test_physics_not_run_by_default(tmp_path, monkeypatch):
    """Off by default, the physics heuristics are not run at all (they are slow on large files)."""
    import sys

    verify_module = sys.modules["pylamda.verify"]  # pylamda.verify is also the function
    calls = []
    monkeypatch.setattr(verify_module, "_check_physics", lambda data, add: calls.append(data))
    path = tmp_path / "x.dat"
    path.write_text(make_lamda())
    verify(path)
    verify(path, select=["C", "F"])
    assert calls == []
    verify(path, physics=True)
    verify(path, select=["P04"])
    assert len(calls) == 2


def test_reader_issues_included():
    issues = verify(LAMDA / "s-h.dat", select=["S"])
    assert [(i.code, i.line) for i in issues] == [("S04", 25), ("S04", 27)]


def test_sorted_by_line(tmp_path):
    path = tmp_path / "x.dat"
    path.write_text(make_lamda(nlev="3", partner="1 X + p-H2", molecule="13CO"))
    lines = [i.line for i in verify(path)]
    assert lines == sorted(lines)


def test_data_without_source_map(tmp_path):
    """Data not read from a file (or changed after reading) give issues without line numbers."""
    d = read(LAMDA / "catom.dat")
    d.source_map = None
    issues = verify(d)
    assert issues and all(i.line is None for i in issues)
    d = read(LAMDA / "catom.dat")
    d.collisions[0].up = d.collisions[0].up[:-1]  # no longer fits the source map
    d.collisions[0].low = d.collisions[0].low[:-1]
    d.collisions[0].rates = d.collisions[0].rates[:-1]
    assert all(i.line is None for i in verify(d, select=["F03"]) if "block 1 " in i.message)


def test_unreadable_file(tmp_path):
    path = tmp_path / "x.dat"
    path.write_text("!MOLECULE\nX\n")
    with pytest.raises(pylamda.LamdaReadError, match="S01"):
        verify(path)


def test_every_code_documented():
    """Every code has a severity of SEVERITIES and appears in the docs (docs/checks.rst)."""
    text = (Path(__file__).parents[1] / "docs" / "checks.rst").read_text()
    for code, (severity, description) in pylamda.CODES.items():
        assert severity in pylamda.SEVERITIES and description
        assert f"``{code}``" in text, code


def test_constants():
    assert pylamda.verify is verify
    assert math.isclose(HC_K, 1.438776877, rel_tol=1e-9)


def test_partner_prefix(tmp_path):
    """Collision messages start with the pair as written and the ID (design: no block numbers)."""
    # a partner line without "<molecule> + <partner>": only the ID
    issue = _one(_check(tmp_path, "C", partner=f"1 CO-H2 scaled from He {DOI}"), "C13")
    assert issue.message.startswith("partner ID 1 = H2: ")
    # the "-- Ref:" form
    issue = _one(_check(tmp_path, "C", partner=f"1 X + H2 -- Ref: scaled, {DOI}"), "C13")
    assert issue.message.startswith('partner "X + H2" (ID 1 = H2): ')
    # an ID outside the standard
    issue = _one(_check(tmp_path, "C", partner=f"9 X + H2 ! {DOI}"), "C01")
    assert issue.message.startswith('partner "X + H2" (ID 9): the ID is not in 1..7')
