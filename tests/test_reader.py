"""Reading: real LAMDA and EMAA files, partners by ID, issues, source map."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from helpers import EMAA_CO, LAMDA, make_lamda

import pylamda
from pylamda import LamdaReadError, read


def _write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "x.dat"
    path.write_text(text)
    return path


def _codes(data: pylamda.LamdaData) -> list[str]:
    return [issue.code for issue in data.issues]


# ------------------------------------------------------------------- real files
def test_read_emaa_file():
    """EMAA: references in the file; partners by ID, so the H2O and CO blocks read as p-H2, o-H2, H."""
    d = read(EMAA_CO)
    assert (d.molecule, d.mass, d.n_levels, d.n_lines) == ("CO", 28.0, 14, 13)
    assert d.dataset_doi == "10.17178/EMAA_CO_Rotation_3d0e1a5c" and d.title == ""
    assert d.level_type == "ROTATION" and d.level_reference.startswith("CDMS (Oct. 2000) Müller")
    assert d.level_energy[1] == 3.845033 and d.level_g[1] == 3.0 and d.level_quanta[1] == "1"
    assert d.line_freq[0] == 115.2712018 and d.line_eup[0] == 5.532
    assert [c.partner_id for c in d.collisions] == [2, 3, 2, 3, 5]
    assert d.partners == ["p-H2", "o-H2", "p-H2", "o-H2", "H"]
    assert [c.text_partner for c in d.collisions] == ["para-H2O", "ortho-H2O", "para-H2", "ortho-H2", "CO"]
    p = d.collisions[2]
    assert p.reference == "Dagdigian MNRAS 514 2214 (2022) https://doi.org/10.1093/mnras/stac1430"
    assert p.rates.shape == (p.n_transitions, len(p.temperatures))
    assert d.partner("para-H2") is d.collisions[0]  # the first block with ID 2
    refs = [(r["role"], r["partner"], r["dois"]) for r in d.references()]
    assert refs[0] == ("dataset", None, ["10.17178/EMAA_CO_Rotation_3d0e1a5c"])
    assert refs[1] == ("levels", None, ["10.1016/j.molstruc.2005.01.027"])
    assert ("collisions", "p-H2", ["10.1093/mnras/stac1430"]) in refs
    assert d.issues == []


def test_read_lamda_partners():
    """LAMDA: partner from the ID; ProDiMo-style ' + <partner> ! <ref>' lines give clean refs."""
    d = read(LAMDA / "catom.dat")
    assert d.partners == ["H", "e", "H+", "He", "p-H2", "o-H2"]
    assert d.partner("H").reference == "Launay & Roueff 1977"
    hcn = read(LAMDA / "hcn@hfs.dat")
    assert hcn.molecule == "HCN" and hcn.title.startswith("hyperfine structure up to J=8")
    assert hcn.level_quanta[:2] == ["00_01", "01_01"]
    # "1 H2 = HCN-pH2 collision data ...": the ID (1 = H2) decides; the text names no partner
    assert hcn.partners == ["H2"] and hcn.collisions[0].reference.startswith("H2 = HCN-pH2")
    assert hcn.collisions[0].text_partner is None


def test_read_quirks():
    """Counts decide (as in RADEX/ProDiMo); surplus values and blocks are issues."""
    s = read(LAMDA / "s-h.dat")  # declares 20 temperatures, gives 25
    h = s.partner("H")
    assert len(h.temperatures) == 20 and h.rates.shape == (3, 20) and h.temperatures[-1] == 3000.0
    assert _codes(s) == ["S04", "S04"]
    assert "more than the 20 declared temperatures" in s.issues[0].message
    assert "3 rate rows have more than 20 values" in s.issues[1].message
    assert all(issue.fixable and issue.severity == "warning" for issue in s.issues)
    o = read(LAMDA / "oatom@lique.dat")  # declares 5 partners, holds 6 blocks
    assert len(o.collisions) == 5 and "e" not in o.partners
    assert any("O + e" in line for line in o.notes)
    assert _codes(o) == ["S03"] and "after the last of the 5 declared" in o.issues[0].message
    assert o.source_map.raw[o.issues[0].line - 1].strip() and o.issues[0].line >= o.source_map.notes
    c = read(LAMDA / "c+.dat")  # old rates commented out after the data
    assert c.notes and all(line.lstrip().startswith("!") for line in c.notes if line.strip())
    assert c.issues == []
    oh = read(LAMDA / "oh.dat")
    assert oh.notes[0] == "!NOTES:"


# ------------------------------------------------------------------- partners
@pytest.mark.parametrize(
    ("line", "pid", "partner", "named", "reference"),
    [
        ("2 X + p-H2 ! a reference", 2, "p-H2", "p-H2", "a reference"),
        ("1 X + p-H2 ! a reference", 1, "H2", "p-H2", "a reference"),  # name and ID disagree: ID decides
        ("1 X + o-H2 !", 1, "H2", "o-H2", ""),
        ("1 X + o-H2 !  ", 1, "H2", "o-H2", ""),
        ("1 X + o-H2", 1, "H2", "o-H2", ""),
        ("4 X + e -- Ref: Smith 2001", 4, "e", "e", "Smith 2001"),
        ("2 X + pH2O -- Ref: Smith 2001", 2, "p-H2", "pH2O", "Smith 2001"),  # unknown spelling, kept as written
        ("3 X-oH2 from Yang et al. (2010)", 3, "o-H2", None, "X-oH2 from Yang et al. (2010)"),
        ("9 X + He ! ref", 9, None, "He", "ref"),  # ID outside the standard
    ],
)
def test_partner_by_id(tmp_path, line, pid, partner, named, reference):
    d = read(_write(tmp_path, make_lamda(partner=line)))
    c = d.collisions[0]
    assert (c.partner_id, c.partner, c.text_partner, c.reference) == (pid, partner, named, reference)
    assert c.partner_text == line.rstrip()
    assert d.issues == []
    if partner is None:
        assert c.label == "ID 9"
        with pytest.raises(KeyError, match=r"no collision partner 'He' \(have: ID 9\)"):
            d.partner("He")


def test_partner_name_spellings():
    assert pylamda.partner_name("Electrons") == "e" and pylamda.partner_name(" ortho-H2 ") == "o-H2"
    assert pylamda.partner_name("para-h2o") == "p-H2O"
    assert pylamda.partner_name("CO2") is None


def test_split_partner_line():
    split = pylamda.split_partner_line
    assert split("5 C + H ! Launay & Roueff 1977") == ("H", "Launay & Roueff 1977")
    assert split("2 CO + p-H2 -- Ref: Yang 2010 https://doi.org/10.1/x") == ("p-H2", "Yang 2010 https://doi.org/10.1/x")
    assert split("2 CO -- Ref: Yang 2010") == (None, "Yang 2010")
    assert split("2 CO + unknown ! Yang") == (None, "CO + unknown ! Yang")  # ProDiMo form needs a known name
    assert split("2") == (None, "")


# ------------------------------------------------------------------- issues
def test_issue_file_too_short(tmp_path):
    text = "!MOLECULE\nCO\n!WEIGHT\n28.0\n!NLEV\n3\n1 0.0 1.0\n2 3.8 3.0\n"
    with pytest.raises(LamdaReadError, match="line 8: S01 error: the file ends before level 3") as exc:
        read(_write(tmp_path, text))
    assert [i.code for i in exc.value.issues] == ["S01"]
    assert isinstance(exc.value, ValueError)


def test_issue_numbering(tmp_path):
    d = read(_write(tmp_path, make_lamda(levels="1 0.0 1.0\n3 3.8 3.0\n4 11.5 5.0")))
    assert [(i.code, i.line) for i in d.issues] == [("S02", 9)]
    assert d.issues[0].message == "levels: row 2 has the number 3 (1 more rows)"
    np.testing.assert_array_equal(d.level_g, [1.0, 3.0, 5.0])  # read by position
    d = read(_write(tmp_path, make_lamda(rates="1 2 1 1e-11 2e-11\n5 3 1 3e-11 4e-11\n3 3 2 5e-11 6e-11")))
    assert [(i.code, i.line) for i in d.issues] == [("S02", 28)]
    assert d.issues[0].message.startswith("collision block 1: row 2 has the number 5")


def test_issue_comment_without_mark(tmp_path):
    """A comment line without '!' before a count is skipped (S05)."""
    d = read(_write(tmp_path, make_lamda(ncol_comment="NUMBER OF COLL TRANS")))
    assert d.collisions[0].rates[0, 0] == 1.0e-11
    issue = d.issues[0]
    assert (issue.code, issue.line, issue.severity, issue.fixable) == ("S05", 20, "error", True)
    assert "before the number of collisional transitions of collision block 1" in issue.message
    # a line that is neither a count nor a comment stops the reader
    with pytest.raises(LamdaReadError, match="line 20: S06 error: cannot read the number of collisional"):
        read(_write(tmp_path, make_lamda(ncol_comment="-")))


def test_issue_unreadable_numbers(tmp_path):
    lines = "1 2 1 7,2e-08 115.27 5.53\n2 3 2 6.9D-07 230.54 16.6"
    d = read(_write(tmp_path, make_lamda(lines=lines, temps="10.0 2O.0")))
    assert [(i.code, i.line) for i in d.issues] == [("S06", 14), ("S06", 25)]
    assert "'7,2e-08'" in d.issues[0].message
    assert np.isnan(d.line_a[0]) and d.line_a[1] == 6.9e-07  # D exponent is Fortran
    assert np.isnan(d.collisions[0].temperatures[1])


def test_issue_short_rows(tmp_path):
    d = read(_write(tmp_path, make_lamda(levels="1 0.0 1.0\n2 3.8\n3 11.5 5.0", lines="1 2 1\n2 3 2 6.9e-07")))
    assert [(i.code, i.line) for i in d.issues] == [("S07", 9), ("S07", 14)]
    assert np.isnan(d.level_g[1]) and np.isnan(d.line_a[0]) and d.line_up[0] == 2


def test_issue_partner_without_id(tmp_path):
    """The standard requires the partner ID; a name alone is an issue, the block is still read."""
    d = read(_write(tmp_path, make_lamda(partner="X + e- ! a reference")))
    assert [(i.code, i.line) for i in d.issues] == [("S06", 19)]
    assert "must start with the partner ID" in d.issues[0].message
    assert d.collisions[0].partner_id == 0 and d.collisions[0].partner is None
    assert d.collisions[0].n_transitions == 3


def test_values_across_lines(tmp_path):
    """Temperatures and rate rows continue on the next line, as in Fortran list-directed reads."""
    d = read(
        _write(
            tmp_path, make_lamda(temps="10.0\n20.0", rates="1 2 1\n1e-11 2e-11\n2 3 1 3e-11 4e-11\n3 3 2 5e-11 6e-11")
        )
    )
    c = d.collisions[0]
    np.testing.assert_array_equal(c.temperatures, [10.0, 20.0])
    np.testing.assert_array_equal(c.rates[0], [1e-11, 2e-11])
    assert d.issues == []
    np.testing.assert_array_equal(d.source_map.blocks[0].rows, [28, 30, 31])


def test_strict(tmp_path):
    path = _write(tmp_path, make_lamda(temps="10.0 20.0 30.0"))
    assert _codes(read(path)) == ["S04"]
    with pytest.raises(LamdaReadError, match="S04 warning: collision block 1: more than the 2 declared") as exc:
        read(path, strict=True)
    assert exc.value.path == str(path)
    read(_write(tmp_path, make_lamda()), strict=True)  # no issues: no error


def test_issue_as_text_and_dict():
    issue = pylamda.issues.make_issue("S04", "too many values", 12, fixable=True)
    assert str(issue) == "line 12: S04 warning: too many values"
    assert issue.as_dict() == {
        "code": "S04",
        "severity": "warning",
        "message": "too many values",
        "line": 12,
        "fixable": True,
    }
    assert set(pylamda.CODES) >= {"S01", "S02", "S03", "S04", "S05", "S06", "S07"}


# ------------------------------------------------------------------- source map
def test_source_map(tmp_path):
    path = _write(tmp_path, make_lamda(after="!NOTES\n! a note\n"))
    d = read(path)
    m = d.source_map
    assert m.path == str(path) and m.raw == path.read_text().splitlines()
    assert (m.molecule, m.mass, m.nlev, m.nlin, m.npart) == (2, 4, 6, 12, 17)
    np.testing.assert_array_equal(m.levels, [8, 9, 10])
    np.testing.assert_array_equal(m.lines, [14, 15])
    b = m.blocks[0]
    assert (b.partner, b.ncol, b.ntemp, b.temperatures) == (19, 21, 23, 25)
    np.testing.assert_array_equal(b.rows, [27, 28, 29])
    assert m.notes == 30 and m.raw[m.notes - 1] == "!NOTES"
    assert d.notes == ["! a note"]
    # the map points at the lines of the items
    assert m.raw[m.lines[1] - 1].split()[3] == "6.9e-07"
    assert m.raw[b.partner - 1] == d.collisions[0].partner_text
    assert read(_write(tmp_path, make_lamda())).source_map.notes is None
