"""fix() and patch(): one synthetic case per fix, both modes, and the command line."""

from __future__ import annotations

import numpy as np
import pytest
from helpers import LAMDA, collision_block, make_lamda

import pylamda
from pylamda import FixError, add_reference, compare, fix, patch, read, verify
from pylamda.cli import main
from pylamda.fix import format_like, replace_tokens
from pylamda.reader import read_string
from pylamda.verify import HC_K

DOI = "https://doi.org/10.1088/0004-637X/718/2/1062"


def _file(tmp_path, **kw):
    path = tmp_path / "x.dat"
    path.write_text(make_lamda(**kw))
    return path


def _codes(data_or_text) -> set[str]:
    data = read_string(data_or_text) if isinstance(data_or_text, str) else data_or_text
    return {i.code for i in verify(data)}


def _same(path, **options) -> tuple[str, list]:
    """patch() and fix() give the same data; returns the patched text and the log."""
    data, log = fix(path, **options)
    text, log2 = patch(path, **options)
    assert [r.code for r in log] == [r.code for r in log2]
    assert not compare(data, read_string(text)), compare(data, read_string(text)).summary()
    return text, log


# ------------------------------------------------------------------- helpers
def test_replace_tokens():
    line = "    1     2     1  7.2e-08   115.27"
    assert replace_tokens(line, {1: "12"}) == "    1    12     1  7.2e-08   115.27"  # right end kept
    assert replace_tokens(line, {4: "115.2711894"}) == "    1     2     1  7.2e-08 115.2711894"
    assert replace_tokens("1 2 3", {1: "22"}) == "1 22 3"  # at least one blank
    assert replace_tokens("\t1\t2", {1: "3"}) == "\t1\t3"  # tabs kept
    assert replace_tokens("1 2  ! c", {0: "9"}) == "9 2  ! c"


@pytest.mark.parametrize(
    ("token", "value", "expected"),
    [
        ("1.2e-11", 1.6e-11, "1.60e-11"),  # at least 3 significant digits
        ("1.234E-11", 2.0e-10, "2.000E-10"),
        ("3.0D-10", 4.5e-10, "4.50D-10"),
        ("1.0e+00", 12.0, "1.20e+01"),
        ("5.0e-01", 1.5, "1.50e+00"),
        ("115.27", 115.2711894, "115.27"),
        ("5", 5.6, "6"),
    ],
)
def test_format_like(token, value, expected):
    assert format_like(token, value) == expected


# ------------------------------------------------------------------- structure
def test_nothing_to_fix(tmp_path):
    path = _file(tmp_path)
    text, log = patch(path)
    assert log == [] and text == path.read_text()  # byte for byte
    data, log = fix(path)
    assert log == [] and not compare(read(path), data)


def test_line_endings_kept(tmp_path):
    path = tmp_path / "x.dat"
    path.write_bytes(make_lamda(ncol_comment="NUMBER OF COLL TRANS").replace("\n", "\r\n").encode())
    text, log = patch(path)
    assert [r.code for r in log] == ["S05"]
    assert text.count("\r\n") == text.count("\n") and "\r\n!NUMBER OF COLL TRANS\r\n" in text


def test_s05_comment_mark(tmp_path):
    text, log = _same(_file(tmp_path, ncol_comment="NUMBER OF COLL TRANS"))
    assert [r.code for r in log] == ["S05"] and "\n!NUMBER OF COLL TRANS\n" in text
    assert "S05" not in _codes(text)
    assert text.rstrip().endswith("fix: S05 line 20: '!' added to the comment line")


def test_s04_surplus_values(tmp_path):
    rates = "1 2 1 1.0e-11 2.0e-11 9.9e-11\n2 3 1 3.0e-11 4.0e-11\n3 3 2 5.0e-11 6.0e-11"
    text, log = _same(_file(tmp_path, temps="10.0 20.0 30.0", rates=rates))
    assert [(r.code, r.count) for r in log] == [("S04", 2)]
    lines = text.splitlines()
    assert lines[24] == "10.0 20.0" and lines[26] == "1 2 1 1.0e-11 2.0e-11"
    assert "S04" not in _codes(text)


def test_s02_row_numbers(tmp_path):
    levels = "1 0.0 1.0 J=0\n5 3.845033 3.0 J=1\n3 11.534919 5.0 J=2"
    rates = "1 2 1 1.0e-11 2.0e-11\n2 3 1 3.0e-11 4.0e-11\n   30 3 2 5.0e-11 6.0e-11"
    text, log = _same(
        _file(
            tmp_path,
            levels=levels,
            rates=rates,
            lines="1 2 1 7.203e-08 115.2711894 5.532\n x 3 2 6.910e-07 230.5369826 16.596",
        )
    )
    assert [(r.code, r.message) for r in log] == [
        ("S02", "levels: 1 row number corrected"),
        ("S02", "radiative transitions: 1 row number corrected"),
        ("S02", "collision block 1: 1 row number corrected"),
    ]
    out = text.splitlines()
    assert out[8] == "2 3.845033 3.0 J=1" and out[14].startswith(" 2 3 2") and out[28] == "    3 3 2 5.0e-11 6.0e-11"
    assert not {"S02", "S06"} & _codes(text)


def test_s03_extra_blocks(tmp_path):
    after = collision_block(f"6 X + He ! Smith 2001 {DOI}")
    path = _file(tmp_path, after=after)
    assert "S03" in _codes(read(path))
    text, log = _same(path, extra_blocks="include")
    assert [r.message for r in log] == ["NPART 1 -> 2: the blocks after the last one counted in"]
    assert read_string(text).partners == ["p-H2", "He"] and text.splitlines()[16] == "2"
    text, log = _same(path, extra_blocks="comment")
    assert log[0].code == "S03" and read_string(text).partners == ["p-H2"]
    assert "\n!6 X + He ! Smith 2001" in text and "S03" not in _codes(text)
    text, log = patch(path)  # no option: not fixed
    assert log == [] and "S03" in _codes(text)
    bad = _file(tmp_path, after="1 2 3\n")
    with pytest.raises(FixError, match="do not form complete collision blocks"):
        patch(bad, extra_blocks="include")


# ------------------------------------------------------------------- partners
@pytest.mark.parametrize(
    ("trust", "line", "pid"),
    [("id", f"1 X + H2 -- Ref: Smith 2001 {DOI}", 1), ("name", f"2 X + p-H2 -- Ref: Smith 2001 {DOI}", 2)],
)
def test_c02_trust(tmp_path, trust, line, pid):
    path = _file(tmp_path, partner=f"1 X + para-H2 ! Smith 2001 {DOI}")
    text, log = _same(path, trust=trust)
    assert text.splitlines()[18] == line and read_string(text).collisions[0].partner_id == pid
    assert log[0].code == "C02" and "C02" not in _codes(text)
    text, log = patch(path)  # without --trust: no safe default
    assert log == [] and "C02" in _codes(text)


def test_c12_trust_id(tmp_path):
    path = _file(tmp_path, partner=f"1 X + He ! Smith 2001, scaled by 1.37 {DOI}")
    text, log = _same(path, trust="id")
    assert [(r.code, r.message) for r in log] == [("C12", "collision block 1: name He -> H2 (trusted ID 1)")]
    assert text.splitlines()[18] == f"1 X + H2 -- Ref: Smith 2001, scaled by 1.37 {DOI}"


# ------------------------------------------------------------------- levels
def test_l01_sort_levels(tmp_path):
    levels = "1 0.0 1.0 J=0\n2 11.534919 5.0 J=2\n3 3.845033 3.0 J=1"
    lines = "1 3 1 7.203e-08 115.2711894 5.532\n2 2 3 6.910e-07 230.5369826 16.596"
    rates = "1 3 1 1.0e-11 2.0e-11\n2 2 1 3.0e-11 4.0e-11\n3 2 3 5.0e-11 6.0e-11"
    path = _file(tmp_path, levels=levels, lines=lines, rates=rates)
    text, log = _same(path)
    assert [(r.code, r.count) for r in log] == [("L01", 2)]
    out = text.splitlines()
    assert out[7:10] == ["1 0.0 1.0 J=0", "2 3.845033 3.0 J=1", "3 11.534919 5.0 J=2"]
    assert out[13:15] == ["1 2 1 7.203e-08 115.2711894 5.532", "2 3 2 6.910e-07 230.5369826 16.596"]
    assert out[26:29] == ["1 2 1 1.0e-11 2.0e-11", "2 3 1 3.0e-11 4.0e-11", "3 3 2 5.0e-11 6.0e-11"]
    assert verify(read_string(text)) == []  # sorted, renumbered, nothing else changed
    assert patch(path, sort_levels=False)[1] == []


# ------------------------------------------------------------------- inverted rows
@pytest.mark.parametrize("method", ["swap", "chianti", "balance"])
def test_c08_inverted(tmp_path, method):
    rates = "1 1 2 1.000e-11 2.000e-11\n2 3 1 3.0e-11 4.0e-11\n3 3 2 5.0e-11 6.0e-11"
    path = _file(tmp_path, rates=rates)
    data, log = fix(path, inverted=method)
    assert [(r.code, r.count) for r in log] == [("C08", 1)]
    c = data.collisions[0]
    assert (c.up[0], c.low[0]) == (2, 1)
    de = 3.845033 * HC_K
    expected = {
        "swap": [1e-11, 2e-11],
        "chianti": [1e-11 * 1.0 / 3.0, 2e-11 * 1.0 / 3.0],
        "balance": [1e-11 * (1.0 / 3.0) * np.exp(de / 10.0), 2e-11 * (1.0 / 3.0) * np.exp(de / 20.0)],
    }[method]
    np.testing.assert_allclose(c.rates[0], expected, rtol=1e-12)
    text, log2 = patch(path, inverted=method)
    assert text.splitlines()[26].split()[1:3] == ["2", "1"] and "C08" not in _codes(text)
    np.testing.assert_allclose(read_string(text).collisions[0].rates[0], expected, rtol=5e-4)  # 4 digits as written
    assert patch(path)[1] == []  # no option: not fixed


def test_c08_balance_overflow(tmp_path):
    levels = "1 0.0 1.0 J=0\n2 3.845033 3.0 J=1\n3 5000.0 5.0 J=2"
    rates = "1 2 1 1.0e-11 2.0e-11\n2 1 3 3.0e-11 4.0e-11\n3 3 2 5.0e-11 6.0e-11"
    path = _file(tmp_path, levels=levels, rates=rates, lines="1 2 1 7.203e-08 115.2711894 5.532\n2 2 1 7e-08")
    with pytest.raises(FixError, match="detailed balance overflows"):
        fix(path, inverted="balance")


# ------------------------------------------------------------------- line columns
def test_r04_r05_missing_columns(tmp_path):
    path = _file(tmp_path, lines="1 2 1 7.203e-08\n2 3 2 6.910e-07 230.5369826")
    text, log = _same(path)
    assert [(r.code, r.count) for r in log] == [("R04", 1), ("R05", 2)]
    out = text.splitlines()
    assert out[13].split() == ["1", "2", "1", "7.203e-08", "115.2711894", "5.53214"]
    assert out[14].split() == ["2", "3", "2", "6.910e-07", "230.5369826", "16.5962"]
    assert not {"R04", "R05"} & _codes(text)
    assert patch(path, missing_columns=False)[1] == []


def test_r04_recompute(tmp_path):
    path = _file(tmp_path, lines="1 2 1 7.203e-08 116.0 5.532\n2 3 2 6.910e-07 230.5369826 16.596")
    assert patch(path)[1] == []  # a value that is there needs --recompute
    text, log = _same(path, recompute=True)
    assert [(r.code, r.count) for r in log] == [("R04", 1)]
    assert text.splitlines()[13].split()[4] == "115.2711894"


# ------------------------------------------------------------------- references
def test_refs(tmp_path):
    path = _file(tmp_path, nlev="3", partner="2 X + p-H2 ! Yang et al. (2010)")
    refs = ["levels=10.1016/j.molstruc.2005.01.027", "p-H2=2010ApJ...718.1062Y", "dataset=10.17178/EMAA_X_1"]
    text, log = _same(path, refs=refs)
    out = text.splitlines()
    assert out[1] == "X -- EMAA_DATASET_DOI: 10.17178/EMAA_X_1"
    assert out[5] == "3 - Ref: doi:10.1016/j.molstruc.2005.01.027"
    assert out[18] == "2 X + p-H2 ! Yang et al. (2010) https://scixplorer.org/abs/2010ApJ...718.1062Y"
    assert [r.code for r in log] == ["F03", "F03", "F01"]
    assert not {"F01", "F02", "F03"} & _codes(text)
    again, log = patch(read_string(text), refs=refs)  # already there: not added again
    assert log == [] and again == text


def test_ref_to_text_with_ref(tmp_path):
    path = _file(tmp_path, nlev="3 - ROTATION - Ref: CDMS", partner="2 X + p-H2")
    text, log = patch(path, refs={"levels": "arXiv:2004.11230", "pH2": "doi:10.1000/x"})
    assert text.splitlines()[5] == "3 - ROTATION - Ref: CDMS arXiv:2004.11230"
    assert text.splitlines()[18] == "2 X + p-H2 doi:10.1000/x"
    assert [r.code for r in log] == ["F03", "F02"]


def test_ref_errors(tmp_path):
    path = _file(tmp_path)
    with pytest.raises(FixError, match="not one DOI, bibcode or arXiv ID"):
        fix(path, refs=["levels=Yang 2010"])
    with pytest.raises(FixError, match="no collision partner 'He'"):
        fix(path, refs=["He=10.1000/x"])
    with pytest.raises(FixError, match="use TARGET=ID"):
        fix(path, refs=["10.1000/x"])
    with pytest.raises(FixError, match="must be a DOI"):
        fix(path, refs=["dataset=2010ApJ...718.1062Y"])
    two = _file(tmp_path, npart="2", after=collision_block("2 X + p-H2 ! a"))
    with pytest.raises(FixError, match="ambiguous"):
        fix(two, refs=["p-H2=10.1000/x"])


def test_add_reference():
    d = read(LAMDA / "catom.dat")
    out, log = add_reference(d, "He", bibcode="1991JPhB...24.2343S")
    assert out.partner("He").reference.endswith(" https://scixplorer.org/abs/1991JPhB...24.2343S")
    assert out.notes[-1].endswith(
        "fix: F03 collision block 4 (He): https://scixplorer.org/abs/1991JPhB...24.2343S added"
    )
    assert "scixplorer" not in d.partner("He").reference  # the input is not changed
    with pytest.raises(FixError, match="exactly one"):
        add_reference(d, "He")


# ------------------------------------------------------------------- general
def test_options_and_sources(tmp_path):
    path = _file(tmp_path)
    for kw, msg in (({"trust": "both"}, "trust"), ({"inverted": "x"}, "inverted"), ({"extra_blocks": "x"}, "extra")):
        with pytest.raises(FixError, match=msg):
            fix(path, **kw)
    d = read(path)
    d.source_map = None
    with pytest.raises(FixError, match="patch mode needs data read from a file"):
        patch(d)
    data, log = fix(d)  # data without a file: the data fixes still work
    assert log == []


def test_notes_in_normalised_output(tmp_path):
    path = _file(tmp_path, ncol_comment="NUMBER OF COLL TRANS", partner=f"1 X + para-H2 ! Smith 2001 {DOI}")
    data, log = fix(path, trust="name")
    assert [r.code for r in log] == ["S05", "C02"]
    assert data.notes[-2:] == [r.note() for r in log]
    pylamda.write(data, tmp_path / "out.dat")
    assert "! pylamda" in (tmp_path / "out.dat").read_text() and read(tmp_path / "out.dat").issues == []


# ------------------------------------------------------------------- command line
def test_fix_command(tmp_path, capsys):
    path = _file(tmp_path, partner=f"1 X + para-H2 ! Smith 2001 {DOI}", ncol_comment="NUMBER OF COLL TRANS")
    original = path.read_text()
    assert main(["fix", str(path), "--dry-run"]) == 0
    captured = capsys.readouterr()
    assert "-NUMBER OF COLL TRANS\n+!NUMBER OF COLL TRANS" in captured.out
    assert "C02 remains; it needs --trust id|name" in captured.err and path.read_text() == original
    out = tmp_path / "out.dat"
    assert main(["fix", str(path), "-o", str(out), "--trust", "name"]) == 0  # normalised
    assert out.read_text().startswith("!MOLECULE\nX\n!MOLECULAR WEIGHT") and "\n2 X + p-H2 -- Ref:" in out.read_text()
    assert path.read_text() == original  # input not touched
    assert main(["fix", str(path), "-o", str(out), "--patch"]) == 0
    assert out.read_text().count("\n") == original.count("\n") + 1  # one note line
    capsys.readouterr()
    assert main(["fix", str(path), "--in-place", "--trust", "id"]) == 0  # patch mode
    assert "\n1 X + H2 -- Ref: Smith 2001" in path.read_text() and "remains" not in capsys.readouterr().err
    assert main(["fix", str(path), "--in-place"]) == 0
    assert "nothing to fix" in capsys.readouterr().err
    assert main(["fix", str(tmp_path / "missing.dat"), "-o", str(out)]) == 2
    assert main(["fix", str(path), "-o", str(out), "--ref", "levels=nonsense"]) == 2
    with pytest.raises(SystemExit):
        main(["fix", str(path)])  # -o, --in-place or --dry-run is required
