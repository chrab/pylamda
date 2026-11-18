"""The reference catalogue at run time: lookup, bib and check."""

from __future__ import annotations

import pytest
from helpers import LAMDA, make_lamda

import pylamda
from pylamda import bibtex, verify
from pylamda.catalogue import Catalogue, load_catalogue, normalise, split_segments
from pylamda.cli import main

ROWS = [
    {"database": "lamda", "molecule": "CO", "file": "co.dat", "role": "collisions", "partner": "p-H2",
     "text": "from yang et al 2010", "key": "Yang2010", "status": "found"},
    {"database": "prodimo", "molecule": "HCN", "file": "hcn.dat", "role": "collisions", "partner": "He",
     "text": "green 1974", "key": "Green1974", "status": "found"},
    {"database": "prodimo", "molecule": "SiO", "file": "sio.dat", "role": "collisions", "partner": "He",
     "text": "green 1974", "key": "Green1974b", "status": "found"},
    {"database": "prodimo", "molecule": "CO", "file": "co.dat", "role": "collisions", "partner": "e",
     "text": "private communication 2001", "key": None, "status": "not in SciX"},
]  # fmt: skip
ENTRIES = {
    "Yang2010": {"bibcode": "2010ApJ...718.1062Y", "doi": "10.1088/0004-637X/718/2/1062", "reviewed": True},
    "Green1974": {"bibcode": "1974ApJ...191..653G", "doi": None, "reviewed": True},
    "Green1974b": {"bibcode": "1974ApJ...999....1G", "doi": None, "reviewed": True},
}


@pytest.fixture
def cat():
    return Catalogue(ENTRIES, ROWS, {"Yang2010": "@ARTICLE{Yang2010,\n  author = {{Yang}, B.},\n  year = 2010,\n}"})


def test_normalise_and_split():
    assert normalise("from Yang et al. (2010).") == "from yang et al 2010"
    assert split_segments("from Flower (2001) & Wernli et al. (2006)") == ["from Flower (2001)", "Wernli et al. (2006)"]


def test_find(cat):
    hit = cat.find("from Yang et al. (2010)", molecule="CO", file="co.dat", role="collisions", partner="p-H2")
    assert (hit.key, hit.bibcode, hit.match) == ("Yang2010", "2010ApJ...718.1062Y", "place")
    assert cat.find("from Yang et al (2010).", molecule="CO", file="other.dat").match == "molecule"
    assert cat.find("from Yang et al 2010", molecule="HCO+").match == "text"  # unique in the catalogue
    assert cat.find("Green 1974", molecule="HCN").key == "Green1974"  # by molecule
    assert cat.find("Green 1974", molecule="CS") is None  # two papers with this text: no guess
    hit = cat.find("private communication 2001", molecule="CO", file="co.dat", role="collisions", partner="e")
    assert hit.key is None and hit.match == "place"  # reviewed, not in SciX
    assert cat.find("private communication 2001", molecule="H2O") is None
    assert cat.find("unknown text", molecule="CO") is None


def test_shipped_catalogue():
    cat = load_catalogue()
    assert len(cat.entries) >= 70 and len(cat.rows) >= 150
    assert set(cat.bibtex) == set(cat.entries)
    assert all(r["key"] is None or r["key"] in cat.entries for r in cat.rows)
    hit = cat.find("Launay & Roueff 1977", molecule="C", file="catom.dat", role="collisions", partner="H")
    assert hit.bibcode == "1977A&A....56..289L"  # neutral C: the A&A paper, not the C+ one


def test_bib_uses_the_catalogue():
    result = bibtex([LAMDA / "catom.dat"], database_refs=False)
    assert result.free_text == []
    assert {e.source for e in result.entries} == {"catalogue"}
    keys = [e.key for e in result.entries]
    assert keys[0] == "Launay1977" and "Schroder1991" in keys
    assert "% C: collision rates with H (pylamda reference catalogue)\n@ARTICLE{Launay1977," in result.text()
    assert all(not w.startswith("catom.dat") for w in result.warnings)


def test_bib_partly_known(tmp_path):
    path = tmp_path / "x.dat"
    path.write_text(
        make_lamda(molecule="CO", nlev="3", partner="2 CO + p-H2 ! from Yang et al. (2010) & Nobody (1999)")
    )
    result = bibtex([path], database_refs=False)
    assert [e.key for e in result.entries] == ["Yang2010"]
    assert [r.text for r in result.free_text] == ["Nobody (1999)"]  # only the unknown part stays


def test_bib_text_only_match_warns(tmp_path):
    path = tmp_path / "x.dat"
    path.write_text(make_lamda(molecule="X", nlev="3", partner="2 X + p-H2 ! Launay & Roueff 1977"))
    result = bibtex([path], database_refs=False)
    # the text is in the catalogue only for C: found by its text alone, with a warning
    assert [e.key for e in result.entries] == ["Launay1977"] and result.free_text == []
    assert "found in the reference catalogue by its text only (Launay1977" in result.warnings[0]


def test_check_f03_uses_the_catalogue():
    def partner_f03(**kw):
        return [i for i in verify(LAMDA / "catom.dat", select=["F03"], **kw) if "collision block" in i.message]

    assert partner_f03() == []  # all six partner references are in the catalogue
    assert len(partner_f03(catalogue=False)) == 6


def test_cli_no_catalogue(capsys):
    assert main(["check", "--select", "F03", str(LAMDA / "catom.dat")]) == 0
    assert "collision block" not in capsys.readouterr().out
    main(["check", "--select", "F03", "--no-catalogue", str(LAMDA / "catom.dat")])
    assert "collision block 1 (H)" in capsys.readouterr().out
    assert main(["bib", str(LAMDA / "catom.dat"), "--no-database-refs"]) == 0
    captured = capsys.readouterr()
    assert "@ARTICLE{Launay1977," in captured.out and "(5 catalogue)" in captured.err  # p-H2 and o-H2 share one paper
    assert pylamda.bibtex is bibtex
