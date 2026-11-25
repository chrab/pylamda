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
        return [i for i in verify(LAMDA / "catom.dat", select=["F03"], **kw) if i.message.startswith("partner ")]

    assert partner_f03() == []  # all six partner references are in the catalogue
    assert len(partner_f03(catalogue=False)) == 6


def test_cli_no_catalogue(capsys):
    assert main(["check", "--select", "F03", str(LAMDA / "catom.dat")]) == 0
    assert ": partner " not in capsys.readouterr().out
    main(["check", "--select", "F03", "--no-catalogue", str(LAMDA / "catom.dat")])
    assert 'F03 info: partner "C + H" (ID 5 = H)' in capsys.readouterr().out
    assert main(["bib", str(LAMDA / "catom.dat"), "--no-database-refs"]) == 0
    captured = capsys.readouterr()
    assert "@ARTICLE{Launay1977," in captured.out and "(5 catalogue)" in captured.err  # p-H2 and o-H2 share one paper
    assert pylamda.bibtex is bibtex


# ------------------------------------------------------------------- identifiers from the catalogue
def test_bib_uses_catalogue_for_known_identifiers(tmp_path):
    """A bibcode or DOI the catalogue knows gets its reviewed BibTeX, offline."""
    path = tmp_path / "x.dat"
    path.write_text(make_lamda(molecule="CO", nlev="3", partner="2 CO + p-H2 ! doi:10.1088/0004-637X/718/2/1062"))
    result = bibtex([path], database_refs=False)
    assert [(e.key, e.source) for e in result.entries] == [("Yang2010", "catalogue")]
    assert bibtex([path], database_refs=False, catalogue=False).entries[0].source == "unresolved"


def test_links_with_ampersand():
    from pylamda.references import find_bibcodes

    assert find_bibcodes("x https://scixplorer.org/abs/1977A&A....56..289L y") == ["1977A&A....56..289L"]
    assert find_bibcodes("https://scixplorer.org/abs/1977A%26A....56..289L") == ["1977A&A....56..289L"]


# ------------------------------------------------------------------- fix --add-refs
def test_fix_add_refs_from_the_catalogue(tmp_path):
    text, log = pylamda.patch(LAMDA / "catom.dat", add_refs=True)
    assert [r.code for r in log] == ["F03"] * 6
    assert "\n5 C + H  ! Launay & Roueff 1977 https://scixplorer.org/abs/1977A%26A....56..289L\n" in text
    assert log[0].message.endswith("added (reference catalogue)")
    fixed = tmp_path / "catom.dat"
    fixed.write_text(text)
    issues = verify(fixed, select=["F"], catalogue=False)
    assert [i.message for i in issues if i.message.startswith("partner ")] == []  # the links are identifiers now
    assert [i.code for i in issues] == ["F03"]  # the levels: still free text, the fix notes do not count
    assert pylamda.patch(fixed, add_refs=True)[1] == []  # a second run adds nothing
    assert {e.source for e in bibtex([fixed], database_refs=False).entries} == {"catalogue"}


def test_fix_add_refs_only_place_or_molecule(tmp_path):
    """A text found only by its text (another molecule) is not written into a file."""
    path = tmp_path / "x.dat"
    path.write_text(make_lamda(molecule="X", nlev="3", partner="2 X + p-H2 ! Launay & Roueff 1977"))
    assert pylamda.patch(path, add_refs=True)[1] == []


def test_fix_add_refs_from_a_review_file(tmp_path):
    from pylamda.refsearch import Group, Use, save

    path = tmp_path / "x.dat"
    path.write_text(make_lamda(molecule="X", nlev="3", partner="2 X + p-H2 ! Smith 1999 & Jones 2001"))
    use = Use("user/x.dat/collisions/p-H2#1", "user", "x.dat", "collisions", "p-H2", None, "Smith 1999 & Jones 2001")
    use2 = Use("user/x.dat/collisions/p-H2#2", "user", "x.dat", "collisions", "p-H2", None, "Smith 1999 & Jones 2001")
    review = tmp_path / "review.json"
    save(
        [
            Group("X|smith 1999", "X", "Smith 1999", "smith 1999", [use], decision="1999ApJ...500..100S"),
            Group("X|jones 2001", "X", "Jones 2001", "jones 2001", [use2], proposal="2001ApJ...600..200J"),
        ],
        review,
    )
    text, log = pylamda.patch(path, add_refs=review)
    assert [r.message for r in log] == [
        "collision block 1 (p-H2): https://scixplorer.org/abs/1999ApJ...500..100S added (review file)"
    ]  # Jones 2001 has only a proposal, no decision: not added
    data, log2 = pylamda.fix(path, add_refs=review)
    assert data.collisions[0].reference.endswith("https://scixplorer.org/abs/1999ApJ...500..100S")
    # the same from the command line, with the option before FILE
    out = tmp_path / "out.dat"
    assert main(["fix", "--add-refs-from", str(review), "-o", str(out), str(path)]) == 0
    assert pylamda.read(out).collisions[0].reference.endswith("https://scixplorer.org/abs/1999ApJ...500..100S")


def test_cli_add_refs_before_the_file(tmp_path, capsys):
    """--add-refs takes no value: a FILE after it is the file to fix (not a review file)."""
    out = tmp_path / "catom.dat"
    assert main(["fix", "--trust", "name", "--add-refs", str(LAMDA / "catom.dat"), "-o", str(out)]) == 0
    assert "scixplorer.org/abs/1977A%26A....56..289L" in out.read_text()
    capsys.readouterr()
    with pytest.raises(SystemExit):  # only one of the two
        main(["fix", "--add-refs", "--add-refs-from", "r.json", str(LAMDA / "catom.dat"), "--dry-run"])
    assert "not allowed with argument" in capsys.readouterr().err


# ------------------------------------------------------------------- --search
YANG_DOC = {
    "bibcode": "2010ApJ...718.1062Y",
    "title": ["Rotational Quenching of CO due to H<SUB>2</SUB> Collisions"],
    "author": ["Yang, Benhui"],
    "year": "2010",
    "volume": "718",
    "page": ["1062"],
    "doctype": "article",
}
SMITH_DOC = {
    "bibcode": "1999ApJ...500..100S",
    "title": ["Excitation of CO by H2"],
    "author": ["Smith, A."],
    "year": "1999",
    "volume": "500",
    "page": ["100"],
    "doctype": "article",
}


class SearchClient:
    has_token = True

    def __init__(self):
        self.queries = []

    def scix_search(self, query, fields):
        self.queries.append(query)
        if 'author:"^Smith"' in query:
            return [SMITH_DOC]
        if 'author:"^Yang"' in query:
            return [YANG_DOC]
        return []


def test_search_references(tmp_path):
    from pylamda.refsearch import load, merge_review, search_references

    path = tmp_path / "co.dat"
    path.write_text(make_lamda(molecule="CO", nlev="3", partner="2 CO + p-H2 ! from Yang et al. (2010) & Smith 1999"))
    client = SearchClient()
    groups = search_references([path], client=client, hints=False, catalogue=load_catalogue())
    assert [g.segment for g in groups] == ["Smith 1999"]  # Yang 2010 is in the catalogue
    assert groups[0].candidates[0]["bibcode"] == "1999ApJ...500..100S" and groups[0].proposal
    assert all('"^Yang"' not in q for q in client.queries)
    review = tmp_path / "review.json"
    merge_review(groups, review)
    first = load(review)
    first[0].decision = "ok"
    from pylamda.refsearch import save

    save(first, review)
    again = search_references([path], client=SearchClient(), hints=False, catalogue=load_catalogue())
    merged = merge_review(again, review)
    assert len(merged) == 1 and merged[0].decision == "ok"  # the decision is kept


def test_bib_search_command(tmp_path, monkeypatch, capsys):
    import pylamda.bib as bib
    from pylamda import refsearch

    path = tmp_path / "co.dat"
    path.write_text(make_lamda(molecule="CO", nlev="3", partner="2 CO + p-H2 ! from Yang et al. (2010) & Smith 1999"))
    review = tmp_path / "review.json"
    monkeypatch.setenv("SCIX_API_TOKEN", "test")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setattr(bib, "HttpClient", lambda token: SearchClient())
    monkeypatch.setattr(refsearch, "_get", lambda url, cache: None)
    assert main(["bib", str(path), "--no-database-refs", "--search", "--review", str(review)]) == 0
    captured = capsys.readouterr()
    assert captured.out.count("@") == 1  # only Yang2010 from the catalogue: the search adds nothing
    assert "--search: 1 free-text citations searched" in captured.err
    assert "-> proposal 1999ApJ...500..100S" in captured.err
    assert refsearch.load(review)[0].segment == "Smith 1999"
    monkeypatch.delenv("SCIX_API_TOKEN")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    assert main(["bib", str(path), "--search", "--review", str(review)]) == 1  # no token
