"""The maintainer tool tools/catalogue.py (reference catalogue), without network."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
from helpers import make_lamda

sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
import catalogue as cat  # noqa: E402

PAGE = """<html><body>
<p>The collisional rate coefficients are adopted from
<a href="https://ui.adsabs.harvard.edu/abs/2010ApJ...718.1062Y/abstract">Yang et al. (2010)</a>
and the potential of <a href="https://ui.adsabs.harvard.edu/abs/2005JChPh.123j4301J/abstract">Jankowski
and Szalewicz (2005)</a>. See also <a href="https://doi.org/10.1000/xyz">Smith et al. (1999)</a>.</p>
<p>CDMS (Oct 2000) M&uuml;ller et al. J. Mol. Struct. 742 215 (2005)
<a href="https://emaa.osug.fr/details-doi/1">DOI</a></p>
</body></html>"""

YANG = {
    "bibcode": "2010ApJ...718.1062Y",
    "title": ["Rotational Quenching of CO due to H<SUB>2</SUB> Collisions"],
    "author": ["Yang, Benhui", "Stancil, P. C."],
    "year": "2010",
    "volume": "718",
    "page": ["1062"],
    "pub": "The Astrophysical Journal",
    "doctype": "article",
}
OTHER = {"bibcode": "2010NatNa...5..143Y", "title": ["Above-bandgap voltages"], "author": ["Yang, S. Y."],
         "year": "2010", "volume": "5", "page": ["143"], "doctype": "article"}  # fmt: skip


class FakeClient:
    def __init__(self):
        self.queries = []

    def scix_search(self, query, fields):
        self.queries.append(query)
        if query.startswith('author:"^Yang"'):
            return [OTHER, YANG]
        if query.startswith('author:"^Smith"'):
            return [{**OTHER, "bibcode": "1999ApJ...500..100S", "author": ["Smith, A."], "volume": "500",
                     "page": ["100"]}]  # fmt: skip
        return []


# ------------------------------------------------------------------- parsing
@pytest.mark.parametrize(
    ("text", "parsed"),
    [
        ("from Yang et al. (2010).", ("Yang", [], 2010, None, None, None)),
        ("Roueff & Le Bourlot 1990, A&A, 236, 515", ("Roueff", ["Le Bourlot"], 1990, "A&A", "236", "515")),
        ("K. Schroeder et al. 1991, J. Phys. B, 24, 2487", ("Schroeder", [], 1991, "JPhB", "24", "2487")),
        (
            "Balakrishnan N., Yan M., Dalgarno A.. Quantum-Mechanical Study ... Ap. J. vol568 : p443--447, 2002",
            ("Balakrishnan", ["Yan", "Dalgarno"], 2002, "ApJ", "568", "443"),
        ),
        (
            "Staemmler, V., Flower, D.R. 1991, J. Phys. B, 24, 2343",
            ("Staemmler", ["Flower"], 1991, "JPhB", "24", "2343"),
        ),
        ("Lique et al. 2010, J. Chem. Phys., 132, 024303", ("Lique", [], 2010, "JChPh", "132", "24303")),
    ],
)
def test_parse(text, parsed):
    p = cat.parse(text)
    assert (p["first_author"], p["authors"], p["year"], p["bibstem"], p["volume"], p["page"]) == parsed


@pytest.mark.parametrize(
    ("text", "first", "others"),
    [
        ("C2H-electron Born approximation + IOS recoupling (Faure 2014)", "Faure", []),
        ("Storzer Hollenbach 2000", "Storzer", ["Hollenbach"]),
        ("He+OH Klos, Lique, Alexander 2007", "Klos", ["Lique", "Alexander"]),
        ("Van der Tak et al 2013)", "Van der Tak", []),
        ("OH-pH2 (van den Heuvel et al 2025)", "van den Heuvel", []),
    ],
)
def test_first_author(text, first, others):
    p = cat.parse(text)
    assert (p["first_author"], p["authors"]) == (first, others)


def test_journal_after_the_year():
    p = cat.parse("rates for He scaled by 1.39 Hammami et al. 2009 A&A 507, 1083 & Turpin et al. A&A 511, A28")
    assert (p["bibstem"], p["volume"], p["page"]) == ("A&A", "507", "1083")
    p = cat.parse("Spirko et al. J. Phys B 36, 1645, 2003")
    assert (p["bibstem"], p["volume"], p["page"]) == ("JPhB", "36", "1645")


def test_split_segments():
    assert cat.split_segments("from Flower (2001) & Wernli et al. (2006) + extrapolation") == [
        "from Flower (2001)",
        "Wernli et al. (2006) + extrapolation",
    ]
    assert cat.split_segments("Yang et al. (2010)") == ["Yang et al. (2010)"]
    assert cat.split_segments("Abrahamsson, Krems, Dalgarno (2007, ApJ, 654, 1171) + Storzer Hollenbach 2000") == [
        "Abrahamsson, Krems, Dalgarno (2007, ApJ, 654, 1171)",
        "Storzer Hollenbach 2000",
    ]
    assert cat.split_segments("Lique et al. (2013), Wiesenfeld & Goldsmith (2014)") == [
        "Lique et al. (2013)",
        "Wiesenfeld & Goldsmith (2014)",
    ]
    # a part without an author is not a citation of its own
    assert cat.split_segments("Spirko et al. J. Phys B 36, 1645, 2003 by E. Roueff, oct 2019 at same T") == [
        "Spirko et al. J. Phys B 36, 1645, 2003 by E. Roueff, oct 2019 at same T"
    ]
    assert (
        cat.normalise("from Yang et al. (2010).") == cat.normalise("from  Yang et al (2010)") == "from yang et al 2010"
    )


# ------------------------------------------------------------------- titles and pages
@pytest.mark.parametrize(
    ("molecule", "title", "expected"),
    [
        ("C", "Fine-structure excitation of ground-state C/+/ ions by hydrogen atoms", -1),
        ("C", "Fine structure excitation of carbon and oxygen by atomic hydrogen impact", 1),
        ("C", "Excitation of forbidden C I fine structure transitions by protons", 1),
        ("C+", "Fine-structure excitation of ground-state C/+/ ions by hydrogen atoms", 1),
        ("C+", "Excitation of C II by electrons", 1),
        ("CO", "Rotational Quenching of CO due to H<SUB>2</SUB> Collisions", 1),
        ("CO", "Excitation of CO<SUP>+</SUP>", -1),
        ("o-H2O", "Rotational excitation of H<SUB>2</SUB>O by H<SUB>2</SUB>", 1),
        ("HCO+", "Collisional rates for HCO<SUP>+</SUP>", 1),
        ("CO", "Carbon monoxide in comets", 0),
    ],
)
def test_title_match(molecule, title, expected):
    assert cat.title_match(molecule, title) == expected


def test_page_names():
    assert cat.page_names("e-CH3OH") == ["e-CH3OH", "CH3OH"]
    assert cat.page_names("C18O") == ["C18O", "CO"]
    assert cat.page_names("H13CO+") == ["H13CO+", "HCO+"]
    assert cat.page_names("c-C3H2") == ["c-C3H2", "C3H2"]
    assert cat.page_names("CO") == ["CO"]


def test_page_hints():
    hints = cat.page_hints(PAGE, "LAMDA page CO")
    assert [(h["text"], h["bibcode"], h["doi"]) for h in hints[:3]] == [
        ("Yang et al. (2010)", "2010ApJ...718.1062Y", None),
        ("Jankowski and Szalewicz (2005)", "2005JChPh.123j4301J", None),
        ("Smith et al. (1999)", None, "10.1000/xyz"),
    ]
    emaa = hints[3]  # the citation is the text before an EMAA "DOI" link
    assert emaa["url"] == "https://emaa.osug.fr/details-doi/1" and "Müller et al. J. Mol. Struct. 742" in emaa["text"]


def test_hint_matches():
    g = cat.Group(key="k", molecule="CO", segment="", normalised="", first_author="Yang", year=2010)
    assert cat.hint_matches(g, {"text": "Yang et al. (2010)", "bibcode": None})
    assert not cat.hint_matches(g, {"text": "Yang et al. (2011)", "bibcode": None})
    assert cat.hint_matches(g, {"text": "the CO-H2 paper", "bibcode": None}) is False
    assert cat.hint_matches(g, {"text": "Yang, B.", "bibcode": "2010ApJ...718.1062Y"})  # year from the bibcode
    g2 = cat.Group(key="k", molecule="C", segment="", normalised="", first_author="Schroeder", year=1991)
    assert cat.hint_matches(g2, {"text": "Schröder 1991", "bibcode": None})  # ö -> oe
    g4 = cat.Group(key="k", molecule="C", segment="", normalised="", first_author="Schröder", year=1991)
    assert cat.hint_matches(g4, {"text": "Schroder et al. 1991", "bibcode": None})  # ö -> o
    g3 = cat.Group(key="k", molecule="C", segment="", normalised="", first_author="Muller", year=2005)
    assert cat.hint_matches(g3, {"text": "Müller et al. 2005", "bibcode": None})


# ------------------------------------------------------------------- collecting, searching, proposals
def _files(tmp_path):
    a = tmp_path / "a.dat"
    a.write_text(make_lamda(molecule="CO", nlev="3", partner="2 CO + p-H2 ! from Yang et al. (2010)"))
    b = tmp_path / "b.dat"
    b.write_text(make_lamda(molecule="CO", nlev="3", partner="1 CO + p-H2 ! from Yang et al. (2010)."))
    return a, b


def test_collect_groups(tmp_path):
    a, b = _files(tmp_path)
    groups, datas = cat.collect([str(a), str(b)], "prodimo")
    assert len(groups) == 1  # the same text for the same molecule: one decision
    g = groups[0]
    assert g.key == "CO|from yang et al 2010" and (g.first_author, g.year) == ("Yang", 2010)
    assert [(u.file, u.partner, u.named) for u in g.uses] == [("a.dat", "p-H2", None), ("b.dat", "H2", "p-H2")]


def test_search_and_proposals(tmp_path):
    a, _ = _files(tmp_path)
    groups, datas = cat.collect([str(a)], "prodimo")
    g = groups[0]
    g.hints = [h for h in cat.page_hints(PAGE, "LAMDA page CO") if cat.hint_matches(g, h)]
    cat.search(groups, FakeClient())
    assert g.candidates[0]["bibcode"] == "2010ApJ...718.1062Y" and "named by a hint" in g.candidates[0]["why"]
    assert (g.proposal, g.proposal_reason) == ("2010ApJ...718.1062Y", "named by: LAMDA page CO")
    # without the hint: only the weaker "best candidate" proposal (CO in the title, clear lead)
    g.hints = []
    cat.search(groups, FakeClient())
    assert (g.proposal, g.proposal_reason) == ("2010ApJ...718.1062Y", cat.BEST_REASON)


def test_proposal_from_volume_and_page():
    g = cat.Group(key="k", molecule="X", segment="", normalised="", first_author="Smith", year=1999,
                  volume="500", page="100")  # fmt: skip
    cat.search([g], FakeClient())
    assert (g.proposal, g.proposal_reason) == ("1999ApJ...500..100S", "volume and page match")


def test_is_first_author():
    assert cat.is_first_author("Yang", "Yang, Benhui")
    assert cat.is_first_author("Hernandez", "Hernández Vera, Mario")
    assert cat.is_first_author("Van der Tak", "van der Tak, F. F. S.")
    assert cat.is_first_author("Schroeder", "Schröder, K.")
    assert not cat.is_first_author("Yang", "Yangtze, A.")


def test_best_candidate_needs_a_lead():
    g = cat.Group(key="k", molecule="CO", segment="", normalised="")
    g.candidates = [{"bibcode": "a", "score": 0.4, "why": "CO in title"}, {"bibcode": "b", "score": 0.3, "why": ""}]
    cat.propose(g)
    assert g.proposal is None  # lead 0.1 < 0.15
    g.candidates[1]["score"] = 0.2
    cat.propose(g)
    assert (g.proposal, g.proposal_reason) == ("a", cat.BEST_REASON)


def test_isotopologue_title():
    g = cat.Group(key="k", molecule="13CO", segment="", normalised="")
    assert "13CO in title" in cat.score(g, {"bibcode": "x", "title": ["Excitation of CO by H2"]})[1]


def test_decisions():
    g = cat.Group(key="k", molecule="X", segment="", normalised="", proposal="2010ApJ...718.1062Y")
    assert g.bibcode is None
    g.decision = "ok"
    assert g.bibcode == "2010ApJ...718.1062Y"
    g.decision = "1999ApJ...500..100S"
    assert g.bibcode == "1999ApJ...500..100S"
    g.decision = "none"
    assert g.bibcode is None


def test_main(tmp_path, monkeypatch, capsys):
    a, b = _files(tmp_path)
    review = tmp_path / "review.json"
    monkeypatch.setenv("SCIX_API_TOKEN", "test")
    monkeypatch.setattr(cat, "HttpClient", lambda token: FakeClient())
    monkeypatch.setattr(cat, "_get", lambda url, cache: PAGE if url.endswith("/CO.html") else None)
    assert cat.main(["search", "--database", "prodimo", str(a), "-o", str(review), "--cache-dir", str(tmp_path)]) == 0
    data = json.loads(review.read_text())
    assert len(data) == 1 and data[0]["proposal"] == "2010ApJ...718.1062Y" and data[0]["decision"] is None
    # the reviewer decides; a second run with more files keeps the decision and adds the use
    data[0]["decision"] = "ok"
    review.write_text(json.dumps(data))
    assert cat.main(["search", "--database", "prodimo", str(b), "-o", str(review), "--cache-dir", str(tmp_path)]) == 0
    data = json.loads(review.read_text())
    assert data[0]["decision"] == "ok" and [u["file"] for u in data[0]["uses"]] == ["b.dat", "a.dat"]
    capsys.readouterr()
    assert cat.main(["show", str(review), "--open"]) == 0
    assert capsys.readouterr().out == ""
    data[0]["decision"] = None
    review.write_text(json.dumps(data))
    assert cat.main(["confirm", str(review)]) == 0
    assert json.loads(review.read_text())[0]["decision"] == "ok"  # named by the LAMDA page: a strong proposal
    data[0]["proposal_reason"] = cat.BEST_REASON
    review.write_text(json.dumps(data))
    assert cat.main(["confirm", str(review)]) == 0
    assert json.loads(review.read_text())[0]["decision"] is None  # weak: only with --include-best
    assert cat.main(["confirm", str(review), "--include-best"]) == 0
    assert json.loads(review.read_text())[0]["decision"] == "ok"
    monkeypatch.delenv("SCIX_API_TOKEN")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    assert cat.main(["search", "--database", "x", str(a), "-o", str(review)]) == 1  # no token


# ------------------------------------------------------------------- build
class ExportClient:
    def scix_export(self, bibcodes):
        entries = {
            "2010ApJ...718.1062Y": "@ARTICLE{2010ApJ...718.1062Y,\n  author = {{Yang}, Benhui},\n  year = 2010,\n"
            "  doi = {10.1088/0004-637X/718/2/1062},\n  adsurl = {https://ui.adsabs.harvard.edu/abs/x}\n}",
            "2010ApJ...999....1Y": "@ARTICLE{2010ApJ...999....1Y,\n  author = {{Yang}, A.},\n  year = 2010,\n}",
        }
        return "\n\n".join(entries[b] for b in bibcodes if b in entries)


def _group(key, bibcode, decision, file="a.dat"):
    g = cat.Group(key=key, molecule="CO", segment="", normalised=key.split("|")[1], proposal=bibcode, decision=decision)
    g.uses = [cat.Use(id=f"prodimo/{file}/collisions/p-H2", database="prodimo", file=file, role="collisions",
                      partner="p-H2", named=None, text="from Yang et al. (2010)")]  # fmt: skip
    return g


def test_build(tmp_path):
    groups = [
        _group("CO|from yang et al 2010", "2010ApJ...718.1062Y", "ok"),
        _group("CO|yang 2010 other", None, "2010ApJ...999....1Y", file="b.dat"),
        _group("CO|flower 2001", None, "none", file="c.dat"),
        _group("CO|iron project", None, "skip", file="d.dat"),
        _group("CO|open one", "2010ApJ...718.1062Y", None, file="e.dat"),
    ]
    assert cat.build(groups, ExportClient(), tmp_path) == (2, 3)
    catalogue = json.loads((tmp_path / "references.json").read_text())
    assert catalogue["entries"] == {
        "Yang2010": {"bibcode": "2010ApJ...718.1062Y", "doi": "10.1088/0004-637X/718/2/1062", "reviewed": True},
        "Yang2010a": {"bibcode": "2010ApJ...999....1Y", "doi": None, "reviewed": True},
    }
    rows = {r["file"]: (r["key"], r["status"], r["text"]) for r in catalogue["lookup"]}
    assert rows == {
        "a.dat": ("Yang2010", "found", "from yang et al 2010"),
        "b.dat": ("Yang2010a", "found", "yang 2010 other"),
        "c.dat": (None, "not in SciX", "flower 2001"),
    }  # "skip" and undecided groups are not in the catalogue
    bib = (tmp_path / "references.bib").read_text()
    assert "@ARTICLE{Yang2010," in bib and "adsurl = {https://scixplorer.org/abs/2010ApJ...718.1062Y}" in bib
    # a rebuild keeps the keys, also when the order of the bibcodes changes
    groups[1].decision = "none"
    groups.append(_group("CO|yang again", None, "2010ApJ...999....1Y", file="f.dat"))
    cat.build(groups, ExportClient(), tmp_path)
    assert json.loads((tmp_path / "references.json").read_text())["entries"]["Yang2010a"]["bibcode"] == (
        "2010ApJ...999....1Y"
    )


def test_build_unknown_bibcode(tmp_path):
    with pytest.raises(SystemExit, match="not found in SciX: 1999XXX"):
        cat.build([_group("CO|x", None, "1999XXX.........1Y")], ExportClient(), tmp_path)
