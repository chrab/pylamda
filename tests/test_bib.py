"""BibTeX export: identifiers, collection, offline entries, online resolution (fake and real)."""

from __future__ import annotations

import json
import os
import warnings

import pytest
from helpers import EMAA_CO, LAMDA, make_lamda

import pylamda
from pylamda import bibtex, collect_references, read
from pylamda.bib import (
    BibError,
    database_entries,
    entry_field,
    entry_key,
    identifier_context,
    is_emaa,
    key_from_entry,
    key_from_text,
    scix_token,
    split_entries,
)
from pylamda.cli import main
from pylamda.references import find_arxiv, find_bibcodes, find_dois, find_identifiers

YANG_DOI = "10.1088/0004-637X/718/2/1062"
YANG_BIBCODE = "2010ApJ...718.1062Y"
YANG_ADS = f"""@ARTICLE{{{YANG_BIBCODE},
       author = {{{{Yang}}, Benhui and {{Stancil}}, P.~C. and {{Balakrishnan}}, N. and {{Forrey}}, R.~C.}},
        title = "{{Rotational Quenching of CO due to H$_{{2}}$ Collisions}}",
      journal = {{\\apj}},
         year = 2010,
       volume = {{718}},
        pages = {{1062-1069}},
          doi = {{{YANG_DOI}}},
       eprint = {{1004.3923}},
       adsurl = {{https://ui.adsabs.harvard.edu/abs/{YANG_BIBCODE}}},
      adsnote = {{Provided by the SAO/NASA Astrophysics Data System}}
}}"""
CROSSREF = (
    "@article{Dagdigian_2022, title={Collisional excitation}, volume={514}, DOI={10.1093/mnras/stac1430}, "
    "journal={MNRAS}, author={Dagdigian, Paul J}, year={2022}, pages={2214–2219} }"
)


# ------------------------------------------------------------------- identifiers
@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("Dagdigian (2022) https://doi.org/10.1093/mnras/stac1430", ["doi:10.1093/mnras/stac1430"]),
        ("DOI:10.1063/1.555928", ["doi:10.1063/1.555928"]),
        ("doi: 10.1051/0004-6361:20041729.", ["doi:10.1051/0004-6361:20041729"]),
        ("(see 10.3390/atoms8020015)", ["doi:10.3390/atoms8020015"]),
        (f"Yang 2010 https://scixplorer.org/abs/{YANG_BIBCODE}", [f"bibcode:{YANG_BIBCODE}"]),
        ("https://ui.adsabs.harvard.edu/abs/2005A%26A...432..369S/abstract", ["bibcode:2005A&A...432..369S"]),
        (f"Yang 2010, {YANG_BIBCODE}", [f"bibcode:{YANG_BIBCODE}"]),
        ("arXiv:2004.11230v2 and arxiv.org/abs/astro-ph/0411110", ["arxiv:2004.11230", "arxiv:astro-ph/0411110"]),
        ("10.48550/arXiv.2004.11230", ["arxiv:2004.11230"]),  # the arXiv DOI is given as arXiv ID
        ("Lim et al 1999 MNRAS 306, 473", []),
    ],
)
def test_find_identifiers(text, expected):
    assert find_identifiers(text) == expected


def test_find_functions():
    assert find_dois("doi:10.1000/a and https://doi.org/10.1000/A") == ["10.1000/a"]  # DOIs are case-insensitive
    assert find_bibcodes(None) == [] and find_arxiv("") == []


# ------------------------------------------------------------------- collecting
def test_collect_references_emaa():
    d = read(EMAA_CO)
    refs = collect_references(d)
    assert [(r.role, r.partner) for r in refs][:7] == [
        ("dataset", None),
        ("levels", None),
        ("collisions", "p-H2"),
        ("collisions", "o-H2"),
        ("collisions", "p-H2"),
        ("collisions", "o-H2"),
        ("collisions", "H"),
    ]
    assert refs[0].identifiers == ["doi:10.17178/EMAA_CO_Rotation_3d0e1a5c"]
    assert refs[2].use == "collision rates with p-H2" and refs[2].identifiers == ["doi:10.1093/mnras/staa242"]
    assert refs[0].path == str(EMAA_CO) and is_emaa(d)
    assert not is_emaa(read(LAMDA / "catom.dat"))


def test_collect_references_comments(tmp_path):
    path = tmp_path / "x.dat"
    path.write_text(make_lamda(after="! notes: rates also in https://doi.org/10.1000/notes\n").replace(
        "!LEVEL + ENERGIES", "! levels from Smith 2001, doi:10.1000/levels\n!LEVEL + ENERGIES"
    ))  # fmt: skip
    refs = collect_references(read(path))
    comments = [r for r in refs if r.role == "comments"]
    assert [r.identifiers for r in comments] == [["doi:10.1000/notes"], ["doi:10.1000/levels"]]
    assert comments[1].text == "levels from Smith 2001, doi:10.1000/levels" and comments[1].use == "named in a comment"


# ------------------------------------------------------------------- BibTeX helpers
def test_bibtex_helpers():
    entries = split_entries(f"% comment\n{YANG_ADS}\n\n{CROSSREF}\n@comment{{x}}")
    assert [entry_key(e) for e in entries] == [YANG_BIBCODE, "Dagdigian_2022"]
    assert entry_field(entries[0], "doi") == YANG_DOI and entry_field(entries[0], "year") == "2010"
    assert entry_field(entries[0], "title") == "{Rotational Quenching of CO due to H$_{2}$ Collisions}"
    assert entry_field(entries[1], "author") == "Dagdigian, Paul J"
    assert key_from_entry(entries[0]) == "Yang2010" and key_from_entry(entries[1]) == "Dagdigian2022"
    assert key_from_entry("@misc{x, author = {{EMAA}}, year = {2023}}") == "EMAA2023"
    assert key_from_entry("@misc{x, title = {no author}}") is None
    assert key_from_text("Faure et al. MNRAS 493 776 (2020)") == "Faure2020"
    assert key_from_text("Launay & Roueff 1977") == "Launay1977"
    assert key_from_text("no year") is None


def test_identifier_context():
    """The key of each identifier is guessed from the text just before it."""
    text = "Flower 2001 doi:10.1000/a & Wernli et al. 2006 doi:10.1000/b; unknown doi:10.9999/c"
    assert key_from_text(identifier_context(text, "doi:10.1000/a")) == "Flower2001"
    assert key_from_text(identifier_context(text, "doi:10.1000/b")) == "Wernli2006"
    assert key_from_text(identifier_context(text, "doi:10.9999/c")) is None
    assert identifier_context("no id here", "doi:10.1000/x") == "no id here"


def test_database_entries():
    db = database_entries()
    assert set(db) == {"Schoier2005", "vanderTak2020", "Faure2025", "EMAA"}
    assert entry_field(db["Schoier2005"], "doi") == "10.1051/0004-6361:20041729"
    assert entry_field(db["EMAA"], "doi") == "10.17178/EMAA"


# ------------------------------------------------------------------- offline
def test_offline_export():
    result = bibtex([EMAA_CO, LAMDA / "catom.dat"])
    keys = [e.key for e in result.entries]
    # database papers first; Faure2025 and EMAA because of the EMAA file
    assert keys[:4] == ["Schoier2005", "vanderTak2020", "Faure2025", "EMAA"]
    by_key = {e.key: e for e in result.entries}
    faure = by_key["Faure2020"]  # key guessed from the text "Faure et al. MNRAS 493 776 (2020)"
    assert faure.source == "unresolved" and faure.identifiers == ["doi:10.1093/mnras/staa242"]
    assert faure.note == "CO: collision rates with p-H2, collision rates with o-H2"
    assert faure.bibtex == (
        "@misc{Faure2020,\n  note = {CO: collision rates with p-H2, collision rates with o-H2},\n"
        "  doi = {10.1093/mnras/staa242},\n}"
    )
    assert [r.partner for r in result.free_text] == ["H", "e", "H+", "He", "p-H2", "o-H2"]  # catom: free text
    text = result.text()
    assert text.startswith("% References of 2 LAMDA files, by pylamda")
    assert "% EMAA asks for this acknowledgement: This research has made use" in text
    assert "% CO: collision rates with p-H2, collision rates with o-H2\n@misc{Faure2020," in text
    assert result.keys["doi:10.1093/mnras/staa242"] == "Faure2020" and result.warnings == []


def test_offline_options(tmp_path):
    assert [e.key for e in bibtex([LAMDA / "catom.dat"]).entries] == ["Schoier2005", "vanderTak2020"]
    assert bibtex([LAMDA / "catom.dat"], database_refs=False).entries == []
    # bibcode and arXiv IDs give @misc entries with a SciX link / eprint
    path = tmp_path / "x.dat"
    path.write_text(
        make_lamda(partner=f"2 X + p-H2 ! Yang 2010 https://scixplorer.org/abs/{YANG_BIBCODE} arXiv:1004.3923")
    )
    entries = {e.identifiers[0]: e.bibtex for e in bibtex([path], database_refs=False).entries}
    assert (
        f"howpublished = {{\\url{{https://scixplorer.org/abs/{YANG_BIBCODE}}}}}" in entries[f"bibcode:{YANG_BIBCODE}"]
    )
    assert "eprint = {1004.3923}" in entries["arxiv:1004.3923"]


def test_database_doi_cited_by_a_file(tmp_path):
    """A file that cites a database paper gets no second entry; its use goes to the database entry."""
    path = tmp_path / "x.dat"
    path.write_text(make_lamda(partner="2 X + p-H2 ! as in Schoier et al. 2005, doi:10.1051/0004-6361:20041729"))
    result = bibtex([path])
    assert [e.key for e in result.entries] == ["Schoier2005", "vanderTak2020"]
    assert result.entries[0].note == "LAMDA database; X: collision rates with p-H2"


def test_keys_unique_and_stable(tmp_path):
    a, b = tmp_path / "a.dat", tmp_path / "b.dat"
    a.write_text(make_lamda(partner="2 X + p-H2 ! Yang et al. 2010, doi:10.1000/one"))
    b.write_text(make_lamda(molecule="Y", partner="2 Y + p-H2 ! Yang et al. 2010, doi:10.1000/two"))
    first = bibtex([a, b], database_refs=False)
    assert [e.key for e in first.entries] == ["Yang2010", "Yang2010a"]
    # the keys of an earlier run are kept, also if the order changes
    again = bibtex([b, a], database_refs=False, keys=first.keys)
    assert {e.identifiers[0]: e.key for e in again.entries} == {
        "doi:10.1000/two": "Yang2010a",
        "doi:10.1000/one": "Yang2010",
    }


# ------------------------------------------------------------------- token
def test_scix_token(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    monkeypatch.delenv("SCIX_API_TOKEN", raising=False)
    assert scix_token() is None
    path = tmp_path / "pylamda" / "scix_token"
    path.parent.mkdir()
    path.write_text("from-file\n")
    path.chmod(0o600)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert scix_token() == "from-file"
    path.chmod(0o644)
    with pytest.warns(UserWarning, match="can be read by other users"):
        assert scix_token() == "from-file"
    monkeypatch.setenv("SCIX_API_TOKEN", " from-env ")
    assert scix_token() == "from-env"
    assert scix_token("given") == "given"


# ------------------------------------------------------------------- online, with a fake client
class FakeClient:
    """SciX and doi.org without network: knows Yang (2010) and Dagdigian (2022)."""

    def __init__(self, has_token: bool = True):
        self.has_token = has_token
        self.calls: list[tuple] = []

    def scix_search(self, query: str, fields: str) -> list[dict]:
        assert self.has_token
        self.calls.append(("search", query))
        docs = []
        if YANG_DOI in query or "arXiv:1004.3923" in query:
            docs.append({"bibcode": YANG_BIBCODE, "doi": [YANG_DOI], "identifier": ["arXiv:1004.3923"]})
        return docs

    def scix_export(self, bibcodes: list[str]) -> str:
        assert self.has_token
        self.calls.append(("export", tuple(bibcodes)))
        return YANG_ADS if YANG_BIBCODE in bibcodes else ""

    def doi_bibtex(self, doi: str) -> str | None:
        self.calls.append(("doi.org", doi))
        return CROSSREF if doi == "10.1093/mnras/stac1430" else None


def _file(tmp_path, partner: str):
    path = tmp_path / "x.dat"
    path.write_text(make_lamda(partner=partner))
    return path


def test_online_resolution(tmp_path):
    path = _file(
        tmp_path,
        f"2 X + p-H2 ! Yang 2010 doi:{YANG_DOI}, arXiv:1004.3923; Dagdigian 2022 doi:10.1093/mnras/stac1430; "
        "unknown doi:10.9999/none",
    )
    client = FakeClient()
    result = bibtex([path], online=True, database_refs=False, cache=tmp_path / "cache.json", client=client)
    by_key = {e.key: e for e in result.entries}
    assert set(by_key) == {"Yang2010", "Dagdigian2022", "doi_10.9999_none"}
    yang = by_key["Yang2010"]
    assert yang.source == "scix" and set(yang.identifiers) == {f"doi:{YANG_DOI}", "arxiv:1004.3923"}  # merged
    assert yang.bibtex.startswith("@ARTICLE{Yang2010,")
    assert f"adsurl = {{https://scixplorer.org/abs/{YANG_BIBCODE}}}" in yang.bibtex
    assert by_key["Dagdigian2022"].source == "doi.org" and by_key["Dagdigian2022"].bibtex.startswith(
        "@article{Dagdigian2022,"
    )
    unknown = [e for e in result.entries if e.identifiers == ["doi:10.9999/none"]][0]
    assert unknown.source == "unresolved" and unknown.bibtex.startswith("@misc{")
    assert result.warnings == ["doi:10.9999/none: not found in SciX or at doi.org"]
    assert sum(1 for c in client.calls if c[0] == "export") == 1  # one export request for all bibcodes
    # the cache: a second run needs no request; the token is never in it
    cache = json.loads((tmp_path / "cache.json").read_text())
    assert (
        cache["bibtex"][YANG_BIBCODE].startswith("@ARTICLE") and cache["doi_bibcode"][YANG_DOI.lower()] == YANG_BIBCODE
    )
    assert cache["doi_bibcode"]["10.1093/mnras/stac1430"] == ""  # not in SciX
    again = FakeClient()
    result2 = bibtex([path], online=True, database_refs=False, cache=tmp_path / "cache.json", client=again)
    assert [e.key for e in result2.entries] == [e.key for e in result.entries]
    assert [c[0] for c in again.calls] == ["doi.org"]  # only the DOI that nobody knows is asked again


def test_online_without_token(tmp_path):
    path = _file(tmp_path, f"2 X + p-H2 ! Dagdigian 2022 doi:10.1093/mnras/stac1430, {YANG_BIBCODE}")
    client = FakeClient(has_token=False)
    result = bibtex([path], online=True, database_refs=False, cache=False, client=client)
    assert [e.source for e in result.entries] == ["doi.org", "unresolved"]
    assert result.warnings[0].startswith("1 bibcodes/arXiv IDs not resolved: SciX needs an API token")
    assert all(c[0] == "doi.org" for c in client.calls)


def test_online_error_is_raised(tmp_path):
    class Rejecting(FakeClient):
        def scix_search(self, query, fields):
            raise BibError("SciX rejected the API token (HTTP 401)")

    path = _file(tmp_path, f"2 X + p-H2 ! doi:{YANG_DOI}")
    with pytest.raises(BibError, match="HTTP 401"):
        bibtex([path], online=True, cache=False, client=Rejecting())


# ------------------------------------------------------------------- command line
def test_bib_command(tmp_path, capsys):
    out = tmp_path / "refs.bib"
    keys = tmp_path / "keys.json"
    assert main(["bib", str(EMAA_CO), "-o", str(out), "--keys", str(keys)]) == 0
    err = capsys.readouterr().err
    assert "entries to" in err and "unresolved" in err
    text = out.read_text()
    assert "@article{Schoier2005," in text and "@misc{Faure2020," in text
    assert json.loads(keys.read_text())["doi:10.1093/mnras/staa242"] == "Faure2020"
    assert main(["bib", str(LAMDA / "catom.dat"), "--no-database-refs"]) == 0
    captured = capsys.readouterr()
    assert captured.out.startswith("% References of 1 LAMDA file,") and "@" not in captured.out
    assert "6 references without DOI, bibcode or arXiv ID" in captured.err
    bad = tmp_path / "bad.dat"
    bad.write_text("!MOLECULE\nX\n")
    assert main(["bib", str(bad), str(LAMDA / "catom.dat")]) == 2


# ------------------------------------------------------------------- online, real (pytest --online)
def _need_token() -> str:
    token = scix_token()
    if not token:
        pytest.skip("needs a SciX API token (SCIX_API_TOKEN)")
    return token


@pytest.mark.online
def test_database_refs_agree_with_scix():
    """The shipped entries of the LAMDA and EMAA papers: DOI, bibcode, year and volume as in SciX."""
    from pylamda.bib import SCIX_ABS, HttpClient

    client = HttpClient(_need_token())
    for key, e in database_entries().items():
        adsurl = entry_field(e, "adsurl")
        if not adsurl:  # the EMAA dataset DOI is not in SciX
            assert client.doi_bibtex(entry_field(e, "doi")) is not None
            continue
        bibcode = adsurl[len(SCIX_ABS) :]
        docs = client.scix_search(f'doi:"{entry_field(e, "doi")}"', "bibcode,year,volume")
        assert [d["bibcode"] for d in docs] == [bibcode], key
        assert docs[0]["year"] == entry_field(e, "year") and docs[0]["volume"] == entry_field(e, "volume"), key


@pytest.mark.online
def test_online_export_emaa(tmp_path):
    _need_token()
    result = bibtex([EMAA_CO], online=True, cache=tmp_path / "cache.json")
    sources = {e.key: e.source for e in result.entries}
    assert sources["Muller2005"] == "scix" and sources["Dagdigian2022"] == "scix"
    assert sources["EMAA2023"] == "doi.org"  # the dataset DOI is not in SciX
    assert result.warnings == [] and os.environ.get("SCIX_API_TOKEN", "x") not in (tmp_path / "cache.json").read_text()


def test_exports():
    assert pylamda.bibtex is bibtex and pylamda.collect_references is collect_references
