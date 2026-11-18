"""The reference catalogue: reviewed free-text references of LAMDA files (design document, section 10.1).

Many files give a reference only as free text ("Yang et al. (2010)"). The catalogue
``pylamda/data/references.json`` (with ``references.bib``) maps such texts to reviewed
SciX records. It is built by the maintainer tool ``tools/catalogue.py``; every match was
confirmed by a person.

A reference text is split into the parts that cite one paper each
(:func:`split_segments`) and normalised (:func:`normalise`); each part is looked up
(:meth:`Catalogue.find`):

1. at the same place: file name, molecule, role (levels/collisions) and partner;
2. else the same text for the same molecule in another file, if it names one paper;
3. else the same text anywhere in the catalogue, if it names one paper.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from functools import lru_cache
from importlib.resources import files

# journal names as written in LAMDA/ProDiMo files -> ADS bibstem
JOURNALS = {
    "mnras": "MNRAS",
    "m.n.r.a.s.": "MNRAS",
    "mon. not. r. astron. soc.": "MNRAS",
    "apj": "ApJ",
    "ap. j.": "ApJ",
    "astrophys. j.": "ApJ",
    "apjs": "ApJS",
    "apjl": "ApJL",
    "a&a": "A&A",
    "astron. astrophys.": "A&A",
    "j. chem. phys.": "JChPh",
    "jchph": "JChPh",
    "j. phys. b": "JPhB",
    "j. phys b": "JPhB",
    "j.phys.b": "JPhB",
    "j phys b": "JPhB",
    "jphysb": "JPhB",
    "jpb": "JPhB",
    "pccp": "PCCP",
    "phys. chem. chem. phys.": "PCCP",
    "phys. rev. a": "PhRvA",
    "atoms": "Atoms",
    "jqsrt": "JQSRT",
    "chem. phys. lett.": "CPL",
    "mol. phys.": "MolPh",
}
_JOURNAL = re.compile("|".join(re.escape(j) for j in sorted(JOURNALS, key=len, reverse=True)), re.I)
YEAR = re.compile(r"\b(19\d\d|20\d\d)[a-z]?\b")
_SURNAME = re.compile(
    r"\b(?:(?:[Vv]an [Dd]e[rn]?|[Vv]on [Dd]er|[Dd]e [Ll]a|[Vv]an|[Vv]on|[Dd]e|[Dd]el|[Dd]i|[Dd]u|[Ll]e|[Ll]a)\s+)?"
    r"[A-Z][a-z]+(?:[-'][A-Z]?[a-z]+)*\b"
)
# capitalised words that are not surnames
_NOT_NAMES = {
    "from", "rates", "rate", "collision", "collisions", "coll", "downward", "data", "scaled", "using", "the",
    "rotational", "vibrational", "excitation", "quantum", "mechanical", "study", "transitions", "induced", "atoms",
    "molecules", "by", "of", "and", "with", "in", "astron", "astrophys", "phys", "chem", "rev", "lett", "mon",
    "not", "soc", "astrophysical", "journal", "vol", "priv", "comm", "private", "communication", "see", "also",
    "table", "extrapolation", "computed", "calculations", "ref", "note", "notes", "hyperfine", "levels",
    "spectroscopy", "corrected", "transmitted", "selected", "points", "same", "temperatures", "strictly", "valid",
}  # fmt: skip


def normalise(text: str) -> str:
    """Lower case; punctuation and blanks as one blank (the lookup key of the catalogue)."""
    return " ".join(re.sub(r"[^\w&+]+", " ", text.lower()).split())


def split_segments(text: str) -> list[str]:
    """Parts of a reference text that cite one paper each (one year per part).

    A part without an author-like name (``"oct 2019 at same temperatures"``) is merged
    into the part before it.
    """
    years = list(YEAR.finditer(text))
    if len(years) <= 1:
        return [text]
    out, start = [], 0
    for i, m in enumerate(years):
        end = m.end()
        if i + 1 < len(years):  # extend over ")" "," journal/volume/page up to the next author
            nxt = years[i + 1].start()
            cut = _cut(text[end:nxt])
            end = end + cut if cut is not None else nxt
        else:
            end = len(text)
        out.append(text[start:end])
        start = end
    merged: list[str] = []
    for part in out:
        head = part[: (YEAR.search(part) or re.search(r"$", part)).start()]
        if merged and not _author_names(_JOURNAL.sub(" ", head)):
            merged[-1] += part
        else:
            merged.append(part)
    return [s.strip(" ,;&+") for s in merged if s.strip(" ,;&+")]


def _cut(between: str) -> int | None:
    """Where the text between two cited years splits: ``;`` or ``+``, else a comma before a
    capitalised word that is not a journal, else ``&`` or "and"."""
    for pattern in (r"[;+]", r",\s(?=[A-Z][a-z])", r"&|\band\b"):
        for m in re.finditer(pattern, between):
            if pattern.startswith(",") and _JOURNAL.match(between, m.end()):
                continue
            return m.start()
    return None


# what follows an author name: et al., &, and, a comma, an initial, the year
_AFTER_AUTHOR = re.compile(r"\s*(?:et\s+al|&|and\b|,|[A-Z]\.|\(?(?:19|20)\d\d|$)")


def _author_names(head: str) -> list[str]:
    """Surnames in the text before the year; author-like ones (followed by et al., &, ...) first."""
    j = _JOURNAL.search(head)
    names_part = head[: j.start()] if j else head
    found = [
        m
        for m in _SURNAME.finditer(names_part)
        if not set(re.split(r"[-' ]", m.group().lower())) & _NOT_NAMES and len(m.group()) >= 3
    ]
    authorish = [m.group() for m in found if _AFTER_AUTHOR.match(names_part, m.end())]
    if not authorish:
        return [m.group() for m in found]
    first = next(i for i, m in enumerate(found) if m.group() == authorish[0])
    # "Storzer Hollenbach 2000": a name directly before the author-like one is the first author
    while first > 0 and not names_part[found[first - 1].end() : found[first].start()].strip():
        first -= 1
    return [m.group() for m in found[first:]]


def parse(segment: str) -> dict:
    """First author, other authors, year, journal (bibstem), volume and page of one citation.

    The first author is the first name followed by "et al.", "&", "and", a comma, an
    initial or the year (so that "Coulomb-Born + recoupling (Faure 2014" gives Faure).
    The journal is the first one after the year, else the last one before it (a title
    can hold a journal word, "4He Atoms").
    """
    year_m = YEAR.search(segment)
    head = segment[: year_m.start()] if year_m else segment
    names = _author_names(head)
    out = {"first_author": names[0] if names else None, "authors": names[1:4]}
    out["year"] = int(year_m.group(1)) if year_m else None
    journals = list(_JOURNAL.finditer(segment))
    after = [j for j in journals if year_m and j.start() >= year_m.end()]
    jm = after[0] if after else (journals[-1] if journals else None)
    out["bibstem"] = JOURNALS[jm.group().lower()] if jm else None
    vp = None
    if jm:
        vp = re.search(r"(?:vol\.?\s*)?(\d{1,4})\s*[:,]?\s*(?:p\.?\s*)?([A-Z]?\d{1,6})", segment[jm.end() :])
    out["volume"], out["page"] = (vp.group(1), vp.group(2).lstrip("0") or "0") if vp else (None, None)
    return out


# ------------------------------------------------------------------- the catalogue
@dataclass(frozen=True)
class CatalogueHit:
    """The result of a lookup.

    Attributes
    ----------
    key : str or None
        The catalogue key, or None if the reference was reviewed and is not in SciX.
    bibcode, doi : str or None
        Identifiers of the paper.
    match : str
        How it was found: ``"place"``, ``"molecule"`` or ``"text"``.
    """

    key: str | None
    bibcode: str | None
    doi: str | None
    match: str


@dataclass
class Catalogue:
    """The shipped reference catalogue (:func:`load_catalogue`).

    Attributes
    ----------
    entries : dict
        Key -> ``{"bibcode", "doi", "reviewed"}``.
    rows : list of dict
        The lookup rows: database, molecule, file, role, partner, text (normalised), key,
        status.
    bibtex : dict of str to str
        Key -> BibTeX entry (from SciX).
    """

    entries: dict[str, dict] = field(default_factory=dict)
    rows: list[dict] = field(default_factory=list)
    bibtex: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self._place: dict[tuple, dict] = {}
        self._molecule: dict[tuple, set] = {}
        self._text: dict[str, set] = {}
        for r in self.rows:
            self._place[(r["file"], r["molecule"], r["role"], r["partner"], r["text"])] = r
            self._molecule.setdefault((r["molecule"], r["text"]), set()).add(r["key"])
            self._text.setdefault(r["text"], set()).add(r["key"])

    def _hit(self, key: str | None, match: str) -> CatalogueHit:
        e = self.entries.get(key, {}) if key else {}
        return CatalogueHit(key, e.get("bibcode"), e.get("doi"), match)

    def find(
        self, text: str, *, molecule: str, file: str | None = None, role: str = "", partner: str | None = None
    ) -> CatalogueHit | None:
        """Look up one cited paper (a part of a reference text).

        Parameters
        ----------
        text : str
            The text of one citation (normalised here with :func:`normalise`).
        molecule : str
            The molecule of the file.
        file : str or None, optional
            The file name (without directories).
        role : str, optional
            ``"levels"`` or ``"collisions"``.
        partner : str or None, optional
            The partner by ID (``"p-H2"``).

        Returns
        -------
        CatalogueHit or None
            The hit, or None if the text is not in the catalogue (or names more than one
            paper there).
        """
        norm = normalise(text)
        row = self._place.get((file, molecule, role, partner, norm)) if file else None
        if row is not None:
            return self._hit(row["key"], "place")
        keys = self._molecule.get((molecule, norm), set())
        if len(keys) == 1:
            return self._hit(next(iter(keys)), "molecule")
        keys = self._text.get(norm, set())
        if len(keys) == 1 and None not in keys:
            return self._hit(next(iter(keys)), "text")
        return None

    def find_all(
        self, text: str, *, molecule: str, file: str | None = None, role: str = "", partner: str | None = None
    ) -> list[tuple[str, CatalogueHit | None]]:
        """The parts of a reference text (:func:`split_segments`) with their hits."""
        return [
            (seg, self.find(seg, molecule=molecule, file=file, role=role, partner=partner))
            for seg in split_segments(text)
        ]


def _entries_by_key(text: str) -> dict[str, str]:
    from .bib import entry_key, split_entries

    return {entry_key(e): e for e in split_entries(text)}


@lru_cache(maxsize=1)
def load_catalogue() -> Catalogue:
    """The catalogue shipped with pylamda (empty if the files are missing).

    Returns
    -------
    Catalogue
        Entries, lookup rows and BibTeX of ``pylamda/data/references.json`` and ``.bib``.
    """
    data = files("pylamda").joinpath("data")
    try:
        cat = json.loads(data.joinpath("references.json").read_text(encoding="utf-8"))
        bib = data.joinpath("references.bib").read_text(encoding="utf-8")
    except (FileNotFoundError, OSError):
        return Catalogue()
    return Catalogue(cat.get("entries", {}), cat.get("lookup", []), _entries_by_key(bib))
