"""Search SciX for free-text references (design document, section 10.1).

Used by the maintainer tool ``tools/catalogue.py`` (to build the reference catalogue)
and by ``pylamda bib --search`` (candidates for references not in the catalogue). A
free-text reference is split into the papers it cites; equal texts of the same
molecule form one :class:`Group` with one decision. For each group:

* hints are collected: the file's comments, the LAMDA molecule page (reference links to
  ADS), the EMAA species page;
* SciX is searched by first author and year, by species, by volume/page/year and by
  compound surname;
* the candidates are scored (volume, page, journal, first author, co-authors, species
  and partner in the title, journal article, named by a hint);
* a ``proposal`` is pre-filled where the match is clear.

Groups are kept in a review file (JSON). A person sets ``decision`` for each group:
``"ok"`` (the proposal), a bibcode, ``"none"`` (not in SciX) or ``"skip"`` (not a
reference). Nothing found here is used without such a decision.
"""

from __future__ import annotations

import html
import json
import re
import sys
import unicodedata
from dataclasses import asdict, dataclass, field
from pathlib import Path

from .bib import HttpClient, collect_references
from .catalogue import YEAR as _YEAR
from .catalogue import Catalogue, normalise, parse, split_segments
from .model import LamdaData, partner_name
from .reader import read
from .references import find_bibcodes, find_dois

LAMDA_URL = "https://home.strw.leidenuniv.nl/~moldata/"
EMAA_URL = "https://emaa.osug.fr/"

ELEMENTS = {
    "H": "hydrogen", "He": "helium", "C": "carbon", "N": "nitrogen", "O": "oxygen", "F": "fluorine",
    "Ne": "neon", "Na": "sodium", "Mg": "magnesium", "Al": "aluminium", "Si": "silicon", "P": "phosphorus",
    "S": "sulfur", "Cl": "chlorine", "Ar": "argon", "K": "potassium", "Ca": "calcium", "Fe": "iron",
}  # fmt: skip
_ROMAN = ["i", "ii", "iii", "iv", "v", "vi"]
# words in a title that name a collision partner
_PARTNER_WORDS = {
    "H2": ("h2", "molecular hydrogen"),
    "p-H2": ("h2", "para-hydrogen", "molecular hydrogen"),
    "o-H2": ("h2", "ortho-hydrogen", "molecular hydrogen"),
    "e": ("electron",),
    "H": ("atomic hydrogen", "hydrogen atom", " h atom", " h impact", " h(", " by h "),
    "He": ("helium", " he ", "4he", "3he", "he atom", "he impact", " he)"),
    "H+": ("proton", "h+"),
}
_SCORE_FOR_PROPOSAL = ("volume", "page")
# a "best candidate" proposal: species in the title, a minimum score and a clear lead
BEST_MIN_SCORE = 0.3
BEST_MARGIN = 0.15
BEST_REASON = "best candidate: species in the title, clear lead (check it)"
# proposals and skips entered by a first pass of an assistant: confirmed only on request
FIRST_PASS = "first pass (Claude)"


# ------------------------------------------------------------------- data
@dataclass
class Use:
    """Where a free-text reference occurs (one lookup row of the catalogue)."""

    id: str
    database: str
    file: str
    role: str
    partner: str | None  # by ID
    named: str | None  # the partner name in the text, if it differs from the ID
    text: str  # the whole reference text


@dataclass
class Group:
    """One cited paper: the same text (normalised) for the same molecule, in one or more files."""

    key: str
    molecule: str
    segment: str
    normalised: str
    uses: list[Use] = field(default_factory=list)
    first_author: str | None = None
    authors: list[str] = field(default_factory=list)
    year: int | None = None
    bibstem: str | None = None
    volume: str | None = None
    page: str | None = None
    hints: list[dict] = field(default_factory=list)
    query: str = ""
    candidates: list[dict] = field(default_factory=list)
    proposal: str | None = None  # bibcode pre-filled by the tool
    proposal_reason: str = ""
    decision: str | None = None  # "ok", a bibcode, "none" or "skip": set by the reviewer
    comment: str = ""

    @classmethod
    def from_dict(cls, d: dict) -> Group:
        d = dict(d)
        d["uses"] = [Use(**u) for u in d.get("uses", [])]
        return cls(**d)

    @property
    def bibcode(self) -> str | None:
        """The decided bibcode (None if not decided, ``none`` or ``skip``)."""
        if self.decision == "ok":
            return self.proposal
        return self.decision if self.decision and self.decision not in ("none", "skip") else None


# ------------------------------------------------------------------- parsing the texts
def ascii_fold(text: str) -> str:
    """Text without accents (``Schröder`` -> ``Schroder``), lower case."""
    return unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode().lower()


def name_forms(text: str) -> set[str]:
    """A text in lower case without accents, in both German spellings (``ö`` -> ``o`` and ``oe``)."""
    umlauts = str.maketrans({"ä": "ae", "ö": "oe", "ü": "ue", "Ä": "Ae", "Ö": "Oe", "Ü": "Ue", "ß": "ss"})
    return {ascii_fold(text), ascii_fold(text.translate(umlauts))}


# ------------------------------------------------------------------- titles
def normalise_title(title: str) -> str:
    """A SciX title in plain lower case: ``C<SUP>+</SUP>`` and ``C/+/`` -> ``c+``, ``H<SUB>2</SUB>`` -> ``h2``."""
    t = re.sub(r"<sup>\s*([+-]+)\s*</sup>|/([+-])/|\$\^\{?([+-])\}?\$|\^\{?([+-])\}?", lambda m: next(
        g for g in m.groups() if g), title, flags=re.I)  # fmt: skip
    t = re.sub(r"<[^>]+>", "", t)
    return " " + " ".join(html.unescape(t).lower().split()) + " "


def species_forms(molecule: str) -> tuple[list[str], list[str]]:
    """Regular expressions for a species in a title, and for the species with another charge.

    For an atom also the element name and the spectroscopic notation (``C II`` = ``C+``).
    """
    name = re.sub(r"^(?:ortho|para|meta|[opmaecl])-", "", molecule, flags=re.I)
    charge = len(name) - len(name.rstrip("+"))
    base = name.rstrip("+-")
    b = re.escape(base.lower())
    # the species: the formula with its charge, not followed by more charge or letters/digits
    own = [rf"(?<![a-z0-9]){b}{re.escape('+' * charge)}(?![a-z0-9+])"]
    other = [rf"(?<![a-z0-9]){b}\+"] if charge == 0 else [rf"(?<![a-z0-9]){b}(?![a-z0-9+])(?! i)"]
    if base in ELEMENTS:  # an atom
        el = ELEMENTS[base]
        own.append(rf"(?<![a-z0-9]){b} {_ROMAN[charge]}(?![a-z])")
        if charge == 0:
            own.append(rf"\b{el}(?! monoxide| dioxide| sulfide| ion| ii)\b")
            other += [rf"(?<![a-z0-9]){b} ii(?![a-z])", rf"\bionized {el}\b", rf"\b{el} ions?\b"]
        else:
            own += [rf"\bionized {el}\b", rf"\b{el} ions?\b"]
    return own, other


def title_match(molecule: str, title: str) -> int:
    """+1 if the title names the species, -1 if it names it only with another charge, else 0."""
    t = normalise_title(title)
    own, other = species_forms(molecule)
    if any(re.search(p, t) for p in own):
        return 1
    if any(re.search(p, t) for p in other):
        return -1
    return 0


# ------------------------------------------------------------------- hints
def page_names(molecule: str) -> list[str]:
    """Guesses for the LAMDA/EMAA page name of a molecule (``e-CH3OH`` -> ``CH3OH``, ``C18O`` -> ``CO``)."""
    names = [molecule]
    despin = re.sub(r"^(?:ortho|para|meta|[opmae])-", "", molecule, flags=re.I)
    noiso = re.sub(r"(?<![A-Za-z])\d{1,3}(?=[A-Z])|(?<=[A-Za-z])1[3-8](?=[A-Z])", "", despin)
    nostruct = re.sub(r"^[cl]-", "", noiso)  # c-C3H2, l-C3H
    for n in (despin, noiso, nostruct):
        if n not in names:
            names.append(n)
    return names


def _get(url: str, cache: Path | None) -> str | None:
    """A web page, from the cache if there; None if not found."""
    path = cache / re.sub(r"[^\w.+-]", "_", url.split("//", 1)[-1]) if cache else None
    if path is not None and path.is_file():
        text = path.read_text(encoding="utf-8")
        return text or None
    import requests

    try:
        r = requests.get(url, timeout=30, headers={"User-Agent": "pylamda-catalogue"})
    except requests.RequestException as exc:
        print(f"catalogue.py: {url}: {type(exc).__name__}", file=sys.stderr)
        return None
    text = r.text if r.status_code == 200 else ""
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")  # "" = not found
    return text or None


def page_hints(page_html: str, source: str) -> list[dict]:
    """Hints from a database page: the text of each link to ADS, SciX, a DOI or an EMAA DOI page."""
    out = []
    for m in re.finditer(r'<a[^>]+href=["\']([^"\']+)["\'][^>]*>(.*?)</a>', page_html, re.S | re.I):
        href, label = m.group(1), " ".join(html.unescape(re.sub(r"<[^>]+>", " ", m.group(2))).split())
        # the citation text: the link text, or the text just before an EMAA "DOI"/"BibTeX" link
        if not re.search(r"\d{4}", label):
            before = re.sub(r"<[^>]+>", " ", page_html[max(0, m.start() - 400) : m.start()])
            label = " ".join(html.unescape(before).split())[-200:]
        bib = find_bibcodes(href)
        doi = find_dois(href)
        emaa = re.search(r"details-doi/\d+", href)
        if bib or doi or emaa:
            out.append({
                "source": source, "text": label, "bibcode": bib[0] if bib else None,
                "doi": doi[0] if doi else None, "url": href if emaa else None,
            })  # fmt: skip
    return out


def file_hints(data: LamdaData, source: str) -> list[dict]:
    """Hints from the comment lines and notes of a file (lines with an identifier)."""
    out = []
    texts = list(data.notes)
    if data.source_map is not None:
        texts += [line for line in data.source_map.raw if line.lstrip().startswith("!")]
    for t in texts:
        bib, doi = find_bibcodes(t), find_dois(t)
        if bib or doi:
            out.append({"source": source, "text": t.strip("! ").strip(), "bibcode": bib[0] if bib else None,
                        "doi": doi[0] if doi else None, "url": None})  # fmt: skip
    return out


def hint_matches(group: Group, hint: dict) -> bool:
    """Whether a hint is about the paper of a group: same first author, and the same year."""
    if not group.first_author or not group.year:
        return False
    texts = name_forms(hint["text"])
    names = name_forms(group.first_author.split()[-1])
    if not any(n in t for n in names for t in texts):
        return False
    years = {int(y) for y in _YEAR.findall(hint["text"])}
    if hint.get("bibcode"):
        years.add(int(hint["bibcode"][:4]))
    return group.year in years


# ------------------------------------------------------------------- collecting
def collect(
    sources: list, database: str, catalogue: Catalogue | None = None
) -> tuple[list[Group], dict[str, LamdaData]]:
    """The free-text references of files, grouped by molecule and normalised text.

    Parameters
    ----------
    sources : list of str, os.PathLike or LamdaData
        Files, or data read from files.
    database : str
        Label of the database (``"lamda"``, ``"prodimo"``, ...) for the review file.
    catalogue : Catalogue or None, optional
        Leave out the papers this catalogue knows (``pylamda bib --search``).

    Returns
    -------
    groups : list of Group
        One group per molecule and normalised citation text.
    datas : dict
        Source label -> data (for the file hints).
    """
    groups: dict[str, Group] = {}
    datas = {}
    for source in sources:
        data = source if isinstance(source, LamdaData) else read(source)
        path = data.source_map.path if data.source_map else str(source)
        name = Path(path).name
        datas[path] = data
        named = {c.label: c.text_partner for c in data.collisions}
        for ref in collect_references(data):
            if ref.identifiers or ref.role not in ("levels", "collisions"):
                continue
            segments = split_segments(ref.text)
            for k, seg in enumerate(segments, 1):
                if catalogue is not None:
                    hit = catalogue.find(seg, molecule=data.molecule, file=name, role=ref.role, partner=ref.partner)
                    if hit is not None:
                        continue
                written = named.get(ref.partner) if ref.partner else None
                written = (partner_name(written) or written) if written else None
                use_id = f"{database}/{name}/{ref.role}" + (f"/{ref.partner}" if ref.partner else "")
                use = Use(
                    id=use_id + (f"#{k}" if len(segments) > 1 else ""),
                    database=database,
                    file=name,
                    role=ref.role,
                    partner=ref.partner,
                    named=written if written != ref.partner else None,
                    text=ref.text,
                )
                key = f"{data.molecule}|{normalise(seg)}"
                g = groups.get(key)
                if g is None:
                    g = groups[key] = Group(key=key, molecule=data.molecule, segment=seg, normalised=normalise(seg))
                    for k2, v in parse(seg).items():
                        setattr(g, k2, v)
                if use.id not in (u.id for u in g.uses):
                    g.uses.append(use)
                g._path = path  # noqa: SLF001  (for the file hints)
    return list(groups.values()), datas


def add_hints(groups: list[Group], datas: dict, database: str, lamda_page: str | None, cache: Path | None) -> None:
    """Hints from the files' comments, the LAMDA molecule page and (EMAA files) the EMAA species page."""
    pages: dict[str, list[dict]] = {}

    def lamda(molecule: str) -> list[dict]:
        names = [lamda_page] if lamda_page else page_names(molecule)
        for n in names:
            if f"lamda:{n}" not in pages:
                text = _get(f"{LAMDA_URL}{n}.html", cache)
                pages[f"lamda:{n}"] = page_hints(text, f"LAMDA page {n}") if text else []
            if pages[f"lamda:{n}"]:
                return pages[f"lamda:{n}"]
        return []

    def emaa(molecule: str) -> list[dict]:
        for n in page_names(molecule):
            if f"emaa:{n}" not in pages:
                text = _get(f"{EMAA_URL}details/{n}", cache)
                pages[f"emaa:{n}"] = page_hints(text, f"EMAA page {n}") if text else []
            if pages[f"emaa:{n}"]:
                return pages[f"emaa:{n}"]
        return []

    for g in groups:
        data = datas[g._path]  # noqa: SLF001
        hints = file_hints(data, "comment in the file") + lamda(g.molecule)
        if database == "emaa":
            hints += emaa(g.molecule)
        g.hints = [h for h in hints if hint_matches(g, h)]


# ------------------------------------------------------------------- SciX
def is_first_author(name: str, scix_author: str) -> bool:
    """Whether a name from a text is the surname of a SciX author (``"Yang, Benhui"``).

    Equal, or the SciX surname starts or ends with it (``Hernandez`` = ``Hernández Vera``,
    ``Tak`` = ``van der Tak``); accents and umlauts in both spellings.
    """
    surname = scix_author.split(",")[0]
    return any(
        s == n or s.startswith(n + " ") or s.endswith(" " + n) for n in name_forms(name) for s in name_forms(surname)
    )


def score(group: Group, doc: dict) -> tuple[float, list[str]]:
    """Score of a SciX document for a group, and why."""
    s, why = 0.0, []
    if group.volume and str(doc.get("volume", "")) == group.volume:
        s, why = s + 0.3, [*why, "volume"]
    pages = [str(p).lstrip("0") for p in (doc.get("page") or [])]
    if group.page and group.page.lstrip("0") in pages:
        s, why = s + 0.3, [*why, "page"]
    if group.bibstem and group.bibstem in doc.get("bibcode", ""):
        s, why = s + 0.1, [*why, "journal"]
    if group.first_author and is_first_author(group.first_author, (doc.get("author") or [""])[0]):
        s, why = s + 0.1, [*why, "first author"]
    authors = ascii_fold(" ".join(doc.get("author", [])))
    hits = [a for a in group.authors if ascii_fold(a) in authors]
    if hits:
        s, why = s + min(0.1 * len(hits), 0.2), [*why, f"co-authors {', '.join(hits)}"]
    title = " ".join(doc.get("title", []))
    # the species or, for an isotopologue or spin isomer, its parent (13CO -> CO)
    matches = [title_match(name, title) for name in page_names(group.molecule)]
    m = 1 if 1 in matches else (-1 if -1 in matches else 0)
    if m > 0:
        s, why = s + 0.15, [*why, f"{group.molecule} in title"]
    elif m < 0:
        s, why = s - 0.15, [*why, f"title about another charge state of {group.molecule}"]
    partners = {u.named or u.partner for u in group.uses} | {u.partner for u in group.uses}
    t = normalise_title(title)
    if any(w in t for p in partners if p for w in _PARTNER_WORDS.get(p, ())):
        s, why = s + 0.1, [*why, "partner in title"]
    if doc.get("doctype") == "article":
        s, why = s + 0.05, [*why, "journal article"]
    if any(h.get("bibcode") == doc.get("bibcode") for h in group.hints):
        s, why = s + 0.5, [*why, "named by a hint"]
    return round(s, 2), why


def _candidate(group: Group, doc: dict) -> dict:
    sc, why = score(group, doc)
    authors = doc.get("author") or []
    return {
        "bibcode": doc["bibcode"],
        "score": sc,
        "why": ", ".join(why),
        "title": (doc.get("title") or [""])[0][:120],
        "authors": "; ".join(authors[:4]) + (" et al." if len(authors) > 4 else ""),
        "pub": f"{doc.get('pub', '')} {doc.get('volume', '')}, {(doc.get('page') or [''])[0]}",
        "doi": (doc.get("doi") or [None])[0],
    }


_FIELDS = "bibcode,title,author,year,volume,page,doi,pub,doctype"


def search(groups: list[Group], client: HttpClient, rows: int = 25, keep: int = 5) -> None:
    """SciX candidates for each group (first author and year), best first, and the proposal."""
    for g in groups:
        # hints with a DOI or an EMAA DOI page: the bibcode from SciX
        for h in g.hints:
            if not h.get("bibcode") and h.get("doi"):
                docs = client.scix_search(f'doi:"{h["doi"]}"', "bibcode")
                h["bibcode"] = docs[0]["bibcode"] if docs else None
        if not g.first_author or not g.year:
            g.comment = "no author or year found in the text: search by hand"
            continue
        g.query = f'author:"^{g.first_author}" year:{g.year}'
        docs = client.scix_search(g.query, _FIELDS)[:rows]
        known = {d["bibcode"] for d in docs}
        # common names: also the papers with the species in title or abstract (Green 1986 for CH3CN)
        terms = page_names(g.molecule)[-1:]
        base = re.sub(r"[+-]+$", "", terms[0])
        if base in ELEMENTS:
            terms.append(ELEMENTS[base])
        species = " OR ".join(f'abs:"{x}" OR title:"{x}"' for x in dict.fromkeys(terms))
        extra = [f"{g.query} ({species})"]
        if g.volume and g.page:  # by volume, page and year: also with a misspelt author name
            extra.append(f"volume:{g.volume} page:{g.page} year:{g.year}")
        for query in extra:
            for d in client.scix_search(query, _FIELDS)[:rows]:
                if d["bibcode"] not in known:
                    docs.append(d)
                    known.add(d["bibcode"])
        # a compound surname ("Hernandez Vera"): the first two names together, if nothing fits well
        best = max((score(g, d)[0] for d in docs), default=0.0)
        if best < BEST_MIN_SCORE and g.authors and f"{g.first_author} {g.authors[0]}" in g.segment:
            for d in client.scix_search(f'author:"^{g.first_author} {g.authors[0]}" year:{g.year}', _FIELDS)[:rows]:
                if d["bibcode"] not in known:
                    docs.append(d)
                    known.add(d["bibcode"])
        for h in g.hints:  # a paper named by a hint, not found by the query (other first-author spelling)
            if h.get("bibcode") and h["bibcode"] not in known:
                docs += client.scix_search(f'bibcode:"{h["bibcode"]}"', _FIELDS)
                known.add(h["bibcode"])
        cands = sorted((_candidate(g, d) for d in docs), key=lambda c: -c["score"])
        g.candidates = cands[:keep]
        if not docs:
            g.comment = "SciX found nothing for first author and year"
        propose(g)


def propose(group: Group) -> None:
    """Pre-fill the proposal where the match is clear.

    Clear: a hint (LAMDA/EMAA page, the file's comments) names the paper for the same
    first author and year; or volume and page of the text match the top candidate, and
    no other candidate matches them.
    """
    group.proposal, group.proposal_reason = None, ""
    named = [c for c in group.candidates if "named by a hint" in c["why"]]
    if len(named) == 1:
        sources = sorted({h["source"] for h in group.hints if h.get("bibcode") == named[0]["bibcode"]})
        group.proposal, group.proposal_reason = named[0]["bibcode"], f"named by: {', '.join(sources)}"
        return
    vp = [c for c in group.candidates if all(w in c["why"].split(", ") for w in _SCORE_FOR_PROPOSAL)]
    if len(vp) > 1:  # the same volume and page in several journals: the right journal or first author
        vp = [c for c in vp if {"journal", "first author"} & set(c["why"].split(", "))]
    if len(vp) == 1:
        group.proposal, group.proposal_reason = vp[0]["bibcode"], "volume and page match"
        return
    # weaker: the best candidate names the species, scores >= 0.3 and leads clearly
    c = group.candidates
    if c and " in title" in c[0]["why"] and c[0]["score"] >= BEST_MIN_SCORE:
        second = c[1]["score"] if len(c) > 1 else 0.0
        if c[0]["score"] - second >= BEST_MARGIN:
            group.proposal, group.proposal_reason = c[0]["bibcode"], BEST_REASON


# ------------------------------------------------------------------- review file
def load(path: Path) -> list[Group]:
    return [Group.from_dict(x) for x in json.loads(path.read_text(encoding="utf-8"))] if path.exists() else []


def save(groups: list[Group], path: Path) -> None:
    out = []
    for g in groups:
        d = asdict(g)
        d.pop("_path", None)
        out.append(d)
    path.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


def show(groups: list[Group], only: str | None = None) -> None:
    """The review file in a readable form (``only``: ``"open"`` or ``"proposals"``)."""
    for g in groups:
        if only == "open" and g.decision:
            continue
        if only == "proposals" and not (g.proposal and not g.decision):
            continue
        if only == "best" and not (g.proposal_reason == BEST_REASON and not g.decision):
            continue
        if only == "no_proposal" and (g.proposal or g.decision):
            continue
        if only == "first_pass" and not (g.proposal_reason.startswith(FIRST_PASS) or g.comment.startswith(FIRST_PASS)):
            continue
        print(f"\n== {g.key}")
        print(f"   text:      {g.segment}")
        for u in g.uses:
            named = f" (text says {u.named})" if u.named else ""
            print(f"   used in:   {u.file}: {u.role} {u.partner or ''}{named}")
        print(f"   parsed:    {g.first_author} {g.year}  {g.bibstem or ''} {g.volume or ''} {g.page or ''}"
              f"  co-authors: {', '.join(g.authors) or '-'}")  # fmt: skip
        for h in g.hints:
            print(f"   hint:      {h['source']}: {h['text'][:80]} -> {h.get('bibcode') or h.get('doi') or ''}")
        if g.comment:
            print(f"   note:      {g.comment}")
        for k, c in enumerate(g.candidates, 1):
            mark = "*" if c["bibcode"] == g.proposal else " "
            print(f"  {mark}{k}. {c['score']:.2f}  {c['bibcode']}  {c['title'][:70]}")
            print(f"            {c['authors'][:60]}; {c['pub']}  [{c['why']}]")
        if g.proposal:
            print(f"   proposal:  {g.proposal} ({g.proposal_reason})")
        print(f"   decision:  {g.decision}")


# ------------------------------------------------------------------- pylamda bib --search
def search_references(
    sources: list,
    *,
    client: HttpClient,
    hints: bool = True,
    cache_dir: Path | None = None,
    catalogue: Catalogue | None = None,
    database: str = "user",
) -> list[Group]:
    """SciX candidates for the free-text references of files that the catalogue does not know.

    The candidates are only proposals: nothing is added to a BibTeX file or a data file
    before a person confirms it in the review file.

    Parameters
    ----------
    sources : list of str, os.PathLike or LamdaData
        Files, or data read from files.
    client : HttpClient
        The HTTP client with the user's SciX token.
    hints : bool, optional
        Read the LAMDA molecule pages (and EMAA species pages) for hints. Default: True.
    cache_dir : Path or None, optional
        Cache of the database pages.
    catalogue : Catalogue or None, optional
        Leave out the papers this catalogue knows.
    database : str, optional
        Label of the files' database in the review file. Default: ``"user"``.

    Returns
    -------
    list of Group
        The groups with candidates and proposals.
    """
    groups, datas = collect(sources, database, catalogue)
    if hints and groups:
        add_hints(groups, datas, database, None, cache_dir)
    else:
        for g in groups:
            g.hints = [h for h in file_hints(datas[g._path], "comment in the file") if hint_matches(g, h)]  # noqa: SLF001
    search(groups, client)
    return groups


def merge_review(groups: list[Group], path: Path) -> list[Group]:
    """Add groups to a review file; decisions and comments already in it are kept.

    Parameters
    ----------
    groups : list of Group
        New or searched-again groups.
    path : Path
        The review file (created if missing).

    Returns
    -------
    list of Group
        All groups now in the file.
    """
    old = {g.key: g for g in load(path)}
    for g in groups:
        before = old.pop(g.key, None)
        if before is not None:
            g.decision, g.comment = before.decision, before.comment or g.comment
            g.uses += [u for u in before.uses if u.id not in {x.id for x in g.uses}]
    merged = [*old.values(), *groups]
    save(merged, path)
    return merged


def main(argv: list[str] | None = None) -> int:
    """``python -m pylamda.refsearch REVIEW [--open]``: show a review file."""
    import argparse

    parser = argparse.ArgumentParser(prog="python -m pylamda.refsearch", description="Show a review file.")
    parser.add_argument("review")
    parser.add_argument("--open", action="store_true", help="only groups without a decision")
    args = parser.parse_args(argv)
    show(load(Path(args.review)), "open" if args.open else None)
    return 0


if __name__ == "__main__":
    sys.exit(main())
