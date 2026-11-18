"""BibTeX export of the references of LAMDA files (design document, section 10).

1. **Collect** the references of each file: dataset DOI, level reference, partner
   references, and identifiers in comment lines, line 2 and the notes
   (:func:`collect_references`).
2. **Extract** the identifiers: DOIs, ADS/SciX bibcodes, arXiv IDs
   (:func:`pylamda.references.find_identifiers`).
3. **Resolve** them (``online=True``). SciX is the source of the bibliographic records:
   DOIs and arXiv IDs are searched in SciX, and all bibcodes are exported as BibTeX in
   one request. DOIs that SciX does not know (datasets) are resolved with doi.org. SciX
   needs the user's API token (:func:`scix_token`); without it, only DOIs are resolved
   (doi.org).
4. Identifiers that are not resolved (or offline) give a ``@misc`` entry with the
   identifier and a note saying what the reference is for.
5. The main papers of LAMDA and EMAA are added from ``pylamda/data/database_refs.bib``.

References given only as free text ("Flower (2001)") are not guessed: they are listed
in :attr:`BibResult.free_text`.
"""

from __future__ import annotations

import json
import os
import re
import stat
import unicodedata
import warnings
from collections.abc import Iterable
from dataclasses import dataclass, field, replace
from importlib.resources import files
from pathlib import Path

from .catalogue import load_catalogue
from .model import LamdaData
from .reader import read
from .references import find_identifiers

#: Base URL of the SciX API (the same API as ADS).
SCIX_API = "https://api.scixplorer.org/v1"
#: SciX abstract page of a bibcode: ``SCIX_ABS + bibcode``.
SCIX_ABS = "https://scixplorer.org/abs/"
#: Keys in ``database_refs.bib`` that every export gets: the LAMDA paper and the 2020 LAMDA
#: status review (decided 2026-11-18, design decision 10).
LAMDA_KEYS: tuple[str, ...] = ("Schoier2005", "vanderTak2020")
#: Keys added if a file comes from EMAA (EMAA's citation policy).
EMAA_KEYS: tuple[str, ...] = ("Faure2025", "EMAA")
#: The acknowledgement EMAA asks for (https://emaa.osug.fr/citation-policy).
EMAA_ACKNOWLEDGEMENT = (
    "This research has made use of spectroscopic and collisional data from the EMAA database "
    "(https://dx.doi.org/10.17178/EMAA). EMAA is supported by the Observatoire des Sciences de l'Univers "
    "de Grenoble (OSUG) (https://ror.org/03vte9x46)."
)
#: What LAMDA asks for besides the citation (https://home.strw.leidenuniv.nl/~moldata/).
LAMDA_REQUEST = (
    "LAMDA asks to mention the date when the database was accessed, and to cite the original papers "
    "of the spectroscopic and collisional data of each molecule."
)

_DATABASE_NOTES = {
    "Schoier2005": "LAMDA database",
    "vanderTak2020": "LAMDA database (2020 status review)",
    "Faure2025": "EMAA database",
    "EMAA": "EMAA database (dataset DOI)",
}
_ROLE_TEXT = {"dataset": "dataset", "levels": "levels and radiative data", "comments": "named in a comment"}
_SEARCH_CHUNK = 900  # characters of a SciX query (the limit is 1000)
_EXPORT_CHUNK = 200  # bibcodes per SciX export request


class BibError(RuntimeError):
    """Online resolution failed (e.g. SciX rejected the token, or requests is missing)."""


# ------------------------------------------------------------------- references in files
@dataclass
class Reference:
    """One reference found in a file.

    Attributes
    ----------
    path : str or None
        The file (None for data not read from a file).
    molecule : str
        The molecule of the file.
    role : str
        ``"dataset"`` (dataset DOI), ``"levels"`` (reference on the NLEV line),
        ``"collisions"`` (partner line) or ``"comments"`` (a comment line, line 2 or
        the notes; only with an identifier).
    partner : str or None
        The partner (by ID) for ``"collisions"``.
    text : str
        The reference text.
    identifiers : list of str
        ``"doi:..."``, ``"bibcode:..."`` and ``"arxiv:..."`` found in the text.
    """

    path: str | None
    molecule: str
    role: str
    partner: str | None
    text: str
    identifiers: list[str]

    @property
    def use(self) -> str:
        """str: What the reference is for (``"collision rates with p-H2"``)."""
        if self.role == "collisions":
            return f"collision rates with {self.partner}"
        return _ROLE_TEXT.get(self.role, self.role)


def collect_references(data: LamdaData) -> list[Reference]:
    """The references of one data set.

    Parameters
    ----------
    data : LamdaData
        The data (comment lines are read from ``data.source_map`` if present).

    Returns
    -------
    list of Reference
        Dataset DOI, level reference and partner references (with or without
        identifiers), then the comment lines, line 2 and notes that hold an identifier.
    """
    path = data.source_map.path if data.source_map else None
    mol = data.molecule
    refs = []
    if data.dataset_doi:
        refs.append(Reference(path, mol, "dataset", None, data.dataset_doi, [f"doi:{data.dataset_doi}"]))
    if data.level_reference.strip():
        text = data.level_reference.strip()
        refs.append(Reference(path, mol, "levels", None, text, find_identifiers(text)))
    for c in data.collisions:
        if c.reference.strip():
            text = c.reference.strip()
            refs.append(Reference(path, mol, "collisions", c.label, text, find_identifiers(text)))
    for text in _comment_texts(data):
        ids = find_identifiers(text)
        if ids:
            refs.append(Reference(path, mol, "comments", None, text.strip().lstrip("!").strip(), ids))
    return refs


def _comment_texts(data: LamdaData) -> list[str]:
    """Line 2 (title), the notes and the comment lines of the file (not the partner lines)."""
    texts = [data.title, *data.notes]
    smap = data.source_map
    if smap is not None:
        notes_from = smap.notes or len(smap.raw) + 1
        texts += [line for i, line in enumerate(smap.raw, 1) if i < notes_from and line.lstrip().startswith("!")]
    fix_note = re.compile(r"^\s*!\s*pylamda \S+ fix:")  # pylamda's own fix notes are not references
    return [t for t in texts if t and t.strip() and not fix_note.match(t)]


def is_emaa(data: LamdaData) -> bool:
    """Whether a data set comes from EMAA (dataset DOI ``10.17178/EMAA...`` or "EMAA" in its texts).

    Parameters
    ----------
    data : LamdaData
        The data.

    Returns
    -------
    bool
        True if there is an EMAA hint.
    """
    if (data.dataset_doi or "").upper().startswith("10.17178/EMAA"):
        return True
    texts = [data.level_reference, data.source, *(c.reference for c in data.collisions), *_comment_texts(data)]
    return any(re.search(r"\bEMAA\b", t or "") for t in texts)


# ------------------------------------------------------------------- BibTeX text
def split_entries(text: str) -> list[str]:
    """The entries of a BibTeX text (``@type{key, ...}``), without comments around them."""
    out = []
    for m in re.finditer(r"@\s*(\w+)\s*([{(])", text):
        if m.group(1).lower() in ("comment", "preamble", "string"):
            continue
        open_, close = m.group(2), "}" if m.group(2) == "{" else ")"
        depth, i = 0, m.end() - 1
        while i < len(text):
            ch = text[i]
            if ch == open_:
                depth += 1
            elif ch == close:
                depth -= 1
                if depth == 0:
                    out.append(text[m.start() : i + 1].strip())
                    break
            i += 1
    return out


def entry_key(entry: str) -> str:
    """The key of a BibTeX entry."""
    m = re.match(r"@\s*\w+\s*[{(]\s*([^,\s]*)\s*,", entry)
    return m.group(1) if m else ""


def _set_key(entry: str, key: str) -> str:
    return re.sub(r"^(@\s*\w+\s*[{(]\s*)[^,]*?(\s*,)", lambda m: m.group(1) + key + m.group(2), entry, count=1)


def _field_span(entry: str, name: str) -> tuple[int, int] | None:
    """(start, end) of the value of field ``name`` (with its braces or quotes)."""
    m = re.search(rf"[,{{\s]{re.escape(name)}\s*=\s*", entry, re.I)
    if not m:
        return None
    i = m.end()
    if i >= len(entry):
        return None
    if entry[i] in "{":
        depth = 0
        for j in range(i, len(entry)):
            depth += {"{": 1, "}": -1}.get(entry[j], 0)
            if depth == 0:
                return i, j + 1
        return None
    if entry[i] == '"':
        j = i + 1
        while j < len(entry) and (entry[j] != '"' or entry[j - 1] == "\\"):
            j += 1
        return i, j + 1
    m2 = re.match(r"[^,}\s]+", entry[i:])
    return (i, i + m2.end()) if m2 else None


def entry_field(entry: str, name: str) -> str | None:
    """The value of a field of a BibTeX entry, without the outer braces or quotes."""
    span = _field_span(entry, name)
    if span is None:
        return None
    value = entry[span[0] : span[1]]
    if value[:1] in '{"':
        value = value[1:-1]
    return value.strip()


def _set_field(entry: str, name: str, value: str) -> str:
    """Set (replace or add) a field, ``name = {value}``."""
    span = _field_span(entry, name)
    if span is not None:
        return entry[: span[0]] + "{" + value + "}" + entry[span[1] :]
    body = entry.rstrip()[:-1].rstrip().rstrip(",")
    return f"{body},\n  {name} = {{{value}}}\n}}"


def _ascii_letters(text: str) -> str:
    """Letters of a name without LaTeX accents, braces and blanks (``Sch{\\"o}ier`` -> ``Schoier``)."""
    text = re.sub(r"\\[\"'`^~=.uvHckbdrt]\s*\{?\s*([A-Za-z])\s*\}?", r"\1", text)
    letters = {"l": "l", "L": "L", "o": "o", "O": "O", "ss": "ss", "ae": "ae", "AE": "AE", "aa": "a", "AA": "A"}
    text = re.sub(r"\\(ss|ae|AE|aa|AA|[lLoOij])(?![A-Za-z])", lambda m: letters.get(m.group(1), m.group(1)), text)
    text = re.sub(r"\\[A-Za-z]+", "", text)
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()
    return re.sub(r"[^A-Za-z]", "", text)


def key_from_entry(entry: str) -> str | None:
    """``<first author surname><year>`` of a BibTeX entry (``Yang2010``), or None."""
    author, year = entry_field(entry, "author"), entry_field(entry, "year")
    if not author or not year or not re.fullmatch(r"\{?\d{4}\}?", year.strip()):
        return None
    first = re.split(r"\s+and\s+", author.strip())[0]
    depth, cut = 0, None
    for i, ch in enumerate(first):  # a comma outside braces separates "Surname, Given"
        depth += {"{": 1, "}": -1}.get(ch, 0)
        if ch == "," and depth == 0:
            cut = i
            break
    surname = first[:cut] if cut is not None else first.split()[-1] if first.split() else ""
    surname = _ascii_letters(surname)
    return f"{surname}{year.strip('{}')}" if surname else None


def key_from_text(text: str) -> str | None:
    """``<author><year>`` guessed from a free reference text (``"Yang et al. (2010)"`` -> ``Yang2010``)."""
    m = re.search(r"([A-Z][A-Za-z'\-]+).{0,60}?\b((?:19|20)\d\d)\b", text or "")
    return f"{_ascii_letters(m.group(1))}{m.group(2)}" if m else None


def identifier_context(text: str, identifier: str) -> str:
    """The part of a reference text that belongs to one identifier (for :func:`key_from_text`).

    The text before the identifier, back to the previous identifier or the last ``;``:
    in ``"Yang 2010 doi:A; Dagdigian 2022 doi:B"``, ``"Dagdigian 2022 "`` for ``doi:B``.

    Parameters
    ----------
    text : str
        The reference text.
    identifier : str
        One of its identifiers (``"doi:..."``, ...).

    Returns
    -------
    str
        The part of the text; the whole text if the identifier is not found in it.
    """
    lower = text.lower()
    value = identifier.split(":", 1)[1].lower()
    idx = lower.find(value)
    if idx < 0:
        return text
    start = lower.rfind(";", 0, idx) + 1
    for other in find_identifiers(text):
        v = other.split(":", 1)[1].lower()
        j = lower.find(v)
        if v != value and 0 <= j and j + len(v) <= idx:
            start = max(start, j + len(v))
    return text[start:idx]


def _escape(text: str) -> str:
    return re.sub(r"([&%#_$])", r"\\\1", text)


def _misc(key: str, identifier: str, note: str) -> str:
    """A ``@misc`` entry for an identifier that is not resolved."""
    kind, value = identifier.split(":", 1)
    if kind == "doi":
        fields = [f"doi = {{{value}}}"]
    elif kind == "bibcode":
        fields = [f"howpublished = {{\\url{{{SCIX_ABS}{value}}}}}", f"adsurl = {{{SCIX_ABS}{value}}}"]
    else:
        fields = [f"eprint = {{{value}}}", "archivePrefix = {arXiv}"]
    return "@misc{" + key + ",\n  " + ",\n  ".join([f"note = {{{_escape(note)}}}", *fields]) + ",\n}"


# ------------------------------------------------------------------- token, cache, HTTP
def scix_token(token: str | None = None) -> str | None:
    """The user's SciX API token.

    Looked up in this order:

    1. the ``token`` argument;
    2. the environment variable ``SCIX_API_TOKEN``;
    3. the file ``$XDG_CONFIG_HOME/pylamda/scix_token`` (default
       ``~/.config/pylamda/scix_token``), first line. A warning is given if other users
       can read the file.

    Parameters
    ----------
    token : str or None, optional
        A token given directly.

    Returns
    -------
    str or None
        The token, or None if there is none. pylamda never writes it to files, logs,
        the cache or messages.
    """
    if token:
        return token.strip()
    if os.environ.get("SCIX_API_TOKEN", "").strip():
        return os.environ["SCIX_API_TOKEN"].strip()
    path = token_file()
    if path.is_file():
        if path.stat().st_mode & (stat.S_IRWXG | stat.S_IRWXO):
            warnings.warn(f"{path} can be read by other users; run: chmod 600 {path}", stacklevel=2)
        lines = path.read_text(encoding="utf-8").splitlines()
        return lines[0].strip() if lines and lines[0].strip() else None
    return None


def token_file() -> Path:
    """The SciX token file, ``$XDG_CONFIG_HOME/pylamda/scix_token`` (default ``~/.config/...``)."""
    return Path(os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config") / "pylamda" / "scix_token"


def default_cache() -> Path:
    """The cache of resolved references, ``$XDG_CACHE_HOME/pylamda/bib_cache.json`` (default ``~/.cache/...``)."""
    return Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "pylamda" / "bib_cache.json"


class _Cache:
    """Resolved references: bibcode -> BibTeX, DOI/arXiv -> bibcode ("" = not in SciX), DOI -> BibTeX."""

    SECTIONS = ("bibtex", "doi_bibcode", "arxiv_bibcode", "doi_bibtex")

    def __init__(self, path: Path | None):
        self.path = path
        self.data: dict[str, dict[str, str]] = {s: {} for s in self.SECTIONS}
        self.changed = False
        if path is not None and path.is_file():
            try:
                stored = json.loads(path.read_text(encoding="utf-8"))
                for s in self.SECTIONS:
                    self.data[s].update(stored.get(s, {}))
            except (OSError, ValueError):
                warnings.warn(f"cannot read the cache {path}; starting a new one", stacklevel=2)

    def get(self, section: str, key: str) -> str | None:
        return self.data[section].get(key)

    def set(self, section: str, key: str, value: str) -> None:
        if self.data[section].get(key) != value:
            self.data[section][key] = value
            self.changed = True

    def save(self) -> None:
        if self.path is None or not self.changed:
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


class HttpClient:
    """SciX and doi.org over HTTPS (needs ``requests``, the ``bib`` extra).

    Parameters
    ----------
    token : str or None
        The SciX API token (None: SciX is not used).
    timeout : float, optional
        Timeout of each request [s]. Default: 30.
    """

    def __init__(self, token: str | None, timeout: float = 30.0):
        try:
            import requests
        except ImportError:  # pragma: no cover (requests is in the dev environment)
            raise BibError("resolving references online needs requests: pip install 'pylamda[bib]'") from None
        from . import __version__

        self._requests = requests
        self._session = requests.Session()
        self._session.headers["User-Agent"] = f"pylamda/{__version__}"
        self._token = token
        self.timeout = timeout

    @property
    def has_token(self) -> bool:
        """bool: Whether a SciX token is set."""
        return bool(self._token)

    def _scix(self, method: str, path: str, **kw) -> dict:
        headers = {"Authorization": f"Bearer {self._token}"}
        try:
            r = self._session.request(method, f"{SCIX_API}{path}", headers=headers, timeout=self.timeout, **kw)
        except self._requests.RequestException as exc:
            raise BibError(f"SciX request failed: {type(exc).__name__}") from None
        if r.status_code == 401:
            raise BibError("SciX rejected the API token (HTTP 401); check SCIX_API_TOKEN or the token file")
        if r.status_code != 200:
            raise BibError(f"SciX request {path} failed: HTTP {r.status_code}")
        return r.json()

    def scix_search(self, query: str, fields: str) -> list[dict]:
        """The documents of a SciX search (``fl`` = ``fields``, up to 2000)."""
        out = self._scix("GET", "/search/query", params={"q": query, "fl": fields, "rows": 2000})
        return out.get("response", {}).get("docs", [])

    def scix_export(self, bibcodes: list[str]) -> str:
        """BibTeX of the bibcodes (SciX export)."""
        return self._scix("POST", "/export/bibtex", json={"bibcode": bibcodes}).get("export", "")

    def doi_bibtex(self, doi: str) -> str | None:
        """BibTeX of a DOI from doi.org (content negotiation), or None if not found."""
        try:
            r = self._session.get(
                f"https://doi.org/{doi}", headers={"Accept": "application/x-bibtex"}, timeout=self.timeout
            )
        except self._requests.RequestException:
            return None
        text = r.text.strip() if r.status_code == 200 else ""
        return text if text.startswith("@") else None


# ------------------------------------------------------------------- result
@dataclass
class BibEntry:
    """One entry of the BibTeX export.

    Attributes
    ----------
    key : str
        The BibTeX key.
    bibtex : str
        The entry.
    source : str
        ``"database"`` (shipped), ``"scix"``, ``"doi.org"`` or ``"unresolved"`` (``@misc``).
    identifiers : list of str
        The identifiers that led to it (``"doi:..."``, ...).
    note : str
        What it is cited for (``"CO: collision rates with p-H2"``).
    """

    key: str
    bibtex: str
    source: str
    identifiers: list[str] = field(default_factory=list)
    note: str = ""


@dataclass
class BibResult:
    """Result of :func:`bibtex`.

    Attributes
    ----------
    entries : list of BibEntry
        The entries: database papers first, then in the order found.
    free_text : list of Reference
        Level and partner references without an identifier (not in the export).
    warnings : list of str
        What could not be resolved, and why.
    keys : dict of str to str
        Identifier -> key of every entry (for ``keys`` in the next run).
    emaa : bool
        Whether a file comes from EMAA (the export then holds EMAA's acknowledgement).
    n_files : int
        Number of data sets.
    """

    entries: list[BibEntry] = field(default_factory=list)
    free_text: list[Reference] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    keys: dict[str, str] = field(default_factory=dict)
    emaa: bool = False
    n_files: int = 0

    def text(self) -> str:
        """The ``.bib`` file: a header with the database requests, then the entries with notes."""
        from . import __version__

        head = [
            f"% References of {self.n_files} LAMDA file{'s' if self.n_files != 1 else ''}, by pylamda {__version__}."
        ]
        if any(e.key in LAMDA_KEYS for e in self.entries):
            head.append(f"% {LAMDA_REQUEST}")
        if self.emaa:
            head.append(f"% EMAA asks for this acknowledgement: {EMAA_ACKNOWLEDGEMENT}")
        parts = ["\n".join(head)]
        for e in self.entries:
            note = f"{e.note} (pylamda reference catalogue)" if e.source == "catalogue" else e.note
            parts.append((f"% {note}\n" if note else "") + e.bibtex)
        return "\n\n".join(parts) + "\n"

    def write(self, path: str | os.PathLike) -> None:
        """Write :meth:`text` to ``path``."""
        Path(path).write_text(self.text(), encoding="utf-8")


# ------------------------------------------------------------------- export
@dataclass
class _Work:
    """One cited work: its identifiers and what it is used for."""

    identifiers: list[str]
    text: str
    uses: dict[str, list[str]] = field(default_factory=dict)  # molecule -> uses
    bibcode: str | None = None
    bibtex: str | None = None
    source: str = "unresolved"
    catalogue_key: str | None = None

    def add_use(self, molecule: str, use: str) -> None:
        uses = self.uses.setdefault(molecule, [])
        if use not in uses:
            uses.append(use)

    @property
    def note(self) -> str:
        return "; ".join(f"{mol}: {', '.join(uses)}" for mol, uses in self.uses.items())


def _from_catalogue(
    cat, data: LamdaData, fname: str | None, ref: Reference, works: dict, result: BibResult
) -> list[str]:
    """Add the papers of a free-text reference that the catalogue knows; return the parts it does not know."""
    unresolved = []
    for seg, hit in cat.find_all(ref.text, molecule=data.molecule, file=fname, role=ref.role, partner=ref.partner):
        if hit is None or hit.key is None:
            unresolved.append(seg)
            continue
        ident = f"bibcode:{hit.bibcode}"
        w = works.get(ident)
        if w is None:
            ids = [ident] + ([f"doi:{hit.doi}"] if hit.doi else [])
            w = works[ident] = _Work(
                ids, seg, bibcode=hit.bibcode, bibtex=cat.bibtex.get(hit.key), source="catalogue", catalogue_key=hit.key
            )
        w.add_use(data.molecule, ref.use)
        if hit.match == "text":
            result.warnings.append(
                f"{fname or data.molecule}: {seg!r} found in the reference catalogue by its text only "
                f"({hit.key}, cited there for another molecule): check it"
            )
    return unresolved


def database_entries() -> dict[str, str]:
    """The shipped BibTeX entries of the LAMDA and EMAA papers, by key.

    Returns
    -------
    dict of str to str
        Key -> entry, from ``pylamda/data/database_refs.bib``.
    """
    text = files("pylamda").joinpath("data/database_refs.bib").read_text(encoding="utf-8")
    return {entry_key(e): e for e in split_entries(text)}


def bibtex(
    sources: Iterable[LamdaData | str | os.PathLike],
    *,
    online: bool = False,
    token: str | None = None,
    database_refs: bool = True,
    keys: dict[str, str] | None = None,
    cache: bool | str | os.PathLike = True,
    catalogue: bool = True,
    client: object | None = None,
) -> BibResult:
    """The BibTeX entries of the references of LAMDA files.

    Parameters
    ----------
    sources : iterable of LamdaData, str or os.PathLike
        Data sets or files (read with :func:`pylamda.read`).
    online : bool, optional
        Resolve the identifiers with SciX and doi.org. Default: False (``@misc``
        entries with the identifiers).
    token : str or None, optional
        SciX API token; default: :func:`scix_token` (environment, token file).
    database_refs : bool, optional
        Add the main papers of LAMDA (always) and EMAA (if a file comes from EMAA).
        Default: True.
    keys : dict of str to str or None, optional
        Identifier -> key from an earlier run, to keep keys stable.
    cache : bool, str or os.PathLike, optional
        Cache of resolved references: True = :func:`default_cache`, False = none, or a
        file. Default: True.
    catalogue : bool, optional
        Look up free-text references in the shipped reference catalogue
        (:mod:`pylamda.catalogue`); found ones become entries with source ``"catalogue"``,
        without network. Default: True.
    client : object or None, optional
        The HTTP client (default :class:`HttpClient`); tests pass a fake one.

    Returns
    -------
    BibResult
        The entries, the free-text references, warnings and the keys.

    Raises
    ------
    LamdaReadError
        If a file cannot be read.
    BibError
        If SciX rejects the token.
    """
    datas = [s if isinstance(s, LamdaData) else read(s) for s in sources]
    result = BibResult(n_files=len(datas), emaa=any(is_emaa(d) for d in datas))

    works: dict[str, _Work] = {}
    cat = load_catalogue() if catalogue else None
    for d in datas:
        fname = Path(d.source_map.path).name if d.source_map else None
        for ref in collect_references(d):
            if not ref.identifiers and ref.role in ("levels", "collisions"):
                unresolved = _from_catalogue(cat, d, fname, ref, works, result) if cat else [ref.text]
                if unresolved:
                    text = "; ".join(unresolved)
                    result.free_text.append(ref if text == ref.text else replace(ref, text=text))
            for ident in ref.identifiers:
                w = works.setdefault(ident, _Work([ident], identifier_context(ref.text, ident)))
                w.add_use(d.molecule, ref.use)

    if cat is not None:  # identifiers whose paper the catalogue knows: its reviewed BibTeX, offline
        known = {}
        for key, e in cat.entries.items():
            known[f"bibcode:{e['bibcode']}".lower()] = key
            if e.get("doi"):
                known[f"doi:{e['doi']}".lower()] = key
        for w in works.values():
            key = next((known[i.lower()] for i in w.identifiers if i.lower() in known), None)
            if key is not None and w.bibtex is None and key in cat.bibtex:
                w.bibcode, w.bibtex, w.source = cat.entries[key]["bibcode"], cat.bibtex[key], "catalogue"
                w.catalogue_key = key

    if online and works:
        cache_path = default_cache() if cache is True else (Path(cache) if cache else None)
        _resolve(list(works.values()), client or HttpClient(scix_token(token)), _Cache(cache_path), result)

    # merge works that are the same paper (same bibcode or DOI)
    merged: list[_Work] = []
    by_id: dict[str, _Work] = {}
    for w in works.values():
        ids = [i.lower() for i in w.identifiers]
        if w.bibcode:
            ids.append(f"bibcode:{w.bibcode}".lower())
        if w.bibtex and entry_field(w.bibtex, "doi"):
            ids.append(f"doi:{entry_field(w.bibtex, 'doi')}".lower())
        same = next((by_id[i] for i in ids if i in by_id), None)
        if same is None:
            merged.append(w)
            same = w
        else:
            same.identifiers += [i for i in w.identifiers if i not in same.identifiers]
            for mol, uses in w.uses.items():
                for u in uses:
                    same.add_use(mol, u)
            if same.bibtex is None and w.bibtex is not None:
                same.bibtex, same.bibcode, same.source = w.bibtex, w.bibcode, w.source
        for i in ids:
            by_id.setdefault(i, same)

    # the database papers
    used_keys: set[str] = set()
    if database_refs:
        db = database_entries()
        wanted = [*LAMDA_KEYS, *(EMAA_KEYS if result.emaa else ())]
        db_ids: dict[str, BibEntry] = {}
        for key in wanted:
            e = BibEntry(key, db[key], "database", note=_DATABASE_NOTES.get(key, "database"))
            doi = entry_field(db[key], "doi")
            adsurl = entry_field(db[key], "adsurl") or ""
            for ident in (f"doi:{doi}" if doi else None, f"bibcode:{adsurl[len(SCIX_ABS) :]}" if adsurl else None):
                if ident:
                    e.identifiers.append(ident)
                    db_ids[ident.lower()] = e
            result.entries.append(e)
            used_keys.add(key)
        kept = []
        for w in merged:
            ids = {i.lower() for i in w.identifiers} | ({f"bibcode:{w.bibcode}".lower()} if w.bibcode else set())
            hit = next((db_ids[i] for i in ids if i in db_ids), None)
            if hit is None:
                kept.append(w)
            else:  # cited by a file too: one entry, with the file's use
                hit.note = f"{hit.note}; {w.note}"
                hit.identifiers += [i for i in w.identifiers if i not in hit.identifiers]
        merged = kept

    # keys and entries
    keys = dict(keys or {})
    for w in merged:
        key = next((keys[i] for i in w.identifiers if i in keys), None)
        if key is None and w.catalogue_key and w.catalogue_key not in used_keys:
            key = w.catalogue_key
        if key is None or key in used_keys:
            base = (key_from_entry(w.bibtex) if w.bibtex else None) or key_from_text(w.text)
            base = base or re.sub(r"[^A-Za-z0-9_.:-]", "_", w.identifiers[0].replace(":", "_", 1))
            key = base
            for suffix in "abcdefghijklmnopqrstuvwxyz":
                if key not in used_keys:
                    break
                key = base + suffix
        used_keys.add(key)
        entry = _set_key(w.bibtex, key) if w.bibtex else _misc(key, w.identifiers[0], w.note)
        if w.source == "scix" and w.bibcode:
            entry = _set_field(entry, "adsurl", f"{SCIX_ABS}{w.bibcode}")
        result.entries.append(BibEntry(key, entry, w.source, list(w.identifiers), w.note))
    for e in result.entries:
        for i in e.identifiers:
            result.keys[i] = e.key
    return result


def _chunks(terms: list[str], field_name: str) -> list[str]:
    """SciX queries ``field:"a" OR field:"b" ...`` below the length limit."""
    out, current = [], ""
    for t in terms:
        part = f'{field_name}:"{t}"'
        if current and len(current) + 4 + len(part) > _SEARCH_CHUNK:
            out.append(current)
            current = ""
        current = f"{current} OR {part}" if current else part
    if current:
        out.append(current)
    return out


def _resolve(works: list[_Work], client, cache: _Cache, result: BibResult) -> None:
    """Resolve the works in place: SciX for bibcodes, DOIs, arXiv IDs; doi.org for other DOIs."""
    works = [w for w in works if w.source != "catalogue"]  # these come with their BibTeX
    kind = {w.identifiers[0]: w.identifiers[0].split(":", 1) for w in works}
    has_token = getattr(client, "has_token", True)

    # 1. DOIs and arXiv IDs -> bibcodes (SciX search)
    for w in works:
        k, v = kind[w.identifiers[0]]
        if k == "bibcode":
            w.bibcode = v
        elif k in ("doi", "arxiv"):
            hit = cache.get(f"{k}_bibcode", v.lower())
            if hit is not None:
                w.bibcode = hit or None
    if has_token:
        for k, fieldname, prefix in (("doi", "doi", ""), ("arxiv", "identifier", "arXiv:")):
            todo = sorted(
                {
                    kind[w.identifiers[0]][1]
                    for w in works
                    if kind[w.identifiers[0]][0] == k
                    and cache.get(f"{k}_bibcode", kind[w.identifiers[0]][1].lower()) is None
                }
            )
            if not todo:
                continue
            found: dict[str, str] = {}
            for query in _chunks([prefix + v for v in todo], fieldname):
                for doc in client.scix_search(query, "bibcode,doi,identifier"):
                    values = doc.get("doi" if k == "doi" else "identifier", []) or []
                    for v in todo:
                        target = (prefix + v).lower()
                        if any(str(x).lower() == target for x in values):
                            old = found.get(v.lower())
                            # prefer the journal record to an arXiv preprint record
                            if old is None or ("arXiv" in old and "arXiv" not in doc["bibcode"]):
                                found[v.lower()] = doc["bibcode"]
            for v in todo:
                cache.set(f"{k}_bibcode", v.lower(), found.get(v.lower(), ""))
            for w in works:
                kk, v = kind[w.identifiers[0]]
                if kk == k and v.lower() in found:
                    w.bibcode = found[v.lower()]

    # 2. bibcodes -> BibTeX (SciX export, one request per chunk)
    need = sorted({w.bibcode for w in works if w.bibcode})
    if has_token:
        missing = [b for b in need if cache.get("bibtex", b) is None]
        for i in range(0, len(missing), _EXPORT_CHUNK):
            for entry in split_entries(client.scix_export(missing[i : i + _EXPORT_CHUNK])):
                cache.set("bibtex", entry_key(entry), entry)
    for w in works:
        if w.bibcode and cache.get("bibtex", w.bibcode):
            w.bibtex, w.source = cache.get("bibtex", w.bibcode), "scix"

    # 3. DOIs not in SciX -> doi.org
    for w in works:
        k, v = kind[w.identifiers[0]]
        if w.bibtex is None and k == "doi":
            text = cache.get("doi_bibtex", v.lower())
            if text is None:
                text = client.doi_bibtex(v)
                if text:
                    cache.set("doi_bibtex", v.lower(), text)
            entries = split_entries(text or "")
            if entries:
                w.bibtex, w.source = entries[0], "doi.org"
            else:
                result.warnings.append(f"{w.identifiers[0]}: not found in SciX or at doi.org")
    cache.save()

    if not has_token:
        n = sum(1 for w in works if w.bibtex is None and kind[w.identifiers[0]][0] in ("bibcode", "arxiv"))
        if n:
            result.warnings.append(
                f"{n} bibcodes/arXiv IDs not resolved: SciX needs an API token (SCIX_API_TOKEN or {token_file()})"
            )
    for w in works:
        if w.bibcode and w.bibtex is None and has_token:
            result.warnings.append(f"bibcode {w.bibcode}: not found in SciX")
        elif kind[w.identifiers[0]][0] == "arxiv" and w.bibtex is None and has_token:
            result.warnings.append(f"{w.identifiers[0]}: not found in SciX")
