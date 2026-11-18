"""Identifiers in reference texts: DOIs, ADS/SciX bibcodes and arXiv IDs."""

from __future__ import annotations

import re
from urllib.parse import unquote

# a DOI, with or without "https://doi.org/", "doi:" or "DOI " in front
_DOI = re.compile(r"(?:https?://(?:dx\.)?doi\.org/|\bdoi\s*[:=]?\s*|(?<![\w./]))(10\.\d{4,9}/[^\s,;\"'<>{}]+)", re.I)
# bibcodes in SciX or ADS links
_BIBCODE_LINK = re.compile(
    r"(?:scixplorer\.org|ui\.adsabs\.harvard\.edu|adsabs\.harvard\.edu)/(?:#?abs/|cgi-bin/nph-bib_query\?bibcode=)"
    r"([^/\s?#&\"'<>]+(?:%26[^/\s?#\"'<>]+)*)",
    re.I,
)
# a bibcode on its own: YYYYJJJJJVVVVMPPPPA (19 characters)
_BIBCODE = re.compile(r"(?<![\w.&/])(\d{4}[A-Za-z&.]{5}[\w.&]{4}[A-Za-z.][\w.]{4}[A-Z.])(?![\w&])")
# arXiv IDs: new style (2004.11230, with optional version) and old style (astro-ph/0411110)
_ARXIV = re.compile(
    r"(?:arxiv\s*:\s*|arxiv\.org/(?:abs|pdf)/)(\d{4}\.\d{4,5}|[a-z-]+(?:\.[A-Z]{2})?/\d{7})(?:v\d+)?", re.I
)


def find_dois(text: str | None) -> list[str]:
    """DOIs in a reference text (``https://doi.org/10...``, ``doi:10...``, ``DOI 10...``, ``10...``).

    Parameters
    ----------
    text : str or None
        Reference text.

    Returns
    -------
    list of str
        The DOIs (``10.xxxx/...``) in the order found, without trailing ``.``, ``)`` or
        ``]``, each once; empty if there are none.
    """
    out: list[str] = []
    for m in _DOI.findall(text or ""):
        doi = m.rstrip(".)]")
        if doi.lower() not in (d.lower() for d in out):
            out.append(doi)
    return out


def find_bibcodes(text: str | None) -> list[str]:
    """ADS/SciX bibcodes in a reference text: from SciX and ADS links, and on their own.

    Parameters
    ----------
    text : str or None
        Reference text, e.g. ``"Yang (2010) https://scixplorer.org/abs/2010ApJ...718.1062Y"``.

    Returns
    -------
    list of str
        The bibcodes (19 characters, ``%26`` in links decoded to ``&``) in the order
        found, each once.
    """
    out: list[str] = []
    text = text or ""
    for m in _BIBCODE_LINK.findall(text):
        code = unquote(m)
        if len(code) == 19 and code not in out:
            out.append(code)
    for code in _BIBCODE.findall(text):
        if code not in out:
            out.append(code)
    return out


def find_arxiv(text: str | None) -> list[str]:
    """arXiv IDs in a reference text (``arXiv:2004.11230``, ``arxiv.org/abs/astro-ph/0411110``).

    Parameters
    ----------
    text : str or None
        Reference text.

    Returns
    -------
    list of str
        The IDs without ``arXiv:`` and version (``2004.11230``, ``astro-ph/0411110``),
        in the order found, each once.
    """
    out: list[str] = []
    for m in _ARXIV.findall(text or ""):
        if m not in out:
            out.append(m)
    return out


def find_identifiers(text: str | None) -> list[str]:
    """All identifiers in a text, as ``"doi:..."``, ``"bibcode:..."`` and ``"arxiv:..."``.

    A bibcode that is part of a DOI (none in practice) and a DOI of an arXiv preprint
    (``10.48550/arXiv....``, given as arXiv ID instead) are not repeated.

    Parameters
    ----------
    text : str or None
        Reference text.

    Returns
    -------
    list of str
        The identifiers, DOIs first, then bibcodes, then arXiv IDs.
    """
    out = []
    arxiv = find_arxiv(text)
    for doi in find_dois(text):
        m = re.fullmatch(r"10\.48550/arxiv\.(.+)", doi, re.I)
        if m:
            if m.group(1) not in arxiv:
                arxiv.append(m.group(1))
            continue
        out.append(f"doi:{doi}")
    out += [f"bibcode:{b}" for b in find_bibcodes(text)]
    out += [f"arxiv:{a}" for a in arxiv]
    return out
