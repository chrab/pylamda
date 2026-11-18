"""Identifiers in reference texts (DOIs)."""

from __future__ import annotations

import re

_DOI = re.compile(r"(?:https?://(?:dx\.)?doi\.org/|\bdoi\s*[:=]?\s*)(10\.\d{4,9}/[^\s,;]+)", re.I)


def find_dois(text: str | None) -> list[str]:
    """DOIs in a reference text (``https://doi.org/10...``, ``doi:10...``, ``DOI 10...``).

    Parameters
    ----------
    text : str or None
        Reference text.

    Returns
    -------
    list of str
        The DOIs (``10.xxxx/...``) in the order found, without trailing ``.`` or ``)``;
        empty if there are none.
    """
    return [m.rstrip(".)") for m in _DOI.findall(text or "")]
