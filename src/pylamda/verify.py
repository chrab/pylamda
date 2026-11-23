"""Verify a LAMDA file or data set: the checks S, L, R, C, F and P (design document, section 6).

:func:`verify` returns the structure issues of the reader (``S``) and the issues of the
checks of the data:

* ``L``: levels;
* ``R``: radiative transitions;
* ``C``: collisions;
* ``F``: references;
* ``P``: physics heuristics (``info``; they can give false alarms; off by default).

Checks of many rows (rates, transitions) give one issue per block or section, with the
number of rows (:attr:`Issue.count <pylamda.issues.Issue.count>`), the first example and
its line.
"""

from __future__ import annotations

import math
import os
import re
from collections.abc import Callable, Iterable
from pathlib import Path

import numpy as np

from .catalogue import Catalogue, load_catalogue
from .formula import formula_mass
from .issues import CODES, Issue, make_issue, match_codes
from .model import PARTNERS, BlockSource, CollisionSet, LamdaData, SourceMap, partner_name
from .reader import read
from .references import find_dois

#: Speed of light [GHz cm]: frequency [GHz] = energy difference [cm-1] * C_GHZ_CM.
C_GHZ_CM = 29.9792458
#: Second radiation constant hc/k [K cm]: energy [K] = energy [cm-1] * HC_K.
HC_K = 1.438776877

#: Largest plausible downward rate coefficient [cm3 s-1] per partner ID (check P03):
#: 1e-8 for neutral partners, 1e-7 for H+, 1e-4 for electrons.
P03_MAX_RATE: dict[int, float] = {1: 1e-8, 2: 1e-8, 3: 1e-8, 4: 1e-4, 5: 1e-8, 6: 1e-8, 7: 1e-7}
#: Check P01 uses the temperatures below this fraction of ΔE/k.
P01_T_FRACTION = 0.5
#: Check P01 flags a row if its slope d ln k / d(1/T) is -ΔE/k times a factor in this range.
P01_SLOPE_RANGE = (2 / 3, 3 / 2)
#: Check P04 flags a mass that differs from the formula mass by more than this [Da].
P04_TOLERANCE = 0.5

_NUMBER = re.compile(r"[+-]?(\d*)(?:\.(\d*))?(?:[eEdD]([+-]?\d+))?")
# links, DOIs, arXiv IDs and ADS/SciX bibcodes (YYYYJJJJJVVVVMPPPPA)
_IDENTIFIER = re.compile(
    r"https?://|\b10\.\d{4,9}/|\barxiv\s*:?\s*\d{4}\.\d{4,5}|\b\d{4}[A-Za-z&.]{5}[\w.]{4}[A-Za-z.][\w.]{4}[A-Z.]",
    re.I,
)
_REFERENCE_LIKE = re.compile(
    r"\bet\s+al\b|\b(?:19|20)\d\d\b|\b(?:CDMS|JPL|NIST|HITRAN|ExoMol|CHIANTI|BASECOL|LAMDA|EMAA)\b", re.I
)
_FIX_NOTE = re.compile(r"^\s*!\s*pylamda \S+ fix:")  # a note line of pylamda fix (see pylamda.fix.FIX_NOTE)
# C13: the partner text says the rates are scaled
_SCALED = re.compile(r"\bscal(?:ed|ing)\b|\bfrom\s+(?:[op]-|ortho-|para-)?H2\b", re.I)
# C12: "scaled"/"scaling", or a scaling factor such as "*1.39", "x 1.4", "×1.37"
_HE_SCALED = re.compile(r"\bscal(?:ed|ing)\b|(?:\*|×|\bx)\s*1\.\d+", re.I)


class _Collector:
    """Collects the issues of the checks."""

    def __init__(self) -> None:
        self.issues: list[Issue] = []

    def __call__(
        self,
        code: str,
        message: str,
        line: int | None = None,
        *,
        fixable: bool = False,
        count: int = 1,
        severity: str | None = None,
    ) -> None:
        self.issues.append(make_issue(code, message, line, fixable, count, severity))


# ------------------------------------------------------------------- helpers
def _halfunit(token: str | None) -> float:
    """Half a unit of the last digit of a number as written (0 if unknown, and for 0, which is exact)."""
    if not token:
        return 0.0
    m = _NUMBER.fullmatch(token)
    if not m or not (m.group(1) or m.group(2)) or not (m.group(1) + (m.group(2) or "")).strip("0"):
        return 0.0
    return 0.5 * 10.0 ** (int(m.group(3) or 0) - len(m.group(2) or ""))


def _token(smap: SourceMap | None, line: int | None, index: int) -> str | None:
    """The ``index``-th whitespace-separated value of a line of the file, if there is one."""
    if smap is None or line is None or not 0 < line <= len(smap.raw):
        return None
    parts = smap.raw[line - 1].split()
    return parts[index] if index < len(parts) else None


def _lines(smap: SourceMap | None, attr: str, n: int) -> np.ndarray | None:
    """Line numbers of the levels or lines, if the source map still fits the data."""
    if smap is None:
        return None
    lines = getattr(smap, attr)
    return lines if len(lines) == n else None


def _block_source(data: LamdaData, b: int) -> BlockSource | None:
    """The BlockSource of collision block ``b`` (0-based), if the source map still fits."""
    smap = data.source_map
    if smap is None or len(smap.blocks) != len(data.collisions):
        return None
    src = smap.blocks[b]
    return src if len(src.rows) == data.collisions[b].n_transitions else None


def _at(lines: np.ndarray | None, k: int) -> int | None:
    """Line number of item ``k``, or None."""
    return int(lines[k]) if lines is not None else None


def _first(idx: np.ndarray | list[int], describe: Callable[[int], str]) -> str:
    """The first item described, and how many more there are."""
    more = f" and {len(idx) - 1} more" if len(idx) > 1 else ""
    return describe(int(idx[0])) + more


def _duplicates(pairs: Iterable[tuple]) -> list[tuple[int, int]]:
    """(index, index of the first occurrence) of every repeated item."""
    seen: dict[tuple, int] = {}
    out = []
    for k, key in enumerate(pairs):
        if key in seen:
            out.append((k, seen[key]))
        else:
            seen[key] = k
    return out


def _block_name(b: int, c: CollisionSet) -> str:
    return f"collision block {b + 1} ({c.label})"


# ------------------------------------------------------------------- L: levels
def _check_levels(data: LamdaData, add: _Collector) -> None:
    lines = _lines(data.source_map, "levels", data.n_levels)
    e, g = data.level_energy, data.level_g
    if data.n_levels == 0:
        return

    def below(k: int) -> str:
        return f"level {k + 1} ({e[k]:g} cm-1) below level {k} ({e[k - 1]:g} cm-1)"

    def weight(k: int) -> str:
        return f"level {k + 1} (g = {g[k]:g})"

    down = np.flatnonzero(np.diff(e) < 0) + 1
    if down.size:
        add(
            "L01", f"energies not increasing: {_first(down, below)}", _at(lines, down[0]), fixable=True, count=down.size
        )

    bad = np.flatnonzero(~np.isfinite(g) | (g <= 0))
    if bad.size:
        add("L02", f"weight g <= 0 or not finite: {_first(bad, weight)}", _at(lines, bad[0]), count=bad.size)
    frac = np.flatnonzero(np.isfinite(g) & (g > 0) & (g != np.round(g)))
    if frac.size:
        add(
            "L02",
            f"weight g not an integer: {_first(frac, weight)}",
            _at(lines, frac[0]),
            count=frac.size,
            severity="info",
        )

    quanta = [" ".join(q.split()) for q in data.level_quanta]
    dupes = _duplicates((float(e[k]), float(g[k]), quanta[k]) for k in range(data.n_levels))
    if dupes:
        k, j = dupes[0]
        more = f" and {len(dupes) - 1} more" if len(dupes) > 1 else ""
        add(
            "L03", f"level {k + 1} has the same E, g and quanta as level {j + 1}{more}", _at(lines, k), count=len(dupes)
        )

    empty = [k for k, q in enumerate(quanta) if not q]
    if empty:
        add("L04", f"{len(empty)} of {data.n_levels} levels without quanta", _at(lines, empty[0]), count=len(empty))
    groups: dict[str, list[int]] = {}
    for k, q in enumerate(quanta):
        if q:
            groups.setdefault(q, []).append(k)
    repeated = [ks for ks in groups.values() if len(ks) > 1]
    if repeated:
        ks = repeated[0]
        more = f" ({len(repeated) - 1} more quanta repeated)" if len(repeated) > 1 else ""
        add(
            "L04",
            f"quanta not unique: {quanta[ks[0]]!r} at levels {', '.join(str(k + 1) for k in ks)}{more}",
            _at(lines, ks[1]),
            count=sum(len(ks) for ks in repeated),
        )

    finite = e[np.isfinite(e)]
    if finite.size and finite.min() != 0:
        k = int(np.flatnonzero(e == finite.min())[0])
        add("L05", f"the lowest level (level {k + 1}) has E = {finite.min():g} cm-1, not 0", _at(lines, k))


# ------------------------------------------------------------------- R: radiative transitions
def _check_lines(data: LamdaData, add: _Collector, rtol: float) -> None:
    smap = data.source_map
    lines = _lines(smap, "lines", data.n_lines)
    up, low, a = data.line_up, data.line_low, data.line_a
    nlev = data.n_levels
    if data.n_lines == 0:
        return

    def name(k: int) -> str:
        return f"line {k + 1} ({int(up[k])} -> {int(low[k])})"

    def with_a(k: int) -> str:
        return f"{name(k)}, A = {a[k]:g}"

    outside = (up < 1) | (up > nlev) | (low < 1) | (low > nlev) | (up == low)
    bad = np.flatnonzero(outside)
    if bad.size:
        add(
            "R01",
            f"level numbers outside 1..{nlev} or upper = lower: {_first(bad, name)}",
            _at(lines, bad[0]),
            count=bad.size,
        )
    ok = np.flatnonzero(~outside)
    eu = np.full(len(up), np.nan)
    el = np.full(len(up), np.nan)
    eu[ok], el[ok] = data.level_energy[up[ok] - 1], data.level_energy[low[ok] - 1]

    inverted = ok[eu[ok] <= el[ok]]
    if inverted.size:
        add("R02", f"E(up) <= E(low): {_first(inverted, name)}", _at(lines, inverted[0]), count=inverted.size)

    bad = np.flatnonzero(~np.isfinite(a) | (a <= 0))
    if bad.size:
        add("R03", f"Einstein A <= 0 or not finite: {_first(bad, with_a)}", _at(lines, bad[0]), count=bad.size)

    for code, what, values in (("R04", "frequency", data.line_freq), ("R05", "E_up [K]", data.line_eup)):
        missing = np.flatnonzero(np.isnan(values))
        if missing.size:
            add(
                code,
                f"{what} missing: {_first(missing, name)}",
                _at(lines, missing[0]),
                fixable=True,
                count=missing.size,
            )
    good = ok[eu[ok] > el[ok]]
    _compare_column(data, add, "R04", (eu - el) * C_GHZ_CM, good, rtol, name, upper_only=False)
    _compare_column(data, add, "R05", eu * HC_K, ok, rtol, name, upper_only=True)

    dupes = _duplicates((int(u), int(lo)) for u, lo in zip(up, low, strict=True))
    if dupes:
        k, j = dupes[0]
        more = f" and {len(dupes) - 1} more" if len(dupes) > 1 else ""
        add("R06", f"{name(k)} repeats line {j + 1}{more}", _at(lines, k), count=len(dupes))


def _compare_column(
    data: LamdaData,
    add: _Collector,
    code: str,
    computed: np.ndarray,
    rows: np.ndarray,
    rtol: float,
    name: Callable[[int], str],
    *,
    upper_only: bool,
) -> None:
    """R04/R05: report the frequency or E_up values that differ from the level energies."""
    lines = _lines(data.source_map, "lines", data.n_lines)
    given, what, unit = (data.line_freq, "frequency", "GHz") if code == "R04" else (data.line_eup, "E_up", "K")
    bad, worst, worst_k = _mismatch_rows(data, code, computed, rows, rtol)
    if bad:
        source = "E(up)" if upper_only else "the level energies"
        k = worst_k
        add(
            code,
            f"{what} differs from {source} beyond rounding and rtol {rtol:g} in {len(bad)} lines; largest: "
            f"{name(k)}: {given[k]:.10g} {unit}, computed {computed[k]:.10g} (relative {worst:.2g})",
            _at(lines, bad[0]),
            fixable=True,
            count=len(bad),
        )


def _mismatch_rows(
    data: LamdaData, code: str, computed: np.ndarray, rows: np.ndarray, rtol: float
) -> tuple[list[int], float, int]:
    """Rows (0-based) whose frequency (R04) or E_up (R05) differs from ``computed``.

    The tolerance is ``rtol`` (relative) plus the rounding of the numbers as written in
    the file: half a unit of the last digit of the column value and of the energies.
    Returns the rows, the largest relative difference and its row.
    """
    smap = data.source_map
    lines = _lines(smap, "lines", data.n_lines)
    level_lines = _lines(smap, "levels", data.n_levels)
    if code == "R04":
        given, column, factor, upper_only = data.line_freq, 4, C_GHZ_CM, False
    else:
        given, column, factor, upper_only = data.line_eup, 5, HC_K, True
    bad: list[int] = []
    worst, worst_k = -1.0, -1
    for k in (int(k) for k in rows):
        if not (math.isfinite(given[k]) and math.isfinite(computed[k])):
            continue
        tol = rtol * abs(computed[k]) + _halfunit(_token(smap, _at(lines, k), column))
        if level_lines is not None:
            tol += factor * _halfunit(_token(smap, _at(level_lines, int(data.line_up[k]) - 1), 1))
            if not upper_only:
                tol += factor * _halfunit(_token(smap, _at(level_lines, int(data.line_low[k]) - 1), 1))
        diff = abs(given[k] - computed[k])
        if diff > tol:
            bad.append(k)
            rel = diff / abs(computed[k]) if computed[k] else math.inf
            if rel > worst:
                worst, worst_k = rel, k
    return bad, worst, worst_k


def line_columns(data: LamdaData) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Frequency [GHz] and E_up [K] of each line computed from the level energies.

    Parameters
    ----------
    data : LamdaData
        The data.

    Returns
    -------
    freq, eup : numpy.ndarray
        Computed frequency and upper-level energy; NaN for lines with level numbers
        outside 1..NLEV.
    valid : numpy.ndarray
        Indices of the lines with valid level numbers and E(up) > E(low).
    """
    up, low, n = data.line_up, data.line_low, data.n_levels
    ok = np.flatnonzero((up >= 1) & (up <= n) & (low >= 1) & (low <= n) & (up != low))
    eu = np.full(len(up), np.nan)
    el = np.full(len(up), np.nan)
    eu[ok], el[ok] = data.level_energy[up[ok] - 1], data.level_energy[low[ok] - 1]
    valid = ok[eu[ok] > el[ok]]
    return (eu - el) * C_GHZ_CM, eu * HC_K, valid


def column_mismatches(data: LamdaData, rtol: float = 1e-4) -> dict[str, list[int]]:
    """Lines whose frequency (R04) or E_up (R05) differs from the level energies.

    Parameters
    ----------
    data : LamdaData
        The data (rounding of the numbers from ``data.source_map``, if it fits).
    rtol : float, optional
        Relative tolerance on top of the rounding. Default: 1e-4.

    Returns
    -------
    dict
        ``{"R04": [line indices], "R05": [line indices]}`` (0-based).
    """
    freq, eup, valid = line_columns(data)
    return {
        "R04": _mismatch_rows(data, "R04", freq, valid, rtol)[0],
        "R05": _mismatch_rows(data, "R05", eup, valid, rtol)[0],
    }


# ------------------------------------------------------------------- C: collisions
def _check_partner(data: LamdaData, b: int, c: CollisionSet, add: _Collector, line: int | None) -> None:
    """C01, C02, C11, C12: the partner ID and the name in the partner text."""
    if c.partner is None:
        add(
            "C01",
            f"{_block_name(b, c)}: partner ID {c.partner_id} is not in 1..7 (1=H2 2=p-H2 3=o-H2 4=e 5=H 6=He 7=H+)",
            line,
        )
    written = c.text_partner
    if not written:
        return
    canonical = partner_name(written)
    text = c.partner_text.strip()
    if canonical not in PARTNERS.values():
        add(
            "C11",
            f"partner line {text!r}: {written!r} is not a standard partner; the block is read by its "
            f"ID {c.partner_id} as {c.label}",
            line,
        )
    elif canonical == "He" and c.partner in ("H2", "p-H2", "o-H2") and _HE_SCALED.search(text):
        add(
            "C12",
            f"partner line {text!r}: the name says He, the ID {c.partner_id} is {c.partner}, and the text says the "
            f"rates are scaled: probably He rates scaled to {c.partner}, with the name not changed "
            f"(read as {c.partner}; ProDiMo would read He)",
            line,
            fixable=True,
        )
    elif c.partner is not None and canonical != c.partner:
        add(
            "C02",
            f"partner line {text!r}: ID {c.partner_id} is {c.partner}, the name says {canonical} "
            f"(read as {c.partner}; ProDiMo would read {canonical})",
            line,
            fixable=True,
        )


def _check_rows(data: LamdaData, b: int, c: CollisionSet, add: _Collector, rows: np.ndarray | None) -> None:
    """C06-C10: the rate rows of one block."""
    block = _block_name(b, c)
    nlev = data.n_levels

    def row(k: int) -> str:
        return f"row {k + 1} ({int(c.up[k])} -> {int(c.low[k])})"

    r = c.rates
    if r.size:
        bad = np.flatnonzero(np.any(~np.isfinite(r) | (r < 0), axis=1))
        if bad.size:
            add(
                "C06",
                f"{block}: rates < 0 or not finite in {bad.size} rows: {_first(bad, row)}",
                _at(rows, bad[0]),
                count=bad.size,
            )
        zero = np.flatnonzero(np.all(r == 0, axis=1))
        if zero.size:
            add(
                "C07",
                f"{block}: {zero.size} rows are 0 at all temperatures: {_first(zero, row)}",
                _at(rows, zero[0]),
                count=zero.size,
            )

    up, low = c.up, c.low
    outside = (up < 1) | (up > nlev) | (low < 1) | (low > nlev) | (up == low)
    bad = np.flatnonzero(outside)
    if bad.size:
        add(
            "C10",
            f"{block}: level numbers outside 1..{nlev} or upper = lower: {_first(bad, row)}",
            _at(rows, bad[0]),
            count=bad.size,
        )
    ok = np.flatnonzero(~outside)
    e = data.level_energy
    inverted = ok[e[up[ok] - 1] < e[low[ok] - 1]]
    if inverted.size:
        add(
            "C08",
            f"{block}: E(up) < E(low) in {inverted.size} rows (wrong order or excitation rates): "
            f"{_first(inverted, row)}",
            _at(rows, inverted[0]),
            fixable=True,
            count=inverted.size,
        )

    dupes = _duplicates((int(u), int(lo)) for u, lo in zip(up, low, strict=True))
    if dupes:
        k, j = dupes[0]
        more = f" and {len(dupes) - 1} more" if len(dupes) > 1 else ""
        add("C09", f"{block}: {row(k)} repeats row {j + 1}{more}", _at(rows, k), count=len(dupes))


def _check_collisions(data: LamdaData, add: _Collector) -> None:
    first_block: dict[int, int] = {}
    for b, c in enumerate(data.collisions):
        src = _block_source(data, b)
        line = src.partner if src else None
        _check_partner(data, b, c, add, line)
        text = (c.partner_text or c.reference).strip()
        if _SCALED.search(text):
            add(
                "C13",
                f"{_block_name(b, c)}: the partner text says the rates are scaled ({text!r}); codes that scale "
                "rates themselves (e.g. ProDiMo for He) may scale twice",
                line,
            )
        if c.partner_id in first_block:
            add(
                "C03",
                f"{_block_name(b, c)}: partner ID {c.partner_id} also in collision block "
                f"{first_block[c.partner_id] + 1}",
                line,
            )
        else:
            first_block[c.partner_id] = b
        t = c.temperatures
        if t.size and (not np.all(np.isfinite(t)) or np.any(t <= 0) or np.any(np.diff(t) <= 0)):
            add(
                "C05",
                f"{_block_name(b, c)}: temperatures not strictly increasing or <= 0: {' '.join(f'{x:g}' for x in t)}",
                src.temperatures if src else None,
            )
        _check_rows(data, b, c, add, src.rows if src else None)

    partners = data.partners
    if "H2" in partners and ("p-H2" in partners or "o-H2" in partners):
        both = "/".join(p for p in ("p-H2", "o-H2") if p in partners)
        b = partners.index("H2")
        src = _block_source(data, b)
        add(
            "C04",
            f"H2 (collision block {b + 1}) together with {both}: some codes count H2 collisions twice",
            src.partner if src else None,
        )


# ------------------------------------------------------------------- F: references
def _has_identifier(text: str) -> bool:
    """Whether a text holds a DOI, a link, an arXiv ID or a bibcode."""
    return bool(find_dois(text) or _IDENTIFIER.search(text))


def _in_catalogue(cat: Catalogue | None, data: LamdaData, role: str, partner: str | None, text: str) -> bool:
    """Whether every paper of a free-text reference is in the reference catalogue."""
    if cat is None or not text.strip():
        return False
    fname = Path(data.source_map.path).name if data.source_map else None
    hits = cat.find_all(text, molecule=data.molecule, file=fname, role=role, partner=partner)
    return all(hit is not None and hit.key is not None for _, hit in hits)


def _check_references(data: LamdaData, add: _Collector, cat: Catalogue | None = None) -> None:
    smap = data.source_map
    nlev_line = smap.nlev if smap else None
    # texts that can hold the reference of the levels and lines
    texts = [data.level_reference, data.title, *data.notes]
    if data.dataset_doi:
        texts.append(f"doi:{data.dataset_doi}")
    if smap is not None:
        partner_lines = {s.partner for s in smap.blocks}
        texts += [
            line for i, line in enumerate(smap.raw, 1) if line.lstrip().startswith("!") and i not in partner_lines
        ]
    texts = [t for t in texts if t and t.strip() and not _FIX_NOTE.match(t)]  # not pylamda's own fix notes
    if not any(_has_identifier(t) or _REFERENCE_LIKE.search(t) for t in texts):
        add(
            "F01",
            "no reference for the levels and lines (not on the NLEV line, line 2, in comments or notes)",
            nlev_line,
            fixable=True,
        )
    elif not any(_has_identifier(t) for t in texts) and not _in_catalogue(
        cat, data, "levels", None, data.level_reference
    ):
        add(
            "F03",
            "the levels and lines have a reference only as free text (no DOI, bibcode or link, "
            "not in the reference catalogue)",
            nlev_line,
            fixable=True,
        )
    for b, c in enumerate(data.collisions):
        src = _block_source(data, b)
        line = src.partner if src else None
        if not c.reference.strip():
            add("F02", f"{_block_name(b, c)}: no reference on the partner line", line, fixable=True)
        elif not _has_identifier(c.reference) and not _in_catalogue(cat, data, "collisions", c.label, c.reference):
            add("F03", f"{_block_name(b, c)}: reference only as free text: {c.reference.strip()!r}", line, fixable=True)


# ------------------------------------------------------------------- P: physics heuristics
def excitation_slope(temperatures: np.ndarray, rates: np.ndarray, delta_e: float) -> float | None:
    """Slope d ln k / d(1/T) [K] of a rate row at low temperature (check P01).

    An excitation rate falls like exp(-ΔE/kT) at T ≪ ΔE/k, so its slope is about
    -ΔE/k. A de-excitation rate changes slowly with T; its slope is much smaller.

    Parameters
    ----------
    temperatures : numpy.ndarray
        Temperatures [K].
    rates : numpy.ndarray
        The rate coefficients of one row [cm3 s-1].
    delta_e : float
        Energy difference of the two levels, ΔE/k [K].

    Returns
    -------
    float or None
        The least-squares slope over the temperatures below ``P01_T_FRACTION * delta_e``
        with a positive rate, or None if there are fewer than two.
    """
    use = (temperatures < P01_T_FRACTION * delta_e) & np.isfinite(rates) & (rates > 0)
    if use.sum() < 2:
        return None
    x = 1.0 / temperatures[use]
    y = np.log(rates[use])
    x0 = x - x.mean()
    return float((x0 * (y - y.mean())).sum() / (x0 * x0).sum())


def _check_physics(data: LamdaData, add: _Collector) -> None:
    e = data.level_energy
    nlev = data.n_levels
    for b, c in enumerate(data.collisions):
        src = _block_source(data, b)
        rows = src.rows if src else None
        block = _block_name(b, c)

        flagged = []  # P01
        for k in range(c.n_transitions):
            u, lo = int(c.up[k]), int(c.low[k])
            if not (1 <= u <= nlev and 1 <= lo <= nlev):
                continue
            delta_e = (e[u - 1] - e[lo - 1]) * HC_K
            if not delta_e > 0:
                continue
            slope = excitation_slope(c.temperatures, c.rates[k], delta_e)
            if slope is not None and P01_SLOPE_RANGE[0] <= -slope / delta_e <= P01_SLOPE_RANGE[1]:
                flagged.append((k, slope, delta_e))
        if flagged:
            k, slope, delta_e = flagged[0]
            more = f" and {len(flagged) - 1} more" if len(flagged) > 1 else ""
            add(
                "P01",
                f"{block}: {len(flagged)} rows fall like exp(-ΔE/kT) at T < {P01_T_FRACTION:g} ΔE/k, as "
                f"excitation rates do: row {k + 1} ({int(c.up[k])} -> {int(c.low[k])}): d ln k/d(1/T) = "
                f"{slope:.3g} K, ΔE/k = {delta_e:.3g} K{more}",
                _at(rows, k),
                count=len(flagged),
            )

        limit = P03_MAX_RATE.get(c.partner_id)
        if limit is not None and c.rates.size:
            big = np.flatnonzero(np.any(np.isfinite(c.rates) & (c.rates > limit), axis=1))
            if big.size:
                k = int(big[0])
                add(
                    "P03",
                    f"{block}: rates > {limit:g} cm3 s-1 in {big.size} rows: row {k + 1} "
                    f"({int(c.up[k])} -> {int(c.low[k])}), max {np.nanmax(c.rates[k]):.3g}",
                    _at(rows, k),
                    count=big.size,
                )

    expected = formula_mass(data.molecule)
    if expected is not None and math.isfinite(data.mass) and abs(data.mass - expected) > P04_TOLERANCE:
        smap = data.source_map
        add(
            "P04",
            f"mass {data.mass:g} differs from the formula mass {expected:.3f} of {data.molecule!r}",
            smap.mass if smap else None,
        )


# ------------------------------------------------------------------- verify
def check_data(data: LamdaData, rtol: float = 1e-4, catalogue: bool = True, physics: bool = True) -> list[Issue]:
    """The issues of the checks L, R, C, F and P of ``data`` (not the reader's S issues).

    Parameters
    ----------
    data : LamdaData
        The data. Line numbers come from ``data.source_map`` where it still fits the data.
    rtol : float, optional
        Relative tolerance of the checks R04 and R05. Default: 1e-4.
    catalogue : bool, optional
        F03: a free-text reference found in the reference catalogue is not reported.
        Default: True.
    physics : bool, optional
        Run the physics heuristics (P checks). Default: True.

    Returns
    -------
    list of Issue
        The issues, in the order of the checks.
    """
    add = _Collector()
    _check_levels(data, add)
    _check_lines(data, add, rtol)
    _check_collisions(data, add)
    _check_references(data, add, load_catalogue() if catalogue else None)
    if physics:
        _check_physics(data, add)
    return add.issues


def verify(
    source: str | os.PathLike | LamdaData,
    *,
    physics: bool | None = None,
    select: Iterable[str] | None = None,
    ignore: Iterable[str] = (),
    error_on: Iterable[str] = (),
    rtol: float = 1e-4,
    catalogue: bool = True,
) -> list[Issue]:
    """Check a LAMDA file or data set against the standard.

    Parameters
    ----------
    source : str, os.PathLike or LamdaData
        A file (read with :func:`pylamda.read`) or data already read.
    physics : bool or None, optional
        Run the physics heuristics (P checks): True runs them, False does not (the same
        as ``ignore=["P"]``). Default: None, they run only if ``select`` names P codes
        (``select=["P"]``, ``["P01"]``). They are slow on large files (P01 fits every rate
        row) and can give false alarms.
    select : iterable of str or None, optional
        Report only these codes (prefixes, as in ruff: ``"C"``, ``"C0"``, ``"C02"``).
        Default: all.
    ignore : iterable of str, optional
        Do not report these codes (prefixes).
    error_on : iterable of str, optional
        Raise the severity of these codes (prefixes) to ``"error"``, e.g. ``["C02"]`` for
        files that must have no C02 issue.
    rtol : float, optional
        Relative tolerance of the checks R04 (frequency) and R05 (E_up), on top of the
        rounding of the numbers as written. Default: 1e-4.
    catalogue : bool, optional
        Use the reference catalogue for F03 (:mod:`pylamda.catalogue`). Default: True.

    Returns
    -------
    list of Issue
        The issues, sorted by line (issues without a line first).

    Raises
    ------
    LamdaReadError
        If ``source`` is a file that cannot be read (see :func:`pylamda.read`).
    OSError
        If the file cannot be opened.
    """
    select = None if select is None else list(select)
    if physics is None:
        physics = select is not None and any(match_codes(code, select) for code in CODES if code.startswith("P"))
    data = source if isinstance(source, LamdaData) else read(source)
    issues = [*data.issues, *check_data(data, rtol=rtol, catalogue=catalogue, physics=physics)]
    # a partner line without the ID (S06) is read as ID 0: do not report it again as C01
    s06_lines = {i.line for i in issues if i.code == "S06"}
    issues = [i for i in issues if not (i.code == "C01" and i.line in s06_lines and i.line is not None)]
    ignore = [*ignore, *([] if physics else ["P"])]
    error_on = list(error_on)
    out = []
    for issue in issues:
        if select is not None and not match_codes(issue.code, select):
            continue
        if match_codes(issue.code, ignore):
            continue
        if match_codes(issue.code, error_on) and issue.severity != "error":
            issue = make_issue(issue.code, issue.message, issue.line, issue.fixable, issue.count, "error")
        out.append(issue)
    return sorted(out, key=lambda i: i.line or 0)
