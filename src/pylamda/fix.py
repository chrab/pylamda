"""Fix LAMDA files: what is unambiguous by default, the rest only with an explicit option.

Design document, sections 7 and 10.2. Two output modes:

* **normalised** (:func:`fix`, ``pylamda fix in.dat -o out.dat``): the fixed data, written
  in pylamda's layout (:func:`pylamda.write`);
* **patch** (:func:`patch`, ``pylamda fix in.dat --in-place``): only the lines that
  change are edited, token by token, keeping the column widths, the number style, the
  comments and the line endings. Small git diffs for hand-maintained files.

The fixes, in the order applied:

=========  =============================================================  ====================
code       fix                                                            needs
=========  =============================================================  ====================
S05        add ``!`` to a comment line without it                         default
S04        drop surplus values (what all readers do)                      default
S02        correct the running numbers of the rows (no reader uses them)  default
S03        count extra blocks in (``include``) or comment them out        ``extra_blocks``
C02, C12   correct the partner name (``trust="id"``) or the ID            ``trust``
           (``trust="name"``)
L01        sort the levels by energy and renumber all transitions         default
C08        inverted rate rows: ``swap``, ``chianti`` or ``balance``       ``inverted``
R04, R05   compute missing frequencies and E_up; differing ones only      default /
           with ``recompute``                                             ``recompute``
F01-F03    append identifiers to reference texts                          ``refs``
=========  =============================================================  ====================

Every fix is recorded as a :class:`FixRecord` and as a note line
``! pylamda <version> fix: <code> <what>`` at the end of the file.
"""

from __future__ import annotations

import math
import os
import re
from collections.abc import Iterable
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from .catalogue import load_catalogue, normalise, split_segments
from .issues import LamdaReadError
from .model import PARTNER_IDS, LamdaData, SourceMap, partner_name
from .reader import read_string
from .references import find_identifiers
from .verify import HC_K, column_mismatches, line_columns
from .writer import write

#: The options of ``inverted`` (check C08).
INVERTED_METHODS: tuple[str, ...] = ("swap", "chianti", "balance")
#: The options of ``extra_blocks`` (check S03).
EXTRA_BLOCK_METHODS: tuple[str, ...] = ("include", "comment")
#: The options of ``trust`` (checks C02, C12).
TRUST_OPTIONS: tuple[str, ...] = ("id", "name")


#: A note line written by :func:`fix` (``! pylamda 0.1.0 fix: ...``); not a reference text.
FIX_NOTE = re.compile(r"^\s*!\s*pylamda \S+ fix:")


class FixError(ValueError):
    """A fix cannot be applied (e.g. the extra data do not form a collision block)."""


@dataclass(frozen=True)
class FixRecord:
    """One fix that was applied.

    Attributes
    ----------
    code : str
        The check code of the issue fixed (``"C02"``, ``"L01"``, ...).
    message : str
        What was changed.
    count : int
        Number of items changed.
    """

    code: str
    message: str
    count: int = 1

    def note(self) -> str:
        """The note line written to the file (``! pylamda 0.1.0 fix: C02 ...``)."""
        from . import __version__

        return f"! pylamda {__version__} fix: {self.code} {self.message}"


# ------------------------------------------------------------------- editing lines
_TOKEN = re.compile(r"\S+")
_NUMBER = re.compile(r"([+-]?)(\d*)(?:\.(\d*))?(?:([eEdD])([+-]?)(\d+))?")


def replace_tokens(line: str, replacements: dict[int, str]) -> str:
    """Replace whitespace-separated tokens of a line, keeping the layout.

    A replaced token keeps its right end where the blanks before it allow (right-aligned
    columns stay aligned); at least one blank stays between tokens.

    Parameters
    ----------
    line : str
        The line.
    replacements : dict of int to str
        Token index (0-based) -> new text.

    Returns
    -------
    str
        The edited line.
    """
    out, pos = [], 0
    for i, m in enumerate(_TOKEN.finditer(line)):
        gap, token = line[pos : m.start()], m.group()
        if i in replacements:
            new = replacements[i]
            diff = len(new) - len(token)
            if diff > 0 and gap:
                keep = 1 if m.start() > 0 else 0
                cut = min(diff, len(gap) - keep) if gap.strip(" ") == "" else 0
                gap = gap[: len(gap) - cut]
            elif diff < 0 and gap and gap.strip(" ") == "":
                gap += " " * (-diff)
            token = new
        out += [gap, token]
        pos = m.end()
    out.append(line[pos:])
    return "".join(out)


def format_like(token: str, value: float) -> str:
    """``value`` written in the style of ``token`` (decimals, exponent letter and digits).

    At least three significant digits are kept in exponential notation.

    Parameters
    ----------
    token : str
        A number as written in the file (``"1.2e-11"``, ``"3.0D-10"``, ``"115.27"``).
    value : float
        The new value.

    Returns
    -------
    str
        The new text.
    """
    m = _NUMBER.fullmatch(token)
    if not math.isfinite(value) or not m or not (m.group(2) or m.group(3)):
        return repr(float(value)) if math.isfinite(value) else "NaN"
    decimals = len(m.group(3) or "")
    if m.group(4):
        text = f"{value:.{max(decimals, 2)}e}"
        mant, exp = text.split("e")
        sign = "-" if int(exp) < 0 else ("+" if m.group(5) else "")
        digits = str(abs(int(exp))).rjust(len(m.group(6)), "0")
        return f"{mant}{m.group(4)}{sign}{digits}"
    if "." in token:
        return f"{value:.{decimals}f}"
    return str(round(value))


class _Doc:
    """The lines of a file being patched, with the source map of its items."""

    def __init__(self, text: str, data: LamdaData):
        self.newline = "\r\n" if "\r\n" in text else "\n"
        self.final_newline = text.endswith(("\n", "\r"))
        self.lines = text.splitlines()
        self.smap: SourceMap = data.source_map

    def text(self) -> str:
        return self.newline.join(self.lines) + (self.newline if self.final_newline or not self.lines else "")

    def tokens(self, line: int) -> list[str]:
        return self.lines[line - 1].split()

    def replace(self, line: int, replacements: dict[int, str]) -> None:
        self.lines[line - 1] = replace_tokens(self.lines[line - 1], replacements)

    def group(self, first: int, n: int) -> list[tuple[int, int]]:
        """(line, token index) of ``n`` values that start on line ``first`` and may continue."""
        out: list[tuple[int, int]] = []
        line = first
        while len(out) < n and line <= len(self.lines):
            text = self.lines[line - 1]
            if text.strip() and not text.lstrip().startswith("!"):
                for t in range(len(text.split())):
                    if len(out) < n:
                        out.append((line, t))
            line += 1
        if len(out) < n:
            raise FixError(f"line {first}: fewer than {n} values")
        return out

    def append(self, text: str) -> None:
        self.lines.append(text)


# ------------------------------------------------------------------- stage 1: structure
def _fix_structure(text: str, data: LamdaData, extra_blocks: str | None, log: list[FixRecord]) -> str:
    """S05, S03, S04, S02 on the raw lines; returns the new text.

    S03 comes before S04 and S02, so that blocks counted in are fixed as well.
    """
    smap = data.source_map
    doc = _Doc(text, data)
    for issue in smap.issues:
        if issue.code == "S05" and issue.line:
            doc.lines[issue.line - 1] = "!" + doc.lines[issue.line - 1]
            log.append(FixRecord("S05", f"line {issue.line}: '!' added to the comment line"))
    has_s03 = any(i.code == "S03" for i in smap.issues)
    if has_s03 and extra_blocks == "comment":
        start = smap.notes or len(doc.lines) + 1
        n = 0
        for i in range(start, len(doc.lines) + 1):
            line = doc.lines[i - 1]
            if line.strip() and not line.lstrip().startswith("!"):
                doc.lines[i - 1] = "!" + line
                n += 1
        log.append(FixRecord("S03", f"{n} data lines after the last declared block commented out", n))
    elif has_s03 and extra_blocks == "include":
        npart = len(data.collisions)
        for k in range(1, 21):
            trial = list(doc.lines)
            trial[smap.npart - 1] = replace_tokens(trial[smap.npart - 1], {0: str(npart + k)})
            try:
                back = read_string("\n".join(trial))
            except LamdaReadError:
                continue
            if not any(i.code in ("S03", "S01") for i in back.issues):
                doc.lines = trial
                log.append(
                    FixRecord("S03", f"NPART {npart} -> {npart + k}: the blocks after the last one counted in", k)
                )
                break
        else:
            raise FixError("the data after the last declared block do not form complete collision blocks")
    # S04, on the blocks as they are now (with the ones counted in)
    current = read_string("\n".join(doc.lines))
    n_surplus = 0
    for b, src in enumerate(current.source_map.blocks):
        ntemp = len(current.collisions[b].temperatures)
        for first, n in [(src.temperatures, ntemp)] + [(int(r), 3 + ntemp) for r in src.rows]:
            line, t = doc.group(first, n)[-1]
            parts = list(_TOKEN.finditer(doc.lines[line - 1]))
            if len(parts) > t + 1:
                doc.lines[line - 1] = doc.lines[line - 1][: parts[t].end()]
                n_surplus += len(parts) - t - 1
    if n_surplus:
        log.append(FixRecord("S04", f"{n_surplus} surplus values dropped", n_surplus))
    # S02: the running numbers of the level, line and rate rows (no reader uses them)
    smap = current.source_map
    sections = [("levels", smap.levels), ("radiative transitions", smap.lines)] + [
        (f"collision block {b + 1}", src.rows) for b, src in enumerate(smap.blocks)
    ]
    for what, rows in sections:
        n = 0
        for k, line in enumerate(rows):
            line = int(line)
            first = doc.tokens(line)[0] if doc.tokens(line) else ""
            if first != str(k + 1) and not (first.lstrip("+").isdigit() and int(first) == k + 1):
                doc.replace(line, {0: str(k + 1)})
                n += 1
        if n:
            log.append(FixRecord("S02", f"{what}: {n} row number{'s' if n != 1 else ''} corrected", n))
    return doc.text()


# ------------------------------------------------------------------- stage 2: data
def _partner_line(molecule: str, pid: int, name: str, reference: str) -> str:
    return f"{pid} {molecule} + {name}" + (f" -- Ref: {reference}" if reference else "")


def _fix_partners(data: LamdaData, doc: _Doc | None, trust: str | None, log: list[FixRecord]) -> None:
    """C02 and C12: make the name in the partner text and the ID agree."""
    if trust is None:
        return
    for b, c in enumerate(data.collisions):
        written = c.text_partner
        canonical = partner_name(written) if written else None
        if canonical is None or canonical not in PARTNER_IDS or c.partner is None or canonical == c.partner:
            continue
        code = "C12" if canonical == "He" and c.partner in ("H2", "p-H2", "o-H2") else "C02"
        if trust == "id":
            new_id, new_name = c.partner_id, c.partner
            what = f"name {written} -> {new_name} (trusted ID {c.partner_id})"
        else:
            new_id, new_name = PARTNER_IDS[canonical], canonical
            what = f"ID {c.partner_id} -> {new_id} (trusted name {written})"
        line = _partner_line(data.molecule, new_id, new_name, c.reference)
        data.collisions[b] = replace(c, partner_id=new_id, partner_text=line)
        if doc is not None:
            doc.lines[doc.smap.blocks[b].partner - 1] = line
        log.append(FixRecord(code, f"collision block {b + 1}: {what}"))


def _fix_level_order(data: LamdaData, doc: _Doc | None, log: list[FixRecord]) -> None:
    """L01: sort the levels by energy (stable) and renumber all transitions."""
    order = np.argsort(data.level_energy, kind="stable")
    if np.array_equal(order, np.arange(data.n_levels)):
        return
    new_number = np.empty(data.n_levels, dtype=int)
    new_number[order] = np.arange(1, data.n_levels + 1)

    def renumber(x: np.ndarray) -> np.ndarray:
        ok = (x >= 1) & (x <= data.n_levels)
        return np.where(ok, new_number[np.clip(x, 1, data.n_levels) - 1], x)

    if doc is not None:
        smap = doc.smap
        old_lines = [doc.lines[int(n) - 1] for n in smap.levels]
        for i, k in enumerate(order):
            doc.lines[int(smap.levels[i]) - 1] = replace_tokens(old_lines[k], {0: str(i + 1)})
        for k, line in enumerate(smap.lines):
            up, low = int(data.line_up[k]), int(data.line_low[k])
            nu, nl = int(renumber(np.array([up]))[0]), int(renumber(np.array([low]))[0])
            if (nu, nl) != (up, low):
                doc.replace(int(line), {1: str(nu), 2: str(nl)})
        for b, c in enumerate(data.collisions):
            for k, line in enumerate(smap.blocks[b].rows):
                up, low = int(c.up[k]), int(c.low[k])
                nu, nl = int(renumber(np.array([up]))[0]), int(renumber(np.array([low]))[0])
                if (nu, nl) != (up, low):
                    pos = doc.group(int(line), 3)
                    doc.replace(pos[1][0], {pos[1][1]: str(nu)})
                    doc.replace(pos[2][0], {pos[2][1]: str(nl)})
    n_moved = int((order != np.arange(data.n_levels)).sum())
    data.level_energy = data.level_energy[order]
    data.level_g = data.level_g[order]
    data.level_quanta = [data.level_quanta[k] for k in order]
    data.line_up, data.line_low = renumber(data.line_up), renumber(data.line_low)
    data.collisions = [replace(c, up=renumber(c.up), low=renumber(c.low)) for c in data.collisions]
    log.append(FixRecord("L01", f"levels sorted by energy ({n_moved} moved), transitions renumbered", n_moved))


def _fix_inverted(data: LamdaData, doc: _Doc | None, method: str | None, log: list[FixRecord]) -> None:
    """C08: rate rows with E(up) < E(low)."""
    if method is None:
        return
    e, g, n = data.level_energy, data.level_g, data.n_levels
    for b, c in enumerate(data.collisions):
        ok = (c.up >= 1) & (c.up <= n) & (c.low >= 1) & (c.low <= n)
        rows = np.flatnonzero(ok & (e[np.clip(c.up, 1, n) - 1] < e[np.clip(c.low, 1, n) - 1]))
        if not rows.size:
            continue
        up, low, rates = c.up.copy(), c.low.copy(), c.rates.copy()
        for k in rows:
            a, z = int(c.up[k]), int(c.low[k])  # E(a) < E(z): new upper level z, new lower level a
            if method == "chianti":
                rates[k] = c.rates[k] * g[a - 1] / g[z - 1]
            elif method == "balance":
                delta = (e[z - 1] - e[a - 1]) * HC_K
                with np.errstate(over="ignore"):
                    factor = g[a - 1] / g[z - 1] * np.exp(delta / c.temperatures)
                if not np.all(np.isfinite(factor)):
                    raise FixError(
                        f"collision block {b + 1} row {k + 1}: detailed balance overflows at the lowest "
                        f"temperatures (ΔE/k = {delta:.3g} K); use --inverted swap or chianti"
                    )
                rates[k] = c.rates[k] * factor
            up[k], low[k] = z, a
            if doc is not None:
                ntemp = len(c.temperatures)
                pos = doc.group(int(doc.smap.blocks[b].rows[k]), 3 + ntemp)
                doc.replace(pos[1][0], {pos[1][1]: str(z)})
                doc.replace(pos[2][0], {pos[2][1]: str(a)})
                if method != "swap":
                    for j, (line, t) in enumerate(pos[3:]):
                        old = doc.tokens(line)[t]
                        doc.replace(line, {t: format_like(old, float(rates[k, j]))})
        data.collisions[b] = replace(c, up=up, low=low, rates=rates)
        what = {
            "swap": "up and low swapped",
            "chianti": "swapped, rates times g(low)/g(up)",
            "balance": "excitation rates turned into de-excitation rates (detailed balance)",
        }[method]
        log.append(
            FixRecord("C08", f"collision block {b + 1} ({c.label}): {rows.size} inverted rows: {what}", rows.size)
        )


def _fix_columns(data: LamdaData, doc: _Doc | None, missing: bool, recompute: bool, log: list[FixRecord]) -> None:
    """R04, R05: missing frequencies and E_up computed; differing ones with ``recompute``."""
    freq, eup, valid = line_columns(data)
    differ = column_mismatches(data) if recompute else {"R04": [], "R05": []}
    valid_set = set(int(k) for k in valid)
    counts = {"R04": 0, "R05": 0}
    for k in range(data.n_lines):
        if k not in valid_set:
            continue
        set_freq = (missing and math.isnan(data.line_freq[k])) or k in differ["R04"]
        set_eup = (
            (missing and math.isnan(data.line_eup[k]))
            or k in differ["R05"]
            or (set_freq and math.isnan(data.line_eup[k]))
        )
        if not (set_freq or set_eup):
            continue
        if set_freq:
            data.line_freq[k] = float(f"{freq[k]:.10g}")
            counts["R04"] += 1
        if set_eup:
            data.line_eup[k] = float(f"{eup[k]:.6g}")
            counts["R05"] += 1
        if doc is not None:
            line = int(doc.smap.lines[k])
            tokens = doc.tokens(line)
            texts = {4: f"{data.line_freq[k]:.10g}", 5: f"{data.line_eup[k]:.6g}"}
            flags = {4: set_freq, 5: set_eup}
            # values that are there: replaced (recomputed values with their full precision)
            doc.replace(line, {t: texts[t] for t in (4, 5) if flags[t] and t < len(tokens)})
            # missing columns: appended (E_up only after a frequency)
            for t in (4, 5):
                if flags[t] and len(doc.tokens(line)) == t:
                    doc.lines[line - 1] = doc.lines[line - 1].rstrip() + f" {texts[t]:>12s}"
    for code, what in (("R04", "frequencies"), ("R05", "E_up values")):
        if counts[code]:
            log.append(FixRecord(code, f"{counts[code]} {what} computed from the level energies", counts[code]))


# ------------------------------------------------------------------- references
def _identifier_text(identifier: str) -> tuple[str, str]:
    """(normalised identifier ``kind:value``, text to append) of a DOI, bibcode or arXiv ID."""
    ident = identifier.strip()
    found = find_identifiers(ident)
    if not found and re.fullmatch(r"\d{4}\.\d{4,5}(v\d+)?|[a-z-]+/\d{7}", ident):
        found = [f"arxiv:{ident}"]
    if len(found) != 1:
        raise FixError(f"{identifier!r} is not one DOI, bibcode or arXiv ID")
    kind, value = found[0].split(":", 1)
    link = f"https://scixplorer.org/abs/{value.replace('&', '%26')}"  # "&" encoded, as ADS does
    text = {"doi": f"doi:{value}", "bibcode": link, "arxiv": f"arXiv:{value}"}[kind]
    return found[0], text


def _add_reference(data: LamdaData, doc: _Doc | None, target: str, identifier: str, log: list[FixRecord]) -> None:
    ident, text = _identifier_text(identifier)
    key = ident.lower()
    if target == "dataset":
        if not ident.startswith("doi:"):
            raise FixError("the dataset reference must be a DOI")
        doi = ident[4:]
        if data.dataset_doi:
            if data.dataset_doi.lower() == doi.lower():
                return
            raise FixError(f"the data already have the dataset DOI {data.dataset_doi}")
        data.dataset_doi = doi
        if doc is not None:
            label = "EMAA_DATASET_DOI" if doi.startswith("10.17178/EMAA") else "DATASET_DOI"
            line = doc.smap.molecule
            doc.lines[line - 1] = doc.lines[line - 1].rstrip() + f" -- {label}: {doi}"
        log.append(FixRecord("F01", f"dataset DOI {doi} added"))
        return
    if target == "levels":
        _append_to_levels(data, doc, key, text, log)
        return
    name = partner_name(target) or target
    blocks = [b for b, c in enumerate(data.collisions) if c.partner == name]
    if not blocks:
        raise FixError(f"no collision partner {target!r} (targets: levels, dataset, or a partner name)")
    if len(blocks) > 1:
        raise FixError(f"several blocks have partner {name}; the reference target is ambiguous")
    _append_to_block(data, doc, blocks[0], key, text, log)


def _append_to_levels(
    data: LamdaData, doc: _Doc | None, key: str, text: str, log: list[FixRecord], how: str = ""
) -> None:
    """Append an identifier (as ``text``) to the level reference, unless it is there."""
    if key in (i.lower() for i in find_identifiers(data.level_reference)):
        return
    data.level_reference = f"{data.level_reference} {text}".strip()
    if doc is not None:
        line = doc.smap.nlev
        old = doc.lines[line - 1].rstrip()
        doc.lines[line - 1] = f"{old} {text}" if "Ref:" in old else f"{old} - Ref: {text}"
    log.append(FixRecord("F03", f"levels: {text} added to the reference{how}"))


def _append_to_block(
    data: LamdaData, doc: _Doc | None, b: int, key: str, text: str, log: list[FixRecord], how: str = ""
) -> None:
    """Append an identifier (as ``text``) to the reference of collision block ``b``, unless it is there."""
    c = data.collisions[b]
    if key in (i.lower() for i in find_identifiers(c.reference)):
        return
    new_line = f"{c.partner_text.rstrip()} {text}" if c.partner_text else ""
    data.collisions[b] = replace(c, reference=f"{c.reference} {text}".strip(), partner_text=new_line)
    if doc is not None:
        line = doc.smap.blocks[b].partner
        doc.lines[line - 1] = doc.lines[line - 1].rstrip() + f" {text}"
    code = "F03" if c.reference.strip() else "F02"
    log.append(FixRecord(code, f"collision block {b + 1} ({c.label}): {text} added{how}"))


# ------------------------------------------------------------------- identifiers from the catalogue
def _free_text_references(data: LamdaData) -> list[tuple[int | None, str, str | None, str]]:
    """(block or None for the levels, role, partner, text) of references without any identifier."""
    out: list[tuple[int | None, str, str | None, str]] = []
    if data.level_reference.strip() and not find_identifiers(data.level_reference):
        out.append((None, "levels", None, data.level_reference))
    for b, c in enumerate(data.collisions):
        if c.reference.strip() and not find_identifiers(c.reference):
            out.append((b, "collisions", c.label, c.reference))
    return out


def _catalogue_refs(data: LamdaData, fname: str | None) -> list[tuple[int | None, str]]:
    """(block or None, bibcode) of the papers the shipped catalogue knows (by place or molecule)."""
    cat = load_catalogue()
    out = []
    for b, role, partner, text in _free_text_references(data):
        for _, hit in cat.find_all(text, molecule=data.molecule, file=fname, role=role, partner=partner):
            if hit is not None and hit.key is not None and hit.match in ("place", "molecule"):
                out.append((b, hit.bibcode))
    return out


def _review_refs(data: LamdaData, fname: str | None, review: Path) -> list[tuple[int | None, str]]:
    """(block or None, bibcode) of the papers decided in a review file for this file."""
    from .refsearch import load

    decided = {}
    for g in load(review):
        if g.bibcode and g.molecule == data.molecule:
            for u in g.uses:
                if u.file == fname:
                    decided[(u.role, u.partner, g.normalised)] = g.bibcode
    out = []
    for b, role, partner, text in _free_text_references(data):
        for seg in split_segments(text):
            bib = decided.get((role, partner, normalise(seg)))
            if bib:
                out.append((b, bib))
    return out


def add_reference(
    data: LamdaData, target: str, *, doi: str | None = None, bibcode: str | None = None, arxiv: str | None = None
) -> tuple[LamdaData, list[FixRecord]]:
    """Append an identifier to a reference text (design document, section 10.2).

    The existing text stays; the identifier is appended: a bibcode as SciX link
    (``https://scixplorer.org/abs/<bibcode>``), a DOI as ``doi:<doi>``, an arXiv ID as
    ``arXiv:<id>``. An identifier already in the text is not added again.

    Parameters
    ----------
    data : LamdaData
        The data (not changed).
    target : str
        ``"levels"`` (the level reference), ``"dataset"`` (the dataset DOI, only a DOI) or
        a partner name (``"p-H2"``, any spelling).
    doi, bibcode, arxiv : str or None
        The identifier; exactly one must be given.

    Returns
    -------
    data : LamdaData
        A copy with the identifier and a note line.
    log : list of FixRecord
        What was changed (empty if the identifier was there already).

    Raises
    ------
    FixError
        If the target is unknown or ambiguous, or not exactly one identifier is given.
    """
    candidates = (doi and f"doi:{doi}", bibcode and f"https://scixplorer.org/abs/{bibcode}", arxiv and f"arXiv:{arxiv}")
    given = [x for x in candidates if x]
    if len(given) != 1:
        raise FixError("give exactly one of doi, bibcode, arxiv")
    out = data.copy()
    log: list[FixRecord] = []
    _add_reference(out, None, target, given[0], log)
    out.notes = [*out.notes, *(r.note() for r in log)]
    return out, log


# ------------------------------------------------------------------- the fix
def _parse_refs(refs: dict[str, str] | Iterable[tuple[str, str]] | Iterable[str] | None) -> list[tuple[str, str]]:
    if not refs:
        return []
    if isinstance(refs, dict):
        return list(refs.items())
    out = []
    for r in refs:
        if isinstance(r, str):
            target, sep, ident = r.partition("=")
            if not sep:
                raise FixError(f"--ref {r!r}: use TARGET=ID, e.g. p-H2=2010ApJ...718.1062Y")
            out.append((target.strip(), ident.strip()))
        else:
            out.append((r[0], r[1]))
    return out


def _run(
    source: LamdaData | str | os.PathLike,
    *,
    keep_layout: bool,
    trust: str | None,
    sort_levels: bool,
    inverted: str | None,
    missing_columns: bool,
    recompute: bool,
    extra_blocks: str | None,
    refs,
    add_refs: bool | str | os.PathLike = False,
) -> tuple[LamdaData, _Doc | None, list[FixRecord]]:
    if trust not in (None, *TRUST_OPTIONS):
        raise FixError(f"trust must be 'id' or 'name', not {trust!r}")
    if inverted not in (None, *INVERTED_METHODS):
        raise FixError(f"inverted must be one of {', '.join(INVERTED_METHODS)}, not {inverted!r}")
    if extra_blocks not in (None, *EXTRA_BLOCK_METHODS):
        raise FixError(f"extra_blocks must be 'include' or 'comment', not {extra_blocks!r}")
    ref_list = _parse_refs(refs)
    if isinstance(source, LamdaData):
        data = source.copy()
        smap = data.source_map
        text = "\n".join(smap.raw) + "\n" if smap is not None else None
        name = smap.path if smap is not None else "<data>"
    else:
        name = str(source)
        text = Path(source).read_bytes().decode("utf-8", errors="replace")  # keeps "\r\n"
        data = read_string(text, name)
    if keep_layout and text is None:
        raise FixError("patch mode needs data read from a file")
    log: list[FixRecord] = []
    if text is not None and data.source_map is not None:
        text = _fix_structure(text, data, extra_blocks, log)
        data = read_string(text, name)
    elif extra_blocks == "include":
        raise FixError("extra_blocks='include' needs data read from a file")
    doc = _Doc(text, data) if keep_layout else None
    _fix_partners(data, doc, trust, log)
    if sort_levels:
        _fix_level_order(data, doc, log)
    _fix_inverted(data, doc, inverted, log)
    _fix_columns(data, doc, missing_columns, recompute, log)
    for target, identifier in ref_list:
        _add_reference(data, doc, target, identifier, log)
    if add_refs:
        fname = Path(name).name if name != "<data>" else None
        review = None if add_refs is True else Path(add_refs)
        found = _review_refs(data, fname, review) if review else _catalogue_refs(data, fname)
        how = " (review file)" if review else " (reference catalogue)"
        for b, bibcode in found:
            ident, text = _identifier_text(bibcode)
            if b is None:
                _append_to_levels(data, doc, ident.lower(), text, log, how)
            else:
                _append_to_block(data, doc, b, ident.lower(), text, log, how)
    notes = [r.note() for r in log]
    data.notes = [*data.notes, *notes]
    if doc is not None:
        for line in notes:
            doc.append(line)
    data.source_map = None
    return data, doc, log


def fix(
    source: LamdaData | str | os.PathLike,
    *,
    trust: str | None = None,
    sort_levels: bool = True,
    inverted: str | None = None,
    missing_columns: bool = True,
    recompute: bool = False,
    extra_blocks: str | None = None,
    refs: dict[str, str] | Iterable[tuple[str, str]] | Iterable[str] | None = None,
    add_refs: bool | str | os.PathLike = False,
) -> tuple[LamdaData, list[FixRecord]]:
    """Fix a LAMDA data set: the unambiguous issues by default, the others with an option.

    Parameters
    ----------
    source : LamdaData, str or os.PathLike
        The data or file (not changed). Structure fixes (S03-S05) need the file's lines
        (a file, or data read from one).
    trust : {"id", "name"} or None, optional
        C02/C12, the partner name and ID disagree: ``"id"`` corrects the name,
        ``"name"`` the ID. Default: None (not fixed; there is no safe default).
    sort_levels : bool, optional
        L01: sort the levels by energy and renumber all transitions. Default: True.
    inverted : {"swap", "chianti", "balance"} or None, optional
        C08, rate rows with E(up) < E(low): ``"swap"`` exchanges up and low and keeps the
        rate (wrong order only); ``"chianti"`` also multiplies the rate by g(new low) /
        g(new up) (rates divided by the wrong g); ``"balance"`` turns excitation rates
        into de-excitation rates by detailed balance. Default: None (not fixed).
    missing_columns : bool, optional
        R04/R05: compute missing frequencies and E_up values from the level energies.
        Default: True.
    recompute : bool, optional
        R04/R05: also replace values that differ from the level energies. Default: False.
    extra_blocks : {"include", "comment"} or None, optional
        S03, data after the last declared block: count them in (NPART raised) or comment
        them out. Default: None (they stay as notes).
    refs : dict, iterable of (target, id) or of "TARGET=ID", optional
        Identifiers to append to reference texts (:func:`add_reference`): target
        ``"levels"``, ``"dataset"`` or a partner name; id a DOI, bibcode or arXiv ID.
    add_refs : bool, str or os.PathLike, optional
        Append the SciX links of free-text references (references without any
        identifier): True takes them from the shipped reference catalogue (only papers
        found at the same place or for the same molecule); a path takes them from the
        decisions in a review file (``pylamda bib --search``, ``tools/catalogue.py``).
        Default: False.

    Returns
    -------
    data : LamdaData
        The fixed data (a copy), with one note line per fix.
    log : list of FixRecord
        The fixes applied.

    Raises
    ------
    FixError
        If an option is unknown or a fix cannot be applied.
    LamdaReadError
        If the file cannot be read.
    """
    data, _, log = _run(
        source,
        keep_layout=False,
        trust=trust,
        sort_levels=sort_levels,
        inverted=inverted,
        missing_columns=missing_columns,
        recompute=recompute,
        extra_blocks=extra_blocks,
        refs=refs,
        add_refs=add_refs,
    )
    return data, log


def patch(
    source: LamdaData | str | os.PathLike,
    *,
    trust: str | None = None,
    sort_levels: bool = True,
    inverted: str | None = None,
    missing_columns: bool = True,
    recompute: bool = False,
    extra_blocks: str | None = None,
    refs: dict[str, str] | Iterable[tuple[str, str]] | Iterable[str] | None = None,
    add_refs: bool | str | os.PathLike = False,
) -> tuple[str, list[FixRecord]]:
    """Fix a LAMDA file and return its text with only the changed lines edited.

    The same fixes and options as :func:`fix`. Tokens are replaced in place: column
    widths, number styles, comments and line endings stay. The note lines are appended
    at the end of the file.

    Parameters
    ----------
    source : LamdaData, str or os.PathLike
        The file, or data read from a file.
    trust, sort_levels, inverted, missing_columns, recompute, extra_blocks, refs, add_refs
        As in :func:`fix`.

    Returns
    -------
    text : str
        The patched file.
    log : list of FixRecord
        The fixes applied (empty: ``text`` is the file as it was).

    Raises
    ------
    FixError
        If an option is unknown, a fix cannot be applied, or ``source`` is data not read
        from a file.
    LamdaReadError
        If the file cannot be read.
    """
    _, doc, log = _run(
        source,
        keep_layout=True,
        trust=trust,
        sort_levels=sort_levels,
        inverted=inverted,
        missing_columns=missing_columns,
        recompute=recompute,
        extra_blocks=extra_blocks,
        refs=refs,
        add_refs=add_refs,
    )
    return doc.text(), log


def fix_file(path: str | os.PathLike, output: str | os.PathLike, *, keep_layout: bool, **options) -> list[FixRecord]:
    """Fix a file and write the result (normalised or patched).

    Parameters
    ----------
    path : str or os.PathLike
        The input file.
    output : str or os.PathLike
        The output file (may be ``path``).
    keep_layout : bool
        Patch mode (True) or the normalised layout (False).
    **options
        The options of :func:`fix`.

    Returns
    -------
    list of FixRecord
        The fixes applied.
    """
    if keep_layout:
        text, log = patch(path, **options)
        Path(output).write_bytes(text.encode("utf-8"))
        return log
    data, log = fix(path, **options)
    write(data, output)
    return log


# what each fixable issue needs, for hints
OPTION_HINTS: dict[str, str] = {
    "C02": "--trust id|name",
    "C12": "--trust id|name",
    "C08": "--inverted swap|chianti|balance",
    "S03": "--extra-blocks include|comment",
    "R04": "--recompute",
    "R05": "--recompute",
    "F01": "--ref levels=ID",
    "F02": "--ref PARTNER=ID",
    "F03": "--ref TARGET=ID or --add-refs",
}
