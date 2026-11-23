"""Read LAMDA files by position, as RADEX and ProDiMo do (Fortran list-directed reads).

Rules (design document, section 5):

* Lines whose first non-blank character is ``!`` are comments; blank lines are skipped.
* The declared counts decide; header text is not used.
* Surplus values on a line are ignored and give an issue (S04).
* Everything after the last declared block is kept as notes (data there: S03).
* Collision partners are identified by their ID only.
* Deviations are collected as :class:`~pylamda.issues.Issue` objects, not logged.
* The header (molecule name, molecular weight, number of levels, level rows, number of
  radiative transitions) must read; else the file is not a LAMDA file
  (:class:`~pylamda.issues.NotLamdaError`). Lines are read only as needed, so the reader
  stops there without reading the rest of a large file.
"""

from __future__ import annotations

import math
import os
import re
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import NoReturn, TextIO

import numpy as np

from .issues import Issue, LamdaReadError, NotLamdaError, make_issue
from .model import BlockSource, CollisionSet, LamdaData, SourceMap, split_partner_line

# numbers as Fortran list-directed input reads them (D exponents allowed, no commas)
_REAL = re.compile(r"[+-]?(?:(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?|nan|inf|infinity)", re.I)
_INT = re.compile(r"[+-]?\d+")
_DATASET_DOI = re.compile(r"\s*--\s*(?:\w+_)?DATASET_DOI:\s*(\S+)\s*$")
# a number without "nan"/"inf": a molecule name that reads as one is not a molecule name
_NUMBER = re.compile(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?")

#: A line longer than this (characters) means the file is not a LAMDA file (e.g. a
#: binary file without line ends); longer lines are never read whole.
MAX_LINE_LENGTH = 1_000_000


def _file_lines(f: TextIO) -> Iterator[str]:
    """The lines of an open file, one at a time, split as :meth:`str.splitlines` splits them.

    A line longer than :data:`MAX_LINE_LENGTH` is yielded as its first
    ``MAX_LINE_LENGTH + 1`` characters (the reader stops there).
    """
    while True:
        piece = f.readline(MAX_LINE_LENGTH + 1)
        if not piece:
            return
        if len(piece) > MAX_LINE_LENGTH and not piece.endswith("\n"):  # text mode: all line ends are "\n"
            yield piece
            return
        yield from piece.splitlines()


class _Reader:
    """The non-comment lines of a file, read one value or line at a time, with the issues found.

    Lines are taken from ``lines`` only when needed. While ``header`` is True, a failure
    means that the file is not a LAMDA file (:class:`NotLamdaError`).
    """

    def __init__(self, path: str, lines: Iterable[str]):
        self.path = path
        self.raw: list[str] = []  # the lines taken so far (all of them after rest())
        self.data: list[tuple[int, str]] = []  # (line number, text) of the non-comment lines taken so far
        self._lines = iter(lines)
        self.pos = 0
        self.issues: list[Issue] = []
        self.header = True

    def issue(self, code: str, message: str, line: int | None = None, fixable: bool = False) -> None:
        """Record an issue."""
        self.issues.append(make_issue(code, message, line, fixable))

    def not_lamda(self, reason: str, line: int | None = None) -> NoReturn:
        """Stop reading: the file is not a LAMDA file."""
        raise NotLamdaError(self.path, self.issues, f"line {line}: {reason}" if line else reason)

    def fail(self, code: str, message: str, line: int | None = None) -> NoReturn:
        """Record an issue and stop reading (in the header: not a LAMDA file)."""
        self.issue(code, message, line)
        if self.header:
            self.not_lamda(message, line)
        raise LamdaReadError(self.path, self.issues)

    def _take(self) -> bool:
        """Take lines up to the next non-comment line; False at the end of the file."""
        for text in self._lines:
            self.raw.append(text)
            lineno = len(self.raw)
            if len(text) > MAX_LINE_LENGTH:
                self.not_lamda(f"longer than {MAX_LINE_LENGTH} characters", lineno)
            if self.header and "\x00" in text:
                self.not_lamda("binary data (NUL characters)", lineno)
            if text.strip() and not text.lstrip().startswith("!"):
                self.data.append((lineno, text))
                return True
        return False

    def line(self, what: str) -> tuple[int, str]:
        """The next non-comment line: (line number, text)."""
        if self.pos >= len(self.data) and not self._take():
            self.fail("S01", f"the file ends before {what}", len(self.raw) or None)
        self.pos += 1
        return self.data[self.pos - 1]

    def values(self, n: int, what: str) -> tuple[list[str], int, bool]:
        """``n`` whitespace-separated values, continuing on the next lines if needed.

        Returns the values, the number of the first line, and whether the last line held
        surplus values (ignored).
        """
        out: list[str] = []
        first = 0
        surplus = False
        while len(out) < n:
            lineno, text = self.line(what)
            first = first or lineno
            tokens = text.split()
            surplus = len(out) + len(tokens) > n
            out += tokens
        return out[:n], first, surplus

    def real(self, token: str, lineno: int, what: str) -> float:
        """A real number; NaN and an S06 issue if Fortran cannot read it."""
        if _REAL.fullmatch(token):
            return float(token.replace("D", "E").replace("d", "e"))
        self.issue("S06", f"cannot read {what} from {token!r}", lineno)
        return math.nan

    def integer(self, token: str, lineno: int, what: str) -> int:
        """An integer; 0 and an S06 issue if Fortran cannot read it."""
        if _INT.fullmatch(token):
            return int(token)
        self.issue("S06", f"cannot read {what} from {token!r}", lineno)
        return 0

    def count(self, what: str) -> tuple[int, int, str]:
        """A count line: (count, line number, text).

        The standard puts a comment line before every count; RADEX and ProDiMo skip it by
        position, whatever it holds. This reader skips comment lines by their ``!``, so a
        comment line without ``!`` would be taken for the count. A line that starts with
        a letter where a count is expected is therefore skipped as such a comment (S05).
        """
        lineno, text = self.line(what)
        if text.lstrip()[:1].isalpha():
            self.issue("S05", f"comment line without '!' before {what}: {text.strip()!r}", lineno, fixable=True)
            lineno, text = self.line(what)
        token = text.split()[0]
        if not _INT.fullmatch(token) or int(token) < 0:
            self.fail("S06", f"cannot read {what} from {text.strip()!r}", lineno)
        return int(token), lineno, text

    def rest(self) -> tuple[int | None, list[str]]:
        """The first line number after the last line read, and all lines from there (reads the file to its end)."""
        self.raw.extend(self._lines)
        if any(len(text) > MAX_LINE_LENGTH for text in self.raw[-1:]):
            self.not_lamda(f"longer than {MAX_LINE_LENGTH} characters", len(self.raw))
        last = self.data[self.pos - 1][0] if self.pos else 0
        return (last + 1 if last < len(self.raw) else None), self.raw[last:]


def _numbering(reader: _Reader, numbers: list[tuple[int, int]], what: str) -> None:
    """One S02 issue if the row numbers are not 1..N in order: ``numbers`` = (number, line)."""
    bad = [(k, n, line) for k, (n, line) in enumerate(numbers, 1) if n != k]
    if bad:
        k, n, line = bad[0]
        more = f" ({len(bad) - 1} more rows)" if len(bad) > 1 else ""
        reader.issue("S02", f"{what}: row {k} has the number {n}{more}", line, fixable=True)


def read(path: str | os.PathLike, strict: bool = False) -> LamdaData:
    """Read a LAMDA file (LAMDA, EMAA, ProDiMo or pylamda's own layout).

    Parameters
    ----------
    path : str or os.PathLike
        The file.
    strict : bool, optional
        Raise :class:`~pylamda.issues.LamdaReadError` if the reader finds any structure
        issue (for CI). Default: False, the issues are kept in
        :attr:`LamdaData.issues <pylamda.model.LamdaData.issues>`.

    Returns
    -------
    LamdaData
        The data, with a :class:`~pylamda.model.SourceMap` (line numbers, raw lines and
        the issues found) in ``source_map``.

    Raises
    ------
    NotLamdaError
        If the header does not read (molecule name, molecular weight, number of levels,
        level rows, number of radiative transitions): the file is not a LAMDA file. The
        rest of the file is not read. A subclass of LamdaReadError.
    LamdaReadError
        If the file ends before its counts are satisfied (S01), a count cannot be read
        (S06), or ``strict`` is True and there is any structure issue. A subclass of
        ValueError.
    OSError
        If the file cannot be opened.
    """
    path = Path(path)
    with path.open(encoding="utf-8", errors="replace") as f:
        return _parse(_file_lines(f), str(path), strict)


def read_string(text: str, name: str = "<string>", strict: bool = False) -> LamdaData:
    """Read LAMDA data from a string (as :func:`read` reads a file).

    Parameters
    ----------
    text : str
        The contents of a LAMDA file.
    name : str, optional
        Name used as the path in messages and the source map. Default: ``"<string>"``.
    strict : bool, optional
        As in :func:`read`.

    Returns
    -------
    LamdaData
        The data, with a source map of the text.

    Raises
    ------
    NotLamdaError, LamdaReadError
        As :func:`read`.
    """
    return _parse(text.splitlines(), name, strict)


def _parse(lines: Iterable[str], path: str, strict: bool) -> LamdaData:
    """Read the lines of a LAMDA file (see :func:`read`)."""
    reader = _Reader(path, lines)
    smap = SourceMap(path=path, raw=reader.raw)  # the same list, complete after reader.rest()

    smap.molecule, name_line = reader.line("the molecule name")
    m = _DATASET_DOI.search(name_line)
    dataset_doi = m.group(1) if m else None
    if m:
        name_line = name_line[: m.start()]
    molecule, _, title = name_line.strip().partition(" ")
    if _NUMBER.fullmatch(molecule):
        reader.not_lamda(f"the molecule name {molecule!r} is a number", smap.molecule)
    smap.mass, text = reader.line("the molecular weight")
    mass = reader.real(text.split()[0], smap.mass, "the molecular weight")
    if math.isnan(mass):
        reader.not_lamda(f"cannot read the molecular weight from {text.split()[0]!r}", smap.mass)
    if not (math.isfinite(mass) and mass > 0):
        reader.not_lamda(f"the molecular weight {mass:g} is not > 0", smap.mass)

    # levels
    nlev, smap.nlev, text = reader.count("the number of levels")
    if nlev < 1:
        reader.not_lamda("the number of levels is 0", smap.nlev)
    level_type, level_reference = "", ""
    after = text.strip().split(None, 1)[1] if len(text.split()) > 1 else ""
    if after:  # EMAA: "14 - ROTATION - Ref: ..."
        head, _, ref = after.partition("Ref:")
        level_reference = ref.strip()
        level_type = head.strip(" -")
    energy = np.full(nlev, np.nan)
    g = np.full(nlev, np.nan)
    quanta: list[str] = []
    smap.levels = np.zeros(nlev, dtype=int)
    numbers = []
    for i in range(nlev):
        lineno, text = reader.line(f"level {i + 1}")
        smap.levels[i] = lineno
        parts = text.split(None, 3)
        if len(parts) < 3:
            reader.issue("S07", f"level row {i + 1} needs i, E, g: {text.strip()!r}", lineno)
            parts = [*parts, "", "", ""][:3]
        numbers.append((reader.integer(parts[0], lineno, "a level number"), lineno))
        if parts[1]:
            energy[i] = reader.real(parts[1], lineno, f"the energy of level {i + 1}")
        if parts[2]:
            g[i] = reader.real(parts[2], lineno, f"the weight of level {i + 1}")
        quanta.append(parts[3].strip() if len(parts) > 3 else "")
    _numbering(reader, numbers, "levels")

    # radiative transitions
    nlin, smap.nlin, _ = reader.count("the number of radiative transitions")
    reader.header = False  # the header has been read: a LAMDA file from here on
    up = np.zeros(nlin, dtype=int)
    low = np.zeros(nlin, dtype=int)
    a = np.full(nlin, np.nan)
    freq = np.full(nlin, np.nan)
    eup = np.full(nlin, np.nan)
    smap.lines = np.zeros(nlin, dtype=int)
    numbers = []
    for i in range(nlin):
        lineno, text = reader.line(f"radiative transition {i + 1}")
        smap.lines[i] = lineno
        parts = text.split()
        if len(parts) < 4:
            reader.issue("S07", f"radiative transition {i + 1} needs i, up, low, A: {text.strip()!r}", lineno)
            parts = [*parts, "", "", "", ""][:4]
        what = f"radiative transition {i + 1}"
        numbers.append((reader.integer(parts[0], lineno, f"the number of {what}"), lineno))
        if parts[1]:
            up[i] = reader.integer(parts[1], lineno, f"the upper level of {what}")
        if parts[2]:
            low[i] = reader.integer(parts[2], lineno, f"the lower level of {what}")
        if parts[3]:
            a[i] = reader.real(parts[3], lineno, f"the Einstein A of {what}")
        if len(parts) > 4:
            freq[i] = reader.real(parts[4], lineno, f"the frequency of {what}")
        if len(parts) > 5:
            eup[i] = reader.real(parts[5], lineno, f"the upper-level energy of {what}")
    _numbering(reader, numbers, "radiative transitions")

    # collision blocks
    npart, smap.npart, _ = reader.count("the number of collision partners")
    collisions = []
    for b in range(1, npart + 1):
        lineno_p, pline = reader.line(f"the partner line of collision block {b}")
        token = pline.split()[0]
        if _INT.fullmatch(token):
            pid = int(token)
        else:
            reader.issue(
                "S06",
                f"the partner line must start with the partner ID (1=H2 2=p-H2 3=o-H2 4=e 5=H 6=He 7=H+): "
                f"{pline.strip()!r}",
                lineno_p,
            )
            pid = 0
        _, reference = split_partner_line(pline)
        label = f"collision block {b}"
        ncol, lineno_ncol, _ = reader.count(f"the number of collisional transitions of {label}")
        ntemp, lineno_ntemp, _ = reader.count(f"the number of temperatures of {label}")
        temps, lineno_t, surplus = reader.values(ntemp, f"the temperatures of {label}")
        if surplus:
            reader.issue(
                "S04",
                f"{label}: more than the {ntemp} declared temperatures; the rest is ignored",
                reader.data[reader.pos - 1][0],
                fixable=True,
            )
        temperatures = np.array([reader.real(t, lineno_t, f"a temperature of {label}") for t in temps])
        cup = np.zeros(ncol, dtype=int)
        clow = np.zeros(ncol, dtype=int)
        rates = np.full((ncol, ntemp), np.nan)
        rows = np.zeros(ncol, dtype=int)
        numbers = []
        first_surplus, n_surplus = 0, 0
        for k in range(ncol):
            vals, lineno, surplus = reader.values(3 + ntemp, f"collisional transition {k + 1} of {label}")
            rows[k] = lineno
            if surplus:
                n_surplus += 1
                first_surplus = first_surplus or reader.data[reader.pos - 1][0]
            what = f"collisional transition {k + 1} of {label}"
            numbers.append((reader.integer(vals[0], lineno, f"the number of {what}"), lineno))
            cup[k] = reader.integer(vals[1], lineno, f"the upper level of {what}")
            clow[k] = reader.integer(vals[2], lineno, f"the lower level of {what}")
            rates[k] = [reader.real(v, lineno, f"a rate of {what}") for v in vals[3:]]
        _numbering(reader, numbers, label)
        if n_surplus:
            reader.issue(
                "S04",
                f"{label}: {n_surplus} rate rows have more than {ntemp} values; the rest is ignored",
                first_surplus,
                fixable=True,
            )
        smap.blocks.append(BlockSource(lineno_p, lineno_ncol, lineno_ntemp, lineno_t, rows))
        collisions.append(CollisionSet(pid, temperatures, cup, clow, rates, reference, pline.rstrip()))

    # notes
    smap.notes, rest = reader.rest()
    notes = [line.rstrip() for line in rest]
    while notes and not notes[-1].strip():
        notes.pop()
    data_lines = [k for k, line in enumerate(notes) if line.strip() and not line.lstrip().startswith("!")]
    if data_lines:
        reader.issue(
            "S03",
            f"{len(data_lines)} data lines after the last of the {npart} declared collision blocks "
            "(ignored, kept as notes)",
            (smap.notes or 0) + data_lines[0],
            fixable=True,
        )
    source = ""
    for i, line in enumerate(notes):
        if line.startswith("!SOURCE:"):
            source = line[len("!SOURCE:") :].strip()
            del notes[i]
            break
    if notes and notes[0].strip() == "!NOTES":
        notes = notes[1:]

    smap.issues = reader.issues
    if strict and reader.issues:
        raise LamdaReadError(str(path), reader.issues)
    return LamdaData(
        molecule,
        mass,
        energy,
        g,
        quanta,
        up,
        low,
        a,
        freq,
        eup,
        collisions,
        title=title.strip(),
        dataset_doi=dataset_doi,
        level_type=level_type,
        level_reference=level_reference,
        source=source,
        notes=notes,
        source_map=smap,
    )
