"""Issues: deviations from the LAMDA standard found in a file or a data set.

The reader (:func:`pylamda.read`) collects the structure issues (codes ``S..``) it finds
while reading; :func:`pylamda.verify` adds the checks of the data (``L``, ``R``, ``C``,
``F``, ``P``). Each issue has a code, a severity and a message, as in ruff. The codes,
their default severities and short descriptions are in :data:`CODES`.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import asdict, dataclass

#: Severities, from the most to the least serious.
SEVERITIES: tuple[str, ...] = ("error", "warning", "info")

#: Check code -> (default severity, short description).
CODES: dict[str, tuple[str, str]] = {
    # structure (found by the reader)
    "S01": ("error", "file shorter than its counts require"),
    "S02": ("error", "level or transition numbers not 1..N in order"),
    "S03": ("warning", "data after the last declared block (more blocks than NPART)"),
    "S04": ("warning", "more values on a line than declared (NTEMP)"),
    "S05": ("error", "text line without '!' between blocks (read as data by Fortran)"),
    "S06": ("error", "number that Fortran cannot read, or a partner line without the ID"),
    "S07": ("error", "level or radiative-transition row with fewer values than required"),
    # levels
    "L01": ("warning", "energies not in increasing order"),
    "L02": ("error", "weight g <= 0 or not finite (error), or not an integer (info)"),
    "L03": ("warning", "two levels with the same E, g and quanta"),
    "L04": ("info", "quanta empty or not unique"),
    "L05": ("info", "lowest level energy is not 0"),
    # radiative transitions
    "R01": ("error", "upper or lower level outside 1..NLEV, or upper = lower"),
    "R02": ("error", "E(up) <= E(low)"),
    "R03": ("error", "Einstein A <= 0 or not finite"),
    "R04": ("warning", "frequency missing or not (E_up - E_low) * c"),
    "R05": ("warning", "E_up [K] missing or not E(up) * hc/k"),
    "R06": ("error", "the same (up, low) twice"),
    # collisions
    "C01": ("error", "partner ID not in 1..7"),
    "C02": ("warning", "the partner text names another standard partner than the ID"),
    "C03": ("error", "the same partner ID twice"),
    "C04": ("warning", "H2 together with p-H2 or o-H2 (counted twice by some codes)"),
    "C05": ("error", "temperatures not strictly increasing, or <= 0"),
    "C06": ("error", "rate < 0 or not finite"),
    "C07": ("info", "rate row that is 0 at all temperatures"),
    "C08": ("warning", "E(up) < E(low) in a rate row (wrong order or an excitation rate)"),
    "C09": ("error", "the same (up, low) twice in one block"),
    "C10": ("error", "upper or lower level outside 1..NLEV, or upper = lower"),
    "C11": ("warning", "the partner text names a partner outside the standard (read by its ID)"),
    "C12": ("warning", "the partner text names He, the ID is H2, and the text says scaled (He rates scaled to H2)"),
    "C13": ("info", "partner text says the rates are scaled (codes that scale may scale twice)"),
    # references
    "F01": ("warning", "no reference for the levels and lines"),
    "F02": ("warning", "partner without a reference"),
    "F03": ("info", "reference without a DOI, bibcode or link (free text only)"),
    # physics heuristics
    "P01": ("info", "rate row that behaves like an excitation rate"),
    "P03": ("info", "rate outside the plausible range"),
    "P04": ("info", "mass differs from the formula mass of the molecule"),
}


@dataclass(frozen=True)
class Issue:
    """One deviation from the LAMDA standard.

    Attributes
    ----------
    code : str
        Check code (:data:`CODES`), e.g. ``"S04"``.
    severity : str
        ``"error"``, ``"warning"`` or ``"info"``.
    message : str
        What was found, in one line.
    line : int or None
        Line number in the file (1-based), or None if the issue has no single line. For
        an issue about several items, the line of the first one.
    fixable : bool
        Whether ``pylamda fix`` can remove the issue.
    count : int
        Number of items affected (rows, levels, ...): checks report one issue per
        block or section, not one per row.
    """

    code: str
    severity: str
    message: str
    line: int | None = None
    fixable: bool = False
    count: int = 1

    def __str__(self) -> str:
        where = f"line {self.line}: " if self.line is not None else ""
        return f"{where}{self.code} {self.severity}: {self.message}"

    def as_dict(self) -> dict:
        """The issue as a JSON-compatible dict.

        Returns
        -------
        dict
            The fields of the dataclass.
        """
        return asdict(self)


def make_issue(
    code: str,
    message: str,
    line: int | None = None,
    fixable: bool = False,
    count: int = 1,
    severity: str | None = None,
) -> Issue:
    """An :class:`Issue` with the default severity of ``code``.

    Parameters
    ----------
    code : str
        Check code, a key of :data:`CODES`.
    message : str
        What was found.
    line : int or None, optional
        Line number in the file (1-based).
    fixable : bool, optional
        Whether a fix exists. Default: False.
    count : int, optional
        Number of items affected. Default: 1.
    severity : str or None, optional
        Severity if not the default of ``code``.

    Returns
    -------
    Issue
        The issue.

    Raises
    ------
    KeyError
        If ``code`` is not in :data:`CODES`.
    """
    default = CODES[code][0]
    return Issue(code, severity or default, message, None if line is None else int(line), fixable, int(count))


def match_codes(code: str, selectors: Iterable[str]) -> bool:
    """Whether ``code`` matches one of ``selectors`` (code prefixes, as in ruff).

    Parameters
    ----------
    code : str
        A check code, e.g. ``"C02"``.
    selectors : iterable of str
        Code prefixes: ``"C"`` (all C codes), ``"C0"`` (C01-C09) or ``"C02"``; case is
        ignored.

    Returns
    -------
    bool
        True if a selector is a prefix of ``code``.
    """
    return any(code.upper().startswith(s.strip().upper()) for s in selectors if s.strip())


class LamdaReadError(ValueError):
    """A file cannot be read as a LAMDA file, or (with ``strict=True``) has structure issues.

    Attributes
    ----------
    path : str
        The file.
    issues : list of Issue
        The issues found; the last one stopped the reader unless ``strict=True`` was given.
    """

    def __init__(self, path: str, issues: list[Issue]):
        self.path = path
        self.issues = list(issues)
        lines = [f"{path}: {issue}" for issue in self.issues]
        super().__init__("\n".join(lines) if lines else f"{path}: cannot be read")


class NotLamdaError(LamdaReadError):
    """The header of a file does not read as a LAMDA header: the file is not a LAMDA file.

    The header is the molecule name, the molecular weight, the number of levels, the
    level rows and the number of radiative transitions. The reader stops at the first
    of these that fails, so it does not read the rest of a large file. A subclass of
    :class:`LamdaReadError`.

    Attributes
    ----------
    path : str
        The file.
    issues : list of Issue
        The issues found up to that point (may be empty).
    reason : str
        Why the file is not a LAMDA file, in one line.
    """

    def __init__(self, path: str, issues: list[Issue], reason: str):
        self.path = path
        self.issues = list(issues)
        self.reason = reason
        lines = [f"{path}: not a LAMDA file: {reason}", *(f"{path}: {issue}" for issue in self.issues)]
        ValueError.__init__(self, "\n".join(lines))
