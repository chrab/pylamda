"""Issues: deviations from the LAMDA standard found in a file or a data set.

The reader (:func:`pylamda.read`) collects the structure issues (codes ``S..``) it finds
while reading. Each issue has a code, a severity and a message, as in ruff. The codes
and their default severities are in :data:`CODES`; the design document lists all
planned checks (section 6.1).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

#: Severities, from the most to the least serious.
SEVERITIES: tuple[str, ...] = ("error", "warning", "info")

#: Check code -> (default severity, short description).
CODES: dict[str, tuple[str, str]] = {
    "S01": ("error", "file shorter than its counts require"),
    "S02": ("error", "level or transition numbers not 1..N in order"),
    "S03": ("warning", "data after the last declared block (more blocks than NPART)"),
    "S04": ("warning", "more values on a line than declared (NTEMP)"),
    "S05": ("error", "text line without '!' between blocks (read as data by Fortran)"),
    "S06": ("error", "number that Fortran cannot read"),
    "S07": ("error", "level or radiative-transition row with fewer values than required"),
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
        Line number in the file (1-based), or None if the issue has no single line.
    fixable : bool
        Whether ``pylamda fix`` can remove the issue.
    """

    code: str
    severity: str
    message: str
    line: int | None = None
    fixable: bool = False

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


def make_issue(code: str, message: str, line: int | None = None, fixable: bool = False) -> Issue:
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

    Returns
    -------
    Issue
        The issue.

    Raises
    ------
    KeyError
        If ``code`` is not in :data:`CODES`.
    """
    return Issue(code, CODES[code][0], message, line, fixable)


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
