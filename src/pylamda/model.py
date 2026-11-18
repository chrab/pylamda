"""The data model: one LAMDA file as :class:`LamdaData` with one :class:`CollisionSet` per partner.

The format is that of the Leiden Atomic and Molecular Database (LAMDA,
https://home.strw.leidenuniv.nl/~moldata/molformat.html), also used by EMAA
(https://emaa.osug.fr/file_format). Levels and lines are stored as one numpy array per
column, with 1-based level numbers as in the file.

Collision partners are identified by their ID only (:data:`PARTNERS`), as RADEX does.
The name in the free text of the partner line is kept (:attr:`CollisionSet.partner_text`)
so that checks can compare it with the ID.
"""

from __future__ import annotations

import copy
import re
from collections.abc import Iterable
from dataclasses import dataclass, field, replace

import numpy as np

from .issues import Issue
from .references import find_dois

#: The collision partners of the LAMDA standard: ID -> name.
PARTNERS: dict[int, str] = {1: "H2", 2: "p-H2", 3: "o-H2", 4: "e", 5: "H", 6: "He", 7: "H+"}
#: Name -> ID of the standard partners.
PARTNER_IDS: dict[str, int] = {name: pid for pid, name in PARTNERS.items()}
#: The standard partners, in ID order.
STANDARD_PARTNERS: tuple[str, ...] = tuple(PARTNERS.values())

# spellings found in LAMDA/EMAA/ProDiMo files -> canonical name (keys case-folded).
# p-H2O, o-H2O and CO are EMAA partners outside the standard; they are recognised so
# that checks can report them (they reuse the IDs 2, 3 and 5).
_PARTNER_ALIASES = {
    "h2": "H2",
    "p-h2": "p-H2",
    "ph2": "p-H2",
    "para-h2": "p-H2",
    "p_h2": "p-H2",
    "para h2": "p-H2",
    "o-h2": "o-H2",
    "oh2": "o-H2",
    "ortho-h2": "o-H2",
    "o_h2": "o-H2",
    "ortho h2": "o-H2",
    "e": "e",
    "e-": "e",
    "elec": "e",
    "electron": "e",
    "electrons": "e",
    "h": "H",
    "he": "He",
    "h+": "H+",
    "p-h2o": "p-H2O",
    "para-h2o": "p-H2O",
    "o-h2o": "o-H2O",
    "ortho-h2o": "o-H2O",
    "co": "CO",
}

# " ! " between partner and reference on a ProDiMo-style partner line; also a final " !"
_COMMENT_MARK = re.compile(r"\s!(?=\s|$)")
_REF = re.compile(r"\s*-*\s*Ref:\s*", re.I)


def partner_name(text: str) -> str | None:
    """Canonical partner name of a spelling (``"para-H2"`` -> ``"p-H2"``).

    Parameters
    ----------
    text : str
        A partner name as written in LAMDA or EMAA files (``"pH2"``, ``"para-H2"``,
        ``"electron"``, ``"H+"``, ...); case and surrounding blanks are ignored.

    Returns
    -------
    str or None
        The canonical name (a :data:`STANDARD_PARTNERS` name, or ``"p-H2O"``,
        ``"o-H2O"``, ``"CO"`` for EMAA's partners outside the standard), or None if the
        spelling is unknown.
    """
    return _PARTNER_ALIASES.get(text.strip().casefold())


def split_partner_line(line: str) -> tuple[str | None, str]:
    """The partner name and the reference in the free text of a partner line.

    The text after the ID is free in the LAMDA standard. Two forms name the partner:

    * EMAA and pylamda: ``<id> <molecule> + <partner> -- Ref: <reference>``;
    * ProDiMo: ``<id> <molecule> + <partner> ! <reference>`` (the ``!`` part may be
      missing); only recognised if ``<partner>`` is a known spelling (:func:`partner_name`).

    Parameters
    ----------
    line : str
        The partner line, with the ID.

    Returns
    -------
    name : str or None
        The partner name as written (not normalised), or None if the text does not name
        the partner in one of these forms.
    reference : str
        The reference text: after ``Ref:`` or ``!`` in these forms, else all text after
        the ID.
    """
    parts = line.strip().split(None, 1)
    rest = parts[1].strip() if len(parts) > 1 else ""
    left, sep, right = rest.partition("--")
    if sep and _REF.match(right):
        reference = _REF.sub("", right, count=1).strip()
        name = left.rsplit(" + ", 1)[1].strip() if " + " in left else None
        return name or None, reference
    if " + " in rest:
        candidate, after = ([*_COMMENT_MARK.split(rest.split(" + ", 1)[1], maxsplit=1), ""])[:2]
        if partner_name(candidate):
            return candidate.strip(), after.strip()
    return None, rest


@dataclass
class CollisionSet:
    """Collision rate coefficients for one collision partner.

    Attributes
    ----------
    partner_id : int
        Partner ID as written in the file (1-7 in the standard, :data:`PARTNERS`).
    temperatures : numpy.ndarray
        Kinetic temperatures [K], shape (ntemp,).
    up, low : numpy.ndarray
        Upper and lower level numbers (1-based, as in the file), shape (ncol,).
    rates : numpy.ndarray
        Downward rate coefficients [cm3 s-1], shape (ncol, ntemp).
    reference : str
        Reference text of the partner line (without the ID and the partner name where
        these could be separated, see :func:`split_partner_line`).
    partner_text : str
        The partner line as written in the file (``""`` if the data do not come from a
        file).
    """

    partner_id: int
    temperatures: np.ndarray
    up: np.ndarray
    low: np.ndarray
    rates: np.ndarray
    reference: str = ""
    partner_text: str = ""

    @property
    def partner(self) -> str | None:
        """str or None: Standard name of :attr:`partner_id` (None if the ID is not 1-7)."""
        return PARTNERS.get(self.partner_id)

    @property
    def text_partner(self) -> str | None:
        """str or None: The partner name in :attr:`partner_text`, as written (None if none)."""
        return split_partner_line(self.partner_text)[0] if self.partner_text else None

    @property
    def label(self) -> str:
        """str: The partner name, or ``"ID <n>"`` for an ID outside the standard."""
        return self.partner or f"ID {self.partner_id}"

    @property
    def n_transitions(self) -> int:
        """int: Number of collisional transitions."""
        return len(self.up)

    def rate(self, temperature: float, clamp: bool = False) -> np.ndarray:
        """Downward rate coefficients at ``temperature``, linear in T between grid points.

        Parameters
        ----------
        temperature : float
            Kinetic temperature [K].
        clamp : bool, optional
            Outside the grid, use the rates at the nearest grid temperature (as RADEX
            does). Default: False, raise ValueError.

        Returns
        -------
        numpy.ndarray
            Rate coefficients [cm3 s-1], shape (ncol,), in the order of ``up``/``low``.

        Raises
        ------
        ValueError
            If ``temperature`` is outside the grid and ``clamp`` is False.
        """
        t = self.temperatures
        if not clamp and not t[0] <= temperature <= t[-1]:
            raise ValueError(
                f"T = {temperature:g} K outside the {self.label} grid {t[0]:g}-{t[-1]:g} K "
                "(use clamp=True to take the nearest grid values)"
            )
        if len(t) == 1:
            return self.rates[:, 0].copy()
        temperature = min(max(temperature, t[0]), t[-1])
        i = int(np.clip(np.searchsorted(t, temperature), 1, len(t) - 1))
        w = (temperature - t[i - 1]) / (t[i] - t[i - 1])
        return (1 - w) * self.rates[:, i - 1] + w * self.rates[:, i]


@dataclass
class BlockSource:
    """Line numbers of one collision block in the file (see :class:`SourceMap`).

    Attributes
    ----------
    partner : int
        The partner line.
    ncol : int
        The line with the number of collisional transitions.
    ntemp : int
        The line with the number of temperatures.
    temperatures : int
        The first line of the temperatures.
    rows : numpy.ndarray
        The first line of each rate row, shape (ncol,).
    """

    partner: int
    ncol: int
    ntemp: int
    temperatures: int
    rows: np.ndarray


@dataclass
class SourceMap:
    """Where each item of a :class:`LamdaData` is in the file it was read from.

    All line numbers are 1-based. A value that continues on the next lines (temperatures,
    rate rows) has the number of its first line.

    Attributes
    ----------
    path : str
        The file.
    raw : list of str
        All lines of the file, without line endings.
    molecule, mass, nlev, nlin, npart : int
        The lines with the molecule name, the mass and the three counts.
    levels : numpy.ndarray
        The line of each level row, shape (nlev,).
    lines : numpy.ndarray
        The line of each radiative-transition row, shape (nlin,).
    blocks : list of BlockSource
        One entry per collision block, in file order.
    notes : int or None
        The first line after the last declared block (None if the file ends there).
    issues : list of Issue
        The structure issues the reader found.
    """

    path: str
    raw: list[str]
    molecule: int = 0
    mass: int = 0
    nlev: int = 0
    levels: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    nlin: int = 0
    lines: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    npart: int = 0
    blocks: list[BlockSource] = field(default_factory=list)
    notes: int | None = None
    issues: list[Issue] = field(default_factory=list)


@dataclass
class LamdaData:
    """Contents of one LAMDA file.

    Levels and lines are stored as one array per column (level ``i`` of the file is
    index ``i - 1``). Created by :func:`pylamda.read`, written by :func:`pylamda.write`.

    Attributes
    ----------
    molecule : str
        Molecule name (first word of line 2), e.g. ``"CO"``, ``"o-H2O"``.
    mass : float
        Molecular weight [amu].
    level_energy : numpy.ndarray
        Level energies [cm-1].
    level_g : numpy.ndarray
        Statistical weights.
    level_quanta : list of str
        Quantum-number text after E and g (``""`` if none).
    line_up, line_low : numpy.ndarray
        Upper and lower level numbers (1-based) of the radiative transitions.
    line_a : numpy.ndarray
        Einstein A coefficients [s-1].
    line_freq : numpy.ndarray
        Frequencies [GHz]; NaN where the file gives none.
    line_eup : numpy.ndarray
        Upper-level energies [K]; NaN where the file gives none.
    collisions : list of CollisionSet
        One entry per collision block, in file order.
    title : str
        Rest of line 2 after the molecule name (description), without a dataset DOI.
    dataset_doi : str or None
        Dataset DOI (EMAA: ``-- EMAA_DATASET_DOI: ...`` on line 2).
    level_type : str
        Type of the level set given on the NLEV line (EMAA: ``"ROTATION"``, ...).
    level_reference : str
        Reference of the levels and radiative data (EMAA: ``Ref:`` on the NLEV line).
    source : str
        Where the data come from (written as ``!SOURCE:`` comment line).
    notes : list of str
        Lines after the last declared block (comments, notes, surplus data), unchanged.
    source_map : SourceMap or None
        Line numbers and raw lines of the file, and the reader's issues; None if the
        data were not read from a file. Not updated when the data change.
    """

    molecule: str
    mass: float
    level_energy: np.ndarray
    level_g: np.ndarray
    level_quanta: list[str]
    line_up: np.ndarray
    line_low: np.ndarray
    line_a: np.ndarray
    line_freq: np.ndarray
    line_eup: np.ndarray
    collisions: list[CollisionSet] = field(default_factory=list)
    title: str = ""
    dataset_doi: str | None = None
    level_type: str = ""
    level_reference: str = ""
    source: str = ""
    notes: list[str] = field(default_factory=list)
    source_map: SourceMap | None = field(default=None, repr=False, compare=False)

    @property
    def n_levels(self) -> int:
        """int: Number of energy levels."""
        return len(self.level_energy)

    @property
    def n_lines(self) -> int:
        """int: Number of radiative transitions."""
        return len(self.line_a)

    @property
    def partners(self) -> list[str | None]:
        """list of str or None: Collision partner names by ID, in file order (None for an ID outside 1-7)."""
        return [c.partner for c in self.collisions]

    @property
    def issues(self) -> list[Issue]:
        """list of Issue: The issues the reader found (empty if not read from a file)."""
        return list(self.source_map.issues) if self.source_map else []

    def partner(self, name: str) -> CollisionSet:
        """The collision data of partner ``name`` (any known spelling, e.g. ``"pH2"``).

        Partners are identified by their ID: ``"p-H2"`` finds the block with ID 2. If
        the file has several blocks with that ID, the first one is returned.

        Parameters
        ----------
        name : str
            Partner name (:func:`partner_name` spellings).

        Returns
        -------
        CollisionSet
            The data of that partner.

        Raises
        ------
        KeyError
            If the data has no such partner.
        """
        key = partner_name(name) or name
        for c in self.collisions:
            if c.partner == key:
                return c
        have = ", ".join(c.label for c in self.collisions)
        raise KeyError(f"{self.molecule}: no collision partner {name!r} (have: {have})")

    def select_partners(self, names: Iterable[str]) -> LamdaData:
        """A copy with only some collision partners.

        Parameters
        ----------
        names : iterable of str
            Partner names (any spelling of :func:`partner_name`); names the data does
            not have are ignored.

        Returns
        -------
        LamdaData
            A shallow copy with the collision sets of these partners, in file order.
        """
        keep = {partner_name(n) or n for n in names}
        return replace(self, collisions=[c for c in self.collisions if c.partner in keep])

    def references(self) -> list[dict]:
        """References in the data.

        Returns
        -------
        list of dict
            One dict per reference with ``role``: ``"dataset"`` (dataset DOI),
            ``"levels"`` (levels and radiative data) or ``"collisions"`` (one entry per
            partner with a reference text); ``partner`` (name or None); ``text`` (the
            reference text); ``dois`` (list of DOIs found in it).
        """
        refs = []
        if self.dataset_doi:
            refs.append({"role": "dataset", "partner": None, "text": "", "dois": [self.dataset_doi]})
        if self.level_reference:
            refs.append(
                {
                    "role": "levels",
                    "partner": None,
                    "text": self.level_reference,
                    "dois": find_dois(self.level_reference),
                }
            )
        for c in self.collisions:
            if c.reference:
                refs.append(
                    {"role": "collisions", "partner": c.partner, "text": c.reference, "dois": find_dois(c.reference)}
                )
        return refs

    def copy(self) -> LamdaData:
        """A deep copy (arrays, collision sets, notes and source map).

        Returns
        -------
        LamdaData
            The copy.
        """
        return copy.deepcopy(self)
