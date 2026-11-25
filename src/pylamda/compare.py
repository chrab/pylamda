"""Compare two LAMDA data sets: levels, lines, collision rates and critical densities.

Merged from astrolinedata's ``compare_lamda`` and ``compare_rates`` of TEST_LAMDA_FILES
(design document, section 8):

* Levels are matched first (:func:`match_levels`); lines and rate rows are matched
  through the level map, so that inserted levels do not show up as changes of all later
  entries.
* Rates are compared on the common temperatures exactly (real changes against rounding,
  :attr:`PartnerDiff.n_rates_changed`) and on the union grid of both tables in the
  common temperature range, interpolated linearly in log T and log k (trends also when
  the grids differ, :attr:`PartnerDiff.rates_interpolated`).
* Per upper level and partner, the total de-excitation rate sum_l k_ul(T) is
  compared. It uses all rows of each data set, so added or removed transitions count
  too. (The critical density, :func:`critical_density`, is not compared: it changes
  only with the rates or the A values, which are compared themselves.)
* Changes are also summarised per E_up range (:data:`EUP_BINS`), so that a change of
  the high levels does not hide behind many unchanged low levels.
"""

from __future__ import annotations

import math
import os
import re
from dataclasses import asdict, dataclass, field

import numpy as np

from .model import CollisionSet, LamdaData
from .reader import read
from .verify import HC_K

#: Upper-level energy ranges [K] of the summaries per E_up range (bin edges).
EUP_BINS: tuple[float, ...] = (0.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0, math.inf)

# rates below this are treated as 0 in ratios and logarithms [cm3 s-1]
_RATE_FLOOR = 1e-90


# ------------------------------------------------------------------- level matching
@dataclass
class LevelMatch:
    """How the levels of two data sets correspond (result of :func:`match_levels`).

    Attributes
    ----------
    mapping : dict of int to int
        Level number in ``a`` -> level number in ``b`` (1-based, as in the files).
    method : str
        ``"quanta"``, ``"quanta+energy"`` (by quanta, the rest by energy and g),
        ``"energy"`` or ``"index"``.
    max_energy_difference : float
        Largest energy difference of matched levels [cm-1] (0 if none).
    """

    mapping: dict[int, int]
    method: str
    max_energy_difference: float = 0.0

    @property
    def n_matched(self) -> int:
        """int: Number of matched levels."""
        return len(self.mapping)


def normalise_quanta(text: str) -> str:
    """Quantum-number text for matching: lower case; blanks, ``_``, ``,`` and ``;`` as one ``_``.

    Parameters
    ----------
    text : str
        Quantum-number text of a level (``"J=1 F=2"``, ``"1_1_0"``, ...).

    Returns
    -------
    str
        The normalised text (``"j=1_f=2"``, ``"1_1_0"``).
    """
    return re.sub(r"[\s_,;]+", "_", text.strip().lower()).strip("_")


def _match_energy(
    a: LamdaData, b: LamdaData, atol: float, rtol: float, skip_a: set[int], used_b: set[int]
) -> dict[int, int]:
    """Match levels (0-based) by energy within ``atol + rtol * E`` and equal g, closest first."""
    out: dict[int, int] = {}
    used = set(used_b)
    order = np.argsort(b.level_energy, kind="stable")
    eb = b.level_energy[order]
    for i in range(a.n_levels):
        if i in skip_a:
            continue
        e = a.level_energy[i]
        tol = atol + rtol * abs(e)
        lo, hi = np.searchsorted(eb, e - tol), np.searchsorted(eb, e + tol, side="right")
        cands = [
            int(order[k]) for k in range(lo, hi) if int(order[k]) not in used and b.level_g[order[k]] == a.level_g[i]
        ]
        if cands:
            j = min(cands, key=lambda k: abs(b.level_energy[k] - e))
            out[i] = j
            used.add(j)
    return out


def match_levels(a: LamdaData, b: LamdaData, match: str = "auto", atol: float = 1e-3, rtol: float = 1e-5) -> LevelMatch:
    """Match the levels of ``a`` to those of ``b``.

    Parameters
    ----------
    a, b : LamdaData
        The two data sets.
    match : str, optional
        How levels are matched:

        * ``"quanta"``: by the quantum-number text (:func:`normalise_quanta`), which must
          be present and unique in both;
        * ``"energy"``: by energy (within ``atol + rtol * E``) and equal weight g, the
          closest level first;
        * ``"index"``: by level number;
        * ``"auto"`` (default): by quanta if they are unique in both and at least half
          of the levels match that way (same notation); the levels left over are then
          matched by energy and g (method ``"quanta+energy"``). Otherwise by energy.
    atol : float, optional
        Energy tolerance [cm-1] for matching by energy. Default: 1e-3.
    rtol : float, optional
        Relative energy tolerance for matching by energy. Default: 1e-5.

    Returns
    -------
    LevelMatch
        The level map (1-based) and the method used.

    Raises
    ------
    ValueError
        If ``match`` is unknown, or ``"quanta"`` and the quanta are missing or not unique.
    """
    if match not in ("auto", "quanta", "energy", "index"):
        raise ValueError(f"match must be 'auto', 'quanta', 'energy' or 'index', not {match!r}")
    mapping: dict[int, int] | None = None
    method = ""
    if match in ("auto", "quanta"):
        ka = [normalise_quanta(q) for q in a.level_quanta]
        kb = [normalise_quanta(q) for q in b.level_quanta]
        unique = all(ka) and all(kb) and len(set(ka)) == len(ka) and len(set(kb)) == len(kb)
        if unique:
            index_b = {k: j for j, k in enumerate(kb)}
            by_quanta = {i: index_b[k] for i, k in enumerate(ka) if k in index_b}
            if match == "quanta":
                mapping, method = by_quanta, "quanta"
            elif len(by_quanta) >= 0.5 * min(len(ka), len(kb)):
                rest = _match_energy(a, b, atol, rtol, set(by_quanta), set(by_quanta.values()))
                mapping, method = {**by_quanta, **rest}, "quanta+energy" if rest else "quanta"
        elif match == "quanta":
            raise ValueError("level quanta are missing or not unique; use match='energy' or 'index'")
    if mapping is None and match in ("auto", "energy"):
        mapping, method = _match_energy(a, b, atol, rtol, set(), set()), "energy"
    if mapping is None:  # index
        mapping, method = {i: i for i in range(min(a.n_levels, b.n_levels))}, "index"
    de = [abs(float(b.level_energy[j] - a.level_energy[i])) for i, j in mapping.items()]
    return LevelMatch({i + 1: j + 1 for i, j in sorted(mapping.items())}, method, max(de, default=0.0))


# ------------------------------------------------------------------- rates on a grid
def interpolate_rates(temperatures: np.ndarray, rates: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Rate coefficients on another temperature grid, linear in log T and log k.

    Outside the table, the rates at the nearest table temperature are used (as RADEX
    does). A rate of 0 stays (about) 0.

    Parameters
    ----------
    temperatures : numpy.ndarray
        Temperatures of the table [K], increasing, shape (ntemp,).
    rates : numpy.ndarray
        Rate coefficients [cm3 s-1], shape (nrow, ntemp).
    grid : numpy.ndarray
        Temperatures to interpolate to [K], shape (ngrid,).

    Returns
    -------
    numpy.ndarray
        Rate coefficients [cm3 s-1], shape (nrow, ngrid); values below 1e-90 are set to 0.
    """
    rates = np.atleast_2d(np.asarray(rates, dtype=float))
    grid = np.asarray(grid, dtype=float)
    out = np.zeros((rates.shape[0], len(grid)))
    if rates.size == 0 or len(grid) == 0:
        return out
    lt, lg = np.log(temperatures), np.log(grid)
    lr = np.log(np.maximum(np.nan_to_num(rates, nan=0.0), _RATE_FLOOR * 1e-10))
    for k in range(rates.shape[0]):
        out[k] = np.exp(np.interp(lg, lt, lr[k]))
    out[out < _RATE_FLOOR] = 0.0
    # table values exactly at the table temperatures and outside the table (no round-off)
    t = np.asarray(temperatures, dtype=float)
    idx = np.clip(np.searchsorted(t, grid), 0, len(t) - 1)
    exact = (t[idx] == grid) | (grid <= t[0]) | (grid >= t[-1])
    idx = np.where(grid <= t[0], 0, np.where(grid >= t[-1], len(t) - 1, idx))
    out[:, exact] = np.nan_to_num(rates[:, idx[exact]], nan=0.0)
    return out


def total_rates(data: LamdaData, collisions: CollisionSet, temperatures: np.ndarray | None = None) -> np.ndarray:
    """Total de-excitation rate coefficient of each level, sum_l k_ul(T).

    Parameters
    ----------
    data : LamdaData
        The data set (for the number of levels).
    collisions : CollisionSet
        The collision data of one partner.
    temperatures : numpy.ndarray or None, optional
        Temperatures [K]. Default: the partner's grid. Rates between grid temperatures
        are interpolated in log T and log k (:func:`interpolate_rates`).

    Returns
    -------
    numpy.ndarray
        Shape (n_levels, ntemp) [cm3 s-1], row ``u - 1`` for upper level ``u``; 0 for
        levels without rows. Rows with level numbers outside 1..NLEV are left out.
    """
    t = collisions.temperatures if temperatures is None else np.asarray(temperatures, dtype=float)
    k = interpolate_rates(collisions.temperatures, collisions.rates, t)
    out = np.zeros((data.n_levels, len(t)))
    up = collisions.up
    ok = (up >= 1) & (up <= data.n_levels) & (collisions.low >= 1) & (collisions.low <= data.n_levels)
    np.add.at(out, up[ok] - 1, k[ok])
    return out


def critical_density(
    data: LamdaData, partner: str | CollisionSet, temperatures: np.ndarray | None = None
) -> np.ndarray:
    """Critical density of each level, n_crit,u(T) = sum_l A_ul / sum_l k_ul(T).

    Parameters
    ----------
    data : LamdaData
        The data set.
    partner : str or CollisionSet
        The collision partner (a name, see :meth:`LamdaData.partner
        <pylamda.model.LamdaData.partner>`, or its collision data).
    temperatures : numpy.ndarray or None, optional
        Temperatures [K]. Default: the partner's grid.

    Returns
    -------
    numpy.ndarray
        Critical densities [cm-3], shape (n_levels, ntemp), row ``u - 1`` for level
        ``u``; NaN where a level has no radiative decay or no collisional de-excitation.
    """
    c = data.partner(partner) if isinstance(partner, str) else partner
    k = total_rates(data, c, temperatures)
    a = np.zeros(data.n_levels)
    up = data.line_up
    ok = (up >= 1) & (up <= data.n_levels) & np.isfinite(data.line_a)
    np.add.at(a, up[ok] - 1, data.line_a[ok])
    with np.errstate(divide="ignore", invalid="ignore"):
        n = a[:, None] / k
    n[~((a[:, None] > 0) & (k > 0))] = np.nan
    return n


# ------------------------------------------------------------------- ratio statistics
@dataclass
class RatioSummary:
    """Statistics of ratios b/a of one quantity (rates, critical densities, A values).

    Attributes
    ----------
    n : int
        Number of finite positive ratios (items x temperatures).
    n_items : int
        Number of items (transitions or levels) with at least one ratio.
    n_changed : int
        Items whose ratio differs from 1 by more than ``limit`` at some temperature.
    limit : float
        The relative change that counts in ``n_changed``.
    median, p05, p95 : float
        Median, 5th and 95th percentile of all ratios.
    max_factor : float
        Largest deviation as a factor ``max(b/a, a/b)``.
    max_ratio : float
        The ratio b/a of that deviation (below 1 if b is smaller).
    max_at : str
        Item (and temperature) of that deviation.
    by_temperature : list of dict
        Per temperature: ``temperature``, ``n``, ``median``, ``p05``, ``p95``,
        ``max_factor``, ``max_ratio`` (empty for a quantity without temperatures).
    by_eup : list of dict
        Per E_up range (:data:`EUP_BINS`) with items: ``eup_min``, ``eup_max`` [K],
        ``n_items``, ``n_changed``, ``median``, ``max_factor``, ``max_ratio``.
    """

    n: int = 0
    n_items: int = 0
    n_changed: int = 0
    limit: float = 0.1
    median: float = 1.0
    p05: float = 1.0
    p95: float = 1.0
    max_factor: float = 1.0
    max_ratio: float = 1.0
    max_at: str = ""
    by_temperature: list[dict] = field(default_factory=list)
    by_eup: list[dict] = field(default_factory=list)

    def eup_text(self) -> str:
        """The changes per E_up range in one line.

        ``"0-100 K 0 of 12, 100-300 K 3 of 40 (max change B/A 2.1)"``: per range, the items
        changed by more than :attr:`limit` of all items, and the largest change as the
        ratio B/A, as all ratios in the summary.
        """
        parts = []
        for x in self.by_eup:
            hi = "" if math.isinf(x["eup_max"]) else f"-{x['eup_max']:g}"
            text = f"{x['eup_min']:g}{hi if hi else '+'} K {x['n_changed']} of {x['n_items']}"
            if x["n_changed"]:
                text += f" (max change B/A {x['max_ratio']:.3g})"
            parts.append(text)
        return ", ".join(parts)


def summarise_ratios(
    ratio: np.ndarray,
    labels: list[str],
    eup: np.ndarray,
    temperatures: list[float] | np.ndarray | None,
    limit: float,
) -> RatioSummary | None:
    """Statistics of a ratio table b/a (items x temperatures).

    Parameters
    ----------
    ratio : numpy.ndarray
        Ratios, shape (nitems, ntemp) or (nitems,); NaN or <= 0 where undefined.
    labels : list of str
        One label per item (for ``max_at``).
    eup : numpy.ndarray
        Upper-level energy [K] of each item (for the E_up ranges).
    temperatures : sequence of float or None
        The temperatures of the columns, or None for one column without temperature.
    limit : float
        Relative change that counts as a change.

    Returns
    -------
    RatioSummary or None
        The statistics, or None if there is no finite positive ratio.
    """
    q = np.asarray(ratio, dtype=float)
    if q.ndim == 1:
        q = q[:, None]
    ok = np.isfinite(q) & (q > 0)
    if not ok.any():
        return None
    factor = np.where(ok, np.maximum(q, 1.0 / np.where(ok, q, 1.0)), 1.0)
    changed = ok & (np.abs(q - 1.0) > limit)
    k, j = np.unravel_index(int(np.argmax(factor)), q.shape)
    at = labels[k] + (f" at {temperatures[j]:g} K" if temperatures is not None else "")
    vals = q[ok]
    s = RatioSummary(
        n=int(ok.sum()),
        n_items=int(ok.any(axis=1).sum()),
        n_changed=int(changed.any(axis=1).sum()),
        limit=limit,
        median=float(np.median(vals)),
        p05=float(np.percentile(vals, 5)),
        p95=float(np.percentile(vals, 95)),
        max_factor=float(factor[k, j]),
        max_ratio=float(q[k, j]),
        max_at=at,
    )
    if temperatures is not None:
        for j, t in enumerate(temperatures):
            col = ok[:, j]
            if col.any():
                v = q[col, j]
                s.by_temperature.append(
                    {
                        "temperature": float(t),
                        "n": int(col.sum()),
                        "median": float(np.median(v)),
                        "p05": float(np.percentile(v, 5)),
                        "p95": float(np.percentile(v, 95)),
                        "max_factor": float(factor[col, j].max()),
                        "max_ratio": float(v[np.argmax(factor[col, j])]),
                    }
                )
    for lo, hi in zip(EUP_BINS[:-1], EUP_BINS[1:], strict=True):
        rows = (eup >= lo) & (eup < hi) & ok.any(axis=1)
        if rows.any():
            v = q[rows][ok[rows]]
            s.by_eup.append(
                {
                    "eup_min": lo,
                    "eup_max": hi,
                    "n_items": int(rows.sum()),
                    "n_changed": int(changed[rows].any(axis=1).sum()),
                    "median": float(np.median(v)),
                    "max_factor": float(factor[rows].max()),
                    "max_ratio": float(np.where(ok, q, 1.0)[rows].flat[np.argmax(factor[rows])]),
                }
            )
    return s


# ------------------------------------------------------------------- results
@dataclass
class PartnerDiff:
    """Differences of the collision data of one partner (see :class:`LamdaDiff`).

    Attributes
    ----------
    n_transitions : tuple of int
        Number of collisional transitions in ``a`` and ``b``.
    t_range : tuple of tuple of float
        (first, last) temperature in ``a`` and ``b`` [K].
    temperatures_added, temperatures_removed : list of float
        Grid temperatures only in ``b`` / only in ``a``.
    transitions_added, transitions_removed : list of str
        Transitions only in ``b`` / only in ``a`` (``"upper -> lower"`` level labels).
    n_rates_compared : int
        Number of (transition, temperature) rates compared: common transitions on the
        common temperatures.
    n_temperatures_compared : int
        Number of common temperatures.
    n_rates_changed : int
        Number of compared rates that differ: by more than ``rtol`` relative to ``a``
        where both are positive, and every difference where one of them is 0.
    n_rates_zero_changed : int
        How many of these go from 0 to a value or from a value to 0 (no ratio exists).
    max_rate_factor : float
        Largest deviation of a rate with both values positive, as a factor
        ``max(b/a, a/b)`` (1 if none).
    max_rate_ratio : float
        The ratio b/a of that same rate (1 if none).
    max_rate_change : float
        The relative change ``|b - a| / a`` of that same rate (0 if none).
    max_rate_change_at : str
        Transition and temperature of that rate ("" if none).
    ratio : dict
        Ratio ``b/a`` of all compared rates (both positive): ``median``, ``p05``,
        ``p95``; empty if no compared pair has both rates positive.
    ratio_by_temperature : list of dict
        The same per common temperature: ``temperature``, ``median``, ``p05``, ``p95``,
        ``max_factor``.
    changed_transitions : list of dict
        Transitions with at least one changed rate: ``transition``, ``ratio`` (b/a of its
        largest deviation; None if the change is to or from 0), ``temperature``,
        ``rate_a`` and ``rate_b`` there, ``n_changed`` (changed temperatures). Changes to
        or from 0 come first, then by deviation, largest first.
    reference_changed : bool
        Whether the reference text differs.
    interp_temperatures : list of float
        The union grid of both tables in the common temperature range [K].
    rates_interpolated : RatioSummary or None
        Rates b/a of the common transitions on ``interp_temperatures`` (log-log
        interpolated), with the summary per E_up range.
    total_rate : RatioSummary or None
        Total de-excitation rates b/a of the matched levels on ``interp_temperatures``.
    """

    n_transitions: tuple[int, int]
    t_range: tuple[tuple[float, float], tuple[float, float]]
    temperatures_added: list[float] = field(default_factory=list)
    temperatures_removed: list[float] = field(default_factory=list)
    transitions_added: list[str] = field(default_factory=list)
    transitions_removed: list[str] = field(default_factory=list)
    n_rates_compared: int = 0
    n_temperatures_compared: int = 0
    n_rates_changed: int = 0
    n_rates_zero_changed: int = 0
    max_rate_change: float = 0.0
    max_rate_factor: float = 1.0
    max_rate_ratio: float = 1.0
    max_rate_change_at: str = ""
    ratio: dict = field(default_factory=dict)
    ratio_by_temperature: list[dict] = field(default_factory=list)
    changed_transitions: list[dict] = field(default_factory=list)
    reference_changed: bool = False
    interp_temperatures: list[float] = field(default_factory=list)
    rates_interpolated: RatioSummary | None = None
    total_rate: RatioSummary | None = None

    def __bool__(self) -> bool:
        return bool(
            self.n_transitions[0] != self.n_transitions[1]
            or self.temperatures_added
            or self.temperatures_removed
            or self.transitions_added
            or self.transitions_removed
            or self.n_rates_changed
            or self.reference_changed
        )


@dataclass
class LamdaDiff:
    """Result of :func:`compare`: what differs from ``a`` to ``b``.

    ``bool(diff)`` is True if anything differs; :meth:`summary` gives a readable report,
    :meth:`as_dict` a JSON-compatible one.

    Attributes
    ----------
    header : list of str
        Changes of molecule name, mass, dataset DOI, level type or level reference.
    level_matching : str
        How levels were matched (:attr:`LevelMatch.method`).
    n_levels, n_lines : tuple of int
        Counts in ``a`` and ``b``.
    levels_added, levels_removed : list of str
        Level labels only in ``b`` / only in ``a``.
    levels_changed : list of str
        Matched levels whose energy (beyond the tolerance) or weight differs.
    max_energy_change : float
        Largest energy change of a matched level [cm-1].
    lines_added, lines_removed : list of str
        Radiative transitions only in ``b`` / only in ``a``.
    n_a_changed : int
        Matched transitions whose A differs by more than ``rtol``.
    max_a_change : float
        Largest relative change of A.
    max_a_factor : float
        Largest deviation of A as a factor, ``max(b/a, a/b)``.
    max_a_ratio : float
        The ratio b/a of that A value (1 if none).
    n_freq_changed : int
        Matched transitions whose frequency differs by more than ``rtol``.
    a_ratio : RatioSummary or None
        A values b/a of the matched transitions, with the summary per E_up range.
    partners_added, partners_removed : list of str
        Collision partners only in ``b`` / only in ``a``.
    partners : dict of str to PartnerDiff
        Differences per common partner (only partners with differences). Partners are
        identified by their ID; a second block with the same ID is ``"<name> #2"``.
    rtol : float
        The relative tolerance of the exact comparison (see :func:`compare`).
    limit : float
        The relative change that counts in the summaries on the union grid, per level
        and per E_up range.
    """

    header: list[str] = field(default_factory=list)
    level_matching: str = ""
    n_levels: tuple[int, int] = (0, 0)
    levels_added: list[str] = field(default_factory=list)
    levels_removed: list[str] = field(default_factory=list)
    levels_changed: list[str] = field(default_factory=list)
    max_energy_change: float = 0.0
    n_lines: tuple[int, int] = (0, 0)
    lines_added: list[str] = field(default_factory=list)
    lines_removed: list[str] = field(default_factory=list)
    n_a_changed: int = 0
    max_a_change: float = 0.0
    max_a_factor: float = 1.0
    max_a_ratio: float = 1.0
    n_freq_changed: int = 0
    a_ratio: RatioSummary | None = None
    partners_added: list[str] = field(default_factory=list)
    partners_removed: list[str] = field(default_factory=list)
    partners: dict[str, PartnerDiff] = field(default_factory=dict)
    rtol: float = 1e-6
    limit: float = 0.1

    def __bool__(self) -> bool:
        return bool(
            self.header
            or self.n_levels[0] != self.n_levels[1]
            or self.levels_added
            or self.levels_removed
            or self.levels_changed
            or self.n_lines[0] != self.n_lines[1]
            or self.lines_added
            or self.lines_removed
            or self.n_a_changed
            or self.n_freq_changed
            or self.partners_added
            or self.partners_removed
            or self.partners
        )

    def summary(self, max_items: int = 5) -> list[str]:
        """The differences as lines of text (empty if there are none).

        Parameters
        ----------
        max_items : int, optional
            How many added/removed items to name before "...". Default: 5.

        Returns
        -------
        list of str
            One line per difference (partner details indented by two blanks).
        """

        def items(xs: list[str]) -> str:
            return ", ".join(xs[:max_items]) + (f", ... ({len(xs)} in all)" if len(xs) > max_items else "")

        out = list(self.header)
        if self.n_levels[0] != self.n_levels[1] or self.levels_added or self.levels_removed or self.levels_changed:
            out.append(f"levels: {self.n_levels[0]} -> {self.n_levels[1]} (matched by {self.level_matching})")
            if self.levels_added:
                out.append(f"  added: {items(self.levels_added)}")
            if self.levels_removed:
                out.append(f"  removed: {items(self.levels_removed)}")
            if self.levels_changed:
                out.append(
                    f"  E or g changed: {items(self.levels_changed)} (max |dE| = {self.max_energy_change:.3g} cm-1)"
                )
        if (
            self.n_lines[0] != self.n_lines[1]
            or self.lines_added
            or self.lines_removed
            or self.n_a_changed
            or self.n_freq_changed
        ):
            out.append(f"radiative transitions: {self.n_lines[0]} -> {self.n_lines[1]}")
            if self.lines_added:
                out.append(f"  added: {items(self.lines_added)}")
            if self.lines_removed:
                out.append(f"  removed: {items(self.lines_removed)}")
            if self.n_a_changed:
                out.append(
                    f"  A differs for {self.n_a_changed} by more than {self.rtol:g} (relative), "
                    f"max change B/A {self.max_a_ratio:.3g}"
                )
                if self.a_ratio and self.a_ratio.n_changed:
                    out.append(
                        f"  A by E_up of the upper level (changed by > {_percent(self.limit)}): "
                        f"{self.a_ratio.eup_text()}"
                    )
            if self.n_freq_changed:
                out.append(f"  frequency changed for {self.n_freq_changed}")
        for p in self.partners_added:
            out.append(f"new collision partner {p}")
        for p in self.partners_removed:
            out.append(f"collision partner {p} removed")
        for name, d in self.partners.items():
            out += _partner_summary(name, d, self.rtol, self.limit, items)
        return out

    def as_dict(self) -> dict:
        """The differences as a JSON-compatible dict.

        Returns
        -------
        dict
            All fields of the dataclass (partners and ratio summaries as dicts;
            an infinite E_up bin edge as None).
        """
        return _json_safe(asdict(self))


def _json_safe(x):
    """Replace inf/nan floats by None (JSON has no such numbers)."""
    if isinstance(x, dict):
        return {k: _json_safe(v) for k, v in x.items()}
    if isinstance(x, list | tuple):
        return [_json_safe(v) for v in x]
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def _percent(x: float) -> str:
    """A relative change as percent (``0.1`` -> ``"10%"``, ``0.005`` -> ``"0.5%"``)."""
    return f"{x * 100:g}%"


def _partner_summary(name: str, d: PartnerDiff, rtol: float, limit: float, items) -> list[str]:
    """The summary lines of one partner."""
    (t0a, t1a), (t0b, t1b) = d.t_range
    out = [
        f"{name}: {d.n_transitions[0]} -> {d.n_transitions[1]} transitions, T {t0a:g}-{t1a:g} K -> {t0b:g}-{t1b:g} K"
    ]
    if d.temperatures_added:
        out.append(f"  temperatures added: {items([f'{t:g}' for t in d.temperatures_added])}")
    if d.temperatures_removed:
        out.append(f"  temperatures removed: {items([f'{t:g}' for t in d.temperatures_removed])}")
    if d.transitions_added:
        out.append(f"  transitions added: {items(d.transitions_added)}")
    if d.transitions_removed:
        out.append(f"  transitions removed: {items(d.transitions_removed)}")
    if d.n_rates_changed:
        r = d.ratio
        where = f"{d.n_rates_compared} compared on {d.n_temperatures_compared} common temperatures"
        if r:
            out.append(f"  rates B/A ({where}): median {r['median']:.3g}, 90% within {r['p05']:.3g}-{r['p95']:.3g}")
        else:  # no pair with both rates > 0 (e.g. all old rates 0): no ratios
            out.append(f"  rates B/A ({where}): no pair with both rates > 0")
        if len(d.ratio_by_temperature) > 1:
            # only the temperature whose median is farthest from 1 (the full list is in ratio_by_temperature)
            x = max(d.ratio_by_temperature, key=lambda x: abs(math.log(x["median"])) if x["median"] > 0 else math.inf)
            out.append(
                f"  largest deviation by temperature at {x['temperature']:g} K: median {x['median']:.3g}, "
                f"90% within {x['p05']:.3g}-{x['p95']:.3g}"
            )
        # one line: the counts and only the largest change (all changed transitions are in changed_transitions)
        n, z, nt = d.n_rates_changed, d.n_rates_zero_changed, len(d.changed_transitions)
        line = f"  {n} rate{'s' if n != 1 else ''} in {nt} transition{'s' if nt != 1 else ''} "
        line += "differ" if n != 1 else "differs"
        if n > z:
            line += f" by more than {rtol:g} (relative), max change B/A {d.max_rate_ratio:.3g} ({d.max_rate_change_at})"
        if z:
            line += f"{';' if n > z else ':'} {'all' if z == n else z} to or from a zero rate"
            c = next(c for c in d.changed_transitions if c["ratio"] is None)  # these come first
            line += f" ({c['transition']}: {c['rate_a']:.3g} -> {c['rate_b']:.3g} at {c['temperature']:g} K)"
        out.append(line)
    if d.reference_changed:
        out.append("  reference changed")
    s = d.rates_interpolated
    if s is not None and s.n_changed and (d.temperatures_added or d.temperatures_removed):
        out.append(
            f"  rates B/A on the union grid ({len(d.interp_temperatures)} temperatures, log-log interpolated): "
            f"median {s.median:.3g}, 90% within {s.p05:.3g}-{s.p95:.3g}, max change {s.max_ratio:.3g} "
            f"({s.max_at})"
        )
    if s is not None and s.n_changed:
        out.append(f"  by E_up of the upper level (transitions changed by > {_percent(limit)}): {s.eup_text()}")
    s = d.total_rate
    if s is not None and s.n_changed:
        out.append(
            f"  total de-excitation rate per level B/A: {s.n_changed} of {s.n_items} levels changed by > "
            f"{_percent(limit)}, median {s.median:.3g}, max change {s.max_ratio:.3g} ({s.max_at}); "
            f"by level energy: {s.eup_text()}"
        )
    return out


# ------------------------------------------------------------------- comparing
def _label(data: LamdaData, i: int) -> str:
    """Readable label of level index ``i`` (0-based): quanta, else the energy."""
    q = " ".join(data.level_quanta[i].split())
    return q if q else f"E={data.level_energy[i]:.6g}"


def _check_levels(data: LamdaData, what: str) -> None:
    """Raise ValueError if a transition refers to a level outside 1..NLEV."""
    n = data.n_levels
    for name, up, low in [("radiative transitions", data.line_up, data.line_low)] + [
        (f"{c.label} collisional transitions", c.up, c.low) for c in data.collisions
    ]:
        bad = np.flatnonzero((up < 1) | (up > n) | (low < 1) | (low > n))
        if bad.size:
            raise ValueError(
                f"{what}: {name}: {bad.size} transition(s) refer to levels outside 1..{n} (see pylamda.verify, R01/C10)"
            )


def _keys(data: LamdaData) -> dict[str, CollisionSet]:
    """Partner key -> collision set: the partner name, ``"<name> #2"`` for a repeated ID."""
    seen: dict[int, int] = {}
    out = {}
    for c in data.collisions:
        n = seen.get(c.partner_id, 0) + 1
        seen[c.partner_id] = n
        out[c.label if n == 1 else f"{c.label} #{n}"] = c
    return out


def _key_map(ups: np.ndarray, lows: np.ndarray, mapping: dict[int, int] | None = None) -> dict[tuple, int]:
    """{(u, l, n) in b numbering: row} for transitions; ``n`` counts repeated (u, l) pairs."""
    out = {}
    seen: dict[tuple[int, int], int] = {}
    for k, (u, lo) in enumerate(zip(ups, lows, strict=True)):
        if mapping is None:
            key = (int(u), int(lo))
        elif int(u) in mapping and int(lo) in mapping:
            key = (mapping[int(u)], mapping[int(lo)])
        else:
            continue
        n = seen.get(key, 0)
        seen[key] = n + 1
        out[(*key, n)] = k
    return out


def _exact_rates(pd: PartnerDiff, ra: np.ndarray, rb: np.ndarray, labels: list[str], temps: list, rtol: float) -> None:
    """Fill the exact-comparison fields of ``pd`` from matched rates ``ra``, ``rb`` (transitions x T).

    Ratios and factors are only defined where both rates are positive. A change to or
    from a zero rate is counted separately (``n_rates_zero_changed``) and listed with
    its two values; the largest deviation (``max_rate_factor``, ``max_rate_change``,
    ``max_rate_change_at``) refers to one and the same rate with both values positive.
    """
    pos = (ra > 0) & (rb > 0)
    q = np.where(pos, rb / np.where(pos, ra, 1.0), np.nan)
    factor = np.where(pos, np.maximum(q, 1.0 / np.where(pos, q, 1.0)), 1.0)
    changed = np.where(pos, np.abs(q - 1.0) > rtol, ra != rb)
    zero_changed = changed & ~pos
    pd.n_rates_compared = int(ra.size)
    pd.n_temperatures_compared = len(temps)
    pd.n_rates_changed = int(changed.sum())
    pd.n_rates_zero_changed = int(zero_changed.sum())
    if pos.any():
        pd.ratio = {
            "median": float(np.median(q[pos])),
            "p05": float(np.percentile(q[pos], 5)),
            "p95": float(np.percentile(q[pos], 95)),
        }
    for j, t in enumerate(temps):
        col = pos[:, j]
        if col.any():
            pd.ratio_by_temperature.append(
                {
                    "temperature": float(t),
                    "median": float(np.median(q[col, j])),
                    "p05": float(np.percentile(q[col, j], 5)),
                    "p95": float(np.percentile(q[col, j], 95)),
                    "max_factor": float(factor[col, j].max()),
                }
            )
    if not pd.n_rates_changed:
        return
    pos_changed = changed & pos
    if pos_changed.any():
        k, j = np.unravel_index(int(np.argmax(np.where(pos_changed, factor, 0.0))), q.shape)
        pd.max_rate_factor = float(factor[k, j])
        pd.max_rate_ratio = float(q[k, j])
        pd.max_rate_change = float(abs(q[k, j] - 1.0))
        pd.max_rate_change_at = f"{labels[k]} at {temps[j]:g} K"
    rows = []
    for k in np.flatnonzero(changed.any(axis=1)):
        n_changed = int(changed[k].sum())
        if zero_changed[k].any():  # appearing or vanishing rates first, with their values
            j = int(np.flatnonzero(zero_changed[k])[0])
            ratio, key = None, math.inf
        else:
            j = int(np.argmax(np.where(pos_changed[k], factor[k], 0.0)))
            ratio, key = float(q[k, j]), float(factor[k, j])
        row = {
            "transition": labels[k],
            "ratio": ratio,
            "temperature": float(temps[j]),
            "rate_a": float(ra[k, j]),
            "rate_b": float(rb[k, j]),
            "n_changed": n_changed,
        }
        rows.append((key, row))
    rows.sort(key=lambda r: -r[0])
    pd.changed_transitions = [r for _, r in rows]


def _compare_partner(
    a: LamdaData,
    b: LamdaData,
    ca: CollisionSet,
    cb: CollisionSet,
    levels: LevelMatch,
    rtol: float,
    limit: float,
) -> PartnerDiff:
    """The differences of the collision data of one partner."""

    def tlabel(data: LamdaData, u: int, lo: int) -> str:
        return f"{_label(data, u - 1)} -> {_label(data, lo - 1)}"

    pd = PartnerDiff(
        (ca.n_transitions, cb.n_transitions),
        (
            (float(ca.temperatures[0]), float(ca.temperatures[-1])),
            (float(cb.temperatures[0]), float(cb.temperatures[-1])),
        ),
        reference_changed=ca.reference != cb.reference,
    )
    ta = {round(float(t), 6): i for i, t in enumerate(ca.temperatures)}
    tb = {round(float(t), 6): i for i, t in enumerate(cb.temperatures)}
    pd.temperatures_added = [t for t in tb if t not in ta]
    pd.temperatures_removed = [t for t in ta if t not in tb]
    common_t = [t for t in ta if t in tb]
    xa = _key_map(ca.up, ca.low, levels.mapping)
    xb = _key_map(cb.up, cb.low)
    kept = {v for key, v in xa.items() if key in xb}
    pd.transitions_removed = [
        tlabel(a, int(ca.up[k]), int(ca.low[k])) for k in range(ca.n_transitions) if k not in kept
    ]
    pd.transitions_added = [tlabel(b, *key[:2]) for key in xb if key not in xa]
    pairs = [(xa[key], xb[key], key) for key in xa if key in xb]
    labels = [tlabel(b, *p[2][:2]) for p in pairs]
    ia = np.array([p[0] for p in pairs], dtype=int)
    ib = np.array([p[1] for p in pairs], dtype=int)
    if common_t and pairs:
        cols_a = [ta[t] for t in common_t]
        cols_b = [tb[t] for t in common_t]
        _exact_rates(pd, ca.rates[ia][:, cols_a], cb.rates[ib][:, cols_b], labels, common_t, rtol)

    # union grid in the common temperature range, log-log interpolated
    lo_t = max(ca.temperatures.min(), cb.temperatures.min())
    hi_t = min(ca.temperatures.max(), cb.temperatures.max())
    grid = np.unique(np.concatenate([ca.temperatures, cb.temperatures]))
    grid = grid[(grid >= lo_t) & (grid <= hi_t)]
    pd.interp_temperatures = [float(t) for t in grid]
    if not len(grid):
        return pd
    if pairs:
        ra = interpolate_rates(ca.temperatures, ca.rates[ia], grid)
        rb = interpolate_rates(cb.temperatures, cb.rates[ib], grid)
        with np.errstate(divide="ignore", invalid="ignore"):
            q = np.where((ra > 0) & (rb > 0), rb / ra, np.nan)
        eup = a.level_energy[ca.up[ia] - 1] * HC_K
        pd.rates_interpolated = summarise_ratios(q, labels, eup, pd.interp_temperatures, limit)

    # per level: total de-excitation rate (all rows of each set)
    ka, kb = total_rates(a, ca, grid), total_rates(b, cb, grid)
    lev_a = np.array(list(levels.mapping.keys()), dtype=int) - 1
    lev_b = np.array(list(levels.mapping.values()), dtype=int) - 1
    level_labels = [f"level {_label(a, i)}" for i in lev_a]
    eup = a.level_energy[lev_a] * HC_K
    with np.errstate(divide="ignore", invalid="ignore"):
        k_ratio = np.where((ka[lev_a] > 0) & (kb[lev_b] > 0), kb[lev_b] / ka[lev_a], np.nan)
    pd.total_rate = summarise_ratios(k_ratio, level_labels, eup, pd.interp_temperatures, limit)
    return pd


def compare(
    a: LamdaData | str | os.PathLike,
    b: LamdaData | str | os.PathLike,
    match: str = "auto",
    rtol: float = 1e-6,
    atol: float = 1e-3,
    limit: float = 0.1,
) -> LamdaDiff:
    """Compare two LAMDA data sets: what differs from ``a`` to ``b``.

    Parameters
    ----------
    a, b : LamdaData, str or os.PathLike
        The data sets or files (read with :func:`pylamda.read`).
    match : str, optional
        How levels are matched (see :func:`match_levels`). Default: ``"auto"``.
    rtol : float, optional
        Relative change of an energy, A, a frequency or a rate coefficient that counts
        as a change in the exact comparison. Default: 1e-6 (only real changes, not
        rounding).
    atol : float, optional
        Energy tolerance [cm-1] for matching levels by energy. Default: 1e-3.
    limit : float, optional
        Relative change that counts in the summaries on the union temperature grid, per
        level (total rates) and per E_up range. Default: 0.1 (shown as 10%).

    Returns
    -------
    LamdaDiff
        The differences.

    Raises
    ------
    LamdaReadError
        If a file cannot be read (a subclass of ValueError).
    ValueError
        If a transition refers to a level outside 1..NLEV, or ``match`` is unknown or
        cannot be used.
    """
    a = a if isinstance(a, LamdaData) else read(a)
    b = b if isinstance(b, LamdaData) else read(b)
    _check_levels(a, "a")
    _check_levels(b, "b")
    d = LamdaDiff(n_levels=(a.n_levels, b.n_levels), n_lines=(a.n_lines, b.n_lines), rtol=rtol, limit=limit)
    for attr, what in (
        ("molecule", "molecule"),
        ("mass", "mass"),
        ("dataset_doi", "dataset DOI"),
        ("level_type", "level type"),
        ("level_reference", "level reference"),
    ):
        va, vb = getattr(a, attr), getattr(b, attr)
        if va != vb:
            d.header.append(f"{what}: {va!r} -> {vb!r}")

    levels = match_levels(a, b, match, atol, 1e-5)
    d.level_matching = levels.method
    lev = {i - 1: j - 1 for i, j in levels.mapping.items()}  # 0-based
    matched_b = set(lev.values())
    d.levels_removed = [_label(a, i) for i in range(a.n_levels) if i not in lev]
    d.levels_added = [_label(b, j) for j in range(b.n_levels) if j not in matched_b]
    for i, j in lev.items():
        de = abs(b.level_energy[j] - a.level_energy[i])
        if de > rtol * max(abs(a.level_energy[i]), 1.0) or b.level_g[j] != a.level_g[i]:
            d.levels_changed.append(_label(a, i))
            d.max_energy_change = max(d.max_energy_change, float(de))

    def tlabel(data: LamdaData, u: int, lo: int) -> str:
        return f"{_label(data, u - 1)} -> {_label(data, lo - 1)}"

    la = _key_map(a.line_up, a.line_low, levels.mapping)
    lb = _key_map(b.line_up, b.line_low)
    kept = {v for key, v in la.items() if key in lb}
    d.lines_removed = [tlabel(a, int(a.line_up[k]), int(a.line_low[k])) for k in range(a.n_lines) if k not in kept]
    d.lines_added = [tlabel(b, *key[:2]) for key in lb if key not in la]
    pairs = [(ka, lb[key], key) for key, ka in la.items() if key in lb]
    for ka, kb, _ in pairs:
        change = abs(b.line_a[kb] - a.line_a[ka]) / max(abs(a.line_a[ka]), 1e-300)
        if change > rtol:
            d.n_a_changed += 1
            d.max_a_change = max(d.max_a_change, float(change))
            if a.line_a[ka] > 0 and b.line_a[kb] > 0:
                q = float(b.line_a[kb] / a.line_a[ka])
                if max(q, 1 / q) > d.max_a_factor:
                    d.max_a_factor, d.max_a_ratio = max(q, 1 / q), q
        fa, fb = a.line_freq[ka], b.line_freq[kb]
        if math.isfinite(fa) and math.isfinite(fb) and abs(fb - fa) > rtol * abs(fa):
            d.n_freq_changed += 1
    if pairs:
        ia = np.array([p[0] for p in pairs], dtype=int)
        ib = np.array([p[1] for p in pairs], dtype=int)
        with np.errstate(divide="ignore", invalid="ignore"):
            q = np.where((a.line_a[ia] > 0) & (b.line_a[ib] > 0), b.line_a[ib] / a.line_a[ia], np.nan)
        d.a_ratio = summarise_ratios(
            q, [tlabel(b, *p[2][:2]) for p in pairs], a.level_energy[a.line_up[ia] - 1] * HC_K, None, limit
        )

    pa, pb = _keys(a), _keys(b)
    d.partners_added = [p for p in pb if p not in pa]
    d.partners_removed = [p for p in pa if p not in pb]
    for name, ca in pa.items():
        if name in pb:
            pd = _compare_partner(a, b, ca, pb[name], levels, rtol, limit)
            if pd:
                d.partners[name] = pd
    return d
