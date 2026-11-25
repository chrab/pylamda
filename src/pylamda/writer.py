"""Write a :class:`~pylamda.model.LamdaData` as a LAMDA file in one normalised layout."""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np

from .model import CollisionSet, LamdaData, partner_name


def _num(x: float) -> str:
    """Shortest text that reads back as the same float (Fortran-readable)."""
    x = float(x)
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    return repr(x)


def _num_exp(x: float) -> str:
    """As :func:`_num`, always in exponential notation (``1.94e+00``, ``8.542e-05``)."""
    if not math.isfinite(x):
        return _num(x)
    return np.format_float_scientific(float(x), unique=True, trim="0", exp_digits=2)


def _rate_texts(rates: np.ndarray) -> list[list[str]]:
    """The rates of one block in exponential notation, all with the same number of digits.

    Each value is written as the shortest text that reads back unchanged
    (:func:`_num_exp`); then all are padded with zeros to the number of decimals of the
    most precise one (``2.3e-12`` -> ``2.30e-12`` next to ``2.01e-12``). The padding does
    not change the values.
    """
    texts = [[_num_exp(v) for v in row] for row in rates]
    decimals = [len(t.split("e")[0].partition(".")[2]) for row in texts for t in row if "e" in t]
    width = max(decimals, default=1)

    def pad(t: str) -> str:
        if "e" not in t:  # NaN, Inf
            return t
        mantissa, _, exponent = t.partition("e")
        return f"{mantissa}{'0' * (width - len(mantissa.partition('.')[2]))}e{exponent}"

    return [[pad(t) for t in row] for row in texts]


def _check_level_numbers(nlev: int, up: np.ndarray, low: np.ndarray, what: str) -> None:
    """Raise ValueError if a transition refers to a level outside 1..NLEV."""
    bad = np.flatnonzero((up < 1) | (up > nlev) | (low < 1) | (low > nlev))
    if bad.size:
        k = int(bad[0])
        raise ValueError(
            f"{what}: {bad.size} transition(s) refer to levels outside 1..{nlev}, "
            f"the first is transition {k + 1} ({int(up[k])} -> {int(low[k])})"
        )


def _partner_line(molecule: str, c: CollisionSet) -> str:
    """``<id> <molecule> + <name> -- Ref: <reference>``.

    The name is the one in the original partner text if there is one (so that a name
    that disagrees with the ID stays visible), else the standard name of the ID.
    """
    written = c.text_partner
    name = (partner_name(written) or written) if written else c.partner
    line = f"{c.partner_id} {molecule}" + (f" + {name}" if name else "")
    if c.reference:
        line += f" -- Ref: {c.reference}"
    return line


def write(data: LamdaData, path: str | os.PathLike) -> None:
    """Write ``data`` as a LAMDA file in pylamda's normalised layout.

    The layout is that of LAMDA's format description, with the references in EMAA's
    style:

    * ``<molecule> -- EMAA_DATASET_DOI: <doi>`` on line 2 (if there is a dataset DOI);
    * ``<NLEV> - <TYPE> - Ref: <reference>`` on the level-count line (if known);
    * ``<id> <molecule> + <partner> -- Ref: <reference>`` on each partner line.

    The partner ID is written as it is in ``data``. The partner name is the one in the
    original partner text if that text names one (also if it disagrees with the ID: the
    writer does not decide between the two), else the standard name of the ID.

    The source and the notes follow the data as comment lines (``!SOURCE: ...``,
    ``!NOTES``), which RADEX, ProDiMo and :func:`pylamda.read` do not read as data. All
    values are written as the shortest text that reads back unchanged; the Einstein A
    values and the rate coefficients always in exponential notation (``1.94e+00``,
    ``8.542e-05``, also ``0.0e+00``). All rates of one collision block have the same
    number of digits: those of the most precise rate (``2.30e-12`` next to ``2.01e-12``).

    Parameters
    ----------
    data : LamdaData
        The data.
    path : str or os.PathLike
        Output file (overwritten).

    Raises
    ------
    ValueError
        If a transition refers to a level outside 1..NLEV (nothing is written).
    """
    _check_level_numbers(data.n_levels, data.line_up, data.line_low, "radiative transitions")
    for c in data.collisions:
        _check_level_numbers(data.n_levels, c.up, c.low, f"{c.label} collisional transitions")

    out = ["!MOLECULE"]
    line2 = data.molecule + (f" {data.title}" if data.title else "")
    if data.dataset_doi:
        key = "EMAA_DATASET_DOI" if data.dataset_doi.startswith("10.17178/EMAA") else "DATASET_DOI"
        line2 += f" -- {key}: {data.dataset_doi}"
    out += [line2, "!MOLECULAR WEIGHT", _num(data.mass), "!NUMBER OF ENERGY LEVELS"]
    nlev = str(data.n_levels)
    if data.level_type:
        nlev += f" - {data.level_type}"
    if data.level_reference:
        nlev += f" - Ref: {data.level_reference}"
    out += [nlev, "!LEVEL + ENERGIES(cm^-1) + WEIGHT + QUANTUM NUMBERS"]
    for i in range(data.n_levels):
        row = f"{i + 1:5d} {_num(data.level_energy[i]):>18s} {_num(data.level_g[i]):>8s}"
        out.append(row + (f"   {data.level_quanta[i]}" if data.level_quanta[i] else ""))
    out += [
        "!NUMBER OF RADIATIVE TRANSITIONS",
        str(data.n_lines),
        "!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K)",
    ]
    for i in range(data.n_lines):
        row = f"{i + 1:5d} {data.line_up[i]:5d} {data.line_low[i]:5d} {_num_exp(data.line_a[i]):>12s}"
        if math.isfinite(data.line_freq[i]):
            row += f" {_num(data.line_freq[i]):>18s}"
            if math.isfinite(data.line_eup[i]):
                row += f" {_num(data.line_eup[i]):>12s}"
        out.append(row)
    out += ["!NUMBER OF COLL PARTNERS", str(len(data.collisions))]
    for c in data.collisions:
        out += [
            "!COLLISIONS BETWEEN",
            _partner_line(data.molecule, c),
            "!NUMBER OF COLL TRANS",
            str(c.n_transitions),
            "!NUMBER OF COLL TEMPS",
            str(len(c.temperatures)),
            "!COLL TEMPS",
            " ".join(_num(t) for t in c.temperatures),
            "!TRANS + UP + LOW + COLLRATES(cm^3 s^-1)",
        ]
        for k, texts in enumerate(_rate_texts(c.rates)):
            out.append(f"{k + 1:5d} {c.up[k]:5d} {c.low[k]:5d} " + " ".join(texts))
    if data.source:
        out.append(f"!SOURCE: {data.source}")
    if data.notes:
        out.append("!NOTES")
        out += [line if line.lstrip().startswith("!") or not line.strip() else f"!{line}" for line in data.notes]
    Path(path).write_text("\n".join(out) + "\n", encoding="utf-8")
