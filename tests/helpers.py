"""Helpers shared by the tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from pylamda import LamdaData

DATA = Path(__file__).parent / "data"
LAMDA = DATA / "lamda"
EMAA_CO = DATA / "emaa" / "CO.Rotation.3d0e1a5c.radex"
ALL_FILES = [*sorted(LAMDA.glob("*.dat")), EMAA_CO]


def same_data(a: LamdaData, b: LamdaData) -> bool:
    """All values equal (NaN = NaN); notes, partner text and source map excepted."""
    arrays = ("level_energy", "level_g", "line_up", "line_low", "line_a", "line_freq", "line_eup")
    if any(not np.array_equal(getattr(a, n), getattr(b, n), equal_nan=True) for n in arrays):
        return False
    scalars = ("molecule", "mass", "level_quanta", "title", "dataset_doi", "level_type", "level_reference")
    if any(getattr(a, n) != getattr(b, n) for n in scalars if n != "mass"):
        return False
    if not np.array_equal(a.mass, b.mass, equal_nan=True) or len(a.collisions) != len(b.collisions):
        return False
    for c, d in zip(a.collisions, b.collisions, strict=True):
        if (c.partner_id, c.reference) != (d.partner_id, d.reference):
            return False
        for n in ("temperatures", "up", "low", "rates"):
            if not np.array_equal(getattr(c, n), getattr(d, n), equal_nan=True):
                return False
    return True


def make_lamda(
    *,
    molecule: str = "X",
    mass: str = "28.0",
    levels: str = "1 0.0 1.0 J=0\n2 3.845033 3.0 J=1\n3 11.534919 5.0 J=2",
    nlev: str = "3 - ROTATION - Ref: CDMS https://cdms.astro.uni-koeln.de",
    lines: str = "1 2 1 7.203e-08 115.2711894 5.532\n2 3 2 6.910e-07 230.5369826 16.596",
    nlin: str = "2",
    partner: str = "2 X + p-H2 ! Yang et al. (2010) https://doi.org/10.1088/0004-637X/718/2/1062",
    ncol: str = "3",
    ntemp: str = "2",
    temps: str = "10.0 20.0",
    rates: str = "1 2 1 1.0e-11 2.0e-11\n2 3 1 3.0e-11 4.0e-11\n3 3 2 5.0e-11 6.0e-11",
    ncol_comment: str = "!NUMBER OF COLL TRANS",
    npart: str = "1",
    after: str = "",
) -> str:
    """A small synthetic LAMDA file (3 levels, 2 lines, one partner); arguments replace parts.

    With the defaults, :func:`pylamda.verify` finds no issue.
    """
    text = f"""!MOLECULE
{molecule}
!MOLECULAR WEIGHT
{mass}
!NUMBER OF ENERGY LEVELS
{nlev}
!LEVEL + ENERGIES(cm^-1) + WEIGHT + J
{levels}
!NUMBER OF RADIATIVE TRANSITIONS
{nlin}
!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K)
{lines}
!NUMBER OF COLL PARTNERS
{npart}
!COLLISIONS BETWEEN
{partner}
{ncol_comment}
{ncol}
!NUMBER OF COLL TEMPS
{ntemp}
!COLL TEMPS
{temps}
!TRANS + UP + LOW + COLLRATES(cm^3 s^-1)
{rates}
"""
    return text + after


def collision_block(
    partner: str, rates: str = "1 2 1 1.0e-11 2.0e-11", ntemp: str = "2", temps: str = "10.0 20.0"
) -> str:
    """One more collision block for ``make_lamda(after=...)`` (raise ``npart`` with it)."""
    ncol = len(rates.splitlines())
    return (
        f"!COLLISIONS BETWEEN\n{partner}\n!NUMBER OF COLL TRANS\n{ncol}\n!NUMBER OF COLL TEMPS\n{ntemp}\n"
        f"!COLL TEMPS\n{temps}\n!TRANS + UP + LOW + COLLRATES(cm^3 s^-1)\n{rates}\n"
    )
