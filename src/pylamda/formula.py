"""Formula mass of a molecule name as written in LAMDA files (for check P04)."""

from __future__ import annotations

import re

#: Standard atomic weight [Da] of the elements H to Zn: the abridged values of IUPAC 2021
#: (Prohaska et al. 2022, Pure Appl. Chem. 94, 573, Table 2; for elements with an
#: interval, the conventional value). Copied from astrolinedata (``species.py``).
#:
#: :meta hide-value:
STANDARD_ATOMIC_WEIGHT: dict[str, float] = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305,
    "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.95,
    "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996,
    "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
}  # fmt: skip

# deuterium is written as its own symbol (2H, 2.014 Da)
_EXTRA = {"D": 2.014}
# spin isomers and structural prefixes: o-H2O, p-NH3, A-CH3OH, E-CH3OH, c-C3H2, l-C3H
_PREFIX = re.compile(r"^(?:ortho|para|meta|[opmAEcl])-")
_ATOM = re.compile(r"(\d+)?([A-Z][a-z]?)(\d*)")


def formula_mass(name: str) -> float | None:
    """Mass [Da] of a molecule from its name (``"CO"``, ``"13CO"``, ``"C18O"``, ``"o-H3O+"``).

    The name is read as a chemical formula, after an optional prefix (``o-``, ``p-``,
    ``A-``, ``E-``, ``c-``, ``l-``, ...) and before an optional charge (``+``, ``++``,
    ``-``). A number of 10 or more before an element is a mass number (``13C``, ``18O``,
    also in ``C18O`` and ``H13CN``); a number below 10 after an element is a count
    (``H2``, ``HC3N``). An isotope counts with its mass number, the other elements with
    their standard atomic weight. The charge is ignored (the electron mass is below the
    precision of the check).

    Parameters
    ----------
    name : str
        The molecule name, e.g. the first word of line 2 of a LAMDA file.

    Returns
    -------
    float or None
        The mass [Da], or None if the name does not read as a formula of known elements.
    """
    text = _PREFIX.sub("", name.strip())
    text = re.sub(r"[+-]+$", "", text)
    if not text:
        return None
    mass = 0.0
    pending: int | None = None  # a mass number written after the previous element (C18O)
    pos = 0
    while pos < len(text):
        m = _ATOM.match(text, pos)
        if not m:
            return None
        iso_text, element, count_text = m.groups()
        weight = STANDARD_ATOMIC_WEIGHT.get(element, _EXTRA.get(element))
        if weight is None:
            return None
        if iso_text and pending is not None:
            return None
        iso = int(iso_text) if iso_text else pending
        pending = None
        count = int(count_text) if count_text else 1
        if count >= 10 and m.end() < len(text):  # "C18O": 18 is the mass number of the next atom
            pending, count = count, 1
        if count < 1:
            return None
        if iso is not None and not (0.5 * weight <= iso <= 2 * weight + 2):
            return None
        mass += count * (float(iso) if iso is not None else weight)
        pos = m.end()
    return None if pending is not None else mass
