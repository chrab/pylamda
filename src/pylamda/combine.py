"""Combine LAMDA data sets: add, drop or scale the collision data of one partner (design document, section 9).

:func:`add_partner` comes from astrolinedata (``formats/lamda.py``), adapted to pylamda:
partners are identified by their ID, and warnings go into the record instead of a log.
:func:`scale_partner` multiplies the rates of one partner by a factor (e.g. He rates
given as H2 rates). All combine functions return the new data and a record, and add a
note line.
"""

from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass, field, replace

import numpy as np

from .compare import match_levels
from .model import PARTNER_IDS, CollisionSet, LamdaData, partner_name
from .verify import _SCALED

#: Default tolerances for matching levels in :func:`add_partner`: energy within
#: ``ADD_PARTNER_ATOL`` [cm-1] + ``ADD_PARTNER_RTOL_E`` * E. Looser than in
#: :func:`pylamda.compare`, because different sources give the same level with slightly
#: different energies (O I 3P1: LAMDA 158.2687, CHIANTI 158.265 cm-1).
ADD_PARTNER_ATOL = 0.1
#: See :data:`ADD_PARTNER_ATOL`.
ADD_PARTNER_RTOL_E = 1e-4
#: :func:`scale_partner`: significant digits of the scaled rates if the data have no
#: source map (else the digits of the rates as written, at least
#: :data:`SCALE_MIN_DIGITS`).
SCALE_DEFAULT_DIGITS = 6
#: See :data:`SCALE_DEFAULT_DIGITS`.
SCALE_MIN_DIGITS = 3
#: :func:`scale_partner`: relative tolerance for a factor written in the partner text.
SCALE_TEXT_RTOL = 1e-3

_H2_LIKE = ("H2", "p-H2", "o-H2")
# a scaling factor in a partner text: "*1.37", "(*1.363)", "x 1.4", "×1.37", "by 1.39"
_TEXT_FACTOR = re.compile(r"(?:\*|×|\bx|\bby)\s*\(?\s*(\d+\.\d+)(?![\d.])", re.I)
_MANTISSA = re.compile(r"[+-]?(\d*)(?:\.(\d*))?(?:[eEdD][+-]?\d+)?")


class CombineError(ValueError):
    """Data sets cannot be combined (missing partner, levels that do not fit, ...)."""


@dataclass
class AddedPartner:
    """What :func:`add_partner` did.

    Attributes
    ----------
    partner : str
        Name of the partner added (the standard name of its ID).
    partner_id : int
        Its ID.
    level_matching : str
        How the levels were matched (:attr:`LevelMatch.method <pylamda.compare.LevelMatch.method>`).
    level_map : dict of int to int
        Level number in the source -> level number in the target (1-based), for the
        levels used by the added transitions.
    n_transitions : int
        Number of collisional transitions added.
    n_dropped : int
        Transitions of the source not added, because a level is not in the target.
    max_energy_difference : float
        Largest energy difference [cm-1] between matched levels used.
    reference : str
        Reference of the partner in the source.
    source : str
        Description of the source (``source_label``, the source's ``!SOURCE:`` text, or
        its molecule name).
    warnings : list of str
        Things to check (other molecule name, H2 next to p-H2/o-H2).
    """

    partner: str
    partner_id: int
    level_matching: str
    level_map: dict[int, int]
    n_transitions: int
    n_dropped: int
    max_energy_difference: float
    reference: str
    source: str
    warnings: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        """The record as a JSON-compatible dict (``level_map`` keys as strings)."""
        d = asdict(self)
        d["level_map"] = {str(k): v for k, v in self.level_map.items()}
        return d


@dataclass
class DroppedPartner:
    """What :func:`drop_partner` did.

    Attributes
    ----------
    partner : str
        Name of the partner dropped (the standard name of its ID, or ``"ID <n>"``).
    partner_id : int
        Its ID.
    block : int
        The number of the block dropped (1-based, in the data before).
    n_transitions : int
        Number of collisional transitions dropped.
    reference : str
        Reference of the partner.
    """

    partner: str
    partner_id: int
    block: int
    n_transitions: int
    reference: str

    def as_dict(self) -> dict:
        """The record as a JSON-compatible dict."""
        return asdict(self)


@dataclass
class ScaledPartner:
    """What :func:`scale_partner` did.

    Attributes
    ----------
    partner : str
        Name of the partner scaled (the standard name of its ID).
    partner_id : int
        Its ID.
    block : int
        The number of its block (1-based, in the data before).
    new_partner : str
        Name of the partner of the scaled block (``partner`` if it was scaled in place).
    new_partner_id : int
        Its ID.
    new_block : int
        The number of the scaled block (1-based, in the data after).
    factor : float
        The factor the rates were multiplied by.
    n_transitions : int
        Number of collisional transitions scaled.
    kept_source : bool
        The source block was kept and the scaled block added (``keep=True``).
    replaced_block : int or None
        The block of ``new_partner`` that was replaced (1-based, in the data before), or
        None.
    reference : str
        The reference text of the scaled block.
    text_factor : float or None
        A scaling factor written in the partner text of the source (``*1.37``,
        ``by 1.39``, ...), or None. Only used for a warning.
    digits : int
        Significant digits of the scaled rates.
    warnings : list of str
        Things to check (partner text and ID disagree, another factor in the text, an
        unusual direction or partner pair, H2 next to p-H2/o-H2, a reference text that
        mentions scaling).
    """

    partner: str
    partner_id: int
    block: int
    new_partner: str
    new_partner_id: int
    new_block: int
    factor: float
    n_transitions: int
    kept_source: bool
    replaced_block: int | None
    reference: str
    text_factor: float | None
    digits: int
    warnings: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        """The record as a JSON-compatible dict."""
        return asdict(self)


def _h2_warning(molecule: str, key: str, others: list[str | None]) -> str | None:
    """Check C04 for a new block of ``key``: H2 next to p-H2/o-H2 among the ``others`` partners."""
    if key not in _H2_LIKE:
        return None
    other = ("p-H2", "o-H2") if key == "H2" else ("H2",)
    present = [p for p in other if p in others]
    if not present:
        return None
    return f"{molecule}: {key} next to {'/'.join(present)}; a code may count H2 collisions twice"


def _blocks(data: LamdaData, name: str) -> list[int]:
    """Indices (0-based) of the blocks of partner ``name`` (any spelling; by ID)."""
    key = partner_name(name) or name
    return [b for b, c in enumerate(data.collisions) if c.partner == key]


def _one_block(data: LamdaData, name: str, block: int | None, what: str) -> int:
    """The index (0-based) of the block of ``name``; ``block`` (1-based) chooses among several."""
    blocks = _blocks(data, name)
    if block is not None:
        if not 1 <= block <= len(data.collisions):
            raise CombineError(f"{what} has no collision block {block} (it has {len(data.collisions)})")
        if block - 1 not in blocks:
            raise CombineError(f"collision block {block} of {what} is {data.collisions[block - 1].label}, not {name}")
        return block - 1
    if not blocks:
        have = ", ".join(c.label for c in data.collisions)
        raise CombineError(f"{what} ({data.molecule}): no collision partner {name!r} (have: {have})")
    if len(blocks) > 1:
        numbers = ", ".join(str(b + 1) for b in blocks)
        raise CombineError(
            f"{what} ({data.molecule}) has {len(blocks)} blocks with the ID of {name} (blocks {numbers}; check C03, "
            "C11): choose one with block="
        )
    return blocks[0]


def add_partner(
    target: LamdaData,
    source: LamdaData,
    partner: str,
    *,
    match: str = "auto",
    atol: float = ADD_PARTNER_ATOL,
    rtol_e: float = ADD_PARTNER_RTOL_E,
    replace_existing: bool = False,
    source_label: str | None = None,
    block: int | None = None,
) -> tuple[LamdaData, AddedPartner]:
    """Add the collision rates of one partner from another data set of the same species.

    For example the electron rates of a CHIANTI conversion to a LAMDA file that has
    none. The levels of ``source`` are matched to those of ``target``; the collisional
    transitions between matched levels are added with the level numbers of ``target``.
    Temperatures and rate coefficients are copied unchanged. A note line records where
    the partner comes from.

    Parameters
    ----------
    target : LamdaData
        The data set to extend (not changed).
    source : LamdaData
        The data set with the partner.
    partner : str
        Partner name (any spelling of :func:`pylamda.partner_name`, e.g. ``"e"``,
        ``"pH2"``); partners are identified by their ID.
    match : str, optional
        How levels are matched (:func:`pylamda.match_levels`): ``"auto"`` (default),
        ``"quanta"``, ``"energy"`` (within ``atol + rtol_e * E`` and equal weight g) or
        ``"index"``.
    atol : float, optional
        Energy tolerance [cm-1] (default :data:`ADD_PARTNER_ATOL`).
    rtol_e : float, optional
        Relative energy tolerance (default :data:`ADD_PARTNER_RTOL_E`).
    replace_existing : bool, optional
        Replace the partner if ``target`` has it already (default: error).
    source_label : str or None, optional
        Description of the source for the note and :attr:`AddedPartner.source`.
        Default: the source's ``source`` text, else its molecule name.
    block : int or None, optional
        The collision block of ``source`` (1-based), if it has several blocks with the
        partner's ID (EMAA files).

    Returns
    -------
    data : LamdaData
        A copy of ``target`` with the partner appended (or replaced in place) and a
        note line.
    added : AddedPartner
        What was added.

    Raises
    ------
    CombineError
        If ``source`` has no such partner (or several blocks with its ID and no
        ``block``), ``target`` has it and ``replace_existing`` is False, no transition of
        the partner connects two matched levels, a matched level has another weight g,
        or a transition would point from the lower to the higher energy in ``target``
        (the two data sets order these levels differently). A subclass of ValueError.
    """
    # the target first (as in astrolinedata), then the source
    existing = _blocks(target, partner)
    if existing and not replace_existing:
        key = partner_name(partner) or partner
        raise CombineError(f"{target.molecule} already has collision partner {key} (use replace_existing)")
    cs = source.collisions[_one_block(source, partner, block, "the source")]
    key = cs.partner
    if key is None:
        raise CombineError(f"the source block has the ID {cs.partner_id}, which is not a standard partner")
    existing = _blocks(target, key)
    if len(existing) > 1:
        raise CombineError(f"{target.molecule} has {len(existing)} blocks of {key}; replace one by hand")
    warnings = []
    if source.molecule != target.molecule:
        warnings.append(f"adding {key} rates of {source.molecule!r} to {target.molecule!r}: the molecule names differ")
    h2 = _h2_warning(target.molecule, key, target.partners)
    if h2:
        warnings.append(h2)

    levels = match_levels(source, target, match, atol, rtol_e)
    level_map = levels.mapping  # source -> target, 1-based
    rows, used = [], set()
    for k, (u, lo) in enumerate(zip(cs.up, cs.low, strict=True)):
        iu, il = level_map.get(int(u)), level_map.get(int(lo))
        if iu is None or il is None:
            continue
        for i_src, i_tgt in ((int(u), iu), (int(lo), il)):
            if source.level_g[i_src - 1] != target.level_g[i_tgt - 1]:
                raise CombineError(
                    f"level {i_src} of the source (g={source.level_g[i_src - 1]:g}) matched level {i_tgt} of the "
                    f"target (g={target.level_g[i_tgt - 1]:g})"
                )
        if target.level_energy[iu - 1] < target.level_energy[il - 1]:
            raise CombineError(
                f"{key} transition {u}->{lo} of the source would go from target level {iu} to the higher level {il}; "
                "the data sets order these levels differently"
            )
        rows.append((k, iu, il))
        used.update(((int(u), iu), (int(lo), il)))
    if not rows:
        raise CombineError(
            f"no {key} transition of the source connects two levels of the target "
            f"(levels matched by {levels.method}: {levels.n_matched})"
        )
    idx = [k for k, _, _ in rows]
    pid = PARTNER_IDS[key]
    new = CollisionSet(
        partner_id=pid,
        temperatures=np.array(cs.temperatures, dtype=float),
        up=np.array([u for _, u, _ in rows], dtype=int),
        low=np.array([lo for _, _, lo in rows], dtype=int),
        rates=np.array(cs.rates[idx], dtype=float),
        reference=cs.reference,
        partner_text=f"{pid} {target.molecule} + {key}" + (f" -- Ref: {cs.reference}" if cs.reference else ""),
    )
    label = source_label or source.source or source.molecule
    max_de = max(abs(float(source.level_energy[i - 1] - target.level_energy[j - 1])) for i, j in used)
    added = AddedPartner(
        partner=key,
        partner_id=pid,
        level_matching=levels.method,
        level_map=dict(sorted(used)),
        n_transitions=len(rows),
        n_dropped=cs.n_transitions - len(rows),
        max_energy_difference=max_de,
        reference=cs.reference,
        source=label,
        warnings=warnings,
    )
    if existing:
        collisions = [new if b == existing[0] else c for b, c in enumerate(target.collisions)]
    else:
        collisions = [*target.collisions, new]
    note = (
        f"! Collision partner {key} added from {label}: {added.n_transitions} of {cs.n_transitions} transitions, "
        f"levels matched by {levels.method} (max |dE| {max_de:.3g} cm-1)"
    )
    return replace(target, collisions=collisions, notes=[*target.notes, note], source_map=None), added


def drop_partner(data: LamdaData, partner: str, *, block: int | None = None) -> tuple[LamdaData, DroppedPartner]:
    """Remove the collision data of one partner.

    Parameters
    ----------
    data : LamdaData
        The data (not changed).
    partner : str
        Partner name (any spelling; partners are identified by their ID), or
        ``"ID <n>"`` for an ID outside the standard.
    block : int or None, optional
        The collision block (1-based), if the data have several blocks with the
        partner's ID (EMAA files).

    Returns
    -------
    data : LamdaData
        A copy without the block, with a note line.
    dropped : DroppedPartner
        What was removed.

    Raises
    ------
    CombineError
        If there is no such partner, or several blocks with its ID and no ``block``.
    """
    if partner.startswith("ID ") and block is None:
        blocks = [b for b, c in enumerate(data.collisions) if c.label == partner]
        if len(blocks) != 1:
            raise CombineError(f"{data.molecule}: {len(blocks)} blocks with {partner}; choose one with block=")
        b = blocks[0]
    else:
        b = _one_block(data, partner, block, "the data")
    c = data.collisions[b]
    dropped = DroppedPartner(c.label, c.partner_id, b + 1, c.n_transitions, c.reference)
    note = f"! Collision partner {c.label} (block {b + 1}, {c.n_transitions} transitions) dropped"
    collisions = [x for k, x in enumerate(data.collisions) if k != b]
    return replace(data, collisions=collisions, notes=[*data.notes, note], source_map=None), dropped


def _significant_digits(token: str) -> int:
    """Number of significant digits of a number as written (``"5.47E-10"`` -> 3; 0 if not a number)."""
    m = _MANTISSA.fullmatch(token)
    if not m:
        return 0
    return len(((m.group(1) or "") + (m.group(2) or "")).lstrip("0"))


def _block_digits(data: LamdaData, b: int) -> int | None:
    """The largest number of significant digits of the rates of block ``b`` as written.

    None if the data have no source map, or if it does not fit the block any more.
    """
    smap = data.source_map
    c = data.collisions[b]
    if smap is None or len(smap.blocks) != len(data.collisions) or len(smap.blocks[b].rows) != c.n_transitions:
        return None
    ntemp = len(c.temperatures)
    best = 0
    for first in smap.blocks[b].rows:
        tokens: list[str] = []
        line = int(first)
        while len(tokens) < 3 + ntemp and line <= len(smap.raw):
            text = smap.raw[line - 1]
            if text.strip() and not text.lstrip().startswith("!"):
                tokens += text.split()
            line += 1
        best = max([best, *(_significant_digits(t) for t in tokens[3 : 3 + ntemp])])
    return best or None


def _round_significant(values: np.ndarray, digits: int) -> np.ndarray:
    """``values`` rounded to ``digits`` significant digits (non-finite values unchanged)."""
    flat = [float(f"{v:.{digits - 1}e}") if math.isfinite(v) else v for v in values.ravel().tolist()]
    return np.array(flat, dtype=float).reshape(values.shape)


def _scale_warnings(
    data: LamdaData, cs: CollisionSet, key: str, new_key: str, factor: float, text_factor: float | None
) -> list[str]:
    """The warnings of :func:`scale_partner` about the source block and the partner pair."""
    warnings = []
    written = cs.text_partner
    canonical = partner_name(written) if written else None
    if written and canonical != key:
        warnings.append(
            f"the partner text {cs.partner_text.strip()!r} names {canonical or written!r}; the block was taken by "
            f"its ID {cs.partner_id} as {key} (checks C02, C11, C12)"
        )
    if text_factor is not None:
        matches = any(abs(f / text_factor - 1) <= SCALE_TEXT_RTOL for f in (factor, 1 / factor))
        if not matches:
            warnings.append(
                f"the partner text gives the factor {text_factor:g}; the factor used is {factor:.6g} "
                f"(neither {text_factor:g} nor 1/{text_factor:g})"
            )
    if key == "He" and new_key in _H2_LIKE and factor < 1:
        warnings.append(f"He -> {new_key} with a factor < 1: He rates scaled to H2 are usually larger (about 1.37)")
    elif key in _H2_LIKE and new_key == "He" and factor > 1:
        warnings.append(f"{key} -> He with a factor > 1: He rates are usually smaller than H2 rates (about 1/1.37)")
    elif key != new_key and not {key, new_key} <= {"He", *_H2_LIKE}:
        warnings.append(f"{key} -> {new_key}: scaling is a standard approximation only between He and H2 (p-H2, o-H2)")
    return warnings


def scale_partner(
    data: LamdaData,
    partner: str,
    factor: float,
    *,
    to: str | None = None,
    keep: bool = False,
    block: int | None = None,
    replace_existing: bool = False,
    reference: str | None = None,
    digits: int | None = None,
) -> tuple[LamdaData, ScaledPartner]:
    """Multiply the collision rates of one partner by a factor, optionally as another partner.

    Rates of one partner are often scaled rates of another one: He rates multiplied by
    about 1.37 (the square root of the ratio of the reduced masses) and given as H2 rates.
    Three cases:

    * scale in place: ``scale_partner(data, "H2", 1.37 / 1.39)`` (same ID);
    * change the partner: ``scale_partner(data, "H2", 1 / 1.37, to="He")`` (the block
      gets the ID of ``to``);
    * derive a new partner: ``scale_partner(data, "He", 1.37, to="H2", keep=True)``
      (the source block stays, the scaled block is appended).

    Partners are identified by their ID. The scaled block gets the normalised partner
    line ``<id> <molecule> + <partner> -- Ref: <reference>``; a note line records the
    operation. The factor is never taken from the partner text: a factor written there
    only gives a warning if it differs.

    Parameters
    ----------
    data : LamdaData
        The data (not changed).
    partner : str
        The partner to scale (any spelling of :func:`pylamda.partner_name`).
    factor : float
        The factor the rates are multiplied by (finite, > 0).
    to : str or None, optional
        The partner of the scaled block. Default: ``partner`` (scale in place).
    keep : bool, optional
        Keep the source block and append the scaled block (``to`` must be another
        partner). Default: False, the scaled block takes the place of the source block.
    block : int or None, optional
        The block of ``partner`` (1-based), if the data have several blocks with its ID
        (EMAA files).
    replace_existing : bool, optional
        Replace a block of ``to`` that the data have already. With ``keep`` the scaled
        block takes its place; else the scaled block stays at the place of the source
        block and the block of ``to`` is dropped. Default: error.
    reference : str or None, optional
        The reference text of the scaled block. Default: the source's text, unchanged
        (a warning if it mentions scaling, as it may no longer describe the rates).
    digits : int or None, optional
        Significant digits of the scaled rates. Default: as many as the rates of the
        source as written (at least :data:`SCALE_MIN_DIGITS`), or
        :data:`SCALE_DEFAULT_DIGITS` if the data have no source map.

    Returns
    -------
    data : LamdaData
        A copy with the scaled block and a note line.
    scaled : ScaledPartner
        What was done, with warnings.

    Raises
    ------
    CombineError
        If ``factor`` is not finite and > 0, ``digits`` < 1, ``partner`` or ``to`` is not
        a standard partner, the data have no such partner (or several blocks with its ID
        and no ``block``), ``keep`` is set and ``to`` is the source partner, or ``to``
        exists already and ``replace_existing`` is False. A subclass of ValueError.
    """
    if not (math.isfinite(factor) and factor > 0):
        raise CombineError(f"the factor must be finite and > 0, not {factor!r}")
    if digits is not None and digits < 1:
        raise CombineError(f"digits must be at least 1, not {digits}")
    key = partner_name(partner)
    if key not in PARTNER_IDS:
        raise CombineError(f"{partner!r} is not a standard collision partner ({', '.join(PARTNER_IDS)})")
    new_key = key if to is None else partner_name(to)
    if new_key not in PARTNER_IDS:
        raise CombineError(f"{to!r} is not a standard collision partner ({', '.join(PARTNER_IDS)})")
    if keep and new_key == key:
        raise CombineError(f"keep=True needs another partner than {key} (two {key} blocks would be check C03)")
    b = _one_block(data, key, block, "the data")
    cs = data.collisions[b]
    # scaling in place adds no block of another partner (repeated IDs stay as they are: C03)
    existing = [k for k in _blocks(data, new_key) if k != b] if new_key != key else []
    if len(existing) > 1:
        raise CombineError(f"{data.molecule} has {len(existing)} blocks of {new_key}; drop them first")
    if existing and not replace_existing:
        raise CombineError(f"{data.molecule} already has collision partner {new_key} (use replace_existing)")

    text = (cs.partner_text or cs.reference).strip()
    m = _TEXT_FACTOR.search(text)
    text_factor = float(m.group(1)) if m else None
    warnings = _scale_warnings(data, cs, key, new_key, factor, text_factor)
    ref = cs.reference if reference is None else reference.strip()
    if reference is None and _SCALED.search(ref):
        warnings.append(
            f"the reference text {ref!r} mentions scaling; it may no longer describe the rates (give a new text "
            "with reference=)"
        )
    if digits is None:
        written = _block_digits(data, b)
        digits = max(written, SCALE_MIN_DIGITS) if written else SCALE_DEFAULT_DIGITS

    pid = PARTNER_IDS[new_key]
    new = CollisionSet(
        partner_id=pid,
        temperatures=np.array(cs.temperatures, dtype=float),
        up=np.array(cs.up, dtype=int),
        low=np.array(cs.low, dtype=int),
        rates=_round_significant(np.asarray(cs.rates, dtype=float) * factor, digits),
        reference=ref,
        partner_text=f"{pid} {data.molecule} + {new_key}" + (f" -- Ref: {ref}" if ref else ""),
    )
    replaced = existing[0] if existing else None
    if keep:
        collisions = list(data.collisions)
        if replaced is None:
            collisions.append(new)
            new_index = len(collisions) - 1
        else:
            collisions[replaced] = new
            new_index = replaced
    else:
        collisions = [new if k == b else c for k, c in enumerate(data.collisions) if k != replaced]
        new_index = b - 1 if replaced is not None and replaced < b else b
    h2 = _h2_warning(data.molecule, new_key, [c.partner for k, c in enumerate(collisions) if k != new_index])
    if h2:
        warnings.append(h2)

    scaled = ScaledPartner(
        partner=key,
        partner_id=cs.partner_id,
        block=b + 1,
        new_partner=new_key,
        new_partner_id=pid,
        new_block=new_index + 1,
        factor=float(factor),
        n_transitions=cs.n_transitions,
        kept_source=keep,
        replaced_block=None if replaced is None else replaced + 1,
        reference=ref,
        text_factor=text_factor,
        digits=digits,
        warnings=warnings,
    )
    note = f"! Collision partner {key} (block {b + 1}, {cs.n_transitions} transitions) scaled by {factor:.6g}"
    if keep:
        note += f" and added as {new_key} (block {new_index + 1})"
    elif new_key != key:
        note += f" and relabelled {new_key}"
    if replaced is not None:
        note += f"; the former {new_key} block (block {replaced + 1}) replaced"
    return replace(data, collisions=collisions, notes=[*data.notes, note], source_map=None), scaled
