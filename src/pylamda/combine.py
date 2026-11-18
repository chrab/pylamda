"""Combine LAMDA data sets: add or drop the collision data of one partner (design document, section 9).

:func:`add_partner` comes from astrolinedata (``formats/lamda.py``), adapted to pylamda:
partners are identified by their ID, and warnings go into the record instead of a log.
All combine functions return the new data and a record, and add a note line.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace

import numpy as np

from .compare import match_levels
from .model import PARTNER_IDS, CollisionSet, LamdaData, partner_name

#: Default tolerances for matching levels in :func:`add_partner`: energy within
#: ``ADD_PARTNER_ATOL`` [cm-1] + ``ADD_PARTNER_RTOL_E`` * E. Looser than in
#: :func:`pylamda.compare`, because different sources give the same level with slightly
#: different energies (O I 3P1: LAMDA 158.2687, CHIANTI 158.265 cm-1).
ADD_PARTNER_ATOL = 0.1
#: See :data:`ADD_PARTNER_ATOL`.
ADD_PARTNER_RTOL_E = 1e-4


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
    cs = source.collisions[_one_block(source, partner, block, "the source")]
    key = cs.partner
    if key is None:
        raise CombineError(f"the source block has the ID {cs.partner_id}, which is not a standard partner")
    existing = _blocks(target, key)
    if existing and not replace_existing:
        raise CombineError(f"{target.molecule} already has collision partner {key} (use replace_existing)")
    if len(existing) > 1:
        raise CombineError(f"{target.molecule} has {len(existing)} blocks of {key}; replace one by hand")
    warnings = []
    if source.molecule != target.molecule:
        warnings.append(f"adding {key} rates of {source.molecule!r} to {target.molecule!r}: the molecule names differ")
    if key in ("H2", "p-H2", "o-H2"):
        other = ("p-H2", "o-H2") if key == "H2" else ("H2",)
        present = [p for p in other if p in target.partners]
        if present:
            warnings.append(
                f"{target.molecule}: {key} next to {'/'.join(present)}; a code may count H2 collisions twice"
            )

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
