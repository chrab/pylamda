"""Maintainer tool: build the reference catalogue (design document, section 10.1).

Not part of the installed package. Steps:

1. ``search``: collect the free-text references (level and partner references without
   a DOI, bibcode or arXiv ID) of LAMDA files, one item per cited paper; group equal
   texts of the same molecule (one decision per group); collect hints (the file's
   comments, the LAMDA molecule page, the EMAA species page); search SciX; score the
   candidates; pre-fill a ``proposal`` where the match is clear. All goes into a review
   file (JSON).
2. A person reviews each group and sets ``decision``:

   * ``"ok"``: the proposal is right;
   * a bibcode: this paper (also one that is not among the candidates);
   * ``"none"``: the paper is not in SciX;
   * ``"skip"``: the text is not a reference.

   ``show`` prints the review file; ``confirm`` sets ``"ok"`` for all proposals at once
   (after checking them with ``show --proposals``).
3. ``build``: the decided groups become ``pylamda/data/references.json`` (entries and
   lookup rows) and ``references.bib`` (BibTeX from SciX). Keys stay the same across
   rebuilds.

Usage::

    conda activate PYLAMDA        # SCIX_API_TOKEN set
    python tools/catalogue.py search --database prodimo data/CO_Lambda.dat -o review.json
    python tools/catalogue.py show review.json [--open | --proposals]
    python tools/catalogue.py confirm review.json
    python tools/catalogue.py build review.json
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

from pylamda.bib import (
    SCIX_ABS,
    HttpClient,
    _set_field,
    _set_key,
    entry_field,
    entry_key,
    key_from_entry,
    scix_token,
    split_entries,
)
from pylamda.catalogue import JOURNALS, normalise, parse, split_segments  # noqa: F401 (used by the tests)
from pylamda.refsearch import (  # noqa: F401 (also used by the tests)
    BEST_MARGIN,
    BEST_MIN_SCORE,
    BEST_REASON,
    ELEMENTS,
    EMAA_URL,
    FIRST_PASS,
    LAMDA_URL,
    Group,
    Use,
    add_hints,
    ascii_fold,
    collect,
    file_hints,
    hint_matches,
    is_first_author,
    load,
    merge_review,
    name_forms,
    normalise_title,
    page_hints,
    page_names,
    propose,
    save,
    score,
    search,
    show,
    species_forms,
    title_match,
)

# ------------------------------------------------------------------- build
DATA_DIR = Path(__file__).parents[1] / "src" / "pylamda" / "data"
CATALOGUE_VERSION = 1


def build(groups: list[Group], client, out_dir: Path = DATA_DIR) -> tuple[int, int]:
    """Write ``references.json`` and ``references.bib`` from the decided groups.

    Returns the number of entries and of lookup rows. Keys of an existing
    ``references.json`` are kept for the same bibcode.
    """
    from datetime import date

    json_path, bib_path = out_dir / "references.json", out_dir / "references.bib"
    old_keys = {}
    if json_path.exists():
        old = json.loads(json_path.read_text(encoding="utf-8"))
        old_keys = {e["bibcode"]: k for k, e in old.get("entries", {}).items()}
    decided = [g for g in groups if g.decision and g.decision != "skip"]
    bibcodes = sorted({g.bibcode for g in decided if g.bibcode})
    exported: dict[str, str] = {}
    for i in range(0, len(bibcodes), 200):
        for e in split_entries(client.scix_export(bibcodes[i : i + 200])):
            exported[entry_key(e)] = e
    missing = [b for b in bibcodes if b not in exported]
    if missing:
        raise SystemExit(f"catalogue.py: not found in SciX: {', '.join(missing)} (check the decisions)")
    keys: dict[str, str] = {}
    used = set(old_keys.values()) & {old_keys.get(b) for b in bibcodes}
    for b in bibcodes:
        if b in old_keys:
            keys[b] = old_keys[b]
            continue
        base = key_from_entry(exported[b]) or re.sub(r"[^A-Za-z0-9]", "", b)
        key = base
        for suffix in "abcdefghijklmnopqrstuvwxyz":
            if key not in used:
                break
            key = base + suffix
        keys[b] = key
        used.add(key)
    entries, bib = {}, []
    for b in sorted(bibcodes, key=lambda b: keys[b]):
        e = _set_field(_set_key(exported[b], keys[b]), "adsurl", f"{SCIX_ABS}{b}")
        entries[keys[b]] = {"bibcode": b, "doi": entry_field(e, "doi"), "reviewed": True}
        bib.append(e)
    rows = []
    for g in decided:
        for u in g.uses:
            rows.append({
                "database": u.database, "molecule": g.molecule, "file": u.file, "role": u.role,
                "partner": u.partner, "text": g.normalised, "key": keys.get(g.bibcode) if g.bibcode else None,
                "status": "found" if g.bibcode else "not in SciX",
            })  # fmt: skip
    rows.sort(key=lambda r: (r["database"], r["file"], r["role"], r["partner"] or "", r["text"]))
    catalogue = {
        "description": "Reviewed free-text references of LAMDA files (pylamda design document, section 10.1). "
        "entries: key -> bibcode, DOI. lookup: one row per place a reference text occurs (text normalised: "
        "lower case, punctuation as blanks); key null = reviewed, not in SciX. Built by tools/catalogue.py.",
        "version": CATALOGUE_VERSION,
        "updated": date.today().isoformat(),
        "entries": entries,
        "lookup": rows,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(catalogue, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    head = "% Reviewed references of LAMDA files (pylamda reference catalogue); BibTeX from SciX.\n"
    bib_path.write_text(head + "\n" + "\n\n".join(bib) + "\n", encoding="utf-8")
    return len(entries), len(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="catalogue.py", description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    s = sub.add_parser("search", help="collect the free-text references, hints and SciX candidates")
    s.add_argument("files", nargs="+")
    s.add_argument("--database", required=True, help="lamda, emaa, prodimo, ...")
    s.add_argument("-o", "--output", required=True, help="review file (JSON); decisions in it are kept")
    s.add_argument("--lamda-page", help="LAMDA molecule page to use (default: guessed from the molecule)")
    s.add_argument("--no-hints", action="store_true", help="do not read the LAMDA/EMAA pages")
    s.add_argument(
        "--cache-dir",
        default=str(Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "pylamda" / "catalogue"),
        help="cache of the database pages",
    )
    v = sub.add_parser("show", help="show a review file")
    v.add_argument("review")
    which = v.add_mutually_exclusive_group()
    which.add_argument("--open", action="store_true", help="only groups without a decision")
    which.add_argument("--proposals", action="store_true", help="only proposals without a decision")
    which.add_argument("--best", action="store_true", help='only "best candidate" proposals without a decision')
    which.add_argument("--no-proposal", action="store_true", help="only groups without a proposal and a decision")
    which.add_argument("--first-pass", action="store_true", help="only the proposals and skips of a first pass")
    c = sub.add_parser("confirm", help='set decision "ok" for every proposal without a decision')
    c.add_argument("review")
    c.add_argument("--include-best", action="store_true", help='also the weaker "best candidate" proposals')
    c.add_argument("--include-first-pass", action="store_true", help="also the proposals of an assistant's first pass")
    bld = sub.add_parser("build", help="write references.json and references.bib from the decided groups")
    bld.add_argument("review")
    bld.add_argument("--out-dir", default=str(DATA_DIR), help="default: src/pylamda/data")
    args = parser.parse_args(argv)

    if args.command == "show":
        only = next((k for k in ("open", "proposals", "best", "no_proposal", "first_pass") if getattr(args, k)), None)
        show(load(Path(args.review)), only)
        return 0
    if args.command == "confirm":
        groups = load(Path(args.review))
        n = 0
        for g in groups:
            weak = g.proposal_reason == BEST_REASON
            first_pass = g.proposal_reason.startswith(FIRST_PASS)
            if not g.proposal or g.decision or (weak and not args.include_best):
                continue
            if first_pass and not args.include_first_pass:
                continue
            g.decision, n = "ok", n + 1
        save(groups, Path(args.review))
        print(f"{n} proposals confirmed", file=sys.stderr)
        return 0

    token = scix_token()
    if args.command == "build":
        if not token:
            print("catalogue.py: needs a SciX API token (SCIX_API_TOKEN)", file=sys.stderr)
            return 1
        groups = load(Path(args.review))
        open_ = sum(1 for g in groups if not g.decision)
        n_entries, n_rows = build(groups, HttpClient(token), Path(args.out_dir))
        print(f"{n_entries} entries, {n_rows} lookup rows written to {args.out_dir}; {open_} groups still open",
              file=sys.stderr)  # fmt: skip
        return 0

    if not token:
        print("catalogue.py: needs a SciX API token (SCIX_API_TOKEN)", file=sys.stderr)
        return 1
    out = Path(args.output)
    groups, datas = collect(args.files, args.database)
    if not args.no_hints:
        add_hints(groups, datas, args.database, args.lamda_page, Path(args.cache_dir))
    search(groups, HttpClient(token))
    merged = merge_review(groups, out)
    n_prop = sum(1 for g in groups if g.proposal)
    print(f"{len(groups)} groups searched ({n_prop} with a proposal), {len(merged)} in {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
