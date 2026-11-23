"""Command line: ``pylamda check``, ``compare``, ``bib``, ``fix`` and the combine commands
(``add-partner``, ``drop-partner``, ``scale-partner``).

Exit codes of ``pylamda check``: 0 = no errors, 1 = at least one issue of severity
``error``, 2 = at least one file cannot be read or is not a LAMDA file. Of ``pylamda compare`` (as ``diff``):
0 = no differences, 1 = differences, 2 = a file cannot be read. Of ``pylamda bib``: 0 = done,
1 = online resolution failed, 2 = a file cannot be read. Of ``pylamda fix``: 0 = done,
2 = the file cannot be read or a fix cannot be applied.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from collections.abc import Sequence

from . import __version__
from .bib import BibError, bibtex
from .combine import CombineError, add_partner, drop_partner, scale_partner
from .compare import compare
from .fix import OPTION_HINTS, FixError, fix, patch
from .issues import CODES, SEVERITIES, Issue, LamdaReadError, NotLamdaError
from .model import LamdaData
from .reader import read, read_string
from .verify import verify
from .writer import write

_NOT_LAMDA = "not a LAMDA file"


def _codes(values: list[str] | None) -> list[str]:
    """Code selectors from repeated, comma-separated options (``--select C02,C1``)."""
    return [code.strip() for value in values or [] for code in value.split(",") if code.strip()]


def _factor(text: str) -> float:
    """A scaling factor from the command line: a number or a fraction ``a/b`` (``1/1.37``)."""
    try:
        if "/" in text:
            num, den = text.split("/")
            return float(num) / float(den)
        return float(text)
    except (ValueError, ZeroDivisionError):
        raise argparse.ArgumentTypeError(f"not a number or a fraction a/b: {text!r}") from None


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pylamda",
        description="Toolkit for files in the LAMDA format (Leiden Atomic and Molecular Database).",
    )
    parser.add_argument("--version", action="version", version=f"pylamda {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True

    check = sub.add_parser(
        "check",
        help="check files against the LAMDA standard",
        description="Check LAMDA files against the standard. Exit code: 0 = no errors, 1 = errors, "
        "2 = a file cannot be read.",
    )
    check.add_argument("files", nargs="*", metavar="FILE", help="LAMDA files")
    check.add_argument(
        "--select", action="append", metavar="CODES", help="report only these codes (comma-separated prefixes: C,F01)"
    )
    check.add_argument("--ignore", action="append", metavar="CODES", help="do not report these codes (prefixes)")
    check.add_argument(
        "--error-on", action="append", metavar="CODES", help="report these codes (prefixes) as errors, e.g. C02"
    )
    check.add_argument("--no-physics", action="store_true", help="leave out the physics heuristics (= --ignore P)")
    check.add_argument(
        "--rtol", type=float, default=1e-4, help="relative tolerance of R04/R05, on top of rounding (default 1e-4)"
    )
    check.add_argument("--format", choices=("text", "json"), default="text", help="output format (default text)")
    check.add_argument(
        "--statistics", action="store_true", help="show the number of issues per code instead of each issue"
    )
    check.add_argument("--list-codes", action="store_true", help="list all check codes and exit")
    check.add_argument("--no-catalogue", action="store_true", help="F03: do not use the reference catalogue")

    cmp = sub.add_parser(
        "compare",
        help="compare two LAMDA files",
        description="Compare two LAMDA files: levels, lines, rates, critical densities. Exit code: "
        "0 = no differences, 1 = differences, 2 = a file cannot be read.",
    )
    cmp.add_argument("a", metavar="A", help="the first (old) file")
    cmp.add_argument("b", metavar="B", help="the second (new) file")
    cmp.add_argument(
        "--match",
        choices=("auto", "quanta", "energy", "index"),
        default="auto",
        help="how levels are matched (default auto: quanta, the rest by energy and g)",
    )
    cmp.add_argument(
        "--rtol", type=float, default=1e-6, help="relative change that counts in the exact comparison (default 1e-6)"
    )
    cmp.add_argument("--atol", type=float, default=1e-3, help="energy tolerance [cm-1] for matching (default 1e-3)")
    cmp.add_argument(
        "--limit",
        type=float,
        default=0.1,
        help="relative change that counts on the union grid, per level and per E_up range (default 0.1)",
    )
    cmp.add_argument("--format", choices=("text", "json"), default="text", help="output format (default text)")
    cmp.add_argument("--max-items", type=int, default=5, help="items named per list before '...' (default 5)")

    bib = sub.add_parser(
        "bib",
        help="BibTeX file of the references of LAMDA files",
        description="Write the references of LAMDA files as BibTeX. Offline: @misc entries with the DOIs, "
        "bibcodes and arXiv IDs. --online: resolved with SciX (needs the user's API token) and doi.org.",
    )
    bib.add_argument("files", nargs="+", metavar="FILE", help="LAMDA files")
    bib.add_argument("-o", "--output", metavar="FILE", help="the .bib file (default: standard output)")
    bib.add_argument("--online", action="store_true", help="resolve the identifiers with SciX and doi.org")
    bib.add_argument(
        "--scix-token-file",
        metavar="FILE",
        help="file with the SciX API token (default: SCIX_API_TOKEN, then ~/.config/pylamda/scix_token)",
    )
    bib.add_argument("--no-database-refs", action="store_true", help="leave out the LAMDA and EMAA papers")
    bib.add_argument("--keys", metavar="FILE", help="JSON file identifier -> key; read and updated (stable keys)")
    bib.add_argument("--no-cache", action="store_true", help="do not use the cache of resolved references")
    bib.add_argument("--no-catalogue", action="store_true", help="do not look up free-text references in the catalogue")
    bib.add_argument(
        "--search",
        action="store_true",
        help="search SciX for free-text references not in the catalogue; candidates go to the review file only",
    )
    bib.add_argument(
        "--review",
        metavar="FILE",
        default="references_review.json",
        help="review file of --search (default references_review.json; decisions in it are kept)",
    )
    bib.add_argument("--no-hints", action="store_true", help="--search: do not read the LAMDA/EMAA molecule pages")

    fx = sub.add_parser(
        "fix",
        help="fix a LAMDA file",
        description="Fix a LAMDA file: unambiguous issues by default, the others only with an option. "
        "-o writes the normalised layout; --in-place patches only the changed lines.",
    )
    fx.add_argument("file", metavar="FILE", help="the LAMDA file")
    out = fx.add_mutually_exclusive_group(required=True)
    out.add_argument("-o", "--output", metavar="OUT", help="write the fixed file to OUT (normalised layout)")
    out.add_argument("--in-place", action="store_true", help="change FILE itself (patch mode)")
    out.add_argument("--dry-run", action="store_true", help="show the changes (patch mode) as a diff; write nothing")
    mode = fx.add_mutually_exclusive_group()
    mode.add_argument("--patch", action="store_true", help="patch mode also with -o: only the changed lines")
    mode.add_argument("--normalise", action="store_true", help="normalised layout also with --in-place")
    fx.add_argument("--trust", choices=("id", "name"), help="C02/C12: correct the partner name (id) or the ID (name)")
    fx.add_argument("--no-sort", action="store_true", help="L01: do not sort the levels by energy")
    fx.add_argument(
        "--inverted", choices=("swap", "chianti", "balance"), help="C08: how to fix rate rows with E(up) < E(low)"
    )
    fx.add_argument("--no-missing-columns", action="store_true", help="R04/R05: do not compute missing columns")
    fx.add_argument("--recompute", action="store_true", help="R04/R05: also replace differing frequencies and E_up")
    fx.add_argument(
        "--extra-blocks", choices=("include", "comment"), help="S03: count extra blocks in or comment them out"
    )
    fx.add_argument(
        "--ref",
        action="append",
        metavar="TARGET=ID",
        help="append an identifier to a reference: TARGET levels, dataset or a partner; ID a DOI, bibcode or arXiv ID",
    )
    fx.add_argument(
        "--add-refs",
        nargs="?",
        const=True,
        default=False,
        metavar="REVIEW",
        help="append the SciX links of free-text references: from the reference catalogue, or from the decisions "
        "in a review file (pylamda bib --search)",
    )

    add = sub.add_parser(
        "add-partner",
        help="add the collision rates of one partner from another file",
        description="Add the collision rates of one partner from SOURCE to TARGET (same species): the levels are "
        "matched, the transitions between matched levels are copied unchanged. Writes the normalised layout.",
    )
    add.add_argument("target", metavar="TARGET", help="the file to extend")
    add.add_argument("source", metavar="SOURCE", help="the file with the partner")
    add.add_argument("partner", metavar="PARTNER", help="partner name (e, p-H2, ...)")
    add.add_argument("-o", "--output", required=True, metavar="OUT", help="the new file (not TARGET)")
    add.add_argument("--match", choices=("auto", "quanta", "energy", "index"), default="auto", help="level matching")
    add.add_argument("--atol", type=float, default=0.1, help="energy tolerance [cm-1] (default 0.1)")
    add.add_argument("--rtol-e", type=float, default=1e-4, help="relative energy tolerance (default 1e-4)")
    add.add_argument("--replace", action="store_true", help="replace the partner if TARGET has it")
    add.add_argument("--block", type=int, help="the block of SOURCE, if it has several with the partner's ID")
    add.add_argument("--source-label", help="description of the source for the note line")

    drop = sub.add_parser(
        "drop-partner",
        help="remove the collision data of one partner",
        description="Remove the collision block of one partner. Writes the normalised layout.",
    )
    drop.add_argument("file", metavar="FILE")
    drop.add_argument("partner", metavar="PARTNER", help="partner name (by ID), or 'ID <n>'")
    drop.add_argument("-o", "--output", required=True, metavar="OUT", help="the new file (not FILE)")
    drop.add_argument("--block", type=int, help="the block, if the file has several with the partner's ID")

    scale = sub.add_parser(
        "scale-partner",
        help="multiply the collision rates of one partner by a factor",
        description="Multiply the collision rates of one partner by FACTOR: in place, as another partner (--to), "
        "or as an added partner (--to with --keep). Writes the normalised layout and shows the C checks of the "
        "result.",
    )
    scale.add_argument("file", metavar="FILE")
    scale.add_argument("partner", metavar="PARTNER", help="the partner to scale (by ID)")
    scale.add_argument("factor", metavar="FACTOR", type=_factor, help="a number or a fraction a/b, e.g. 1/1.37")
    scale.add_argument("-o", "--output", required=True, metavar="OUT", help="the new file (not FILE)")
    scale.add_argument("--to", metavar="NAME", help="the partner of the scaled block (default: PARTNER)")
    scale.add_argument("--keep", action="store_true", help="keep the PARTNER block and add the scaled block")
    scale.add_argument("--block", type=int, help="the block, if the file has several with the partner's ID")
    scale.add_argument("--replace", action="store_true", help="replace the block of --to if the file has it")
    scale.add_argument("--reference", metavar="TEXT", help="new reference text of the scaled block")
    scale.add_argument(
        "--digits", type=int, help="significant digits of the scaled rates (default: as written in FILE)"
    )
    return parser


def _same_file(a: str, b: str) -> bool:
    from pathlib import Path

    return Path(a).resolve() == Path(b).resolve()


def _combine(args: argparse.Namespace) -> int:
    inputs = [args.target, args.source] if args.command == "add-partner" else [args.file]
    warnings: list[str] = []
    if any(_same_file(args.output, f) for f in inputs):
        print(f"pylamda {args.command}: the output must be a new file, not {args.output}", file=sys.stderr)
        return 2
    try:
        if args.command == "add-partner":
            data, record = add_partner(
                read(args.target),
                read(args.source),
                args.partner,
                match=args.match,
                atol=args.atol,
                rtol_e=args.rtol_e,
                replace_existing=args.replace,
                source_label=args.source_label,
                block=args.block,
            )
            warnings = record.warnings
        elif args.command == "scale-partner":
            data, record = scale_partner(
                read(args.file),
                args.partner,
                args.factor,
                to=args.to,
                keep=args.keep,
                block=args.block,
                replace_existing=args.replace,
                reference=args.reference,
                digits=args.digits,
            )
            warnings = record.warnings
        else:
            data, record = drop_partner(read(args.file), args.partner, block=args.block)
    except (LamdaReadError, CombineError) as exc:
        print(f"pylamda {args.command}: error: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"pylamda {args.command}: cannot open: {exc.strerror or exc}", file=sys.stderr)
        return 2
    for w in warnings:
        print(f"pylamda {args.command}: warning: {w}", file=sys.stderr)
    write(data, args.output)
    print(data.notes[-1].lstrip("! "), file=sys.stderr)
    if args.command == "scale-partner":
        for issue in verify(args.output, select=["C"]):
            print(_issue_line(args.output, issue), file=sys.stderr)
    return 0


def _list_codes() -> None:
    for code, (severity, description) in CODES.items():
        print(f"{code}  {severity:<8s} {description}")


def _issue_line(path: str, issue: Issue) -> str:
    where = f"{path}:{issue.line}" if issue.line is not None else path
    return f"{where}: {issue.code} {issue.severity}: {issue.message}"


def _check(args: argparse.Namespace) -> int:
    if args.list_codes:
        _list_codes()
        return 0
    if not args.files:
        print("pylamda check: no files given", file=sys.stderr)
        return 2
    results = []  # (path, readable, error text, issues)
    for path in args.files:
        try:
            issues = verify(
                path,
                physics=not args.no_physics,
                select=_codes(args.select) or None,
                ignore=_codes(args.ignore),
                error_on=_codes(args.error_on),
                rtol=args.rtol,
                catalogue=not args.no_catalogue,
            )
            results.append((path, True, None, issues))
        except NotLamdaError as exc:  # one statement, not the structure issues of a file that is not LAMDA
            results.append((path, False, f"{_NOT_LAMDA}: {exc.reason}", []))
        except LamdaReadError as exc:
            results.append((path, False, "not readable", exc.issues))
        except OSError as exc:
            results.append((path, False, f"cannot open: {exc.strerror or exc}", []))

    if args.format == "json":
        out = [
            {"path": p, "readable": ok, "error": err, "issues": [i.as_dict() for i in issues]}
            for p, ok, err, issues in results
        ]
        print(json.dumps(out, indent=2, ensure_ascii=False))
    elif args.statistics:
        _print_statistics(results)
    else:
        for path, _, err, issues in results:
            for issue in issues:
                print(_issue_line(path, issue))
            if err:
                print(f"{path}: {err}")
        _print_summary(results)

    if any(not ok for _, ok, _, _ in results):
        return 2
    return 1 if any(i.severity == "error" for *_, issues in results for i in issues) else 0


def _print_summary(results: list) -> None:
    counts = Counter(i.severity for *_, issues in results for i in issues)
    n = sum(counts.values())
    parts = ", ".join(f"{counts[s]} {s}{'s' if counts[s] != 1 else ''}" for s in SEVERITIES if counts[s])
    not_lamda = sum(1 for _, ok, err, _ in results if not ok and err.startswith(_NOT_LAMDA))
    unreadable = sum(1 for _, ok, _, _ in results if not ok) - not_lamda
    files = f"{len(results)} file{'s' if len(results) != 1 else ''}"
    text = f"Found {n} issue{'s' if n != 1 else ''}" + (f" ({parts})" if parts else "") + f" in {files}"
    if not_lamda:
        text += f"; {not_lamda} not {'a LAMDA file' if not_lamda == 1 else 'LAMDA files'}"
    if unreadable:
        text += f"; {unreadable} not readable"
    print(text + ".")


def _print_statistics(results: list) -> None:
    issues_per_code: Counter = Counter()
    files_per_code: dict[str, set] = {}
    severities: dict[str, set] = {}
    for path, _, _, issues in results:
        for i in issues:
            issues_per_code[i.code] += 1
            files_per_code.setdefault(i.code, set()).add(path)
            severities.setdefault(i.code, set()).add(i.severity)
    for code in sorted(issues_per_code):
        severity = "/".join(s for s in SEVERITIES if s in severities[code])
        print(
            f"{issues_per_code[code]:6d}  {code}  {severity:<13s} in {len(files_per_code[code]):4d} files  "
            f"{CODES[code][1]}"
        )
    _print_summary(results)


def _compare(args: argparse.Namespace) -> int:
    try:
        diff = compare(args.a, args.b, match=args.match, rtol=args.rtol, atol=args.atol, limit=args.limit)
    except LamdaReadError as exc:
        print(f"pylamda compare: error: {exc}", file=sys.stderr)
        return 2
    except (OSError, ValueError) as exc:
        print(f"pylamda compare: error: {exc}", file=sys.stderr)
        return 2
    if args.format == "json":
        print(json.dumps(diff.as_dict(), indent=2, ensure_ascii=False))
    elif diff:
        print("\n".join(diff.summary(max_items=args.max_items)))
    else:
        print("no differences")
    return 1 if diff else 0


def _bib_search(args: argparse.Namespace, datas: list, token: str | None) -> int:
    """``pylamda bib --search``: SciX candidates into the review file; returns an exit status."""
    import os
    from pathlib import Path

    from .bib import HttpClient, scix_token
    from .catalogue import load_catalogue
    from .refsearch import merge_review, search_references

    token = scix_token(token)
    if not token:
        print("pylamda bib: --search needs a SciX API token (SCIX_API_TOKEN)", file=sys.stderr)
        return 1
    cache = Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "pylamda" / "catalogue"
    try:
        groups = search_references(
            datas,
            client=HttpClient(token),
            hints=not args.no_hints,
            cache_dir=cache,
            catalogue=None if args.no_catalogue else load_catalogue(),
        )
    except BibError as exc:
        print(f"pylamda bib: error: {exc}", file=sys.stderr)
        return 1
    merged = merge_review(groups, Path(args.review))
    print(
        f"pylamda bib: --search: {len(groups)} free-text citations searched; candidates in {args.review} "
        f"({len(merged)} groups); nothing is added before you decide there",
        file=sys.stderr,
    )
    for g in groups:
        if g.proposal:
            what = f"proposal {g.proposal} ({g.proposal_reason})"
        elif g.candidates:
            c = g.candidates[0]
            weak = " (weak)" if c["score"] < 0.3 else ""
            what = f"best candidate {c['bibcode']}, score {c['score']:.2f}{weak}"
        else:
            what = "no candidate"
        print(f"  {g.molecule}: {g.segment[:60]!r} -> {what}", file=sys.stderr)
    if groups:
        print(f"  show them: python -m pylamda.refsearch {args.review}", file=sys.stderr)
    return 0


def _bib(args: argparse.Namespace) -> int:
    datas: list[LamdaData] = []
    status = 0
    for path in args.files:
        try:
            datas.append(read(path))
        except LamdaReadError as exc:
            print(f"pylamda bib: {path}: not readable:\n{exc}", file=sys.stderr)
            status = 2
        except OSError as exc:
            print(f"pylamda bib: {path}: cannot open: {exc.strerror or exc}", file=sys.stderr)
            status = 2
    token = None
    if args.scix_token_file:
        from pathlib import Path

        lines = Path(args.scix_token_file).read_text(encoding="utf-8").splitlines()
        token = lines[0].strip() if lines else None
    keys_path = None
    keys = None
    if args.keys:
        from pathlib import Path

        keys_path = Path(args.keys)
        keys = json.loads(keys_path.read_text(encoding="utf-8")) if keys_path.is_file() else {}
    try:
        result = bibtex(
            datas,
            online=args.online,
            token=token,
            database_refs=not args.no_database_refs,
            keys=keys,
            cache=not args.no_cache,
            catalogue=not args.no_catalogue,
        )
    except BibError as exc:
        print(f"pylamda bib: error: {exc}", file=sys.stderr)
        return 1
    if args.output:
        result.write(args.output)
    else:
        print(result.text(), end="")
    if keys_path is not None:
        keys_path.write_text(json.dumps(result.keys, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    for w in result.warnings:
        print(f"pylamda bib: warning: {w}", file=sys.stderr)
    if result.free_text:
        print(
            f"pylamda bib: {len(result.free_text)} references without DOI, bibcode or arXiv ID (not in the "
            "BibTeX file):",
            file=sys.stderr,
        )
        for r in result.free_text:
            print(f"  {r.path}: {r.use}: {r.text}", file=sys.stderr)
    if args.search:
        status = max(status, _bib_search(args, datas, token))
    counts = {}
    for e in result.entries:
        counts[e.source] = counts.get(e.source, 0) + 1
    where = f" to {args.output}" if args.output else ""
    parts = ", ".join(f"{n} {s}" for s, n in counts.items())
    print(f"pylamda bib: {len(result.entries)} entries{where} ({parts})", file=sys.stderr)
    return status


def _fix(args: argparse.Namespace) -> int:
    from pathlib import Path

    options = {
        "trust": args.trust,
        "sort_levels": not args.no_sort,
        "inverted": args.inverted,
        "missing_columns": not args.no_missing_columns,
        "recompute": args.recompute,
        "extra_blocks": args.extra_blocks,
        "refs": args.ref,
        "add_refs": args.add_refs,
    }
    keep_layout = args.dry_run or args.patch or (args.in_place and not args.normalise)
    try:
        if keep_layout:
            text, log = patch(args.file, **options)
        else:
            data, log = fix(args.file, **options)
    except (LamdaReadError, FixError) as exc:
        print(f"pylamda fix: error: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"pylamda fix: {args.file}: cannot open: {exc.strerror or exc}", file=sys.stderr)
        return 2
    if args.dry_run:
        import difflib

        old = Path(args.file).read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
        sys.stdout.writelines(difflib.unified_diff(old, text.splitlines(keepends=True), args.file, "fixed"))
    else:
        target = args.file if args.in_place else args.output
        if keep_layout:
            Path(target).write_bytes(text.encode("utf-8"))
        else:
            write(data, target)
    for r in log:
        print(f"fixed: {r.code} {r.message}", file=sys.stderr)
    if not log:
        print("pylamda fix: nothing to fix with these options", file=sys.stderr)
    # fixable issues that need an option
    try:
        remaining = verify(read_string(text) if keep_layout else data)
    except LamdaReadError:  # pragma: no cover (the fixed data read back)
        remaining = []
    hints = sorted({i.code for i in remaining if i.code in OPTION_HINTS})
    for code in hints:
        print(f"pylamda fix: {code} remains; it needs {OPTION_HINTS[code]}", file=sys.stderr)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command line.

    Parameters
    ----------
    argv : sequence of str or None, optional
        The arguments without the program name. Default: ``sys.argv[1:]``.

    Returns
    -------
    int
        The exit code.
    """
    args = _parser().parse_args(argv)
    if args.command == "check":
        return _check(args)
    if args.command == "compare":
        return _compare(args)
    if args.command == "bib":
        return _bib(args)
    if args.command == "fix":
        return _fix(args)
    if args.command in ("add-partner", "drop-partner", "scale-partner"):
        return _combine(args)
    return 2  # pragma: no cover (argparse requires a command)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
