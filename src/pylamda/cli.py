"""Command line: ``pylamda check FILE ...`` and ``pylamda compare A B``.

Exit codes of ``pylamda check``: 0 = no errors, 1 = at least one issue of severity
``error``, 2 = at least one file cannot be read. Of ``pylamda compare`` (as ``diff``):
0 = no differences, 1 = differences, 2 = a file cannot be read.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from collections.abc import Sequence

from . import __version__
from .compare import compare
from .issues import CODES, SEVERITIES, Issue, LamdaReadError
from .verify import verify


def _codes(values: list[str] | None) -> list[str]:
    """Code selectors from repeated, comma-separated options (``--select C02,C1``)."""
    return [code.strip() for value in values or [] for code in value.split(",") if code.strip()]


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
    return parser


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
            )
            results.append((path, True, None, issues))
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
    unreadable = sum(1 for _, ok, _, _ in results if not ok)
    files = f"{len(results)} file{'s' if len(results) != 1 else ''}"
    text = f"Found {n} issue{'s' if n != 1 else ''}" + (f" ({parts})" if parts else "") + f" in {files}"
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
    return 2  # pragma: no cover (argparse requires a command)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
