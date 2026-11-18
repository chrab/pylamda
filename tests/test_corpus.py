"""Corpus tests on real files outside the repository (``pytest --corpus PATH``).

``--corpus-snapshot FILE`` compares the number of issues per file and code with a JSON
snapshot (design document, section 12), so that changes of the checks are visible. The
snapshot is written if the file does not exist, or with ``--update-snapshot``.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest
from conftest import corpus_files
from helpers import same_data

from pylamda import CODES, SEVERITIES, LamdaReadError, bibtex, compare, read, verify, write
from pylamda.bib import split_entries

pytestmark = pytest.mark.corpus


def test_corpus_round_trip(corpus_file, tmp_path):
    """Every file reads; write -> read gives the same values; writing again gives the same text."""
    try:
        d = read(corpus_file)
    except LamdaReadError as exc:
        pytest.fail(f"not readable:\n{exc}")
    write(d, tmp_path / "a.dat")
    d2 = read(tmp_path / "a.dat")
    assert same_data(d, d2)
    assert d2.issues == []
    write(d2, tmp_path / "b.dat")
    assert (tmp_path / "a.dat").read_text() == (tmp_path / "b.dat").read_text()
    diff = compare(d, d2)
    assert not diff, diff.summary()


def test_corpus_verify(corpus_file):
    """verify() runs on every file; codes, severities and line numbers are valid."""
    n_lines = len(Path(corpus_file).read_text(encoding="utf-8", errors="replace").splitlines())
    for issue in verify(corpus_file):
        assert issue.code in CODES and issue.severity in SEVERITIES and issue.count >= 1
        assert issue.line is None or 1 <= issue.line <= n_lines, issue


def test_corpus_bib(corpus_file):
    """The offline BibTeX export works on every file and gives one parsable entry per BibEntry."""
    result = bibtex([corpus_file])
    assert len(split_entries(result.text())) == len(result.entries)
    assert len({e.key for e in result.entries}) == len(result.entries)


def _counts(path: Path) -> dict[str, int]:
    counts: Counter = Counter()
    for issue in verify(path):
        counts[issue.code] += issue.count
    return dict(sorted(counts.items()))


def test_corpus_snapshot(request):
    """The issue counts per file and code agree with the snapshot (--corpus-snapshot)."""
    snapshot = request.config.getoption("--corpus-snapshot")
    files = corpus_files(request.config)
    if not snapshot or not files:
        pytest.skip("needs --corpus and --corpus-snapshot")
    current = {str(path): _counts(path) for path in files}
    snapshot = Path(snapshot)
    if request.config.getoption("--update-snapshot") or not snapshot.exists():
        snapshot.write_text(json.dumps(current, indent=1) + "\n")
        return
    old = json.loads(snapshot.read_text())
    changed = []
    for path in sorted(set(old) | set(current)):
        a, b = old.get(path), current.get(path)
        if a != b:
            codes = sorted(set(a or {}) | set(b or {}))
            diff = ", ".join(
                f"{c} {(a or {}).get(c, 0)}->{(b or {}).get(c, 0)}"
                for c in codes
                if (a or {}).get(c) != (b or {}).get(c)
            )
            changed.append(f"{path}: {diff if a and b else 'only in snapshot' if a else 'new file'}")
    assert not changed, "issue counts changed (rerun with --update-snapshot to accept):\n" + "\n".join(changed)
