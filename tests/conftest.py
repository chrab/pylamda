"""pytest options: ``--corpus`` runs the corpus tests on real files outside the repository."""

from __future__ import annotations

import glob
from pathlib import Path

import pytest

_SUFFIXES = (".dat", ".radex", ".txt")


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--corpus",
        action="append",
        default=[],
        metavar="PATH",
        help="directory, file or quoted glob pattern of LAMDA files for the corpus tests (repeatable); "
        f"a directory gives all its {', '.join(_SUFFIXES)} files",
    )
    parser.addoption(
        "--corpus-snapshot",
        metavar="FILE",
        help="JSON file with the issue counts per corpus file and code; compared, or written if missing",
    )
    parser.addoption("--update-snapshot", action="store_true", help="write the --corpus-snapshot file anew")


def corpus_files(config: pytest.Config) -> list[Path]:
    """The files given with ``--corpus``, sorted, without duplicates."""
    files: set[Path] = set()
    for spec in config.getoption("--corpus"):
        path = Path(spec).expanduser()
        if path.is_dir():
            files.update(p for p in path.iterdir() if p.is_file() and p.suffix in _SUFFIXES)
        elif path.is_file():
            files.add(path)
        else:
            matches = [Path(p) for p in glob.glob(str(path))]
            if not matches:
                raise pytest.UsageError(f"--corpus {spec}: no such file, directory or matching files")
            files.update(p for p in matches if p.is_file())
    return sorted(files)


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    if "corpus_file" in metafunc.fixturenames:
        files = corpus_files(metafunc.config)
        if files:
            metafunc.parametrize("corpus_file", files, ids=[str(p) for p in files])
        else:
            metafunc.parametrize(
                "corpus_file", [pytest.param(None, marks=pytest.mark.skip(reason="no --corpus given"))]
            )
