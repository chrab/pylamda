"""Corpus tests on real files outside the repository (``pytest --corpus PATH``)."""

from __future__ import annotations

import pytest
from helpers import same_data

from pylamda import LamdaReadError, read, write

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
