"""Writing: lossless round trip, layout, partner lines, and a RADEX-style Fortran read."""

from __future__ import annotations

import re
import shutil
import subprocess

import numpy as np
import pytest
from helpers import ALL_FILES, EMAA_CO, LAMDA, make_lamda, same_data

from pylamda import partner_name, read, write


@pytest.mark.parametrize("path", ALL_FILES, ids=lambda p: p.name)
def test_round_trip(path, tmp_path):
    """write -> read gives the same values; writing again gives the same text."""
    d = read(path)
    write(d, tmp_path / "a.dat")
    d2 = read(tmp_path / "a.dat")
    assert same_data(d, d2)
    names = [partner_name(c.text_partner) or c.text_partner if c.text_partner else c.partner for c in d.collisions]
    assert [c.text_partner for c in d2.collisions] == names
    assert d2.issues == []
    write(d2, tmp_path / "b.dat")
    assert (tmp_path / "a.dat").read_text() == (tmp_path / "b.dat").read_text()


def test_written_layout(tmp_path):
    """EMAA-style references on line 2, the NLEV line and the partner lines."""
    d = read(EMAA_CO).select_partners(["p-H2", "o-H2"])
    d.collisions = d.collisions[2:4]  # the real p-H2 and o-H2 blocks (the first two are H2O)
    d.source = "https://emaa.osug.fr/details/CO, Rotation 2023-03-25, retrieved 2026-11-15"
    write(d, tmp_path / "co.dat")
    lines = (tmp_path / "co.dat").read_text().splitlines()
    assert lines[1] == "CO -- EMAA_DATASET_DOI: 10.17178/EMAA_CO_Rotation_3d0e1a5c"
    assert lines[5].startswith("14 - ROTATION - Ref: CDMS (Oct. 2000)")
    partner_lines = [lines[i + 1] for i, line in enumerate(lines) if line == "!COLLISIONS BETWEEN"]
    assert partner_lines == [
        "2 CO + p-H2 -- Ref: Dagdigian MNRAS 514 2214 (2022) https://doi.org/10.1093/mnras/stac1430",
        "3 CO + o-H2 -- Ref: Dagdigian MNRAS 514 2214 (2022) https://doi.org/10.1093/mnras/stac1430",
    ]
    assert lines[-1].startswith("!SOURCE: https://emaa.osug.fr/details/CO")
    back = read(tmp_path / "co.dat")
    assert back.source == d.source and back.partners == ["p-H2", "o-H2"]
    # Einstein A always in exponential notation, and still lossless
    first = lines.index("!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K)") + 1
    a_texts = [line.split()[3] for line in lines[first : first + d.n_lines]]
    assert all(re.fullmatch(r"\d\.\d+e[+-]\d\d", t) for t in a_texts), a_texts[:3]
    np.testing.assert_array_equal(back.line_a, d.line_a)
    d.line_a[:2] = [1.94, 2.0]
    write(d, tmp_path / "a.dat")
    rows = (tmp_path / "a.dat").read_text().splitlines()[first : first + 2]
    assert [r.split()[3] for r in rows] == ["1.94e+00", "2.0e+00"]


def test_partner_lines(tmp_path):
    """The ID is written as it is; the name is the one in the text, else the standard name of the ID."""
    cat = read(LAMDA / "catom.dat")
    write(cat, tmp_path / "c.dat")
    text = (tmp_path / "c.dat").read_text()
    assert "\n5 C + H -- Ref: Launay & Roueff 1977\n" in text
    # EMAA's partners outside the standard keep their name and ID: nothing is decided silently
    write(read(EMAA_CO), tmp_path / "co.dat")
    text = (tmp_path / "co.dat").read_text()
    assert "\n2 CO + p-H2O -- Ref:" in text and "\n5 CO + CO -- Ref:" in text
    # name and ID disagree: both are kept as they are
    path = tmp_path / "x.dat"
    path.write_text(make_lamda(partner="1 X + para-H2 ! a reference"))
    write(read(path), tmp_path / "y.dat")
    back = read(tmp_path / "y.dat")
    assert back.collisions[0].partner_text == "1 X + p-H2 -- Ref: a reference"
    assert (back.collisions[0].partner, back.collisions[0].text_partner) == ("H2", "p-H2")
    # no name in the text: the standard name of the ID
    path.write_text(make_lamda(partner="3 X-oH2 from Yang et al. (2010)"))
    write(read(path), tmp_path / "y.dat")
    assert "\n3 X + o-H2 -- Ref: X-oH2 from Yang et al. (2010)\n" in (tmp_path / "y.dat").read_text()
    # an ID outside the standard and no name: no name written
    path.write_text(make_lamda(partner="9 something"))
    write(read(path), tmp_path / "y.dat")
    assert "\n9 X -- Ref: something\n" in (tmp_path / "y.dat").read_text()
    assert read(tmp_path / "y.dat").collisions[0].reference == "something"


def test_level_numbers_checked(tmp_path):
    """Transitions must refer to levels 1..NLEV; nothing is written otherwise."""
    co = read(LAMDA / "catom.dat")
    c = co.collisions[0]
    c.low[-1] = co.n_levels + 5
    with pytest.raises(ValueError, match="H collisional transitions: 1 transition"):
        write(co, tmp_path / "out.dat")
    assert not (tmp_path / "out.dat").exists()
    co = read(LAMDA / "catom.dat")
    co.line_up[0] = 0
    with pytest.raises(ValueError, match=r"radiative transitions: .* the first is transition 1 \(0 -> "):
        write(co, tmp_path / "out.dat")


def test_non_finite_values(tmp_path):
    path = tmp_path / "x.dat"
    path.write_text(make_lamda(rates="1 2 1 NaN 2e-11\n2 3 1 Inf 4e-11\n3 3 2 -inf 6e-11"))
    d = read(path)
    write(d, tmp_path / "y.dat")
    back = read(tmp_path / "y.dat")
    assert back.issues == []
    np.testing.assert_array_equal(back.collisions[0].rates[:, 0], [np.nan, np.inf, -np.inf])


def test_notes_written_as_comments(tmp_path):
    d = read(LAMDA / "oatom@lique.dat")  # the 6th block is kept as notes
    write(d, tmp_path / "o.dat")
    back = read(tmp_path / "o.dat")
    assert back.issues == [] and len(back.notes) == len(d.notes)
    assert all(line.startswith("!") or not line.strip() for line in back.notes)


RADEX_STYLE_READER = """
program radex_style
  ! reads a LAMDA file the way RADEX does (list-directed reads, fixed positions)
  implicit none
  character(len=200) :: dummy, name
  integer :: nlev, nlin, npart, ncol, ntemp, i, j, k, idx, iu, il, pid
  real(8) :: amass, e, g, a, f, eu, asum, rsum
  real(8), allocatable :: t(:), r(:)
  open(11, file='data.dat', status='old')
  read(11,*) dummy; read(11,*) name; read(11,*) dummy; read(11,*) amass
  read(11,*) dummy; read(11,*) nlev; read(11,*) dummy
  do i = 1, nlev
    read(11,*) idx, e, g
  end do
  read(11,*) dummy; read(11,*) nlin; read(11,*) dummy
  asum = 0
  do i = 1, nlin
    read(11,*) idx, iu, il, a, f, eu
    asum = asum + a
  end do
  read(11,*) dummy; read(11,*) npart
  rsum = 0
  do j = 1, npart
    read(11,*) dummy; read(11,*) pid; read(11,*) dummy; read(11,*) ncol
    read(11,*) dummy; read(11,*) ntemp; read(11,*) dummy
    allocate(t(ntemp), r(ntemp))
    read(11,*) t
    read(11,*) dummy
    do k = 1, ncol
      read(11,*) idx, iu, il, r
      rsum = rsum + sum(r)
    end do
    write(*,'(a,i0,1x,i0,1x,i0)') 'partner ', pid, ncol, ntemp
    deallocate(t, r)
  end do
  write(*,'(a,a)') 'name ', trim(name)
  write(*,'(a,i0,1x,i0,1x,es24.16,1x,es24.16)') 'sums ', nlev, nlin, asum, rsum
end program radex_style
"""


@pytest.mark.fortran
@pytest.mark.skipif(shutil.which("gfortran") is None, reason="needs gfortran on the PATH")
def test_written_file_reads_like_radex(tmp_path):
    """RADEX-style list-directed reads cope with the '-- Ref:' text and the notes."""
    d = read(EMAA_CO)
    d.collisions = d.collisions[2:4]
    d.notes = ["!NOTES", "!a note"]
    write(d, tmp_path / "data.dat")
    (tmp_path / "radex_style.f90").write_text(RADEX_STYLE_READER)
    subprocess.run(["gfortran", "-o", "radex_style", "radex_style.f90"], cwd=tmp_path, check=True, capture_output=True)
    out = subprocess.run(["./radex_style"], cwd=tmp_path, check=True, capture_output=True, text=True).stdout
    lines = out.splitlines()
    assert lines[:2] == [
        f"partner 2 {d.collisions[0].n_transitions} {len(d.collisions[0].temperatures)}",
        f"partner 3 {d.collisions[1].n_transitions} {len(d.collisions[1].temperatures)}",
    ]
    assert lines[2] == "name CO"
    _, nlev, nlin, asum, rsum = lines[3].split()
    assert (int(nlev), int(nlin)) == (14, 13)
    assert float(asum) == pytest.approx(d.line_a.sum(), rel=1e-14)
    assert float(rsum) == pytest.approx(sum(c.rates.sum() for c in d.collisions), rel=1e-14)
