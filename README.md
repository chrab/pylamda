# pylamda

[![Documentation](https://readthedocs.org/projects/pylamda/badge/?version=latest)](https://pylamda.readthedocs.io/)

Python toolkit for files in the LAMDA format (the
[Leiden Atomic and Molecular Database](https://home.strw.leidenuniv.nl/~moldata/)):
energy levels, Einstein A coefficients and collision rates, as used by LAMDA, EMAA,
RADEX and ProDiMo.

pylamda follows the [LAMDA standard](https://home.strw.leidenuniv.nl/~moldata/molformat.html)
strictly: collision partners are identified by their ID only, and deviations are
reported as issues, never silently interpreted.

**Status:** early development. Reading, writing, checking, comparing, fixing, combining
and the BibTeX export with a reference catalogue work. `cut_levels`, plots and an
online check of DOIs are planned for later (see `docs/design.md`, sections 14 and 15).

**Documentation:** <https://pylamda.readthedocs.io/>

```python
import pylamda

data = pylamda.read("co.dat")
for issue in pylamda.verify(data):  # all checks (see docs/checks.rst)
    print(issue)
pylamda.write(data, "co_normalised.dat")
```

```bash
pylamda check data/*.dat            # exit code 1 if there are errors
pylamda check --statistics data/*.dat
pylamda compare old.dat new.dat     # exit code 1 if they differ
pylamda bib data/*.dat -o refs.bib --online   # SciX needs SCIX_API_TOKEN
pylamda fix in.dat --in-place --trust id      # only the changed lines
pylamda add-partner o.dat o_chianti.dat e -o o_with_e.dat
```

## Database papers

`src/pylamda/data/database_refs.bib` holds the BibTeX entries of the main papers of the
LAMDA and EMAA databases (Schöier et al. 2005; van der Tak et al. 2020; Faure et al.
2025; the EMAA dataset DOI). `pylamda bib` adds them to every export, as the databases
ask. The file is plain BibTeX and can be used on its own.

## Installation

```bash
pip install .                 # numpy only
pip install ".[plot,bib]"     # with matplotlib and requests
pip install -e ".[dev]"       # development: tests, ruff, Sphinx
```

Python ≥ 3.11. License: MIT.
