# pylamda

Python toolkit for files in the LAMDA format (the
[Leiden Atomic and Molecular Database](https://home.strw.leidenuniv.nl/~moldata/)):
energy levels, Einstein A coefficients and collision rates, as used by LAMDA, EMAA,
RADEX and ProDiMo.

pylamda follows the [LAMDA standard](https://home.strw.leidenuniv.nl/~moldata/molformat.html)
strictly: collision partners are identified by their ID only, and deviations are
reported as issues, never silently interpreted.

**Status:** early development. Reading and writing work. Verify, compare, combine, fix
and the BibTeX export follow (see `docs/design.md`, section 14).

```python
import pylamda

data = pylamda.read("co.dat")
for issue in data.issues:  # structure issues found by the reader
    print(issue)
pylamda.write(data, "co_normalised.dat")
```

## Installation

```bash
pip install .                 # numpy only
pip install ".[plot,bib]"     # with matplotlib and requests
pip install -e ".[dev]"       # development: tests, ruff, Sphinx
```

Python ≥ 3.11. License: MIT.
