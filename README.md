# pylamda

[![PyPI](https://img.shields.io/pypi/v/pylamda)](https://pypi.org/project/pylamda/)
[![Documentation](https://img.shields.io/readthedocs/pylamda)](https://pylamda.readthedocs.io/)
[![Coverage](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/chrab/pylamda/python-coverage-comment-action-data/endpoint.json)](https://github.com/chrab/pylamda/tree/python-coverage-comment-action-data)

Python toolkit for working with and on files in the LAMDA format (the
[Leiden Atomic and Molecular Database](https://home.strw.leidenuniv.nl/~moldata/)).

**Documentation:** <https://pylamda.readthedocs.io/>

## Main features

- reading, writing and manipulating (i.e. adding collisional partner) LAMDA files
- validating LAMDA files (i.e. do the follow the LAMDA standard)
- bibliography support (i.e. get all references for a file in bib format)

pylamda follows the [LAMDA standard](https://home.strw.leidenuniv.nl/~moldata/molformat.html)
strictly: collision partners are identified by their ID only, and deviations are
reported as issues, never silently interpreted.

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

## Installation

```bash
pip install .                 # numpy and requests
pip install ".[plot]"         # with matplotlib
pip install -e ".[dev]"       # development: tests, ruff, Sphinx
```

Python ≥ 3.11. License: MIT.
