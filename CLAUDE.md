# pylamda

Python toolkit for files in the LAMDA format (Leiden Atomic and Molecular Database):
read, write, verify, compare, combine, fix, and export references as BibTeX.

**Read `docs/design.md` first.** It holds the scope, the data model, the check codes, the
decisions (section 13) and the implementation steps (section 14).

## Guiding rule

Follow the LAMDA standard strictly (<https://home.strw.leidenuniv.nl/~moldata/molformat.html>):

- collision partners are identified by their ID (1–7) only;
- deviations are reported as issues, never silently interpreted.

## Environment and commands

- conda environment: `PYLAMDA`, defined in `environment.yml` (Python 3.12; numpy,
  matplotlib, requests, pytest, ruff, gfortran, poetry, Sphinx and extensions).
  pylamda itself is installed with `pip install --no-deps -e .`, so that conda keeps all
  dependencies. Never run `pip install -e ".[dev]"` in the conda environment: pip then
  installs a second copy of packages that conda already has (for example ruff).
- Without conda: a venv and `pip install -e ".[dev]"` (also used by CI). Only `gfortran`
  is missing there; the test that needs it is skipped.
- Supported Python: ≥ 3.11.
- Packaging: Poetry. `poetry.toml` sets `virtualenvs.create = false`, so Poetry installs
  into the active environment (conda or venv).
- Dependencies: numpy and requests. Extras: `plot`, `dev`. Keep `environment.yml` in line
  with the dependencies and the `dev` extra in `pyproject.toml`.

```bash
conda activate PYLAMDA
pytest
pytest --cov            # with coverage
ruff format . && ruff check .
sphinx-build -W -b html docs docs/_build/html
```

## Origin of the code

The first code comes from astrolinedata (`/home/rab/git/astrolinedata`,
`src/astrolinedata/formats/lamda.py`, `references.py`, `backends/nlte.py`) and from
`/home/rab/MODELS/TEST_LAMDA_FILES/lamda.py`. The history is not kept. The commit that
adds the code names the astrolinedata commit it comes from.

## Code style

- Type annotations everywhere; NumPy-style docstrings for all public functions.
- ruff for formatting and linting, line length 120 (set in `pyproject.toml`).
- Tests for every new function and every issue code (pytest, synthetic files).
